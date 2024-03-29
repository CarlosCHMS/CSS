#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include"input.h"
#include"mesh.h"

typedef struct{

    double p;
    double T;
    double mach;
    double nx;
    double ny;  
    
    double Uin[4];  

} CONDITION;

typedef struct{

    char* down;
    char* up;
    char* left;
    char* right;

    int Ndown;
    int Nup;
    int Nleft;
    int Nright;

} BOUNDARY;


typedef struct {

    int Nrow;
    int Ncol;
    int pOutFlag;
    int MUSCL;
    int flux;
    int stages;

    double Rgas;
    double gamma;
    double k4;
    double dt;
    double pout; 
    double eFix;       
    double e; 
    double k; 
    double res0[4];
    double res[4];
            
    double ***U;
    double ***R;
    double ***Uaux;        
    
    CONDITION* inlet;
        
    MESH* mesh;
    
    BOUNDARY* bc;

} SOLVER;


int boundaryChoice(char* s)
{
    int ans;

    if(strcmp(s, "wall") == 0)
    {
        ans = 0;
    }
    else if(strcmp(s, "inlet") == 0)
    {
        ans = 1;
    }
    else if(strcmp(s, "outlet") == 0)
    {
        ans = 2;
    }
    else if(strcmp(s, "wallp") == 0)
    {
        ans = 3;
    }

    
    return ans;

}

int fluxChoice(char* s)
{
    int ans;

    if(strcmp(s, "ROE") == 0)
    {
        ans = 0;
    }
    else if(strcmp(s, "AUSMD") == 0)
    {
        ans = 1;
    }
    
    return ans;

}

void boundarySet(BOUNDARY* bc)
{

    bc->Ndown = boundaryChoice(bc->down);
    bc->Nup = boundaryChoice(bc->up);    
    bc->Nleft = boundaryChoice(bc->left);
    bc->Nright = boundaryChoice(bc->right);    
    
}

CONDITION* conditionInit(double p, double T, double mach, double nx, double ny)
{
  
    CONDITION* cond = malloc(sizeof(CONDITION));

    cond->p = p;
    cond->T = T;
    cond->mach = mach;
    cond->nx = nx;
    cond->ny = ny;
    
    return cond;

}

void conditionState(CONDITION* cond, SOLVER* solver)
{

    double r, u, v, E, c;
    
    r = cond->p/(solver->Rgas*cond->T);
    c = sqrt(solver->gamma*solver->Rgas*cond->T);
    u = cond->nx*cond->mach*c;
    v = cond->ny*cond->mach*c;
    E = (solver->Rgas*cond->T)/(solver->gamma-1) + (u*u + v*v)/2;
    
    cond->Uin[0] = r;
    cond->Uin[1] = r*u;
    cond->Uin[2] = r*v;
    cond->Uin[3] = r*E;    
    

}

double conditionVref(CONDITION* cond, SOLVER* solver)
{
      
    double c = sqrt(solver->gamma*solver->Rgas*cond->T);
    double Vref = (cond->mach+1)*c;
        
    return Vref;

}

void solverAllocate(SOLVER* solver)
{
    solver->U = (double***)malloc(4*sizeof(double**));
    for(int jj=0; jj<4; jj++)
    {
        solver->U[jj] = (double**)malloc(solver->Nrow*sizeof(double*));
        for(int ii=0; ii<solver->Nrow; ii++)
        {
            solver->U[jj][ii] = (double*)malloc(solver->Ncol*sizeof(double));
        }
    }
    
    solver->Uaux = (double***)malloc(4*sizeof(double**));
    for(int jj=0; jj<4; jj++)
    {
        solver->Uaux[jj] = (double**)malloc(solver->Nrow*sizeof(double*));
        for(int ii=0; ii<solver->Nrow; ii++)
        {
            solver->Uaux[jj][ii] = (double*)malloc(solver->Ncol*sizeof(double));
        }
    }
    
    solver->R = (double***)malloc(4*sizeof(double**));
    for(int jj=0; jj<4; jj++)
    {
        solver->R[jj] = (double**)malloc(solver->Nrow*sizeof(double*));
        for(int ii=0; ii<solver->Nrow; ii++)
        {
            solver->R[jj][ii] = (double*)malloc(solver->Ncol*sizeof(double));
        }
    }
        
}

void solverInitU(SOLVER* solver, CONDITION* inside)
{

    conditionState(inside, solver);

    for(int kk=0; kk<4; kk++)
    {
        for(int ii=0; ii<solver->Nrow; ii++)
        {
            for(int jj=0; jj<solver->Ncol; jj++)
            {
                solver->U[kk][ii][jj] = inside->Uin[kk];
            }

        }
    }
}

void solverInitUTube(SOLVER* solver, CONDITION* inside1, CONDITION* inside2, double xm)
{

    double x;

    conditionState(inside1, solver);
    conditionState(inside2, solver);

    for(int kk=0; kk<4; kk++)
    {
        for(int ii=0; ii<solver->Nrow; ii++)
        {
            for(int jj=0; jj<solver->Ncol; jj++)
            {
                x = 0.5*(solver->mesh->x[ii+1][jj] + solver->mesh->x[ii][jj]);
                if(x<xm)
                {
                    solver->U[kk][ii][jj] = inside1->Uin[kk];
                }
                else
                {
                    solver->U[kk][ii][jj] = inside2->Uin[kk];
                }
            }

        }
    }
}

void solverResetR(SOLVER* solver)
{
    for(int kk=0; kk<4; kk++)
    {
        for(int ii=0; ii<solver->Nrow; ii++)
        {
            for(int jj=0; jj<solver->Ncol; jj++)
            {
                solver->R[kk][ii][jj] = 0.0;                
            }
        }
    }
}

void solverMatFree(SOLVER* solver, double** M, int Nrow)
{
    for(int ii=0; ii<Nrow; ii++)
    {
        free(M[ii]);
    }
    free(M);
}

void solverFree(SOLVER* solver)
{

    solverMatFree(solver, solver->mesh->x, solver->mesh->Nrow);
    solverMatFree(solver, solver->mesh->y, solver->mesh->Nrow);
    for(int ii=0; ii<4; ii++)
    {
        solverMatFree(solver, solver->U[ii], solver->Nrow);
        solverMatFree(solver, solver->Uaux[ii], solver->Nrow);        
        solverMatFree(solver, solver->R[ii], solver->Nrow);        
    }

}

double solverCalcP(SOLVER* solver, double*** U, int ii, int jj)
{

	double u = U[1][ii][jj]/U[0][ii][jj];
	double v = U[2][ii][jj]/U[0][ii][jj];
    
	return (solver->gamma - 1)*(U[3][ii][jj] - 0.5*(u*u + v*v)*U[0][ii][jj]);
    
}

void solverCalcVel(SOLVER* solver, double*** U, int ii, int jj, double* u, double* v, double* c)
{
    
    double E = U[3][ii][jj]/U[0][ii][jj];
    double aux;

    *u = U[1][ii][jj]/U[0][ii][jj];
    *v = U[2][ii][jj]/U[0][ii][jj];
    
    aux = (E - ((*u)*(*u) + (*v)*(*v))/2);
    aux *= solver->gamma - 1;
    *c = sqrt(aux*solver->gamma);
    
}

void entropyFix(SOLVER* solver, double *l)
{

    // Harten Hyman entropy fix
    if(*l < solver->eFix & *l > -solver->eFix )
    {
        *l = 0.5*(*l * *l/solver->eFix + solver->eFix);
    }

}

void solverFluxRoe(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L, 
               double U0R, double U1R, double U2R, double U3R,
	           double* f)
{

    /*
    Based on: P. L. ROE, Riemann Solvers, Parameter Vectors, and Difference Schemes, (1981)
    */
      
	double uL = U1L/U0L;
	double vL = U2L/U0L;
    double pL = (solver->gamma - 1)*(U3L - (uL*uL + vL*vL)*U0L/2);
    double HL = (U3L + pL)/U0L;

	double uR = U1R/U0R;
	double vR = U2R/U0R;
    double pR = (solver->gamma - 1)*(U3R - (uR*uR + vR*vR)*U0R/2);
    double HR = (U3R + pR)/U0R;

    // Mean values calculation
	double rqL = sqrt(U0L);
    double rqR = sqrt(U0R);

	double ub = (rqL*uL + rqR*uR)/(rqL + rqR);
	double vb = (rqL*vL + rqR*vR)/(rqL + rqR);
	double Hb  = (rqL*HL + rqR*HR)/(rqL + rqR);
	double ab  = sqrt((solver->gamma-1) * (Hb - (ub*ub + vb*vb)/2));

    // Eigenvalues
	double l1 = ub - ab;
	double l2 = ub;
	double l4 = ub;
	double l5 = ub + ab;	

    // Eigenvectors
	double e1[4] = {1.0, ub-ab, vb, Hb-ub*ab};
	double e2v[4] = {0.0, 0.0, 1.0, vb};
	double e4[4] = {1.0, ub, vb, 0.5 * (ub*ub + vb*vb)};
	double e5[4] = {1.0, ub+ab, vb, Hb + ub*ab};

    // Diferences
	double d1 = U0R - U0L;
	double d2 = U1R - U1L;
	double d3 = U2R - U2L;
	double d5 = U3R - U3L;

    // Projections
    double a4 = (Hb - (ub*ub + vb*vb))*d1 + ub*d2 + vb*d3 - d5;
    a4 /= (ab*ab)/(solver->gamma-1);    
    double a2v = d3 - d1*vb;
    double a5 = ((d1 - a4) + (d2 - ub*d1)/ab)*0.5;
    double a1 = ((d1 - a4) - (d2 - ub*d1)/ab)*0.5;

    // Fluxes
	double fL[4] = {U1L, U1L*uL + pL, U1L*vL, uL*(U3L + pL)};
	double fR[4] = {U1R, U1R*uR + pR, U1R*vR, uR*(U3R + pR)};

    // Entropy fix
    entropyFix(solver, &l1);
    entropyFix(solver, &l2);
    entropyFix(solver, &l4);
    entropyFix(solver, &l5);    

	
	for (int ii = 0; ii < 4; ++ii) {
		f[ii] = 0.5 * (fR[ii] + fL[ii] - a1*fabs(l1)*e1[ii] - a2v*fabs(l2)*e2v[ii] - a4*fabs(l4)*e4[ii] - a5*fabs(l5)*e5[ii]);
	}
}

void solverFluxAUSMD(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L, 
               double U0R, double U1R, double U2R, double U3R,
	           double* f)
{

    /*
    Based on: YASUHIRO WADA † AND MENG-SING LIOU, AN ACCURATE AND ROBUST FLUX 
    SPLITTING SCHEME FOR SHOCK AND CONTACT DISCONTINUITIES, (1997)
    */
    
	double uL = U1L/U0L;
	double vL = U2L/U0L;
    double pL = (solver->gamma - 1)*(U3L - (uL*uL + vL*vL)*U0L/2);
    double HL = (U3L + pL)/U0L;

	double uR = U1R/U0R;
	double vR = U2R/U0R;
    double pR = (solver->gamma - 1)*(U3R - (uR*uR + vR*vR)*U0R/2);
    double HR = (U3R + pR)/U0R;

	double cm = fmax(sqrt(solver->gamma*pL/U0L), sqrt(solver->gamma*pR/U0R));

	double alphaL = (2.0*pL/U0L)/(pL/U0L + pR/U0R);
	double alphaR = (2.0*pR/U0R)/(pL/U0L + pR/U0R);

	double uPlus, pPlus;
	if (fabs(uL) < cm) {
		uPlus = 0.25*alphaL*(uL + cm)*(uL + cm)/cm + 0.5*(1.0 - alphaL)*(uL + fabs(uL));
		pPlus = 0.25*pL*(uL + cm)*(uL + cm)*(2.0 - uL/cm)/(cm*cm);
	} else {
		uPlus = 0.5*(uL + fabs(uL));
		pPlus = 0.5*pL*(uL + fabs(uL))/uL;
	}

	double uMinus, pMinus;
	if (fabs(uR) < cm) {
		uMinus = - 0.25*alphaR*(uR - cm)*(uR - cm)/cm + 0.5*(1.0 - alphaR)*(uR - fabs(uR));
		pMinus = 0.25*pR*(uR - cm)*(uR - cm)*(2.0 + uR / cm)/(cm * cm);
	} else {
		uMinus = 0.5*(uR - fabs(uR));
		pMinus = 0.5*pR*(uR - fabs(uR))/uR;
	}


	double rU = uPlus*U0L + uMinus*U0R;
	f[0] = rU;
	f[1] = 0.5*(rU * (uR + uL) - fabs(rU) * (uR - uL)) + (pPlus + pMinus);
	f[2] = 0.5*(rU * (vR + vL) - fabs(rU) * (vR - vL));
	f[3] = 0.5*(rU * (HR + HL) - fabs(rU) * (HR - HL));
}


void solverFluxFree(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L,
	           double* f)
{

	double uL = U1L/U0L;
	double vL = U2L/U0L;
    double pL = (solver->gamma - 1)*(U3L - (uL*uL + vL*vL)*U0L/2);
    double HL = (U3L + pL)/U0L;

    // Fluxes
	f[0] = U1L;
	f[1] = U1L*uL + pL;
	f[2] = U1L*vL;
	f[3] = uL*(U3L + pL);

}

void solverFlux(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L, 
               double U0R, double U1R, double U2R, double U3R,
	           double* f)
{	           
	if(solver->flux == 0)
	{
        solverFluxRoe(solver, U0L, U1L, U2L, U3L, U0R, U1R, U2R, U3R, f);
    }
    else if(solver->flux == 1)     
    {
        solverFluxAUSMD(solver, U0L, U1L, U2L, U3L, U0R, U1R, U2R, U3R, f);
    }
}

void rotation(double* U, double dSx, double dSy, double dS)
{

    double aux = (U[1]*dSx + U[2]*dSy)/dS;
    U[2] = (-U[1]*dSy + U[2]*dSx)/dS;
    U[1] = aux;
	
}

double interpMUSCL_ii(double **U, int ii, int jj, double e)
{

    /*
    Based on: Blazek J., Computacional Fluid Dynamics, Principles and Applications (2001)
    */

    double a = U[ii+1][jj] - U[ii][jj];
    double b = U[ii][jj] - U[ii-1][jj];

    double delta = (a*(b*b+e) + b*(a*a+e))/(a*a+b*b+2*e);
    //double delta = ((2*a*a+e)*b + (b*b+2*e)*a)/(2*a*a+2*b*b-a*b+3*e);
    
    return delta;

}

double interpMUSCL_jj(double **U, int ii, int jj, double e)
{

    /*
    Based on: Blazek J., Computacional Fluid Dynamics, Principles and Applications (2001)
    */

    double a = U[ii][jj+1] - U[ii][jj];
    double b = U[ii][jj] - U[ii][jj-1];

    double delta = (a*(b*b+e) + b*(a*a+e))/(a*a+b*b+2*e);
    //double delta = ((2*a*a+e)*b + (b*b+2*e)*a)/(2*a*a+2*b*b-a*b+3*e);
    
    return delta;

}

void solverGradient(SOLVER* solver, double** U, int ii, int jj, double* dUx, double* dUy)
{

    double dSx, dSy;
    double aux;

    *dUx = 0.0;
    *dUy = 0.0;

    meshCalcDSright(solver->mesh, ii, jj, &dSx, &dSy);
    if(ii<solver->Nrow-1)
    {
        aux = (U[ii][jj] + U[ii+1][jj])*0.5;
    }
    else
    {
        aux = (3*U[ii][jj] - U[ii-1][jj])*0.5;
    }
    *dUx += aux*dSx;
    *dUy += aux*dSy;

    meshCalcDSup(solver->mesh, ii, jj, &dSx, &dSy);
    if(jj<solver->Ncol-1)
    {
        aux = (U[ii][jj] + U[ii][jj+1])*0.5;
    }
    else
    {
        aux = (3*U[ii][jj] - U[ii][jj-1])*0.5;
    }
    *dUx += aux*dSx;
    *dUy += aux*dSy;

    meshCalcDSleft(solver->mesh, ii, jj, &dSx, &dSy);
    if(ii>0)
    {
        aux = (U[ii][jj] + U[ii-1][jj])*0.5;
    }
    else
    {
        aux = (3*U[ii][jj] - U[ii+1][jj])*0.5;
    }
    *dUx += aux*dSx;
    *dUy += aux*dSy;

    meshCalcDSdown(solver->mesh, ii, jj, &dSx, &dSy);
    if(jj>0)
    {
        aux = (U[ii][jj] + U[ii][jj-1])*0.5;
    }
    else
    {
        aux = (3*U[ii][jj] - U[ii][jj+1])*0.5;
    }
    *dUx += aux*dSx;
    *dUy += aux*dSy;

    aux = meshCalcOmega(solver->mesh, ii, jj);

    *dUx /= aux;
    *dUy /= aux;

}

void inter(SOLVER* solver, double ***U)
{  
    // Vertical Element boundary
    # pragma omp parallel for
    for(int jj=0; jj<solver->Ncol; jj++)
    {
    
		int kk;
        double dSx, dSy, dS;
        double aux;
        double UL[4];
		double UR[4];
        double f[4];
        double dUx, dUy, xc, yc, x, y, delta;

        for(int ii=0; ii<solver->Nrow-1; ii++)
        {
     
            meshCalcDSright(solver->mesh, ii, jj, &dSx, &dSy);
            dS = sqrt(dSx*dSx + dSy*dSy);
            
			if(ii>0 & ii<solver->Nrow-2 & solver->MUSCL)
			{
				for(kk=0; kk<4; kk++)
				{
					solverGradient(solver, U[kk], ii, jj, &dUx, &dUy);				
                    meshCalcC(solver->mesh, ii, jj, &xc, &yc);
                    meshCalcCright(solver->mesh, ii, jj, &x, &y);
                    UL[kk] = U[kk][ii][jj] + (dUx*(x-xc) + dUy*(y-yc));

                    solverGradient(solver, U[kk], ii+1, jj, &dUx, &dUy);				
                    meshCalcC(solver->mesh, ii+1, jj, &xc, &yc);
                    meshCalcCleft(solver->mesh, ii+1, jj, &x, &y);
                    UR[kk] = U[kk][ii+1][jj] + (dUx*(x-xc) + dUy*(y-yc));
                }
			}
			else
			{
				for(kk=0; kk<4; kk++)
				{
                    UL[kk] = U[kk][ii][jj];
					UR[kk] = U[kk][ii+1][jj];
				}
			}

            // Rotation of the velocity vectors
			rotation(UL, dSx, dSy, dS);
                        
            // Rotation of the velocity vectors
			rotation(UR, dSx, dSy, dS);
            
            // Flux calculation
            solverFlux(solver, UL[0], UL[1], UL[2], UL[3], UR[0], UR[1], UR[2], UR[3], f);

            // Rotation of the flux
			rotation(f, dSx, -dSy, dS);
                             
            for(kk=0; kk<4; kk++)
            {
                aux = f[kk]*dS;
                solver->R[kk][ii][jj] += aux;
                solver->R[kk][ii+1][jj] -= aux;
            } 
        }        
    }

    // Horizontal Element boundary
    # pragma omp parallel for
    for(int ii=0; ii<solver->Nrow; ii++)
    {
    
		int kk;
        double dSx, dSy, dS;
        double aux;
        double UL[4];
        double UR[4];
        double f[4];
        double dUx, dUy, xc, yc, x, y, delta;
        
        for(int jj=0; jj<solver->Ncol-1; jj++)
        {
     
            meshCalcDSup(solver->mesh, ii, jj, &dSx, &dSy);
            dS = sqrt(dSx*dSx + dSy*dSy);

			if(jj>0 & jj<solver->Ncol-2 & solver->MUSCL)
			{
				for(kk=0; kk<4; kk++)
				{
 					
					solverGradient(solver, U[kk], ii, jj, &dUx, &dUy);				
                    meshCalcC(solver->mesh, ii, jj, &xc, &yc);
                    meshCalcCup(solver->mesh, ii, jj, &x, &y);
                    UL[kk] = U[kk][ii][jj] + (dUx*(x-xc) + dUy*(y-yc));


                    solverGradient(solver, U[kk], ii, jj+1, &dUx, &dUy);				
                    meshCalcC(solver->mesh, ii, jj+1, &xc, &yc);
                    meshCalcCdown(solver->mesh, ii, jj+1, &x, &y);
                    UR[kk] = U[kk][ii][jj+1] + (dUx*(x-xc) + dUy*(y-yc));
				}
			}
			else
			{
				for(kk=0; kk<4; kk++)
				{
					UL[kk] = U[kk][ii][jj];
					UR[kk] = U[kk][ii][jj+1];
				}
			}

            // Rotation of the velocity vectors
			rotation(UL, dSx, dSy, dS);
                        
            // Rotation of the velocity vectors
			rotation(UR, dSx, dSy, dS);
            
            // Flux calculation
            solverFlux(solver, UL[0], UL[1], UL[2], UL[3], UR[0], UR[1], UR[2], UR[3], f);

            // Rotation of the flux
			rotation(f, dSx, -dSy, dS);
                             
            for(kk=0; kk<4; kk++)
            {
                aux = f[kk]*dS;
                solver->R[kk][ii][jj] += aux;
                solver->R[kk][ii][jj+1] -= aux;   
            }             
        } 
    }
}

void boundaryCalc(SOLVER* solver, double*** U, int ii, int jj, double dSx, double dSy, int flagBC, int flagWall)
{

	int kk;
    double dS;
    double UL[4];
    double UR[4];
    double f[4];

    dS = sqrt(dSx*dSx + dSy*dSy);

	for(kk=0; kk<4; kk++)
	{
		UL[kk] = U[kk][ii][jj];
	}
	
    // Rotation of the velocity vectors
    rotation(UL, dSx, dSy, dS);

    if(flagBC == 0)
    {
    
        // Reflexive
        solverFlux(solver, UL[0], UL[1], UL[2], UL[3], UL[0], -UL[1], UL[2], UL[3], f);

    }
    else if(flagBC == 1)
    {

        // Inlet
		for(kk=0; kk<4; kk++)
		{
			UR[kk] = solver->inlet->Uin[kk];
		}

        // Rotation of the velocity vectors
	    rotation(UR, dSx, dSy, dS);

		solverFlux(solver, UL[0], UL[1], UL[2], UL[3], UR[0], UR[1], UR[2], UR[3], f);
    
    }
    else if(flagBC == 2)
    {
    
        // Outlet
        solverFluxFree(solver, UL[0], UL[1], UL[2], UL[3], f);
    
    }
    else if(flagBC == 3)
    {
    
        /*
        Based on: Blazek J., Computacional Fluid Dynamics, Principles and Applications (2001)
        */        

		double p2 = solverCalcP(solver, U, ii, jj);
		double p3, p4;
		
        if(flagWall == 0)
        {
	        p3 = solverCalcP(solver, U, ii, jj+1);
	        p4 = solverCalcP(solver, U, ii, jj+2);
        }
        else if(flagWall == 1)
        {
	        p3 = solverCalcP(solver, U, ii, jj-1);
	        p4 = solverCalcP(solver, U, ii, jj-2);
        }
        else if(flagWall == 2)
        {
	        p3 = solverCalcP(solver, U, ii+1, jj);
	        p4 = solverCalcP(solver, U, ii+2, jj);
        }
        else if(flagWall == 3)
        {
	        p3 = solverCalcP(solver, U, ii-1, jj);
	        p4 = solverCalcP(solver, U, ii-2, jj);
        }
	
        // Outlet
        f[0] = .0;
        f[2] = .0;
        f[3] = .0;

        //f[1] =  p2;
        //f[1] =  0.5*(3*p2 - p3);
        f[1] =  0.125*(15*p2 - 10*p3 + 3*p4);

    }

    // Rotation reverse of the flux
    rotation(f, dSx, -dSy, dS); 
    
    if(dS>0.0)
    {         
        for(int kk=0; kk<4; kk++)
        {
            solver->R[kk][ii][jj] += f[kk]*dS;
        }
    }
}

void boundary(SOLVER* solver, double*** U)
{

    int ii;
    int jj;    
    double dSx, dSy;
    

    for(ii=0; ii<solver->Nrow; ii++)
    {

        // Boundary down
        jj = 0; 
        meshCalcDSdown(solver->mesh, ii, jj, &dSx, &dSy);
        boundaryCalc(solver, U, ii, jj, dSx, dSy, solver->bc->Ndown, 0);
        
        // Boundary up
        jj = solver->Ncol-1;
        meshCalcDSup(solver->mesh, ii, jj, &dSx, &dSy);
        boundaryCalc(solver, U, ii, jj, dSx, dSy, solver->bc->Nup, 1);

    }

    for(jj=0; jj<solver->Ncol; jj++)
    {

        // Boundary left
        ii = 0;
        meshCalcDSleft(solver->mesh, ii, jj, &dSx, &dSy);
        boundaryCalc(solver, U, ii, jj, dSx, dSy, solver->bc->Nleft, 2);

        // Boundary right
        ii = solver->Nrow-1;
        meshCalcDSright(solver->mesh, ii, jj, &dSx, &dSy);
        boundaryCalc(solver, U, ii, jj, dSx, dSy, solver->bc->Nright, 3);

    }

}

void interAxisPressure(SOLVER* solver, double ***U)
{

    double dS;

    for(int ii=0; ii<solver->Nrow; ii++)
    {        
        for(int jj=0; jj<solver->Ncol; jj++)
        {

            dS = meshCalcDSlateral(solver->mesh, ii, jj);
            solver->R[2][ii][jj] -= solverCalcP(solver, U, ii, jj)*dS;

        } 
    }
}

void calcSpectralRad(SOLVER* solver, double*** U, int ii, int jj, double* LcI, double* LcJ)
{

    double u, v, c, dSx, dSy, dS;

    solverCalcVel(solver, U, ii, jj, &u, &v, &c);

    meshCalcDSI(solver->mesh, ii, jj, &dSx, &dSy);
    dS = sqrt(dSx*dSx + dSy*dSy);
    *LcI = (fabs(u*dSx + v*dSy) + c*dS);

    meshCalcDSJ(solver->mesh, ii, jj, &dSx, &dSy);
    dS = sqrt(dSx*dSx + dSy*dSy);
    *LcJ = (fabs(u*dSx + v*dSy) + c*dS);

}

double solverCalcDt(SOLVER* solver)
{
    double dt, dtMin;
    double LcI, LcJ;
    double omega;

    for(int ii=0; ii<solver->Nrow; ii++)
    {
        for(int jj=0; jj<solver->Ncol; jj++)
        {

            omega = meshCalcOmega(solver->mesh, ii, jj);

            calcSpectralRad(solver, solver->U, ii, jj, &LcI, &LcJ);

			dt = 0.5*solver->stages*omega/(LcI + LcJ);
            
            if(ii==0 & jj==0)
            {
                dtMin = dt;
            }
            else
            {
                if(dt<dtMin)
                {
                    dtMin = dt;
                }
            }
        }
    }

    return dtMin;

}

void solverCalcR(SOLVER* solver, double*** U)
{

    solverResetR(solver);

    inter(solver, U);
    
    boundary(solver, U);  
    
    if(solver->mesh->axi==1)
    {
        interAxisPressure(solver, U);
    }

}

void solverRK(SOLVER* solver, double a)
{

    # pragma omp parallel for
    for(int ii=0; ii<solver->mesh->Nrow-1; ii++)
    {
        for(int jj=0; jj<solver->mesh->Ncol-1; jj++)
        {
            double omega = meshCalcOmega(solver->mesh, ii, jj);
            for(int kk=0; kk<4; kk++)
            {
                
                solver->Uaux[kk][ii][jj] = solver->U[kk][ii][jj] - solver->dt*a*solver->R[kk][ii][jj]/omega;
                
            }
        }
    }

}

void solverUpdateU(SOLVER* solver)
{

    # pragma omp parallel for
    for(int ii=0; ii<solver->Nrow; ii++)
    {
        for(int jj=0; jj<solver->Ncol; jj++)
        {
            for(int kk=0; kk<4; kk++)
            {
                solver->U[kk][ii][jj] = solver->Uaux[kk][ii][jj];
            }
        }
    }
}

void solverStepRK(SOLVER* solver)
{   

    if(solver->stages==3)
    {
        solverCalcR(solver, solver->U);
        solverRK(solver, 0.1481);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.4);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 1.0);
    }
    else if(solver->stages==4)
    {
        solverCalcR(solver, solver->U);
        solverRK(solver, 0.0833);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.2069);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.4265);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 1.0);
    }
    else if(solver->stages==5)
    {
        solverCalcR(solver, solver->U);
        solverRK(solver, 0.0533);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.1263);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.2375);

        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 0.4414);
        
        solverCalcR(solver, solver->Uaux);
        solverRK(solver, 1.0);
    }
    
    solverUpdateU(solver);

}

void solverWriteU(SOLVER* solver, char* fileName)
{

    FILE* ff = fopen(fileName, "w");
    fprintf(ff, "4, \n");
    for(int kk=0; kk<4; kk++)
    {
        fprintf(ff, "%i, %i, %i, \n", kk, solver->Nrow, solver->Ncol);
        for(int ii=0; ii<solver->Nrow; ii++)
        {
            for(int jj=0; jj<solver->Ncol; jj++)
            {
         
                fprintf(ff, " %.4f,", solver->U[kk][ii][jj]);   
                        
            }
            fprintf(ff, " \n");
        }
    }

}

void solverCalcRes(SOLVER* solver)
{
    
    double aux;
    
    for(int kk=0; kk<4; kk++)
    {
        solver->res[kk] = fabs(solver->R[kk][0][0]);
    }
    
    for(int ii=0; ii<solver->mesh->Nrow-1; ii++)
    {
        for(int jj=0; jj<solver->mesh->Ncol-1; jj++)
        {
            for(int kk=0; kk<4; kk++)
            {
                if(solver->res[kk] < fabs(solver->R[kk][ii][jj]))
                {
                    solver->res[kk] = fabs(solver->R[kk][ii][jj]);
                }
            }
        }
    }
    
    for(int kk=0; kk<4; kk++)
    {
        printf(" %+.4e,", solver->res[kk]);        
    }
    printf("\n");    
   
}

double duration(struct timeval start, struct timeval stop){

    double tstart, tstop;
    tstart = (double)start.tv_sec + start.tv_usec/1000000.;
    tstop = (double)stop.tv_sec + stop.tv_usec/1000000.;

    return tstop - tstart;

}

int main(int argc, char **argv)
{

    SOLVER* solver = malloc(sizeof(SOLVER));
    
    struct timeval start;
    struct timeval stop;
    
    char s[50];
       
    // Load input   
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "input.dat");
    INPUT* input = inputInit(s, 50);
    printf("Input data:\n");
    inputPrint(input);
   
    // Set number of threads
    omp_set_num_threads(atoi(inputGetValue(input, "threads")));   
        
    // Load mesh    
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "mesh.csv");
    solver->mesh = meshInit(s);
    solver->Nrow = solver->mesh->Nrow-1;
    solver->Ncol = solver->mesh->Ncol-1;        

    // axisymmetric
    solver->mesh->axi = atoi(inputGetValue(input, "axisymmetric"));

	//meshPrintOmega(solver->mesh);
    //meshPrintDStotal(solver->mesh);

    // Memory allocation
    solverAllocate(solver);
 
    // Boundary conditions
    solver->bc = malloc(sizeof(BOUNDARY));
    solver->bc->down = inputGetValue(input, "BCdown");
    solver->bc->up = inputGetValue(input, "BCup");
    solver->bc->left = inputGetValue(input, "BCleft");
    solver->bc->right = inputGetValue(input, "BCright");

    boundarySet(solver->bc);
  
    // Constants
    solver->Rgas = 287.5;
    solver->gamma = 1.4;  
    solver->eFix = 0.1;
    solver->e = strtod(inputGetValue(input, "interpE"), NULL);
    
    // Seletion of MUSCL and flux
    solver->MUSCL = atoi(inputGetValue(input, "MUSCL"));
    solver->flux = fluxChoice(inputGetValue(input, "flux"));
    solver->stages = atoi(inputGetValue(input, "stages"));

    if(atoi(inputGetValue(input, "tube")) == 0)
    {

        //stead state
        solver->inlet = conditionInit(strtod(inputGetValue(input, "pressure"), NULL), 
                                      strtod(inputGetValue(input, "temperature"), NULL), 
                                      strtod(inputGetValue(input, "mach"), NULL), 
                                      strtod(inputGetValue(input, "nx"), NULL),
                                      strtod(inputGetValue(input, "ny"), NULL));
        
        // Initialization of U
        solverInitU(solver, solver->inlet);
        
        // Calculate time step        
        int Nmax = atoi(inputGetValue(input, "Nmax"));

        // Run the solver
        printf("\nRunning solution:\n");
        gettimeofday(&start, NULL);
        for(int ii=0; ii<Nmax; ii++)
        {
            solver->dt = solverCalcDt(solver);
            solverStepRK(solver);
            
            if(ii%100 == 0)
            {
                printf("%i, ", ii);
                solverCalcRes(solver);
            }
            
        }
        gettimeofday(&stop, NULL);  
        printf("Duration %f s\n", duration(start, stop));
 
    }
    else
    {
    
        CONDITION* inside1 = conditionInit(strtod(inputGetValue(input, "pressure1"), NULL), 
                                           strtod(inputGetValue(input, "temperature1"), NULL), 
                                           strtod(inputGetValue(input, "mach1"), NULL), 
                                           strtod(inputGetValue(input, "nx1"), NULL),
                                           strtod(inputGetValue(input, "ny1"), NULL));

        CONDITION* inside2 = conditionInit(strtod(inputGetValue(input, "pressure2"), NULL), 
                                           strtod(inputGetValue(input, "temperature2"), NULL), 
                                           strtod(inputGetValue(input, "mach2"), NULL), 
                                           strtod(inputGetValue(input, "nx2"), NULL),
                                           strtod(inputGetValue(input, "ny2"), NULL));      
    
        // Initialization of U
        solverInitUTube(solver, inside1, inside2, strtod(inputGetValue(input, "xm"), NULL));
        free(inside1);
        free(inside2);
        
        double tmax = strtod(inputGetValue(input, "tmax"), NULL);                

        // Run the solver
        double t = 0.0;
        printf("\nRunning solution:\n");
        gettimeofday(&start, NULL);
        int stopLoop = 0;
        int ii = 0;
        while(stopLoop == 0)
        {
            
            solver->dt = solverCalcDt(solver);
            
            if(t + solver->dt>tmax)
            {
                solver->dt = (tmax-t);
                stopLoop = 1;
            }

            solverStepRK(solver);
            t += solver->dt;
            ii++;

            if(ii%100 == 0)
            {
                printf("%i, ", ii);
                solverCalcRes(solver);
            }

        }
        gettimeofday(&stop, NULL);  
        printf("Duration %f s\n", duration(start, stop));
        printf("time %f s\n", t);

    }

    // Save solution
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "solution.csv");
    solverWriteU(solver, s);

    // Free memory 
    solverFree(solver);
    inputFree(input);

    return 0;


}
