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

    double p;
    double u = U[1][ii][jj]/U[0][ii][jj];
    double v = U[2][ii][jj]/U[0][ii][jj];
    double E = U[3][ii][jj]/U[0][ii][jj];
    double aux;
    
    aux = E - (u*u + v*v)/2;
    aux *= solver->gamma - 1;
    p = aux*solver->U[0][ii][jj];  
    
    return p;  
    
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


void interpMUSCL_iiL1(double **U, int ii, int jj, double *UL, double e)
{

    /*
    Based on: Blazek J., Computacional Fluid Dynamics, Principles and Applications (2001)
    */

    double a = U[ii+1][jj] - U[ii][jj];
    double b = U[ii][jj] - U[ii-1][jj];

    //double delta = (a*(b*b+e) + b*(a*a+e))/(a*a+b*b+2*e);
    double delta = ((2*a*a+e)*b + (b*b+2*e)*a)/(2*a*a+2*b*b-a*b+3*e);
    
    *UL = U[ii][jj] + 0.5*delta;

}

void interpMUSCL_iiR1(double **U, int ii1, int jj, double *UR, double e)
{

    double a = U[ii1+1][jj] - U[ii1][jj];
    double b = U[ii1][jj] - U[ii1-1][jj];

    //double delta = (a*(b*b+e) + b*(a*a+e))/(a*a+b*b+2*e);
    double delta = ((2*a*a+e)*b + (b*b+2*e)*a)/(2*a*a+2*b*b-a*b+3*e);
    
    *UR = U[ii1][jj] - 0.5*delta;

}

void interpMUSCL_jjL1(double **U, int ii, int jj, double *UL, double e)
{

    double a = U[ii][jj+1] - U[ii][jj];
    double b = U[ii][jj] - U[ii][jj-1];

    //double delta = (a*(b*b+e) + b*(a*a+e))/(a*a+b*b+2*e);
    double delta = ((2*a*a+e)*b + (b*b+2*e)*a)/(2*a*a+2*b*b-a*b+3*e);
    
    *UL = U[ii][jj] + 0.5*delta;

}

void interpMUSCL_jjR1(double **U, int ii, int jj1, double *UR, double e)
{

    double a = U[ii][jj1+1] - U[ii][jj1];
    double b = U[ii][jj1] - U[ii][jj1-1];

    //double delta = (a*(b*b+e) + b*(a*a+e))/(a*a+b*b+2*e);
    double delta = ((2*a*a+e)*b + (b*b+2*e)*a)/(2*a*a+2*b*b-a*b+3*e);
    
    *UR = U[ii][jj1] - 0.5*delta;

}


void inter(SOLVER* solver, double ***U)
{

    
    // Vertical Element boundary
    # pragma omp parallel for
    for(int jj=0; jj<solver->Ncol; jj++)
    {
    
        double dSx, dSy;
        double dS;
        double aux;
        double U0L, U1L, U2L, U3L, U0R, U1R, U2R, U3R;
        double f[4];
        double e;

        for(int ii=0; ii<solver->Nrow-1; ii++)
        {
     
            e = sqrt(meshCalcOmega(solver->mesh, ii, jj))*solver->e;
     
            dSx = solver->mesh->y[ii+1][jj+1] - solver->mesh->y[ii+1][jj];
            dSy = -(solver->mesh->x[ii+1][jj+1] - solver->mesh->x[ii+1][jj]);
            dS = sqrt(dSx*dSx+ dSy*dSy);
            
            if(ii>0 & solver->MUSCL)
            {
                                              
                interpMUSCL_iiL1(U[0], ii, jj, &U0L, e);
                interpMUSCL_iiL1(U[1], ii, jj, &U1L, e);
                interpMUSCL_iiL1(U[2], ii, jj, &U2L, e);
                interpMUSCL_iiL1(U[3], ii, jj, &U3L, e);

            }
            else
            {

                U0L = U[0][ii][jj];
                U1L = U[1][ii][jj];
                U2L = U[2][ii][jj];
                U3L = U[3][ii][jj];                                    


            }


            // boudary perpendicular to ii
            if(ii<solver->Nrow-2 & solver->MUSCL)
            {
                                              
                interpMUSCL_iiR1(U[0], ii+1, jj, &U0R, e);
                interpMUSCL_iiR1(U[1], ii+1, jj, &U1R, e);
                interpMUSCL_iiR1(U[2], ii+1, jj, &U2R, e);
                interpMUSCL_iiR1(U[3], ii+1, jj, &U3R, e);

            }
            else
            {

                U0R = U[0][ii+1][jj];
                U1R = U[1][ii+1][jj];
                U2R = U[2][ii+1][jj];
                U3R = U[3][ii+1][jj];                                    

            }

            // Rotation of the velocity vectors
            aux = (U1L*dSx + U2L*dSy)/dS;
            U2L = (-U1L*dSy + U2L*dSx)/dS;
            U1L = aux;
                        
            // Rotation of the velocity vectors
            aux = (U1R*dSx + U2R*dSy)/dS;
            U2R = (-U1R*dSy + U2R*dSx)/dS;
            U1R = aux;
            
            // Flux calculation
            solverFlux(solver, U0L, U1L, U2L, U3L, U0R, U1R, U2R, U3R, f);

            // Rotation of the flux
            aux = (f[1]*dSx - f[2]*dSy)/dS;
            f[2] =  (f[1]*dSy + f[2]*dSx)/dS;
            f[1] = aux;  
                             
            for(int kk=0; kk<4; kk++)
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
    
        double dSx, dSy;
        double dS;
        double aux;
        double U0L, U1L, U2L, U3L, U0R, U1R, U2R, U3R;
        double f[4];
        double e;
        
        for(int jj=0; jj<solver->Ncol-1; jj++)
        {

            e = sqrt(meshCalcOmega(solver->mesh, ii, jj))*solver->e;
     
            dSx = -(solver->mesh->y[ii+1][jj+1] - solver->mesh->y[ii][jj+1]);
            dSy = solver->mesh->x[ii+1][jj+1] - solver->mesh->x[ii][jj+1];
            dS = sqrt(dSx*dSx+ dSy*dSy);
                                      
            if(jj>0 & solver->MUSCL)
            {
                                              
                interpMUSCL_jjL1(U[0], ii, jj, &U0L, e);
                interpMUSCL_jjL1(U[1], ii, jj, &U1L, e);
                interpMUSCL_jjL1(U[2], ii, jj, &U2L, e);
                interpMUSCL_jjL1(U[3], ii, jj, &U3L, e);
                
            }
            else
            {

                U0L = U[0][ii][jj];
                U1L = U[1][ii][jj];
                U2L = U[2][ii][jj];
                U3L = U[3][ii][jj];

            }

            // boudary perpendicular to jj
            if(jj<solver->Ncol-2 & solver->MUSCL)
            {
                                              
                interpMUSCL_jjR1(U[0], ii, jj+1, &U0R, e);
                interpMUSCL_jjR1(U[1], ii, jj+1, &U1R, e);
                interpMUSCL_jjR1(U[2], ii, jj+1, &U2R, e);
                interpMUSCL_jjR1(U[3], ii, jj+1, &U3R, e);

            }
            else
            {

                U0R = U[0][ii][jj+1];
                U1R = U[1][ii][jj+1];
                U2R = U[2][ii][jj+1];
                U3R = U[3][ii][jj+1];                                    

            }

            // Rotation of the velocity vectors
            aux = (U1L*dSx + U2L*dSy)/dS;
            U2L = (-U1L*dSy + U2L*dSx)/dS;
            U1L = aux;
                
            // Rotation of the velocity vectors
            aux = (U1R*dSx + U2R*dSy)/dS;
            U2R = (-U1R*dSy + U2R*dSx)/dS;
            U1R = aux;
            
            // Flux calculation
            solverFlux(solver, U0L, U1L, U2L, U3L, U0R, U1R, U2R, U3R, f);

            // Rotation of the flux
            aux = (f[1]*dSx - f[2]*dSy)/dS;
            f[2] =  (f[1]*dSy + f[2]*dSx)/dS;
            f[1] = aux;  
                             
            for(int kk=0; kk<4; kk++)
            {

                aux = f[kk]*dS;
                solver->R[kk][ii][jj] += aux;
                solver->R[kk][ii][jj+1] -= aux;
                
            }             

        }            
    }
            
}

void boundaryCalc(SOLVER* solver, double*** U, int ii, int jj, double dSx, double dSy, int flag)
{

    double dS;
    double U0L, U1L, U2L, U3L, U0R, U1R, U2R, U3R;
    double aux;
    double f[4];


    dS = sqrt(dSx*dSx + dSy*dSy);

    U0L = U[0][ii][jj];
    U1L = U[1][ii][jj];
    U2L = U[2][ii][jj];
    U3L = U[3][ii][jj];                                    

    // Rotation of the velocity vectors
    aux = (U1L*dSx + U2L*dSy)/dS;
    U2L = (-U1L*dSy + U2L*dSx)/dS;
    U1L = aux;

    if(flag == 0)
    {
    
        // Reflexive
        solverFlux(solver, U0L, U1L, U2L, U3L, U0L, -U1L, U2L, U3L, f);

    }
    else if(flag == 1)
    {

        // Inlet
        U0R = solver->inlet->Uin[0];
        U1R = solver->inlet->Uin[1];
        U2R = solver->inlet->Uin[2];
        U3R = solver->inlet->Uin[3];

        // Rotation of the velocity vectors
        aux = (U1R*dSx + U2R*dSy)/dS;
        U2R = (-U1R*dSy + U2R*dSx)/dS;
        U1R = aux;

        solverFlux(solver, U0L, U1L, U2L, U3L, U0R, U1R, U2R, U3R, f);
    
    }
    else if(flag == 2)
    {
    
        // Outlet
        solverFluxFree(solver, U0L, U1L, U2L, U3L, f);
    
    }
    else if(flag == 3)
    {
    
        // Outlet
        double p = solverCalcP(solver, U, ii, jj);
        f[0] = .0;
        f[1] = p;
        f[2] = .0;
        f[3] = .0;                        
    
    }


    // Rotation reverse of the flux
    aux = (f[1]*dSx - f[2]*dSy)/dS;
    f[2] =  (f[1]*dSy + f[2]*dSx)/dS;
    f[1] = aux;  
                     
    for(int kk=0; kk<4; kk++)
    {

        solver->R[kk][ii][jj] += f[kk]*dS;

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
        dSx = (solver->mesh->y[ii+1][jj] - solver->mesh->y[ii][jj]);
        dSy = -(solver->mesh->x[ii+1][jj] - solver->mesh->x[ii][jj]);

        boundaryCalc(solver, U, ii, jj, dSx, dSy, solver->bc->Ndown);
        
        // Boundary up
        jj = solver->Ncol-1;
        dSx = -(solver->mesh->y[ii+1][jj+1] - solver->mesh->y[ii][jj+1]);
        dSy = solver->mesh->x[ii+1][jj+1] - solver->mesh->x[ii][jj+1];

        boundaryCalc(solver, U, ii, jj, dSx, dSy, solver->bc->Nup);

    }

    for(jj=0; jj<solver->Ncol; jj++)
    {

        // Boundary left
        ii = 0;
        dSx = -(solver->mesh->y[ii][jj+1] - solver->mesh->y[ii][jj]);
        dSy = solver->mesh->x[ii][jj+1] - solver->mesh->x[ii][jj];
        
        boundaryCalc(solver, U, ii, jj, dSx, dSy, solver->bc->Nleft);

        // Boundary right
        ii = solver->Nrow-1;
        dSx = solver->mesh->y[ii+1][jj+1] - solver->mesh->y[ii+1][jj];
        dSy = -(solver->mesh->x[ii+1][jj+1] - solver->mesh->x[ii+1][jj]);

        boundaryCalc(solver, U, ii, jj, dSx, dSy, solver->bc->Nright);

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
    solver->flux=fluxChoice(inputGetValue(input, "flux"));
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
