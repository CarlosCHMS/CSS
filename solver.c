#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"flux.h"
#include"boundary.h"


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
    
    return delta;

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
        double delta;

        for(int ii=0; ii<solver->Nrow-1; ii++)
        {
     
            meshCalcDSright(solver->mesh, ii, jj, &dSx, &dSy);
            dS = sqrt(dSx*dSx + dSy*dSy);
            
			if(ii>0 & solver->MUSCL)
			{
				for(kk=0; kk<4; kk++)
				{
					delta = interpMUSCL_ii(U[kk], ii, jj, solver->e);
					UL[kk] = U[kk][ii][jj] + 0.5*delta;
				}
			}
			else
			{
				for(kk=0; kk<4; kk++)
				{
					UL[kk] = U[kk][ii][jj];
				}
			}

			if(ii<solver->Nrow-2 & solver->MUSCL)
			{
				for(kk=0; kk<4; kk++)
				{
					delta = interpMUSCL_ii(U[kk], ii+1, jj, solver->e);
					UR[kk] = U[kk][ii+1][jj] - 0.5*delta;
				}
			}
			else
			{
				for(kk=0; kk<4; kk++)
				{
					UR[kk] = U[kk][ii+1][jj];
				}
			}

            // Rotation of the velocity vectors
			rotation(UL, dSx, dSy, dS);
                        
            // Rotation of the velocity vectors
			rotation(UR, dSx, dSy, dS);
            
            // Flux calculation
            flux(solver, UL[0], UL[1], UL[2], UL[3], UR[0], UR[1], UR[2], UR[3], f);

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
        double delta;
        
        for(int jj=0; jj<solver->Ncol-1; jj++)
        {
     
            meshCalcDSup(solver->mesh, ii, jj, &dSx, &dSy);
            dS = sqrt(dSx*dSx + dSy*dSy);
                                      
			if(jj>0 & solver->MUSCL)
			{
				for(kk=0; kk<4; kk++)
				{
					delta = interpMUSCL_jj(U[kk], ii, jj, solver->e);
					UL[kk] = U[kk][ii][jj] + 0.5*delta;
				}
			}
			else
			{
				for(kk=0; kk<4; kk++)
				{
					UL[kk] = U[kk][ii][jj];
				}
			}

			if(jj<solver->Ncol-2 & solver->MUSCL)
			{
				for(kk=0; kk<4; kk++)
				{
					delta = interpMUSCL_jj(U[kk], ii, jj+1, solver->e);
					UR[kk] = U[kk][ii][jj+1] - 0.5*delta;
				}
			}
			else
			{
				for(kk=0; kk<4; kk++)
				{
					UR[kk] = U[kk][ii][jj+1];
				}
			}

            // Rotation of the velocity vectors
			rotation(UL, dSx, dSy, dS);
                        
            // Rotation of the velocity vectors
			rotation(UR, dSx, dSy, dS);
            
            // Flux calculation
            flux(solver, UL[0], UL[1], UL[2], UL[3], UR[0], UR[1], UR[2], UR[3], f);

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

    dtMin *= solver->CFL;    

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


