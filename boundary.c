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

int boundaryChoice(char* s)
{
    int ans;

    if(strcmp(s, "symmetry") == 0)
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
    else if(strcmp(s, "wall") == 0)
    {
        ans = 3;
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

