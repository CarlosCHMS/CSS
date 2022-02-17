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

void boundaryInlet(SOLVER* solver, double* Ua, double* Ud, double* Ub, double nx, double ny)
{

    double rd = Ud[0];
    double ud = Ud[1]/Ud[0];
    double vd = Ud[2]/Ud[0];
    double pd = (solver->gamma - 1)*(Ud[3] - 0.5*(ud*ud + vd*vd)*rd);

    double c0 = sqrt(solver->gamma*pd/rd);    

    double m = sqrt(ud*ud + vd*vd)/c0;
    
    if(m < 1.)
    {

        double ra = Ua[0];
        double ua = Ua[1]/Ua[0];
        double va = Ua[2]/Ua[0];
        double pa = (solver->gamma - 1)*(Ua[3] - 0.5*(ua*ua + va*va)*ra);

        double pb = 0.5*(pa + pd - rd*c0*(nx*(ua-ud) + ny*(va-vd)));
        double rb = ra + (pb - pa)/(c0*c0);
        double ub = ua - nx*(pa - pb)/(rd*c0);
        double vb = va - ny*(pa - pb)/(rd*c0);

        Ub[0] = rb;
        Ub[1] = rb*ub;
        Ub[2] = rb*vb;
        Ub[3] = pb/(solver->gamma-1) + 0.5*(ub*ub + vb*vb)*rb;

    }
    else
    {
        for(int ii=0; ii<4; ii++)
        {
        
            Ub[ii] = Ua[ii];

        }
    }    

}

void boundaryOutlet(SOLVER* solver, double* Ud, double* Ub, double nx, double ny)
{

    double rd = Ud[0];
    double ud = Ud[1]/Ud[0];
    double vd = Ud[2]/Ud[0];
    double pd = (solver->gamma - 1)*(Ud[3] - 0.5*(ud*ud + vd*vd)*rd);

    double c0 = sqrt(solver->gamma*pd/rd);    

    double m = sqrt(ud*ud + vd*vd)/c0;
    
    if(m < 1.)
    {
        double pb = solver->pout;
        double rb = rd + (pb - pd)/(c0*c0);
        double ub = ud + nx*(pd - pb)/(rd*c0);
        double vb = vd + ny*(pd - pb)/(rd*c0);

        Ub[0] = rb;
        Ub[1] = rb*ub;
        Ub[2] = rb*vb;
        Ub[3] = pb/(solver->gamma-1) + 0.5*(ub*ub + vb*vb)*rb;
    }
    else
    {
        for(int ii=0; ii<4; ii++)
        {
        
            Ub[ii] = Ud[ii];

        }
    }    

}


void boundaryCalc(SOLVER* solver, double*** U, int ii, int jj, double dSx, double dSy, int flagBC, int flagWall)
{

	int kk;
    double dS;
    double UL[4];
    double Ub[4];
    double f[4];

    dS = sqrt(dSx*dSx + dSy*dSy);

	for(kk=0; kk<4; kk++)
	{
		UL[kk] = U[kk][ii][jj];
	}

    if(flagBC == 0)
    {

        // Rotation of the velocity vectors
        rotation(UL, dSx, dSy, dS);
    
        // Reflexive
        solverFlux(solver, UL[0], UL[1], UL[2], UL[3], UL[0], -UL[1], UL[2], UL[3], f);

    }
    else if(flagBC == 1)
    {

        boundaryInlet(solver, solver->inlet->Uin, UL, Ub, dSx/dS, dSy/dS);

        // Rotation of the velocity vectors
	    rotation(Ub, dSx, dSy, dS);

		solverFluxFree(solver, Ub[0], Ub[1], Ub[2], Ub[3], f);
    
    }
    else if(flagBC == 2)
    {

        boundaryOutlet(solver, UL, Ub, dSx/dS, dSy/dS);

        // Rotation of the velocity vectors
	    rotation(Ub, dSx, dSy, dS);

		solverFluxFree(solver, Ub[0], Ub[1], Ub[2], Ub[3], f);
    
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

