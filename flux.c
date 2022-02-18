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

void entropyFix(SOLVER* solver, double *l)
{

    // Harten Hyman entropy fix
    if(*l < solver->eFix & *l > -solver->eFix )
    {
        *l = 0.5*(*l * *l/solver->eFix + solver->eFix);
    }

}

void fluxRoe(SOLVER* solver, 
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

void fluxAUSMD(SOLVER* solver, 
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


void fluxFree(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L,
	           double* f)
{

	double uL = U1L/U0L;
	double vL = U2L/U0L;
    double pL = (solver->gamma - 1)*(U3L - (uL*uL + vL*vL)*U0L/2);

    // Fluxes
	f[0] = U1L;
	f[1] = U1L*uL + pL;
	f[2] = U1L*vL;
	f[3] = uL*(U3L + pL);

}

void flux(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L, 
               double U0R, double U1R, double U2R, double U3R,
	           double* f)
{	           
	if(solver->flux == 0)
	{
        fluxRoe(solver, U0L, U1L, U2L, U3L, U0R, U1R, U2R, U3R, f);
    }
    else if(solver->flux == 1)     
    {
        fluxAUSMD(solver, U0L, U1L, U2L, U3L, U0R, U1R, U2R, U3R, f);
    }
}
