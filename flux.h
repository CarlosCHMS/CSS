

int fluxChoice(char* s);

void entropyFix(SOLVER* solver, double *l);

void solverFluxRoe(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L, 
               double U0R, double U1R, double U2R, double U3R,
	           double* f);

void solverFluxAUSMD(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L, 
               double U0R, double U1R, double U2R, double U3R,
	           double* f);

void solverFluxFree(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L,
	           double* f);

void solverFlux(SOLVER* solver, 
               double U0L, double U1L, double U2L, double U3L, 
               double U0R, double U1R, double U2R, double U3R,
	           double* f);
