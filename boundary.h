

int boundaryChoice(char* s);

void boundarySet(BOUNDARY* bc);

void boundaryInlet(SOLVER* solver, double* Ua, double* Ud, double* Ub, double nx, double ny);

void boundaryOutlet(SOLVER* solver, double* Ud, double* Ub, double nx, double ny);

void boundaryCalc(SOLVER* solver, double*** U, int ii, int jj, double dSx, double dSy, int flagBC, int flagWall);

void boundary(SOLVER* solver, double*** U);
