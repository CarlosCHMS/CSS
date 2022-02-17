

int boundaryChoice(char* s);

void boundarySet(BOUNDARY* bc);

void boundaryCalc(SOLVER* solver, double*** U, int ii, int jj, double dSx, double dSy, int flagBC, int flagWall);

void boundary(SOLVER* solver, double*** U);
