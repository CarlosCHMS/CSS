
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
    double CFL;
            
    double ***U;
    double ***R;
    double ***Uaux;        
    
    CONDITION* inlet;
        
    MESH* mesh;
    
    BOUNDARY* bc;

} SOLVER;


CONDITION* conditionInit(double p, double T, double mach, double nx, double ny);

void conditionState(CONDITION* cond, SOLVER* solver);

double conditionVref(CONDITION* cond, SOLVER* solver);

void solverAllocate(SOLVER* solver);

void solverInitU(SOLVER* solver, CONDITION* inside);

void solverInitUTube(SOLVER* solver, CONDITION* inside1, CONDITION* inside2, double xm);

void solverResetR(SOLVER* solver);

void solverMatFree(SOLVER* solver, double** M, int Nrow);

void solverFree(SOLVER* solver);

double solverCalcP(SOLVER* solver, double*** U, int ii, int jj);

void solverCalcVel(SOLVER* solver, double*** U, int ii, int jj, double* u, double* v, double* c);

void rotation(double* U, double dSx, double dSy, double dS);

double interpMUSCL_ii(double **U, int ii, int jj, double e);

double interpMUSCL_jj(double **U, int ii, int jj, double e);

void inter(SOLVER* solver, double ***U);

void interAxisPressure(SOLVER* solver, double ***U);

void calcSpectralRad(SOLVER* solver, double*** U, int ii, int jj, double* LcI, double* LcJ);

double solverCalcDt(SOLVER* solver);

void solverCalcR(SOLVER* solver, double*** U);

void solverRK(SOLVER* solver, double a);

void solverUpdateU(SOLVER* solver);

void solverStepRK(SOLVER* solver);

void solverWriteU(SOLVER* solver, char* fileName);

void solverCalcRes(SOLVER* solver);

double duration(struct timeval start, struct timeval stop);


