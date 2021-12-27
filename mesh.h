
typedef struct {

    int Nrow;
    int Ncol;
    
    double **x;
    double **y;

} MESH;

MESH* meshInit(char* fileName);

double meshDeltaMin(MESH* mesh);

double meshCalcOmega(MESH* mesh, int ii, int jj);

void meshPrintOmega(MESH* mesh);

void meshCalcDSI(MESH* mesh, int ii, int jj, double* dSx, double* dSy);

void meshCalcDSJ(MESH* mesh, int ii, int jj, double* dSx, double* dSy);

