
typedef struct {

    int Nrow;
    int Ncol;
    int axi;
    
    double **x;
    double **y;

} MESH;

MESH* meshInit(char* fileName);

double meshDeltaMin(MESH* mesh);

void meshCalcDSup(MESH* mesh, int ii, int jj, double* dSx, double* dSy);

void meshCalcDSright(MESH* mesh, int ii, int jj, double* dSx, double* dSy);

void meshCalcDSdown(MESH* mesh, int ii, int jj, double* dSx, double* dSy);

void meshCalcDSleft(MESH* mesh, int ii, int jj, double* dSx, double* dSy);

double meshCalcDSlateral(MESH* mesh, int ii, int jj);

double meshCalcOmega(MESH* mesh, int ii, int jj);

void meshPrintOmega(MESH* mesh);

void meshCalcDSI(MESH* mesh, int ii, int jj, double* dSx, double* dSy);

void meshCalcDSJ(MESH* mesh, int ii, int jj, double* dSx, double* dSy);

