#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"readTables.h"
#include"mesh.h"

MESH* meshInit(char* fileName)
{

    TABLELIST* tl;
    MESH* mesh = malloc(sizeof(MESH));

    tl = fReadTables(fileName);  
    
    mesh->Nrow = tl->tables[0]->Nrow;
    mesh->Ncol = tl->tables[0]->Ncol;    
    
    mesh->x = tl->tables[0]->values;
    mesh->y = tl->tables[1]->values;    

}

double meshDeltaMin(MESH* mesh)
{
    double delta, aux, dx, dy;
    
    delta = 1.;
    for(int ii=0; ii<mesh->Nrow-1; ii++)
    {
        for(int jj=0; jj<mesh->Ncol-1; jj++)
        {
            dx = mesh->x[ii+1][jj] - mesh->x[ii][jj];
            dy = mesh->y[ii+1][jj] - mesh->y[ii][jj];
            aux = dx*dx + dy*dy;
            if(aux < delta)
            {
                delta = aux;
            }

            dx = mesh->x[ii][jj+1] - mesh->x[ii][jj];
            dy = mesh->y[ii][jj+1] - mesh->y[ii][jj];
            aux = dx*dx + dy*dy;
            if(aux < delta)
            {
                delta = aux;
            }
                    
        }
    }

    delta = sqrt(delta);

    return delta;

}

double meshCalcOmega(MESH* mesh, int ii, int jj)
{

    double xac, xbd, yac, ybd;
    
    xac = mesh->x[ii][jj+1] - mesh->x[ii+1][jj];
    yac = mesh->y[ii][jj+1] - mesh->y[ii+1][jj];    

    xbd = mesh->x[ii][jj] - mesh->x[ii+1][jj+1];
    ybd = mesh->y[ii][jj] - mesh->y[ii+1][jj+1];    

    return 0.5*(xac*ybd - xbd*yac);

}

void meshPrintOmega(MESH* mesh)
{

    for(int ii=0; ii<mesh->Nrow-1; ii++)
    {
        for(int jj=0; jj<mesh->Ncol-1; jj++)
        {
            printf(" %.4f,", meshCalcOmega(mesh, ii, jj));
        }
        printf("\n");
    }

}

void meshCalcDSI(MESH* mesh, int ii, int jj, double* dSx, double* dSy)
{

    double dSx0 = (mesh->y[ii][jj+1] - mesh->y[ii][jj]);
    double dSy0 = -(mesh->x[ii][jj+1] - mesh->x[ii][jj]);
	
	ii += 1;
    double dSx1 = (mesh->y[ii][jj+1] - mesh->y[ii][jj]);
    double dSy1 = -(mesh->x[ii][jj+1] - mesh->x[ii][jj]);

	*dSx = (dSx0 + dSx1)/2.;
	*dSy = (dSy0 + dSy1)/2.;

}

void meshCalcDSJ(MESH* mesh, int ii, int jj, double* dSx, double* dSy)
{

    double dSx0 = -(mesh->y[ii+1][jj] - mesh->y[ii][jj]);
    double dSy0 = (mesh->x[ii+1][jj] - mesh->x[ii][jj]);
	
	jj += 1;
    double dSx1 = -(mesh->y[ii+1][jj] - mesh->y[ii][jj]);
    double dSy1 = (mesh->x[ii+1][jj] - mesh->x[ii][jj]);

	*dSx = (dSx0 + dSx1)/2.;
	*dSy = (dSy0 + dSy1)/2.;

}



