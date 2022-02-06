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

void meshCalcDSup(MESH* mesh, int ii, int jj, double* dSx, double* dSy)
{

    double y;
    *dSx = -(mesh->y[ii+1][jj+1] - mesh->y[ii][jj+1]);
    *dSy = mesh->x[ii+1][jj+1] - mesh->x[ii][jj+1];
    
    if(mesh->axi == 1)
    {
        y = (mesh->y[ii+1][jj+1] + mesh->y[ii][jj+1])*0.5;
        *dSx *= y;
        *dSy *= y;        
    }

}

void meshCalcDSright(MESH* mesh, int ii, int jj, double* dSx, double* dSy)
{

    double y;
    *dSx = mesh->y[ii+1][jj+1] - mesh->y[ii+1][jj];
    *dSy = -(mesh->x[ii+1][jj+1] - mesh->x[ii+1][jj]);

    if(mesh->axi == 1)
    {
        y = (mesh->y[ii+1][jj+1] + mesh->y[ii+1][jj])*0.5;
        *dSx *= y;
        *dSy *= y;        
    }

}

void meshCalcDSdown(MESH* mesh, int ii, int jj, double* dSx, double* dSy)
{

    double y;
    *dSx = (mesh->y[ii+1][jj] - mesh->y[ii][jj]);
    *dSy = -(mesh->x[ii+1][jj] - mesh->x[ii][jj]);

    if(mesh->axi == 1)
    {
        y = (mesh->y[ii+1][jj] + mesh->y[ii][jj])*0.5;
        *dSx *= y;
        *dSy *= y;        
    }

}

void meshCalcDSleft(MESH* mesh, int ii, int jj, double* dSx, double* dSy)
{

    double y;
    *dSx = -(mesh->y[ii][jj+1] - mesh->y[ii][jj]);
    *dSy = mesh->x[ii][jj+1] - mesh->x[ii][jj];
     
    if(mesh->axi == 1)
    {
        y = (mesh->y[ii][jj+1] + mesh->y[ii][jj])*0.5;
        *dSx *= y;
        *dSy *= y;        
    }   

}

double meshCalcDSlateral(MESH* mesh, int ii, int jj)
{

    double xac, xbd, yac, ybd;
    
    xac = mesh->x[ii][jj+1] - mesh->x[ii+1][jj];
    yac = mesh->y[ii][jj+1] - mesh->y[ii+1][jj];    

    xbd = mesh->x[ii][jj] - mesh->x[ii+1][jj+1];
    ybd = mesh->y[ii][jj] - mesh->y[ii+1][jj+1];    

    return 0.5*(xac*ybd - xbd*yac);

}

void meshPrintDStotal(MESH* mesh)
{

    double dSx, dSy;
    double dSxt, dSyt;

    for(int ii=0; ii<mesh->Nrow-1; ii++)
    {
        for(int jj=0; jj<mesh->Ncol-1; jj++)
        {
            dSxt = 0.0;
            dSyt = 0.0;

            meshCalcDSup(mesh, ii, jj, &dSx, &dSy);
            dSxt += dSx;
            dSyt += dSy;

            meshCalcDSdown(mesh, ii, jj, &dSx, &dSy);
            dSxt += dSx;
            dSyt += dSy;

            meshCalcDSleft(mesh, ii, jj, &dSx, &dSy);
            dSxt += dSx;
            dSyt += dSy;
        
            meshCalcDSright(mesh, ii, jj, &dSx, &dSy);
            dSxt += dSx;
            dSyt += dSy;

            if(mesh->axi==1)
            {
                dSyt -= meshCalcDSlateral(mesh, ii, jj);
            }
   
            printf(" %.4e, %.4e,\n", dSxt, dSyt);
        }
    }
}

double meshCalcOmega(MESH* mesh, int ii, int jj)
{

    double xac, xbd, yac, ybd, ans;
    
    xac = mesh->x[ii][jj+1] - mesh->x[ii+1][jj];
    yac = mesh->y[ii][jj+1] - mesh->y[ii+1][jj];    

    xbd = mesh->x[ii][jj] - mesh->x[ii+1][jj+1];
    ybd = mesh->y[ii][jj] - mesh->y[ii+1][jj+1];    

    ans = 0.5*(xac*ybd - xbd*yac);

    if(mesh->axi == 1)
    {
        ans *= (mesh->y[ii][jj+1] + mesh->y[ii+1][jj] + mesh->y[ii][jj] + mesh->y[ii+1][jj+1])*0.25;
    }

    return ans;

}

void meshPrintOmega(MESH* mesh)
{

    for(int ii=0; ii<mesh->Nrow-1; ii++)
    {
        for(int jj=0; jj<mesh->Ncol-1; jj++)
        {
            printf(" %.4e,", meshCalcOmega(mesh, ii, jj));
        }
        printf("\n");
    }

}

void meshCalcDSI(MESH* mesh, int ii, int jj, double* dSx, double* dSy)
{

    double dSx0, dSy0, dSx1, dSy1;

    meshCalcDSleft(mesh, ii, jj, &dSx0, &dSy0);
    meshCalcDSright(mesh, ii, jj, &dSx1, &dSy1);

	*dSx = (dSx1 - dSx0)/2.;
	*dSy = (dSy1 - dSy0)/2.;

}

void meshCalcDSJ(MESH* mesh, int ii, int jj, double* dSx, double* dSy)
{

    double dSx0, dSy0, dSx1, dSy1;

    meshCalcDSdown(mesh, ii, jj, &dSx0, &dSy0);
    meshCalcDSup(mesh, ii, jj, &dSx1, &dSy1);

	*dSx = (dSx1 - dSx0)/2.;
	*dSy = (dSy1 - dSy0)/2.;

}



