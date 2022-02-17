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


int main(int argc, char **argv)
{

    SOLVER* solver = malloc(sizeof(SOLVER));
    
    struct timeval start;
    struct timeval stop;
    
    char s[50];
       
    // Load input   
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "input.dat");
    INPUT* input = inputInit(s, 50);
    printf("Input data:\n");
    inputPrint(input);
   
    // Set number of threads
    omp_set_num_threads(atoi(inputGetValue(input, "threads")));   
        
    // Load mesh    
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "mesh.csv");
    solver->mesh = meshInit(s);
    solver->Nrow = solver->mesh->Nrow-1;
    solver->Ncol = solver->mesh->Ncol-1;        

    // axisymmetric
    solver->mesh->axi = atoi(inputGetValue(input, "axisymmetric"));

	//meshPrintOmega(solver->mesh);
    //meshPrintDStotal(solver->mesh);

    // Memory allocation
    solverAllocate(solver);
 
    // Boundary conditions
    solver->bc = malloc(sizeof(BOUNDARY));
    solver->bc->down = inputGetValue(input, "BCdown");
    solver->bc->up = inputGetValue(input, "BCup");
    solver->bc->left = inputGetValue(input, "BCleft");
    solver->bc->right = inputGetValue(input, "BCright");

    boundarySet(solver->bc);
  
    if(inputNameIsInput(input, "pout"))
    {
        solver->pout = strtod(inputGetValue(input, "pout"), NULL);     
    }      
  
    // Constants
    solver->Rgas = 287.5;
    solver->gamma = 1.4;  
    solver->eFix = 0.1;
    solver->e = strtod(inputGetValue(input, "interpE"), NULL);
    
    // Seletion of MUSCL and flux
    solver->MUSCL = atoi(inputGetValue(input, "order")) - 1;
    solver->flux = fluxChoice(inputGetValue(input, "flux"));
    solver->stages = atoi(inputGetValue(input, "stages"));

    if(atoi(inputGetValue(input, "tube")) == 0)
    {

        //stead state
        solver->inlet = conditionInit(strtod(inputGetValue(input, "pressure"), NULL), 
                                      strtod(inputGetValue(input, "temperature"), NULL), 
                                      strtod(inputGetValue(input, "mach"), NULL), 
                                      strtod(inputGetValue(input, "nx"), NULL),
                                      strtod(inputGetValue(input, "ny"), NULL));
        
        // Initialization of U
        solverInitU(solver, solver->inlet);
        
        // Calculate time step        
        int Nmax = atoi(inputGetValue(input, "Nmax"));

        // Run the solver
        printf("\nRunning solution:\n");
        gettimeofday(&start, NULL);
        for(int ii=0; ii<Nmax; ii++)
        {
            solver->dt = solverCalcDt(solver);
            solverStepRK(solver);
            
            if(ii%100 == 0)
            {
                printf("%i, ", ii);
                solverCalcRes(solver);
            }
            
        }
        gettimeofday(&stop, NULL);  
        printf("Duration %f s\n", duration(start, stop));
 
    }
    else
    {
    
        CONDITION* inside1 = conditionInit(strtod(inputGetValue(input, "pressure1"), NULL), 
                                           strtod(inputGetValue(input, "temperature1"), NULL), 
                                           strtod(inputGetValue(input, "mach1"), NULL), 
                                           strtod(inputGetValue(input, "nx1"), NULL),
                                           strtod(inputGetValue(input, "ny1"), NULL));

        CONDITION* inside2 = conditionInit(strtod(inputGetValue(input, "pressure2"), NULL), 
                                           strtod(inputGetValue(input, "temperature2"), NULL), 
                                           strtod(inputGetValue(input, "mach2"), NULL), 
                                           strtod(inputGetValue(input, "nx2"), NULL),
                                           strtod(inputGetValue(input, "ny2"), NULL));      
    
        // Initialization of U
        solverInitUTube(solver, inside1, inside2, strtod(inputGetValue(input, "xm"), NULL));
        free(inside1);
        free(inside2);
        
        double tmax = strtod(inputGetValue(input, "tmax"), NULL);                

        // Run the solver
        double t = 0.0;
        printf("\nRunning solution:\n");
        gettimeofday(&start, NULL);
        int stopLoop = 0;
        int ii = 0;
        while(stopLoop == 0)
        {
            
            solver->dt = solverCalcDt(solver);
            
            if(t + solver->dt>tmax)
            {
                solver->dt = (tmax-t);
                stopLoop = 1;
            }

            solverStepRK(solver);
            t += solver->dt;
            ii++;

            if(ii%100 == 0)
            {
                printf("%i, ", ii);
                solverCalcRes(solver);
            }

        }
        gettimeofday(&stop, NULL);  
        printf("Duration %f s\n", duration(start, stop));
        printf("time %f s\n", t);

    }

    // Save solution
    s[0] = '\0';
    strcat(s, argv[1]);
    strcat(s, "solution.csv");
    solverWriteU(solver, s);

    // Free memory 
    solverFree(solver);
    inputFree(input);

    return 0;

}
