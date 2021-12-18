#include<stdio.h>
#include<stdlib.h>
#include"readTables.h"

char getWord(FILE* ff, char* s)
{

    int stop;
    int ii;
    char c;
    
    stop = 0;
    ii = 0;
    while(!stop)
    {
        c = fgetc(ff);
                
        if(c == ',')
        {
            stop = 1;
        }
        else if(c == '\n')
        {
            stop = 1;
        }
        else
        {
            s[ii] = c;
            ii++;
        }
        
    }
    s[ii] = '\0';

    return c;

}

void tableMalloc(TABLE* t)
{

    t->values = malloc(t->Nrow*sizeof(double*));

    for(int ii=0; ii<t->Nrow; ii++)
    {    
        t->values[ii] = malloc(t->Ncol*sizeof(double));
    }

}

TABLE* tableRead(FILE* ff)
{

    int ii;
    int jj;
    char s[100];

    TABLE* t = malloc(sizeof(TABLE));
    
    getWord(ff, s); // Read the number of the table
    t->N = atoi(s);
    getWord(ff, s); // Read the number of rows
    t->Nrow = atoi(s);
    getWord(ff, s); // Read the number of collums
    t->Ncol = atoi(s);

    getWord(ff, s);

    tableMalloc(t); // Allocate the values matrix
    
    for(ii=0; ii < t->Nrow; ii++)
    {
        for(jj=0; jj < t->Ncol; jj++)
        {
            getWord(ff, s);
            t->values[ii][jj] = strtod(s, NULL);
        }
        getWord(ff, s);
    }    
 
    return t;

}

void tableDisplay(TABLE* t)
{

    for(int ii=0; ii < t->Nrow; ii++)
    {
        for(int jj=0; jj < t->Ncol; jj++)
        {            
            printf(" %.4f,", t->values[ii][jj]);
        }
        printf("\n");
    }      

}

TABLELIST* fReadTables(char* fileName)
{

    int ii;
    FILE* ff = fopen(fileName, "r");

    TABLELIST* tl = malloc(sizeof(TABLELIST));
    TABLE* t;
    char c;
    char s[100];
    
    //Read number of tables
    c = getWord(ff, s);
    tl->N = atoi(s);
    c = getWord(ff, s);
    
    //Allocate the tables vector
    tl->tables = malloc(tl->N*sizeof(TABLE**));

    for(ii=0; ii<tl->N; ii++)
    {
        t = tableRead(ff);   
        tl->tables[ii] = t;
    }

    fclose(ff);
    
    return tl;

}

