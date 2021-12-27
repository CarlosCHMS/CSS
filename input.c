#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"input.h"


INPUT* inputInit(char* fileName, int N)
{

    INPUT* input = malloc(sizeof(INPUT));
    FILE* ff = fopen(fileName, "r");
    char c;
    int name, value;
    int ii, jj;
    char s[N];
    
    input->Nmax = N;
    
    input->name = (char**)malloc(N*sizeof(char*));    
    for(ii=0; ii<N; ii++)
    {
        input->name[ii] = (char*)malloc(N*sizeof(char));
        input->name[ii][0] = '\0';
    }
    
    input->value = (char**)malloc(N*sizeof(char*));    
    for(ii=0; ii<N; ii++)
    {
        input->value[ii] = (char*)malloc(N*sizeof(char));
        input->value[ii][0] = '\0';
    }
    
    
    name = 1;
    value = 0;
    ii = 0;
    jj = 0;
    while(c != EOF)
    {
    
        c = fgetc(ff);
        if(c != ' ')
        {
        
            if(c == ',')
            {    
                name = 0;
                value = 1;
                input->name[ii][jj] = '\0';
                jj = 0;
            }
            else if(c == '#')
            {
                name = 0;
                value = 0;
                input->value[ii][jj] = '\0';                            
            }
            else if(c == '\n')
            {
                if(value==1)
                {
                    name = 1;
                    value = 0;
                    input->value[ii][jj] = '\0';                
                    ii += 1;
                    jj = 0;
                }
                if(name==1 && jj>0)
                {
                    input->name[ii][jj] = '\0';                
                    ii += 1;
                    jj = 0;
                }
                
            }            
            else if(name)
            {   
                input->name[ii][jj] = c;
                jj += 1;
            }
            else if(value)
            {
                input->value[ii][jj] = c;
                jj += 1;
            }
                    
        }
    
    }
    
    input->N = ii;
    
    return input;

}

void inputPrint(INPUT* input)
{

    for(int ii=0; ii<input->N; ii++)
    {
    
        printf(" %s, %s\n", input->name[ii], input->value[ii]);
    
    }

}

char* inputGetValue(INPUT* input, char* name)
{
    
    
    char* s;
    int found = 0;
    
    for(int ii=0; ii<input->N; ii++)
    {                
        if(strcmp(name, input->name[ii]) == 0)
        {
            s = input->value[ii];
            found = 1;
        }   
    }
    
    if(found == 0)
    {
        printf("Error: Value not found in the input file.\n");
        exit(0);
    }

    return s;
    
}

void inputFree(INPUT* input)
{
    
    for(int ii=0; ii<input->Nmax; ii++)
    {
        free(input->name[ii]);
        free(input->value[ii]);
    }
    free(input->name);
    free(input->value);
    
} 

/*

int main()
{

    INPUT* input = initInput("./input.csv");

    inputPrint(input);

    printf("\n %s \n", inputGetValue(input, "p"));

    return 0;

}

*/
