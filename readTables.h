
typedef struct{

    int N; // Table number
    int Ncol; // Number of collums
    int Nrow; // Number of lines
    double **values;

} TABLE;

typedef struct{

    int N; // Number of tables
    TABLE** tables; // Tables pointers vector

} TABLELIST;


char getWord(FILE* ff, char* s);

void tableMalloc(TABLE* t);

TABLE* tableRead(FILE* ff);

void tableDisplay(TABLE* t);

TABLELIST* fReadTables(char* fileName);

/*

int main()
{

    TABLELIST* tl;

    tl = fReadTables("./mesh.csv");    

    tableDisplay(tl->tables[1]);

    return 0;

}

*/
