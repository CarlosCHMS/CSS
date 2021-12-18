
typedef struct
{

    int N;
    char** name;
    char** value;

} INPUT;

INPUT* inputInit(char* fileName);

void inputPrint(INPUT* input);

char* inputGetValue(INPUT* input, char* name);

void inputFree(INPUT* input);
