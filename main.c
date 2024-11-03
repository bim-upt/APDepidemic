#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#define DEBUG

typedef enum{
    INFECTED = 0,
    SUSCEPTIBLE,
    IMMUNE
}status_e;

typedef enum{
    N = 0,
    S,
    E,
    W
}movementDirection_e;

typedef struct{
    int id;
    int x;
    int y;
    status_e status;
    movementDirection_e dir;
    int movmentAmplitute;
}person_t;



person_t *getPopulation(FILE *f, int size){
    person_t *persons = malloc(sizeof(person_t)*size);
    if(persons == NULL){
        perror("Could not initialize persons vector");
        exit(-1);
    }

    #ifdef DEBUG
        fprintf(stderr, "Reading persons:\n");
    #endif

    for(int i = 0; i < size; i++){
        //id x y status dir amplitute
        fscanf(f, "%d %d %d %d %d %d", &persons[i].id, &persons[i].x, &persons[i].y, (int*)&persons[i].status, (int*)&persons[i].dir, &persons[i].movmentAmplitute);
        #ifdef DEBUG
            fprintf(stderr, "%d %d %d %d %d %d\n", persons[i].id, persons[i].x, persons[i].y, (int)persons[i].status, (int)persons[i].dir, persons[i].movmentAmplitute);
        #endif
    }

    #ifdef DEBUG
        fprintf(stderr, "Finished reading persons\n");
    #endif

    return persons;
}

void initializeSimulationParameters(FILE *f, int *n, int *m, int *size){
    
    fscanf(f, "%d %d", n, m); //get grid size
    fscanf(f, "%d", size);  //get population size
    #ifdef DEBUG
        fprintf(stderr, "Simulation size: %dx%d, with population of %d\n", *n, *m, *size);
    #endif
}


int main(int argc, char **argv)
{
    int n, m, size;

    if(argc != 2){
        fprintf(stderr, "Required arguments: input_file\n");
        exit(-1);
    }

    //initialization
    FILE *f = fopen(argv[1], "r");
    if(f == NULL){
        perror("Could not open initialization file");
        exit(-1);
    }

    initializeSimulationParameters(f, &n, &m, &size);
    person_t *popluation = getPopulation(f, size);
    fclose(f);




    return 0;
}