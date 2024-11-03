#include <stdio.h>
#include <stdlib.h>
#include <string.h>


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
}cardinalDirections_e;

typedef struct{
    int x;  //column
    int y;  //row
}coords_t;

typedef struct{
    int id;
    coords_t coords;
    status_e status;
    status_e futureStatus;
    cardinalDirections_e dir;
    int movmentAmplitute;
}person_t;

typedef struct{
    int data;
    //anticipating for the future when i test stuff
}paddedInt_t;


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
        fscanf(f, "%d %d %d %d %d %d", &persons[i].id, &persons[i].coords.x, &persons[i].coords.y, (int*)&persons[i].status, (int*)&persons[i].dir, &persons[i].movmentAmplitute);
        #ifdef DEBUG
            fprintf(stderr, "%d %d %d %d %d %d\n", persons[i].id, persons[i].coords.x, persons[i].coords.y, (int)persons[i].status, (int)persons[i].dir, persons[i].movmentAmplitute);
        #endif
    }

    #ifdef DEBUG
        fprintf(stderr, "Finished reading persons\n\n");
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

int coordsOutside(int x, int y, int n, int m){
    return x < 0 || y < 0 || x >= m || y >= n;
}


coords_t nextCoords(coords_t coords, int n, int m, cardinalDirections_e dir, int amplitute, cardinalDirections_e *newDirection){
    *newDirection = dir;

    if(dir == N){
        coords.y -= amplitute;
        if(coords.y < 0){
            coords.y *= -1;
            *newDirection = S;
        }
        return coords;
    }

    if(dir == S){
        coords.y += amplitute;
        if(coords.y >= n){
            coords.y = (n - 1) - (coords.y - (n - 1));
            *newDirection = N;
        }
        return coords;
    }

    if(dir == E){
        coords.x += amplitute;
        if(coords.x >= m){
            coords.x = (m - 1) - (coords.x - (m - 1));
            *newDirection = W;
        }
        return coords;

    }

    if(dir == W){
        coords.x -= amplitute;
        if(coords.x < 0){
            coords.x *= -1;
            *newDirection = E;
        }
        return coords;
    }
  

    //shouldn't be able to reach here 
    return coords;

}

void updatePositions(person_t *persons, int start, int end, int n, int m, int round){

    #ifdef DEBUG
        fprintf(stderr, "Coords for t = %d\n", round);
    #endif
    for(int i = start; i < end; i++){
        persons[i].coords = nextCoords(persons[i].coords, n, m, persons[i].dir, persons[i].movmentAmplitute, &persons[i].dir);
        #ifdef DEBUG
            fprintf(stderr, "%d | %d %d\n", persons[i].id, persons[i].coords.x, persons[i].coords.y);
        #endif
    }

}



void makeContagionZone(paddedInt_t *contagionZone, int x, int y, int m){
    contagionZone[y*m+x].data = 1;
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
    person_t *persons = getPopulation(f, size);
    fclose(f);
    
    //1 = someone infected is there, 0 = nope
    paddedInt_t *contagionZone = malloc(sizeof(paddedInt_t)*n*m);
    if(contagionZone == NULL){
        perror("Could not initiate contagious zones");
        exit(-1);
    }
    memset(contagionZone, 0, sizeof(paddedInt_t)*n*m);





    free(contagionZone);
    free(persons);
    return 0;
}