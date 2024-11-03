#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define DEBUG
#define INFECTION_RATE 1
#define IMMUNE_DURATION 0
#define INFECTED_DURATION 50

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
    status_e status;
    int immuneTime;
    int infectionTime;
    int infectionCounter;
}status_t;

typedef struct{
    int id;
    coords_t coords;
    status_t status;
    status_t futureStatus;
    cardinalDirections_e dir;
    int movementAmplitute;
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
        fscanf(f, "%d %d %d %d %d %d", &persons[i].id, &persons[i].coords.x, &persons[i].coords.y, (int*)&persons[i].status, (int*)&persons[i].dir, &persons[i].movementAmplitute);
        
        persons[i].status.immuneTime = 0;
        persons[i].status.infectionTime = 0;
        persons[i].status.infectionCounter = 0;
        if(persons[i].status.status == INFECTED){
            persons[i].status.infectionCounter = 1;
            persons[i].status.infectionTime = INFECTED_DURATION;
        }

        #ifdef DEBUG
            fprintf(stderr, "%d %d %d %d %d %d\n", persons[i].id, persons[i].coords.x, persons[i].coords.y, (int)persons[i].status.status, (int)persons[i].dir, persons[i].movementAmplitute);
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

    for(int i = start; i < end; i++){
        persons[i].coords = nextCoords(persons[i].coords, n, m, persons[i].dir, persons[i].movementAmplitute, &persons[i].dir);
    }

}

void makeContagionZone(paddedInt_t *contagionZone, coords_t coords, int m){
    contagionZone[coords.y*m+coords.x].data = 1;
}

void updateContagionZone(person_t *persons, paddedInt_t *contagionZone, int start, int end, int m){
    for(int i = start; i < end; i++){
        if(persons[i].status.status == INFECTED){
            makeContagionZone(contagionZone, persons[i].coords, m);
        }
    }
}


int coordIsContagious(paddedInt_t *contagionZone, coords_t coords, int m){
    return contagionZone[coords.y*m + coords.x].data;
}




void updateFutureStatus(paddedInt_t *contagionZone, person_t *persons, int start, int end, int m){
    for(int i = start; i < end; i++){
        persons[i].futureStatus = persons[i].status;
        if(persons[i].status.status == SUSCEPTIBLE) {
            if(coordIsContagious(contagionZone, persons[i].coords, m) && rand()/(float)RAND_MAX <= INFECTION_RATE){
                persons[i].futureStatus.status = INFECTED;
                persons[i].futureStatus.infectionCounter++;
                persons[i].futureStatus.infectionTime = INFECTED_DURATION;
            }
        }

        if(persons[i].status.status == IMMUNE){
            if(persons[i].status.immuneTime - 1 <= 0){
                persons[i].futureStatus.status = SUSCEPTIBLE;
            }
            persons[i].futureStatus.immuneTime = persons[i].status.immuneTime > 0 ? persons[i].status.immuneTime - 1 : 0;
        }

        if(persons[i].status.status == INFECTED){
            if(persons[i].status.infectionTime - 1 <= 0){
                persons[i].futureStatus.status = IMMUNE;
                persons[i].futureStatus.immuneTime = IMMUNE_DURATION;
            }
            persons[i].futureStatus.infectionTime = persons[i].status.infectionTime > 0 ? persons[i].status.infectionTime - 1 : 0;
        }
    }
}

void updateStatus(person_t *persons, int start, int end, int round){
    #ifdef DEBUG
        fprintf(stderr, "\n\nstatus for t = %d\n",round);
    #endif
    for(int i = start; i < end; i++){
        #ifdef DEBUG
            char *status;
            char *nxtStatus;
            if(persons[i].status.status == SUSCEPTIBLE){
                status = "SUSCEPTIBLE";
            }else if(persons[i].status.status == IMMUNE){
                status = "IMMUNE";
            }else{
                status = "INFECTED";
            }

            if(persons[i].futureStatus.status == SUSCEPTIBLE){
                nxtStatus = "SUSCEPTIBLE";
            }else if(persons[i].futureStatus.status == IMMUNE){
                nxtStatus = "IMMUNE";
            }else{
                nxtStatus = "INFECTED";
            }
            fprintf(stderr, "id: %d | (%d,%d) | %s -> %s | imn time: %d | inf time: %d | inf count: %d\n", persons[i].id, persons[i].coords.x, persons[i].coords.y, status, nxtStatus, persons[i].status.immuneTime, persons[i].status.infectionTime, persons[i].status.infectionCounter);
        #endif
        persons[i].status = persons[i].futureStatus;
    }
}

void resetContagionZone(paddedInt_t *contagionZone, int start, int nmem){
    memset(contagionZone + start, 0, nmem*sizeof(paddedInt_t));
}

void printVector(paddedInt_t *v, int n, int m){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            fprintf(stderr, "%d ", v[i*m+j].data);
        }
        fprintf(stderr, "\n");
    }
}

void serialEpidemic(int n, int m, int size, person_t *persons, int rounds){
    //1 = someone infected is there, 0 = nope
    paddedInt_t *contagionZone = malloc(sizeof(paddedInt_t)*n*m);
    if(contagionZone == NULL){
        perror("Could not initiate contagious zones");
        exit(-1);
    }
    memset(contagionZone, 0, sizeof(paddedInt_t)*n*m);

    for(int i = 0; i < rounds; i++){
        updatePositions(persons, 0, size, n, m, i);
        updateContagionZone(persons, contagionZone, 0, size, m);
        //in paralel a barrier here probs
        updateFutureStatus(contagionZone, persons, 0, size, m);
        updateStatus(persons, 0, size, i);
        //fprintf(stderr, "%d\n", (int)persons[0].status.status);
        #ifdef DEBUG
            //printVector(contagionZone,n,m);
        #endif
        resetContagionZone(contagionZone, 0, m*n);
    }


    free(contagionZone);
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
    
    serialEpidemic(n, m, size, persons, 2);





    free(persons);
    return 0;
}