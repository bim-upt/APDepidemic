#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <omp.h>


//#define DEBUG
#define IMMUNE_DURATION 2
#define INFECTED_DURATION 7
#define MAX_PATH 256
#define PATH_EXTENSION_SIZE 4
#define PADDING_SIZE 60

#define CHUNKSIZE 800
#define POLICY dynamic


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
    char padding[PADDING_SIZE];
}paddedInt_t;








person_t *getPopulation(FILE *f, int populationSize){
    person_t *persons = malloc(sizeof(person_t)*populationSize);
    if(persons == NULL){
        perror("Could not initialize persons vector");
        exit(-1);
    }

    #ifdef DEBUG
        fprintf(stdout, "Reading persons:\n");
    #endif

    for(int i = 0; i < populationSize; i++){
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
            fprintf(stdout, "%d %d %d %d %d %d\n", persons[i].id, persons[i].coords.x, persons[i].coords.y, (int)persons[i].status.status, (int)persons[i].dir, persons[i].movementAmplitute);
        #endif
    }

    #ifdef DEBUG
        fprintf(stdout, "Finished reading persons\n\n");
    #endif

    return persons;
}

void initializeSimulationParameters(FILE *f, int *n, int *m, int *populationSize){
    
    fscanf(f, "%d %d", m, n); //get grid size
    fscanf(f, "%d", populationSize);  //get populationSize
    #ifdef DEBUG
        fprintf(stdout, "Simulation populationSize: %dx%d, with population of %d\n", *n, *m, *populationSize);
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

void makeContagionZone(paddedInt_t *contagionZone, coords_t coords, int m){
    contagionZone[coords.y*m+coords.x].data = 1;
}


void updatePositionsAndContagionZone(person_t *persons, int start, int end, int n, int m, int round, paddedInt_t *contagionZone){

    for(int i = start; i < end; i++){
        persons[i].coords = nextCoords(persons[i].coords, n, m, persons[i].dir, persons[i].movementAmplitute, &persons[i].dir);
        if(persons[i].status.status == INFECTED && !contagionZone[persons[i].coords.y*m + persons[i].coords.x].data){
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
            if(coordIsContagious(contagionZone, persons[i].coords, m)){
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

void updateStatus(person_t *persons, int start, int end, int rank){
    
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
            fprintf(stdout, "thread %d - id: %d | (%d,%d) | %s -> %s | imn time: %d | inf time: %d | inf count: %d\n", rank, persons[i].id, persons[i].coords.x, persons[i].coords.y, status, nxtStatus, persons[i].status.immuneTime, persons[i].status.infectionTime, persons[i].status.infectionCounter);
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
            fprintf(stdout, "%d ", v[i*m+j].data);
        }
        fprintf(stdout, "\n");
    }
}



void serialEpidemic(int n, int m, int populationSize, person_t *persons, int rounds){

    //1 = someone infected is there, 0 = nope
    paddedInt_t *contagionZone = malloc(sizeof(paddedInt_t)*n*m);
    if(contagionZone == NULL){
        perror("Could not initiate contagious zones");
        exit(-1);
    }
    memset(contagionZone, 0, sizeof(paddedInt_t)*n*m);

    for(int i = 0; i < rounds; i++){
        #ifdef DEBUG
            fprintf(stdout, "\nt = %d\n", i);
        #endif
        updatePositionsAndContagionZone(persons, 0, populationSize, n, m, i, contagionZone);
        updateFutureStatus(contagionZone, persons, 0, populationSize, m);
        updateStatus(persons, 0, populationSize, 0);
        #ifdef DEBUG
            printVector(contagionZone,n,m);
        #endif
        resetContagionZone(contagionZone, 0, m*n);
    }
    free(contagionZone);
}










void writePersonsToFile(person_t *persons, FILE *f, int populationSize){
    fprintf(f, "id x y status inf_cout\n");
    for(int i = 0; i < populationSize; i++){
        char *status;
        if(persons[i].status.status == SUSCEPTIBLE){
            status = "SUSCEPTIBLE";
        }else if(persons[i].status.status == IMMUNE){
            status = "IMMUNE";
        }else{
            status = "INFECTED";
        }
        fprintf(f, "%d %d %d %s %d\n", persons[i].id, persons[i].coords.x, persons[i].coords.y, status, persons[i].status.infectionCounter);
    }
}

void writeResult(person_t *persons, int populationSize, char *path, char *extension){
    char newPath[MAX_PATH];
    strcpy(newPath, path);
    newPath[strlen(newPath) - PATH_EXTENSION_SIZE] = 0;
    strcat(newPath,extension);

 
    FILE *f = fopen(newPath, "w");
    if(f == NULL){
        perror("Could not create output file");
        exit(-1);
    }

    writePersonsToFile(persons, f, populationSize);

    fclose(f);
}

person_t *copyPersonsVector(person_t *persons, int populationSize){
    person_t *copy = malloc(sizeof(person_t)*populationSize);
    if(copy == NULL){
        perror("Could not create persons copy");
        exit(-1);
    }
    for(int i = 0; i < populationSize; i++){
        copy[i] = persons[i];
    }
    return copy;
}

int samePerson(person_t a, person_t b){
    return a.coords.x == b.coords.x && a.coords.y == b.coords.y && a.status.status == b.status.status && a.status.infectionCounter == b.status.infectionCounter && a.id == b.id;
}

int personsVectorsAreEqual(person_t *v1, person_t *v2, int populationSize){
    for(int i = 0; i < populationSize; i++){
        if(!samePerson(v1[i], v2[i])){
            return 0;
        }
    }
    return 1;
}


//https://github.com/bim-upt/APDepidemic
int main(int argc, char **argv)
{
    int n, m, populationSize;
    struct timespec start, finish;
    double elapsedSerial, elapsedMPI;


    if(argc != 3){
        fprintf(stdout, "Required arguments: simulation_time input_file, only %d arguments are present\n", argc-1);
        exit(-1);
    }

    int simulationTime;
    sscanf(argv[1], "%d", &simulationTime);

    //initialization
    char *path = argv[2];
    FILE *f = fopen(path, "r");
    if(f == NULL){
        perror("Could not open initialization file");
        exit(-1);
    }

    initializeSimulationParameters(f, &n, &m, &populationSize);
    person_t *personsSerial = getPopulation(f, populationSize);
    fclose(f);
    person_t *personsMPI = copyPersonsVector(personsSerial, populationSize);

        
    
    //serial
    clock_gettime(CLOCK_MONOTONIC, &start);
    serialEpidemic(n, m, populationSize, personsSerial, simulationTime);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsedSerial = (finish.tv_sec - start.tv_sec);
    elapsedSerial += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

   

 
    

    writeResult(personsSerial, populationSize, path, "_serial_out.txt");
    
  
    
    //fprintf(stdout, "populationSize, simulationTime, threads, t_serial, t_parallel, speedup\n");
    //fprintf(stdout, "%d, %d, %d, %f, %f, %f\n", populationSize, simulationTime, threadNum, elapsedSerial, elapsedParallelV2, elapsedSerial/elapsedParallelV2);

    //fprintf(stdout, "%s, %d, %d, %d, %f\n", "dynamic", CHUNKSIZE, populationSize, simulationTime, elapsedParallelV1);

    //fprintf(stdout, "t_serial = %f\n\nt_parallelV1 = %f\nspeedupV1 = %f\n\nt_parallelV2 = %f\nspeedupV2 = %f\n", elapsedSerial, elapsedParallelV1, elapsedSerial/elapsedParallelV1, elapsedParallelV2, elapsedSerial/elapsedParallelV2);

    free(personsMPI);
    free(personsSerial);
    return 0;
}