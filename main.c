#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>


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

typedef struct{
    int n;
    int m;
    int populationSize;
    int rounds;
    int threadNum;
    int rank;
    person_t *persons;
    paddedInt_t *contagionZone;
}epidemic_DTO;


person_t *getPopulation(FILE *f, int populationSize){
    person_t *persons = malloc(sizeof(person_t)*populationSize);
    if(persons == NULL){
        perror("Could not initialize persons vector");
        exit(-1);
    }

    #ifdef DEBUG
        fprintf(stderr, "Reading persons:\n");
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
            fprintf(stderr, "%d %d %d %d %d %d\n", persons[i].id, persons[i].coords.x, persons[i].coords.y, (int)persons[i].status.status, (int)persons[i].dir, persons[i].movementAmplitute);
        #endif
    }

    #ifdef DEBUG
        fprintf(stderr, "Finished reading persons\n\n");
    #endif

    return persons;
}

void initializeSimulationParameters(FILE *f, int *n, int *m, int *populationSize){
    
    fscanf(f, "%d %d", n, m); //get grid populationSize
    fscanf(f, "%d", populationSize);  //get population populationSize
    #ifdef DEBUG
        fprintf(stderr, "Simulation populationSize: %dx%d, with population of %d\n", *n, *m, *populationSize);
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



void serialEpidemic(int n, int m, int populationSize, person_t *persons, int rounds){

    //1 = someone infected is there, 0 = nope
    paddedInt_t *contagionZone = malloc(sizeof(paddedInt_t)*n*m);
    if(contagionZone == NULL){
        perror("Could not initiate contagious zones");
        exit(-1);
    }
    memset(contagionZone, 0, sizeof(paddedInt_t)*n*m);

    for(int i = 0; i < rounds; i++){
        updatePositionsAndContagionZone(persons, 0, populationSize, n, m, i, contagionZone);
        updateFutureStatus(contagionZone, persons, 0, populationSize, m);
        updateStatus(persons, 0, populationSize, i);
        #ifdef DEBUG
            printVector(contagionZone,n,m);
        #endif
        resetContagionZone(contagionZone, 0, m*n);
    }
    free(contagionZone);
}


void *parallelEpidemic_thread(void *payload){
    epidemic_DTO data = *(epidemic_DTO*)payload;
    int rounds = data.rounds;
    person_t *persons = data.persons;
    int populationSize = data.populationSize;
    int n = data.n;
    int m = data.m;
    paddedInt_t *contagionZone = data.contagionZone;
    int rank = data.rank;
    int threadNum = data.threadNum;

    int workloadSize = populationSize/threadNum;
    int rest;
    if (rank == threadNum - 1){
        rest = n % threadNum;
    }else{
        rest = 0;
    }
    int start = rank*workloadSize;
    int end = (rank+1)*workloadSize + rest;

    for(int i = 0; i < rounds; i++){
        updatePositionsAndContagionZone(persons, start, end, n, m, i, contagionZone);
        
        updateFutureStatus(contagionZone, persons, 0, populationSize, m);
        updateStatus(persons, 0, populationSize, i);
        resetContagionZone(contagionZone, 0, m*n);
    }


}

epidemic_DTO dataToEpidemicDTO(int n, int m, int populationSize, person_t *persons, int rounds, int threadNum, int rank, paddedInt_t *contagion){
    epidemic_DTO payload;
    payload.n = n;
    payload.m = m;
    payload.persons = persons;
    payload.rounds = rounds;
    payload.populationSize = populationSize;
    payload.rank = rank;
    payload.threadNum = threadNum;
    payload.contagionZone = contagion;
    return payload;
}

void parallelEpidemic(int n, int m, int populationSize, person_t *persons, int rounds, int threadNum){
    pthread_t *threads = malloc(sizeof(pthread_t)*threadNum);
    epidemic_DTO *payloads = malloc(sizeof(epidemic_DTO)*threadNum);
    paddedInt_t *contagionZone = malloc(sizeof(paddedInt_t)*n*m);
    if(contagionZone == NULL){
        perror("Could not initiate contagious zones");
        exit(-1);
    }
    memset(contagionZone, 0, sizeof(paddedInt_t)*n*m);
    if(threads == NULL || payloads == NULL){
        perror("Could not create thread array");
        exit(-1);
    }

    for(int i = 0; i < threadNum; i++){
        epidemic_DTO auxPayload = dataToEpidemicDTO(n, m, populationSize, persons, rounds, threadNum, i, contagionZone);
        payloads[i] = auxPayload;
        pthread_create(&threads[i], NULL, parallelEpidemic_thread, (void *)&payloads[i]);
    }

    for (int i = 0; i < threadNum; i++)
    {
        pthread_join(threads[i], NULL);
    }

    free(contagionZone);
}



int main(int argc, char **argv)
{
    int n, m, populationSize;
   
    if(argc != 4){
        fprintf(stderr, "Required arguments: simulation_time input_file thread_num\n");
        exit(-1);
    }

    int simulationTime;
    int threadNum;
    sscanf(argv[1], "%d", &simulationTime);
    sscanf(argv[3], "%d", &threadNum);

    //initialization
    FILE *f = fopen(argv[2], "r");
    if(f == NULL){
        perror("Could not open initialization file");
        exit(-1);
    }

    initializeSimulationParameters(f, &n, &m, &populationSize);
    person_t *persons = getPopulation(f, populationSize);
    fclose(f);
        
    
    //serial
    serialEpidemic(n, m, populationSize, persons, simulationTime);





    free(persons);
    return 0;
}