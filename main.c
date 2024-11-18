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

#define CHUNKSIZE 25600
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

typedef struct{
    int n;
    int m;
    int populationSize;
    int rounds;
    int threadNum;
    int rank;
    person_t *persons;
    paddedInt_t *contagionZone;
    pthread_barrier_t *positionBarrier;
    #ifdef DEBUG
        pthread_barrier_t *printBarrier;
    #endif
    pthread_barrier_t *statusBarrier;
    
    pthread_barrier_t *nextIterationBarrier;
}epidemic_DTO;






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
        rest = populationSize % threadNum;
    }else{
        rest = 0;
    }
    int start = rank*workloadSize;
    int end = (rank+1)*workloadSize + rest;

    for(int i = 0; i < rounds; i++){
        #ifdef DEBUG
            if(rank == threadNum - 1){
                fprintf(stdout, "\nt = %d\n", i);
            }
        #endif
        updatePositionsAndContagionZone(persons, start, end, n, m, i, contagionZone);
        pthread_barrier_wait(data.positionBarrier);
        

        updateFutureStatus(contagionZone, persons, start, end, m);
        updateStatus(persons, start, end, rank);
        pthread_barrier_wait(data.statusBarrier);
        
        #ifdef DEBUG
            if(rank == threadNum - 1){
                printVector(contagionZone, n, m);
            }
            pthread_barrier_wait(data.printBarrier);
        #endif


        resetContagionZone(contagionZone, ((n*m)/threadNum)*rank, (n*m)/threadNum + (rank == threadNum - 1 ? (n*m) % threadNum : 0));
        pthread_barrier_wait(data.nextIterationBarrier);
    }

    return NULL;
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



//manual data partition
void parallelEpidemicV2(int n, int m, int populationSize, person_t *persons, int rounds, int threadNum){
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



    epidemic_DTO defaultPayload = dataToEpidemicDTO(n, m, populationSize, persons, rounds, threadNum, 0, contagionZone);
    defaultPayload.positionBarrier = malloc(sizeof(pthread_barrier_t));
    defaultPayload.nextIterationBarrier = malloc(sizeof(pthread_barrier_t));
    defaultPayload.statusBarrier = malloc(sizeof(pthread_barrier_t));
    #ifdef DEBUG
        defaultPayload.printBarrier = malloc(sizeof(pthread_barrier_t));
        if(defaultPayload.printBarrier == NULL){
            perror("Could not initialize barriers");
            exit(-1);
        }
    #endif
    if(defaultPayload.positionBarrier == NULL || defaultPayload.nextIterationBarrier == NULL || defaultPayload.statusBarrier == NULL){
        perror("Could not initialize barriers");
        exit(-1);
    }

    pthread_barrier_init(defaultPayload.positionBarrier, NULL, threadNum);
    pthread_barrier_init(defaultPayload.nextIterationBarrier, NULL, threadNum);
    pthread_barrier_init(defaultPayload.statusBarrier, NULL, threadNum);
    #ifdef DEBUG
        pthread_barrier_init(defaultPayload.printBarrier, NULL, threadNum);
    #endif


    for(int i = 0; i < threadNum; i++){
        
        payloads[i] = defaultPayload;
        payloads[i].rank = i;
        pthread_create(&threads[i], NULL, parallelEpidemic_thread, (void *)&payloads[i]);
    }

    for (int i = 0; i < threadNum; i++)
    {
        pthread_join(threads[i], NULL);
    }
    

    pthread_barrier_destroy(defaultPayload.positionBarrier);
    pthread_barrier_destroy(defaultPayload.nextIterationBarrier);
    pthread_barrier_destroy(defaultPayload.statusBarrier);
    #ifdef DEBUG
        pthread_barrier_destroy(defaultPayload.printBarrier);
        free(defaultPayload.printBarrier);
    #endif
    free(contagionZone);
    free(threads);
    free(defaultPayload.positionBarrier);
    free(defaultPayload.nextIterationBarrier);
    free(defaultPayload.statusBarrier);
    free(payloads);
}


//parallel for
void parallelEpidemicV1(int n, int m, int populationSize, person_t *persons, int rounds){

    //1 = someone infected is there, 0 = nope
    paddedInt_t *contagionZone = malloc(sizeof(paddedInt_t)*n*m);
    if(contagionZone == NULL){
        perror("Could not initiate contagious zones");
        exit(-1);
    }
    memset(contagionZone, 0, sizeof(paddedInt_t)*n*m);


    #pragma omp parallel
    for(int i = 0; i < rounds; i++){
        #ifdef DEBUG
            #pragma omp single
            fprintf(stdout, "\nt = %d\n", i);
            #pragma omp barrier
        #endif

        #pragma omp for schedule(POLICY, CHUNKSIZE)
        for(int j = 0; j < populationSize; j++){
            updatePositionsAndContagionZone(persons, j, j+1, n, m, i, contagionZone);
        }


        #pragma omp for schedule(POLICY, CHUNKSIZE)
        for(int j = 0; j < populationSize; j++){
            updateFutureStatus(contagionZone, persons, j, j+1, m);
            updateStatus(persons, j, j+1, omp_get_thread_num());
        }



        #ifdef DEBUG
            #pragma omp single
            printVector(contagionZone,n,m);
            #pragma omp barrier
        #endif

        
        resetContagionZone(contagionZone, ((n*m)/omp_get_num_threads())*omp_get_thread_num(), (n*m)/omp_get_num_threads() + (omp_get_thread_num() == omp_get_num_threads() - 1 ? (n*m) % omp_get_num_threads() : 0));
        #pragma omp barrier
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
    double elapsedSerial, elapsedParallelV2, elapsedParallelV1;


    if(argc != 4){
        fprintf(stdout, "Required arguments: simulation_time input_file thread_num, only %d arguments are present\n", argc-1);
        exit(-1);
    }

    int simulationTime;
    int threadNum;
    sscanf(argv[1], "%d", &simulationTime);
    sscanf(argv[3], "%d", &threadNum);

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
    person_t *personsParallelV2 = copyPersonsVector(personsSerial, populationSize);
    person_t *personsParallelV1 = copyPersonsVector(personsSerial, populationSize);

        
    
    //serial
    clock_gettime(CLOCK_MONOTONIC, &start);
    //serialEpidemic(n, m, populationSize, personsSerial, simulationTime);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsedSerial = (finish.tv_sec - start.tv_sec);
    elapsedSerial += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    //parallel V1
    clock_gettime(CLOCK_MONOTONIC, &start); 
    parallelEpidemicV1(n, m, populationSize, personsParallelV1, simulationTime);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsedParallelV1 = (finish.tv_sec - start.tv_sec);
    elapsedParallelV1 += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    /*
    //parallel V2 
    clock_gettime(CLOCK_MONOTONIC, &start); 
    parallelEpidemicV2(n, m, populationSize, personsParallelV2, simulationTime, threadNum);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsedParallelV2 = (finish.tv_sec - start.tv_sec);
    elapsedParallelV2 += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    */

    writeResult(personsSerial, populationSize, path, "_serial_out.txt");
    writeResult(personsParallelV1, populationSize, path, "_omp1_out.txt");
    writeResult(personsParallelV2, populationSize, path, "_omp2_out.txt");
    
    /*
    if(!personsVectorsAreEqual(personsParallelV1, personsSerial, populationSize)){
        fprintf(stdout, "WARNING: serial and parallel V1 differ\n");
    }else{
        fprintf(stdout, "serial and parallel V1 are the same\n");
    }
    

    if(!personsVectorsAreEqual(personsParallelV2, personsSerial, populationSize)){
        fprintf(stdout, "WARNING: serial and parallel V2 differ\n");
    }else{
        fprintf(stdout, "serial and parallel V2 are the same\n");
    }
    */
    //fprintf(stdout, "populationSize, simulationTime, threads, t_serial, t_parallel, speedup\n");
    //fprintf(stdout, "%d, %d, %d, %f, %f, %f\n", populationSize, simulationTime, threadNum, elapsedSerial, elapsedParallelV2, elapsedSerial/elapsedParallelV2);

    fprintf(stdout, "%s, %d, %d, %d, %f\n", "dynamic", CHUNKSIZE, populationSize, simulationTime, elapsedParallelV1);

    //fprintf(stdout, "t_serial = %f\n\nt_parallelV1 = %f\nspeedupV1 = %f\n\nt_parallelV2 = %f\nspeedupV2 = %f\n", elapsedSerial, elapsedParallelV1, elapsedSerial/elapsedParallelV1, elapsedParallelV2, elapsedSerial/elapsedParallelV2);

    free(personsParallelV2);
    free(personsParallelV1);
    free(personsSerial);
    return 0;
}