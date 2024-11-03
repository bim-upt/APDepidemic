#include <stdio.h>
#include <stdbool.h>

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



int main()
{
    return 0;
}