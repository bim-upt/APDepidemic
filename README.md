
## Description
Create a simple simulation for an epidemic. Create four versions and find out which one is best: serial execution, parallel execution with manual data partitioning, parallel execution using parallel for, distributed execution using MPI.
## Input
- PersonID
- Initial coordinates x, y. they must be between 0..MAX_X_COORD, 0..MAX_Y_COORD
- Initial status: (infected=0, susceptible = 1) Initially we consider that there are no immune persons. For the initially infected persons, we consider that they got infected at the moment zero of the simulation.
- Movement pattern direction: (N=0, S=1, E=2, W=3)
- Movement pattern amplitude: an integer number, smaller than the area dimension on the movement direction

## Output
For each person:
   -  Final coordinates x, y
   -  Final status (infected, immune, susceptible)
- Infection counter: how many times during the simulation did the person become infected

## Documentation
Documentation with resulting conclusions are in the pdf. Graphs together with the simulation results are in the excel files. 
