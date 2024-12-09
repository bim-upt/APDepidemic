all:	main	generator_epidemics	mpi_main

main: main.c
	gcc -Wall -o main main.c -fopenmp

mpi_main: mpi_main.c
	mpicc -Wall -o mpi_main mpi_main.c

generator_epidemics: generator_epidemics.c
	gcc -Wall -o generator_epidemics generator_epidemics.c

