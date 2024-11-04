all:	main	generator_epidemics

main: main.c
	gcc -Wall -o main main.c


generator_epidemics: generator_epidemics.c
	gcc -Wall -o generator_epidemics generator_epidemics.c

