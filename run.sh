rounds=(50 100 150 200 500)
inputFile=("epidemics10K.txt" "epidemics20K.txt" "epidemics50K.txt" "epidemics100K.txt" "epidemics500K.txt")
threadNum=( 2 3 4 5 6)

> results.txt
for round in "${rounds[@]}"; do
    for input in "${inputFile[@]}"; do
        for thread in "${threadNum[@]}"; do
            mpiexec "-n" "$thread" "mpi_main" "$round" "$input" >> results.txt
        done
    done
done
