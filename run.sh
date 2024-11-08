rounds=(50 100 150 200 500)
inputFile=("epidemics10K.txt" "epidemics20K.txt" "epidemics50K.txt" "epidemics100K.txt" "epidemics500K.txt")
threadNum=(2 3 4 5 6 7 8 9 10)

> results.txt
for round in "${rounds[@]}"; do
    for input in "${inputFile[@]}"; do
        for thread in "${threadNum[@]}"; do
            ./main "$round" "$input" "$thread" >> results.txt
        done
    done
done
