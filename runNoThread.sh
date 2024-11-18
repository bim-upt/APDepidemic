rounds=(50 100 150 200 500)
inputFile=("epidemics10K.txt" "epidemics20K.txt" "epidemics50K.txt" "epidemics100K.txt" "epidemics500K.txt")


start_time=$(date +%s.%N)

for repeat in {1..3}; do
    for round in "${rounds[@]}"; do
        for input in "${inputFile[@]}"; do
            ./main "$round" "$input" "1" >> results.txt
        done
    done
    echo "" >> results.txt
done

end_time=$(date +%s.%N)

elapsed_time=$(echo "$end_time - $start_time" | bc)


echo "$elapsed_time seconds"