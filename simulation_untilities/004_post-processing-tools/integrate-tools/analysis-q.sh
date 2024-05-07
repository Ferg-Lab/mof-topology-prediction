#!/bin/bash

# Get the name of the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
rm ../q_all.dat
lambda=(0.00000001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

# Iterate over the directories in the parent directory
for dir in ../*/; do
    dir_name="$(basename "$dir")"
    # Check if the directory is the same as the script directory
    if [ "$dir_name" != "$(basename "$SCRIPT_DIR")" ]; then
        # Pass the directory name as $1 to your script
        rm "../$dir_name/2_coul/averaged.dat"
        read -r divisor < "../$dir_name/2_coul/fu.dat"
        for n in $(seq 0 10)
        do
            python block-q.py ${lambda[n]} "../$dir_name/2_coul/q_${lambda[n]}.lmp" "$dir_name"
        done
        python2 trapz.py -d "../$dir_name/2_coul/averaged.dat" | sed -e 's/\+\/-//g' -e 's/^/'"$dir_name"' /' | awk -v divisor="$divisor" '{printf("%s %.2f %.2f\n", $1, $2/divisor, $3/divisor)}' >> "../q_all.dat"
    fi
done
