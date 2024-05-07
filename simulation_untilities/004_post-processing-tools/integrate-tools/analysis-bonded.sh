#!/bin/bash

# Get the name of the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
rm ../bonded_all.dat
lambda=(0.00000001 0.0001 0.001 0.01 0.025 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)

# Iterate over the directories in the parent directory
for dir in ../*/; do
    dir_name="$(basename "$dir")"
    # Check if the directory is the same as the script directory
    if [ "$dir_name" != "$(basename "$SCRIPT_DIR")" ]; then
        # Pass the directory name as $1 to your script
        rm "../$dir_name/0_bonded/averaged.dat"
        divisor=1
        for n in $(seq 0 25)
        do
            python block-bonded.py ${lambda[n]} "../$dir_name/0_bonded/bond_${lambda[n]}.lmp" "$dir_name"
        done
        python2 trapz.py -d "../$dir_name/0_bonded/averaged.dat" | sed -e 's/\+\/-//g' -e 's/^/'"$dir_name"' /' | awk -v divisor="$divisor" '{printf("%s %.2f %.2f\n", $1, $2/divisor, $3/divisor)}' >> "../bonded_all.dat"
    fi
done
