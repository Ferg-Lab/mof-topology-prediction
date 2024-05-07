#!/bin/bash

# Get the name of the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
rm ../LJ_all.dat
lambda=(0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 0.45 0.475 0.5 0.525 0.55 0.575 0.6 0.625 0.65 0.675 0.7 0.725 0.75 0.775 0.8 0.825 0.85 0.875 0.9 0.925 0.95 0.975 1.0)
# Iterate over the directories in the parent directory
for dir in ../*/; do
    dir_name="$(basename "$dir")"
    # Check if the directory is the same as the script directory
    if [ "$dir_name" != "$(basename "$SCRIPT_DIR")" ]; then
        # Pass the directory name as $1 to your script
        rm "../$dir_name/1_LJ/averaged.dat"
        read -r divisor < "../$dir_name/1_LJ/fu.dat"
        for n in $(seq 0 40)
        do
            python block-LJ.py ${lambda[n]} "../$dir_name/1_LJ/LJ_${lambda[n]}.lmp" "$dir_name"
        done
        python2 trapz.py -d "../$dir_name/1_LJ/averaged.dat" | sed -e 's/\+\/-//g' -e 's/^/'"$dir_name"' /' | awk -v divisor="$divisor" '{printf("%s %.2f %.2f\n", $1, $2/divisor, $3/divisor)}' >> "../LJ_all.dat"
    fi
done
