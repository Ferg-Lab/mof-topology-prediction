#!/bin/bash

#### Changing the types #######
# The types of atoms, bonds   #
# and others are changed to   #
# add in the ones for cations #
###############################

# Path of the directory containing the subdirectories
dir_path="./data"

# Define the new values
new_atom_types=7
new_bond_types=6
new_angle_types=9
new_dihedral_types=3
new_improper_types=1

for file in "$dir_path"/*; do
    if [ -f "$file" ]; then
        dataname=$(basename "$file" | cut -d '.' -f 2 | cut -d '_' -f 1)
        mkdir -p "./data/$dataname" && cp "$file" "./data/${dataname}/data.$dataname"
	cd ./data/${dataname}

    	gsed -i -E "s/^[[:space:]]*[0-9]+[[:space:]]*atom[[:space:]]*types/$new_atom_types atom types/" data."$dataname"
    	gsed -i -E "s/^[[:space:]]*[0-9]+[[:space:]]*bond[[:space:]]*types/$new_bond_types bond types/" data."$dataname"
    	gsed -i -E "s/^[[:space:]]*[0-9]+[[:space:]]*angle[[:space:]]*types/$new_angle_types angle types/" data."$dataname"
    	gsed -i -E "s/^[[:space:]]*[0-9]+[[:space:]]*dihedral[[:space:]]*types/$new_dihedral_types dihedral types/" data."$dataname"
    	gsed -i -E "s/^[[:space:]]*[0-9]+[[:space:]]*improper[[:space:]]*types/$new_improper_types improper types/" data."$dataname"
    
    	gsed -i '/^Masses /,/^Atoms$/ { /^Masses /d; { /^Atoms$/!d } }' data."$dataname"
	cd ../../
    fi
done

echo "The data files are now ready. Move the folder to 002_deposit_cations for subsequent simulations"