#!/bin/bash


dir_path="./"
lammps="/Users/LAMMPS/src/lmp_mpi"

for subdir in "$dir_path"/*/; do

    dirname=$(basename "$subdir")
    echo "running simulation on $dirname"
    cd "$subdir"
    cp ../in.deposit .
    
    if [ -e "data.${dirname}_TMA" ]; then
    echo "Cations are happy now in ${dirname}"
    else
    ${lammps} -in in.deposit -screen none -var topo ${dirname} -var cation TMA
    fi

    cd ..
done

