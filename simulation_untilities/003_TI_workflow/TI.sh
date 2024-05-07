#!/bin/bash

lammps="/home/mjianming/Softwares/lammps-23Jun2022/src/lmp_mpi"
dir_path="./"

b_lam=(0.00000001 0.0001 0.001 0.01 0.025 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)
l_lam=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)
q_lam=(0.00000001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
cation="TMA"

for subdir in "$dir_path"/*/; do
    dirname=$(basename "$subdir")
    cd "$subdir"
    
    echo "running bonded simulation in $dirname"
    mkdir 0_bonded; cd 0_bonded
    cp ../../in.bonded .
    
    for n in $(seq 0 25)
    do
      if [ -e "./data.restart_${dirname}_${cation}_${b_lam[n]}" ]; then
      echo "Job with lambda=${b_lam[n]} already completed. Skipping."
      else
      mpirun -np 16 ${lammps} -in in.bonded -screen none -var lambda ${b_lam[n]} -var topo ${dirname} -var cation ${cation}
      fi
    done
    cd ..
       
    echo "running LJ simulation in $dirname"
    mkdir 1_LJ; cd 1_LJ
    cp ../../in.LJ .

    for n in $(seq 0 20)
    do
      if [ -e "./data.restart_${dirname}_${cation}_${l_lam[n]}" ]; then
      echo "Job with lambda=${l_lam[n]} already completed. Skipping."
      else
      mpirun -np 16 ${lammps} -in in.LJ -screen none -var lambda ${l_lam[n]} -var topo ${dirname} -var cation ${cation}
      fi
    done
    cd ..
    
    
    echo "running coul simulation in $dirname"
    mkdir 2_coul; cd 2_coul
    cp ../../in.q .

    for m in $(seq 0 10)
    do
      if [ -e "./data.restart_${dirname}_${cation}_${q_lam[m]}" ]; then
      echo "Job with lambda=${q_lam[m]} already completed. Skipping."
      else
      mpirun -np 16 ${lammps} -in in.q -screen none -var lambda ${q_lam[m]} -var topo ${dirname} -var cation ${cation}
      fi
    done
    cd ..
    
    echo "running non-equi hr simulation in $dirname"
    mkdir 3_restraint; cd 3_restraint
    cp ../../in.hr_nonequi .
    if [ -e "./done.dat" ]; then
      echo "Job done for non-equilibration switching"
    else
      mpirun -np 16 ${lammps} -in in.hr_nonequi -screen none -var topo ${dirname} -var cation ${cation}
    fi
    done
    cd ..

done

