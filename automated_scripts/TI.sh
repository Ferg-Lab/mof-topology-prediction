#!/bin/bash

lammps="/Path/to/LAMMPS/src/lmp_machine"

b_lam=(0.00000001) #0.0001 0.001 0.01 0.025 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)
l_lam=(0) # 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)
q_lam=(0.00000001) # 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

${lammps} -in in.Q_pto-4c_In-3c_BTB_DMA -var lambda ${q_lam[0]}
${lammps} -in in.HR_pto-4c_In-3c_BTB_DMA
${lammps} -in in.BONDED_pto-4c_In-3c_BTB_DMA -var lambda ${b_lam[0]}
${lammps} -in in.LJ_pto-4c_In-3c_BTB_DMA -var lambda ${l_lam[0]}
