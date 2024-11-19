
# Automated MOF Generation and Free Energy Calculation

## Prerequisites

### Pormake is used to build a MOF polymorph

pip install git+https://github.com/aniruddha-seal/PORMAKE.git 

This version allows the partial charges assignment on the building nodes.

### LAMMPS-Interface is used to assign force field parameters

pip install lammps-interface

## Compiling LAMMPS

LAMMPS version 22 Jun 2022 ([LAMMPS](https://github.com/lammps/lammps/releases/tag/stable_23Jun2022))

The following additional packages are required:

- `make yes-molecule`
- `make yes-extra-fix`
- `make yes-extra-compute`
- `make yes-fep`
- `make yes-body`
- `make yes-kspace`
- `make yes-manybody`
- `make yes-rigid`
- `make yes-extra-molecule`

Replace the `fix_ti_spring.cpp` in the `LAMMPS/src` folder with the provided one, which is slightly modified for the final NES step calculation. Finally, compile LAMMPS in the `src` directory.

Export the path to LAMMPS for the code to run equilibration automatically for you. Otherwise, you can run simulations manually using the input files.

export LAMMPS_PATH='/Path/to/LAMMPS/src/lmp_machine'

## MOF Generation and Equilibration

The `gen` command is used to generate a MOF polymorph, for which a linker and a metal node in xyz format are required. You can generate bare frameworks or specify a small molecule to add into an ionic MOF. If the name of the small molecule is available in the `small_molecule` folder (DMA, TMA, TEA, TPA, TPP, MPA, MNP), then the molecule files for small molecules will be created along with the generated MOF. Otherwise, you would need to provide a small molecule in cif format with atomic charges.

The `run_equi` command generates inputs for the small molecule deposition and equilibration. After equilibration, an equilibrated structure in cif format is provided along with the potential energy at 300 K and 0 K for preliminary screening.

Example:

`python main.py gen --node 4c_In --linker 3c_BTB --topos pto --mol DMA`
`python main.py run_equi --mof pto-4c_In-3c_BTB --mol DMA --n_mol 6 --nvt --npt --equi_time 50000 --temp 300 --pressure 0.0`

This is a checkpoint where you can stop and use the equilibrated structures and potential energies to rank the thermodynamic stability of the MOFs.

## Free Energy Calculation Workflow

The equilibrated structures can be used for subsequent free energy calculations. The `run_fe` function is used to generate the input files for calculating the free energies of the MOFs. If a small molecule is used, it is required to specify the atom type of the center atom of the molecules from the LAMMPS `emin` data file.

Example:

`python main.py run_fe --mof pto-4c_In-3c_BTB --mol DMA --center 7`

This will create a list of input files for running free energy simulations and create a folder for storing the output data. Run the free energy calculations manually using the input files.

## Free Energy Integration

Once you obtain all the output files, use the initial codes in `simulation_utilities/004_post-processing-tools` folder to integrate the free energies. 
