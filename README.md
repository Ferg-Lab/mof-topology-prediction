
# MOF Polymorph Generation and Free Energy Calculation

## Prerequisites

Create a conda env and install the dependencies

`conda create -n mof-topology`

`conda install python=3.8 -y`

#### Install kernels for jupyter notebook

`conda install -c conda-forge ipykernel -y`

`python -m ipykernel install --user --name=mof-topology`

#### Pormake is used to build a MOF polymorph

`pip install git+https://github.com/Sangwon91/PORMAKE.git`

This version allows the partial charges assignment on the building nodes.

#### LAMMPS-Interface is used to assign force field parameters

`pip install lammps-interface`

`pip install lammpsio`

Originally cif2lammps and ToBaCCo are used for MOF generation (see folder simulation_utilities); for the automation purposes, the combination of pormake and lammps-interface seems to work better.

## Compile LAMMPS

LAMMPS version 22 Jun 2022 ([LAMMPS](https://github.com/lammps/lammps/releases/tag/stable_23Jun2022))

```
wget https://github.com/lammps/lammps/archive/refs/tags/stable_23Jun2022.tar.gz
tar -xvf stable_23Jun2022.tar.gz
```

The following additional packages are required:

```
cd lammps-stable_23Jun2022/src
make yes-molecule
make yes-extra-fix
make yes-extra-compute
make yes-fep
make yes-body
make yes-kspace
make yes-manybody
make yes-rigid
make yes-extra-molecule
```

Replace the `fix_ti_spring.cpp` in the `LAMMPS/src` folder with the provided one, which is slightly modified for the final NES step calculation, and compile. 

`cp ../../fix_ti_spring.cpp .`

`make mpi -j8` or `make seriel -j8`

`export LAMMPS_PATH=\`pwd\`/lmp_mpi`

## Navigating the Repo

There are three tutorial notebook in the example folder showing the procedures for building MOF polymorph and compute free energies. Example1, example2, and example3  include a ZJU-28 series MOF, a ZIF series MOF, and a Fe4S4-BDT series MOF. One can adapt the notebook or use the command lines to run the pipeline shown in the next section. 

The simulation_utilities folder contains all the scripts of running the free energy calculation on the Fe4S4-BDT MOFs with different cations using the ToBaCCo codes and cif2lammps. There are also guidelines in this folder for running the simulations using this set of tools.

The databank_iron-sulfur-MOFs includes all the relaxed structures of the iron-sulfur MOFs studied in the paper, along with computed free energies, potential energies, and geometric properties. The PDF files are stored in the H5 files that can be accessed using h5py. 

## MOF Generation and Equilibration

The `gen` command is used to generate a MOF polymorph, for which a linker and a metal node in xyz format are required. You can generate bare frameworks or specify a small molecule to add into an ionic MOF. If the name of the small molecule is available in the `small_molecule` folder (DMA, TMA, TEA, TPA, TPP, MPA, MNP), then the molecule files for small molecules will be created along with the generated MOF. Otherwise, you would need to provide a small molecule in cif format with atomic charges.

The `run_equi` command generates inputs for the small molecule deposition and equilibration. After equilibration, an equilibrated structure in cif format is provided along with the potential energy at 300 K and 0 K for preliminary screening.

Example:

`python main.py gen --node 4c_In --linker 3c_BTB --topos pto --mol DMA`

`python main.py run_equi --mof pto-4c_In-3c_BTB --mol DMA --n_mol 6 --nvt --npt --equi_time 50000 --temp 300 --pressure 0.0`

This is a checkpoint where you can stop and use the equilibrated structures and potential energies to rank the thermodynamic stability of the MOFs.

## Free Energy Calculation Workflow

The equilibrated structures can be used for subsequent free energy calculations. The `run_fe` function is used to generate the input files for calculating the free energies of the MOFs. If a small molecule is used, it is required to specify the atom type of the molecules from the LAMMPS `emin` data file.

Example:

`python main.py run_fe --mof pto-4c_In-3c_BTB --mol DMA --center 7 --nproc 8`

This will create a list of input files for running free energy simulations and create a folder for storing the output data.

## Free Energy Integration

The codes for integrating the free energy after the calculation is available in example 1. 
