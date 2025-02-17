import argparse
import sys
from build_mof import *
from write_mol import *
from equilibration import *
from free_energy import *
from utils import *
import subprocess
import glob

def gen(args, out_dir):
    node = args.node
    linker = args.linker
    topos = args.topos
    forcefield = args.forcefield
    replication = args.replication
    mol = args.mol
    
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"Generating MOFs with {node} and {linker} using {topos}")
    build_mof(node, linker, topos, forcefield, replication, out_dir)
    if mol is not None:
        print(f"The small molecule is {mol}.")
        create_mol_file(mol, forcefield, out_dir)
    else:
        print(f"No small molecules.")

def run_equi(args, out_dir):
    mof = args.mof
    mol = args.mol
    n_mol = args.n_mol
    nvt = args.nvt
    npt = args.npt
    equi_time = args.equi_time
    temp = args.temp
    pressure = args.pressure
    
    print("Running Equilibration")
    if mol is not None:
        emin_data = f'{out_dir}/data.emin_{mof}_{mol}'
    else:
        emin_data = f'{out_dir}/data.emin_{mof}'
    write_equi_params(mof, out_dir, mol)
    write_equi_input(mof, out_dir, nvt, npt, mol, n_mol, equi_time, temp, pressure)

    lammps_path = os.getenv('LAMMPS_PATH')
    prev_dir = os.getcwd()  # Save the previous directory
    os.chdir(out_dir) 

    try:
        if os.path.exists(emin_data.split('/')[1]):
            print("Minimization already completed.")
        else:
            commands = f'{lammps_path} -in in.emin_{mof}'
            subprocess.run(commands, shell=True, check=True)
            lmp2cif(mof, mol)
    finally:
        os.chdir(prev_dir) 
    

def run_fe(args,out_dir):
    mof = args.mof
    mol = args.mol
    temp = args.temp
    pressure = args.pressure
    center = args.center
    nproc = args.nproc
    
    print("Running Free Energy Calculation")
    write_non_bonded_lj(mof, out_dir, mol)
    write_non_bonded(mof, out_dir, mol)
    write_bonded(mof,out_dir)
    write_in_fe(mof, out_dir, mol, center, temp, pressure)

    lammps_path = os.getenv('LAMMPS_PATH')
    prev_dir = os.getcwd()  # Save the previous directory
    os.chdir(out_dir) 

    ## the lambdas
    b_lam = [0.00000001]#, 0.0001, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,]  
    l_lam = [0]#, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
    q_lam = [0.00000001] #, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,]
    
    try:
        sub_dirs = {
            "bonded": ("in.BONDED", b_lam),
            "lj": ("in.LJ", l_lam),
            "q": ("in.Q", q_lam),
            "hr": ("in.HR", None)
        }

        for sub_dir, (script, lambdas) in sub_dirs.items():
            sub_dir_path = os.path.join(os.path.join(prev_dir, out_dir), sub_dir)
            if os.path.exists(sub_dir_path):
                os.chdir(sub_dir_path)  # Change to the subdirectory
                print(f"Running {script} in {sub_dir_path}")

                if lambdas is not None:
                    for lam in lambdas:
                        if os.path.exists(glob.glob(f'*{lam}*lmp')):
                            print(f'Completed on {lam}!')
                        else:
                            command = f'mpirun -np {nproc} {lammps_path} -in {script} -var lambda {lam} -screen none'
                            print(f"Executing: {command}")
                            subprocess.run(command, shell=True, check=True)
                else:
                    # For "hr" (no lambda needed)
                    if os.path.exists('done.dat'):
                        print(f'Completed!')
                    else:
                        command = f'mpirun -np {nproc} {lammps_path} -in {script} -screen none'
                        print(f"Executing: {command}")
                        subprocess.run(command, shell=True, check=True)

            else:
                print(f"Skipping {sub_dir}: directory not found.")
    finally:
        os.chdir(prev_dir)  # Restore the previous directory
    

def main():    
    parser = argparse.ArgumentParser(description="MOF PolyRank CLI tool")
    subparsers = parser.add_subparsers(dest='command')
    
    # Subparser for the `gen` command
    parser_gen = subparsers.add_parser('gen', help='Generate MOF and Molecule')
    parser_gen.add_argument('--node', help='The name of the metal node in xyz format')
    parser_gen.add_argument('--linker', help='The name of the organic linker in xyz format')
    parser_gen.add_argument('--topos', help='A list of topologies to construct MOFs')
    parser_gen.add_argument('--mol', help='The name of the small molecule (default is None)', default=None)
    parser_gen.add_argument('--forcefield', help='Force field for MOFs, e.g., UFF, UFF4MOF (default), used in lammps-interface', default='UFF4MOF')
    parser_gen.add_argument('--replication', help='The number to replicate the MOF cell. Default is None, and the MOF will be automatically replicated to ensure the simualtion box size is larger than 25 Å in each direction.', default=None)

    # Subparser for the `run_equi` command
    parser_run_equi = subparsers.add_parser('run_equi', help='Run Equilibration')
    parser_run_equi.add_argument('--mof', help='The name of the mof without extension')
    parser_run_equi.add_argument('--mol', help='The name of the small molecule', default=None)
    parser_run_equi.add_argument('--n_mol', help='The number of the small molecule', default=0)
    parser_run_equi.add_argument('--nvt', help='Run NVT Equilibration', action='store_true')
    parser_run_equi.add_argument('--npt', help='Run NPT Equilibration', action='store_true')
    parser_run_equi.add_argument('--equi_time', help='The equilibration time (default 100 ps)', default=100000)
    parser_run_equi.add_argument('--temp', help='The temperature for simulation (default 300.0)', default=300)
    parser_run_equi.add_argument('--pressure', help='The pressure for simulation (default 0.0)', default=0.0)

    # Subparser for the `run_fe` command
    parser_run_fe = subparsers.add_parser('run_fe', help='Run Free Energy')
    parser_run_fe.add_argument('--mof', help='The name of the mof without extension name')
    parser_run_fe.add_argument('--mol', help='The name of the small molecule', default=None)
    parser_run_fe.add_argument('--center', help='The center of the small molecule', default=None)
    parser_run_fe.add_argument('--temp', help='The temperature for simulation (default 300.0)', default=300)
    parser_run_fe.add_argument('--pressure', help='The pressure for simulation (default 0.0)', default=0.0)
    parser_run_fe.add_argument('--nproc', help='Run the simulation using mpi or no', default=1)

    args = parser.parse_args()

    if args.command == 'gen':
        out_dir = f"{args.topos}-{args.node}-{args.linker}"
        gen(args, out_dir)
    elif args.command == 'run_equi':
        out_dir = f"{args.mof}"
        run_equi(args, out_dir)
    elif args.command == 'run_fe':
        out_dir = f"{args.mof}"
        run_fe(args, out_dir)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
