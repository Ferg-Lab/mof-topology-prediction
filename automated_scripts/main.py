import argparse
import sys
from build_mof import *
from write_mol import *
from equilibration import *
from free_energy import *
from utils import *
import subprocess

def gen(args):
    node = args.node
    linker = args.linker
    topos = args.topos
    forcefield = args.forcefield
    replication = args.replication
    mol = args.mol
    
    print(f"Generating MOFs with {node} and {linker} using {topos}")
    build_mof(node, linker, topos, forcefield, replication)
    if mol is not None:
        print(f"The small molecule is {mol}.")
        create_mol_file(mol, forcefield)
    else:
        print(f"No small molecules.")

def run_equi(args):
    mof = args.mof
    mol = args.mol
    n_mol = args.n_mol
    nvt = args.nvt
    npt = args.npt
    equi_time = args.equi_time
    temp = args.temp
    pressure = args.pressure

    print("Running Equilibration")
    emin_data = f'data.emin_{mof}_{mol}'
    write_equi_params(mof, mol)
    write_equi_input(mof, nvt, npt, mol, n_mol, equi_time, temp, pressure)
    lammps_path = os.getenv('LAMMPS_PATH')
    if os.path.exists(emin_data):
        print("Minimization Done")
    else:
        commands = f'{lammps_path} -in in.emin_{mof}'
        subprocess.run(commands, shell=True)
    lmp2cif(mof, mol)

def run_fe(args):
    mof = args.mof
    mol = args.mol
    temp = args.temp
    pressure = args.pressure
    center = args.center
    
    print("Running Free Energy Calculation")
    write_non_bonded_lj(mof, mol)
    write_non_bonded(mof, mol)
    write_bonded(mof)
    write_in_fe(mof, mol, center, temp, pressure)

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

    args = parser.parse_args()

    if args.command == 'gen':
        gen(args)
    elif args.command == 'run_equi':
        run_equi(args)
    elif args.command == 'run_fe':
        run_fe(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
