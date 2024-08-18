import os
import lammpsio
from pathlib import Path
import numpy as np
from utils import *

def write_header(file, mol_data, styles):
    with open(file, 'a') as f:
        num_atom = mol_data.id.size
        f.write(f'# molecule file for lammps\n\n')
        f.write(f'{num_atom} atoms\n')
        
        if mol_data.bonds is not None:
            num_bonds = mol_data.bonds.id.size
            f.write(f'{num_bonds} bonds\n')
        if mol_data.angles is not None:
            num_angles = mol_data.angles.id.size
            f.write(f'{num_angles} angles\n')
        if mol_data.dihedrals is not None:
            num_dihedrals = mol_data.dihedrals.id.size
            f.write(f'{num_dihedrals} dihedrals\n')
        if mol_data.impropers is not None:
            num_impropers = mol_data.impropers.id.size
            f.write(f'{num_impropers} impropers\n')

def write_section(file, title, data, styles_key=None, styles=None):
    with open(file, 'a') as f:
        f.write(f'\n{title}\n\n')
        for row in data:
            f.write('\t'.join(map(str, row)))
            f.write('\n')

def write_coords(file, mol_data):
    coords = np.column_stack((mol_data.id.astype(int), mol_data.position.astype(str)))
    write_section(file, 'Coords', coords)

def write_atom_types(file, mol_data, styles):
    a_types = np.column_stack((mol_data.id.astype(int), mol_data.typeid.astype(int)))
    write_section(file, f'Types # {styles.get("atom_style")}', a_types)

def write_charges(file, mol_data):
    charges = np.column_stack((mol_data.id.astype(int), mol_data.charge.astype(str)))
    write_section(file, 'Charges', charges)

def write_bonds(file, mol_data, styles):
    bonds = np.column_stack((mol_data.bonds.id.astype(int), mol_data.bonds.typeid.astype(int), mol_data.bonds.members))
    write_section(file, f'Bonds # {styles.get("bond_style")}', bonds)

def write_angles(file, mol_data, styles):
    angles = np.column_stack((mol_data.angles.id.astype(int), mol_data.angles.typeid.astype(int), mol_data.angles.members))
    write_section(file, f'Angles # {styles.get("angle_style")}', angles)

def write_dihedrals(file, mol_data, styles):
    dihedrals = np.column_stack((mol_data.dihedrals.id.astype(int), mol_data.dihedrals.typeid.astype(int), mol_data.dihedrals.members))
    write_section(file, f'Dihedrals # {styles.get("dihedral_style")}', dihedrals)

def write_impropers(file, mol_data, styles):
    impropers = np.column_stack((mol_data.impropers.id.astype(int), mol_data.impropers.typeid.astype(int), mol_data.impropers.members))
    write_section(file, f'Impropers # {styles.get("improper_style")}', impropers)

def cif2lmp(mol_name, forcefield):
    mol_path = './small_molecule'
    mol_cif = os.path.join(mol_path, f'{mol_name}.cif')
    
    if os.path.exists(mol_cif) and not os.path.exists(f'{mol_name}.cif'):
        logging.info('Using the cif of small molecule in the database')
    elif not os.path.exists(mol_cif) and os.path.exists(f'{mol_name}.cif'):
        logging.info('Using the small molecule in current directory')
        mol_cif = f'{mol_name}.cif'
    else:
        logging.error('No small molecule found. A cif file for small molecule in cif format is required.')
        return
    
    command = f'lammps-interface {mol_cif} -ff {forcefield} --replication 1x1x1 > tmp; rm tmp'
    os.system(command)

def create_mol_file(mol_name, forcefield):
    cif2lmp(mol_name, forcefield)
    mol_data = read_data(mol_name)
    
    if mol_data is None:
        logging.error(f"Failed to create molecule file for {mol_name}.")
        return
    
    styles = read_styles(mol_name)
    mol_file = f'{mol_name}.txt'
    
    if os.path.exists(mol_file):
        os.remove(mol_file)
    
    write_header(mol_file, mol_data, styles)
    write_coords(mol_file, mol_data)
    write_atom_types(mol_file, mol_data, styles)
    
    if mol_data.charge is not None:
        write_charges(mol_file, mol_data)
    if mol_data.bonds is not None:
        write_bonds(mol_file, mol_data, styles)
    if mol_data.angles is not None:
        write_angles(mol_file, mol_data, styles)
    if mol_data.dihedrals is not None:
        write_dihedrals(mol_file, mol_data, styles)
    if mol_data.impropers is not None:
        write_impropers(mol_file, mol_data, styles)
    
    print('Small molecule file for LAMMPS deposition has been created!')


