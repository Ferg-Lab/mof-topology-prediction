import os
from pathlib import Path
import random
import lammpsio
from utils import *
import re
from templates import *
#from mpi4py import MPI
#from lammps import lammps

def write_mof_params(f, mof_param, mof_styles, merged_styles):
    f.write("### Parameters for MOF ### \n\n")
    ### write MASSES
    for entry in mof_param['Masses']:
        f.write(f"mass {int(entry[0])} {entry[1]}\n")
    f.write("\n")

    for entry in mof_param['Pair Coeffs']:
        f.write(f"pair_coeff {int(entry[0])} {int(entry[0])} {' '.join([f'{value}' for value in entry[1:]])}\n")
    f.write("\n")
    
    ### write BONDS
    if 'Bond Coeffs' in mof_param.keys():
        mof_b_type = mof_styles.get('bond_style')
        if 'hybrid' in merged_styles.get('bond_style') and 'hybrid' in mof_b_type:
            for entry in mof_param['Bond Coeffs']:
                f.write(f"bond_coeff {int(entry[0])} {' '.join([f'{value}' for value in entry[1:]])}\n")
        elif 'hybrid' in merged_styles.get('bond_style') and 'hybrid' not in mof_b_type:
            for entry in mof_param['Bond Coeffs']:
                f.write(f"bond_coeff {int(entry[0])} {mof_b_type} {' '.join([f'{value}' for value in entry[1:]])}\n")
        else:
            for entry in mof_param['Bond Coeffs']:
                f.write(f"bond_coeff {int(entry[0])} {' '.join([f'{value}' for value in entry[1:]])}\n")
    f.write("\n")
    
    if 'Angle Coeffs' in mof_param.keys():
        mof_a_type = mof_styles.get('angle_style')
        if 'hybrid' in merged_styles.get('angle_style') and 'hybrid' in mof_a_type:
            for entry in mof_param['Angle Coeffs']:
                f.write(f"angle_coeff {int(entry[0])} {' '.join([f'{value}' for value in entry[1:]])}\n")
        elif 'hybrid' in merged_styles.get('angle_style') and 'hybrid' not in mof_a_type:
            for entry in mof_param['Angle Coeffs']:
                f.write(f"angle_coeff {int(entry[0])} {mof_a_type} {' '.join([f'{value}' for value in entry[1:]])}\n")
        else:
            for entry in mof_param['Angle Coeffs']:
                f.write(f"angle_coeff {int(entry[0])} {' '.join([f'{value}' for value in entry[1:]])}\n")
    f.write("\n")

    # Dihedrals
    if 'Dihedral Coeffs' in mof_param.keys():
        mof_d_type = mof_styles.get('dihedral_style')
        if 'hybrid' in merged_styles.get('dihedral_style') and 'hybrid' in mof_d_type:
            for entry in mof_param['Dihedral Coeffs']:
                f.write(f"dihedral_coeff {int(entry[0])} {' '.join([f'{value}' for value in entry[1:]])}\n")
        elif 'hybrid' in merged_styles.get('dihedral_style') and 'hybrid' not in mof_d_type:
            for entry in mof_param['Dihedral Coeffs']:
                f.write(f"dihedral_coeff {int(entry[0])} {mof_d_type} {' '.join([f'{value}' for value in entry[1:]])}\n")
        else:
            for entry in mof_param['Dihedral Coeffs']:
                f.write(f"dihedral_coeff {int(entry[0])} {' '.join([f'{value}' for value in entry[1:]])}\n")
    f.write("\n")

    # Impropers
    if 'Improper Coeffs' in mof_param.keys():
        mof_im_type = mof_styles.get('improper_style')
        if 'hybrid' in merged_styles.get('improper_style') and 'hybrid' in mof_im_type:
            for entry in mof_param['Improper Coeffs']:
                f.write(f"improper_coeff {int(entry[0])} {' '.join([f'{value}' for value in entry[1:]])}\n")
        elif 'hybrid' in merged_styles.get('improper_style') and 'hybrid' not in mof_im_type:
            for entry in mof_param['Improper Coeffs']:
                f.write(f"improper_coeff {int(entry[0])} {mof_im_type} {' '.join([f'{value}' for value in entry[1:]])}\n")
        else:
            for entry in mof_param['Improper Coeffs']:
                f.write(f"improper_coeff {int(entry[0])} {' '.join([f'{value}' for value in entry[1:]])}\n")
    f.write("\n")


def write_mol_params(f, mol_param, mof_param, mol_styles, merged_styles):
    mof_atm_num = len(mof_param['Masses'])
    if 'Bond Coeffs' in mof_param.keys():
        mof_b_num = len(mof_param['Bond Coeffs'])
    else:
        mof_b_num = 0
        
    if 'Angle Coeffs' in mof_param.keys():
        mof_a_num = len(mof_param['Angle Coeffs'])
    else:
        mof_a_num = 0
        
    if 'Dihedral Coeffs' in mof_param.keys():
        mof_d_num = len(mof_param['Dihedral Coeffs'])
    else:
        mof_d_num = 0
    if 'Improper Coeffs' in mof_param.keys():
        mof_im_num = len(mof_param['Improper Coeffs'])
    else:
        mof_im_num = 0

    f.write("### Parameters for Molecule ### \n\n")
    
    for entry in mol_param['Masses']:
        f.write(f"mass {int(entry[0])+mof_atm_num} {entry[1]}\n")
    f.write("\n")

    for entry in mol_param['Pair Coeffs']:
        f.write(f"pair_coeff {int(entry[0])+mof_atm_num} {int(entry[0])+mof_atm_num} {' '.join([f'{value}' for value in entry[1:]])}\n")
    f.write("\n")

    ### Write BONDS
    if 'Bond Coeffs' in mol_param.keys():
        mol_b_type = mol_styles.get('bond_style')
        if 'hybrid' in merged_styles.get('bond_style') and 'hybrid' in mol_b_type:
            for entry in mol_param['Bond Coeffs']:
                f.write(f"bond_coeff {int(entry[0])+mof_b_num} {' '.join([f'{value}' for value in entry[1:]])}\n")
        elif 'hybrid' in merged_styles.get('bond_style') and 'hybrid' not in mol_b_type:
            for entry in mol_param['Bond Coeffs']:
                f.write(f"bond_coeff {int(entry[0])+mof_b_num} {mol_b_type} {' '.join([f'{value}' for value in entry[1:]])}\n")
        else:
            for entry in mol_param['Bond Coeffs']:
                f.write(f"bond_coeff {int(entry[0])+mof_b_num} {' '.join([f'{value}' for value in entry[1:]])}\n")
    f.write("\n")

    ### Write ANGLES
    if 'Angle Coeffs' in mol_param.keys():
        mol_a_type = mol_styles.get('angle_style')
        if 'hybrid' in merged_styles.get('angle_style') and 'hybrid' in mol_a_type:
            for entry in mol_param['Angle Coeffs']:
                f.write(f"angle_coeff {int(entry[0])+mof_a_num} {' '.join([f'{value}' for value in entry[1:]])}\n")
        elif 'hybrid' in merged_styles.get('angle_style') and 'hybrid' not in mol_a_type:
            for entry in mol_param['Angle Coeffs']:
                f.write(f"angle_coeff {int(entry[0])+mof_a_num} {mol_a_type} {' '.join([f'{value}' for value in entry[1:]])}\n")
        else:
            for entry in mol_param['Angle Coeffs']:
                f.write(f"angle_coeff {int(entry[0])+mof_a_num} {' '.join([f'{value}' for value in entry[1:]])}\n")
    f.write("\n")

    ### Write DIHEDRALS
    if 'Dihedral Coeffs' in mol_param.keys():
        mol_d_type = mol_styles.get('dihedral_style')
        if 'hybrid' in merged_styles.get('dihedral_style') and 'hybrid' in mol_d_type:
            for entry in mol_param['Dihedral Coeffs']:
                f.write(f"dihedral_coeff {int(entry[0])+mof_d_num} {' '.join([f'{value}' for value in entry[1:]])}\n")
        elif 'hybrid' in merged_styles.get('dihedral_style') and 'hybrid' not in mol_d_type:
            for entry in mol_param['Dihedral Coeffs']:
                f.write(f"dihedral_coeff {int(entry[0])+mof_d_num} {mol_d_type} {' '.join([f'{value}' for value in entry[1:]])}\n")
        else:
            for entry in mol_param['Dihedral Coeffs']:
                f.write(f"dihedral_coeff {int(entry[0])+mof_d_num} {' '.join([f'{value}' for value in entry[1:]])}\n")
    f.write("\n")

    ### Write IMPROPERS
    if 'Improper Coeffs' in mol_param.keys():
        mol_im_type = mol_styles.get('improper_style')
        if 'hybrid' in merged_styles.get('improper_style') and 'hybrid' in mol_im_type:
            for entry in mol_param['Improper Coeffs']:
                f.write(f"improper_coeff {int(entry[0])+mof_im_num} {' '.join([f'{value}' for value in entry[1:]])}\n")
        elif 'hybrid' in merged_styles.get('improper_style') and 'hybrid' not in mol_im_type:
            for entry in mol_param['Improper Coeffs']:
                f.write(f"improper_coeff {int(entry[0])+mof_im_num} {mol_im_type} {' '.join([f'{value}' for value in entry[1:]])}\n")
        else:
            for entry in mol_param['Improper Coeffs']:
                f.write(f"improper_coeff {int(entry[0])+mof_im_num} {' '.join([f'{value}' for value in entry[1:]])}\n")
    f.write("\n")

def write_equi_params(mof, mol_name=None):
    mof_param = read_param(mof)
    mof_styles = read_styles(mof)
    
    if mol_name is not None:
        mol_param = read_param(mol_name)
        mol_styles = read_styles(mol_name)
        merged_styles = merge_styles(mof_styles, mol_styles)
    else:
        merged_styles = mof_styles
    if os.path.exists(f'emin_{mof}.param'):
        os.remove(f'emin_{mof}.param')
    
    with open(f'emin_{mof}.param', 'a') as f:
        f.write("### Parameters for Equilibration ### \n\n")
        write_mof_params(f, mof_param, mof_styles, merged_styles)
        if mol_name is not None:
            write_mol_params(f, mol_param=mol_param, mof_param=mof_param, mol_styles=mol_styles, merged_styles=merged_styles)

def write_equi_input(mof, nvt, npt, mol_name=None, n_cation=None, equi_time=100000, temp=300, pressure=0.0):
    ## read the input file
    mof_styles = read_styles(mof)
    if mol_name is not None:
        mol_styles = read_styles(mol_name)
        merged_styles = merge_styles(mol_styles, mof_styles)
        if n_cation is None:
            raise ValueError("n_cation must be provided if mol_name is not None")
    else:
        merged_styles = mof_styles
    
    mof_param = read_param(mof)
    mof_atm_num = len(mof_param['Masses'])
    if 'Bond Coeffs' in mof_param.keys():
        mof_b_num = len(mof_param['Bond Coeffs'])
    else:
        mof_b_num = 0
        
    if 'Angle Coeffs' in mof_param.keys():
        mof_a_num = len(mof_param['Angle Coeffs'])
    else:
        mof_a_num = 0
        
    if 'Dihedral Coeffs' in mof_param.keys():
        mof_d_num = len(mof_param['Dihedral Coeffs'])
    else:
        mof_d_num = 0
    if 'Improper Coeffs' in mof_param.keys():
        mof_im_num = len(mof_param['Improper Coeffs'])
    else:
        mof_im_num = 0

    mol_param = read_param(mol_name)
    mol_atm_num = len(mol_param['Masses'])
    if 'Bond Coeffs' in mol_param.keys():
        mol_b_num = len(mol_param['Bond Coeffs'])
    else:
        mol_b_num = 0
        
    if 'Angle Coeffs' in mol_param.keys():
        mol_a_num = len(mol_param['Angle Coeffs'])
    else:
        mol_a_num = 0
        
    if 'Dihedral Coeffs' in mol_param.keys():
        mol_d_num = len(mol_param['Dihedral Coeffs'])
    else:
        mol_d_num = 0
    if 'Improper Coeffs' in mol_param.keys():
        mol_im_num = len(mol_param['Improper Coeffs'])
    else:
        mol_im_num = 0

    total_atm_num = mof_atm_num + mol_atm_num
    total_b_num = mof_b_num + mol_b_num
    total_a_num = mof_a_num + mol_a_num
    total_d_num = mof_d_num + mol_d_num
    total_im_num = mof_im_num + mol_im_num

    modify_mof_data(mof, total_atm_num, total_b_num, total_a_num, total_d_num, total_im_num, 'new')
    
    random_seed = random.randint(1, 10e5)

    ## initilization section
    with open(f'in.emin_{mof}', 'w') as f:
        f.write(em_template.format(seed = random_seed, mof = mof, 
                                   temp = temp, pressure = pressure, 
                                   pstyle = merged_styles.get('pair_style'), 
                                   bstyle = merged_styles.get('bond_style'), 
                                   astyle = merged_styles.get('angle_style'), 
                                   dstyle = merged_styles.get('dihedral_style'), 
                                   istyle = merged_styles.get('improper_style'), 
                                   kstyle = merged_styles.get('kspace_style')
                                  ))
        if mol_name is not None:
            f.write(deposit_template.format(n_cation = n_cation,
                                            cation = mol_name,
                                            Toff = mof_atm_num,
                                            Boff = mof_b_num,
                                            Aoff = mof_a_num,
                                            Doff = mof_d_num,
                                            Ioff = mof_im_num,
                                            step = str(2*int(n_cation))
                                           ))
        
        if nvt and not npt:
            f.write(nvt_template.format(equi_time=equi_time))
        elif npt and not nvt:
            f.write(npt_template.format(equi_time=equi_time))
        elif nvt and npt:
            f.write(nvt_npt_template.format(equi_time=equi_time))

        f.write(output_template.format(mof = mof, cation = mol_name))
            
            
    print('Input for Equilibration has been generated.')
    