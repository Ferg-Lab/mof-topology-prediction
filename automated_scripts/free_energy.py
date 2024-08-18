import os
from pathlib import Path
import random
import lammpsio
from utils import *
import re
from templates import *

def write_non_bonded_lj(mof, mol_name):
    pairij = read_param(f'emin_{mof}_{mol_name}')['PairIJ Coeffs']
    pairij_pairs = [entry[:-1] for entry in pairij]
    with open(f'emin_{mof}.param', 'r') as f:
        lines = f.readlines()
    filtered_lines = [line for line in lines if not line.startswith('pair_coeff')]
    
    cross_term_lines=[f"pair_coeff {' '.join([f'{value}' for value in entry])} 1 12.5\n" for entry in pairij_pairs]
    filtered_lines.extend(cross_term_lines)

    with open(f'non_bonded_{mof}_lj.param', 'w') as f:
        f.writelines(filtered_lines)

def write_non_bonded(mof, mol_name):
    pairij = read_param(f'emin_{mof}_{mol_name}')['PairIJ Coeffs']
    pairij_pairs = [entry[:-1] for entry in pairij]
    with open(f'emin_{mof}.param', 'r') as f:
        lines = f.readlines()
    filtered_lines = [line for line in lines if not line.startswith('pair_coeff')]
    
    cross_term_lines=[f"pair_coeff {' '.join([f'{value}' for value in entry])}\n" for entry in pairij_pairs]
    filtered_lines.extend(cross_term_lines)

    with open(f'non_bonded_{mof}.param', 'w') as f:
        f.writelines(filtered_lines)

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False    

def process_coeff_line(line, var_prefix):
    parts = line.split()
    if is_number(parts[2]):
        energy_index = 2
    else:
        energy_index = 3

    energy_term = parts[energy_index]
    var_name = f'{var_prefix}{parts[1]}'
    variable_assignment = f'variable {var_name} equal ${{lambda}}*{energy_term}'
    parts[energy_index] = f'${{{var_name}}}'
    modified_line = ' '.join(parts)
    
    return variable_assignment, modified_line

def write_bonded(mof):
    input_filename = f'emin_{mof}.param'
    output_filename = f'bonded_{mof}.param'

    with open(input_filename, 'r') as file:
        lines = file.readlines()

    output_lines = []
    in_mof_section = False
    var_prefixes = {
        'bond_coeff': 'b',
        'angle_coeff': 'a',
        'dihedral_coeff': 'd',
        'improper_coeff': 'i'
    }

    for line in lines:
        line = line.strip()
        if "MOF" in line:
            in_mof_section = True
        elif "Molecule" in line:
            in_mof_section = False
            
        if line.startswith('pair_coeff'):
            continue
            
        if in_mof_section and line.startswith(('bond_coeff', 'angle_coeff', 'dihedral_coeff', 'improper_coeff')):
            coeff_type = line.split()[0]
            var_prefix = var_prefixes[coeff_type]
            variable_assignment, modified_line = process_coeff_line(line, var_prefix)
            output_lines.append(variable_assignment)
            output_lines.append(modified_line)
        else:
            output_lines.append(line)
    output_lines.append('pair_coeff * *')
    # Write the output file
    with open(output_filename, 'w') as file:
        for line in output_lines:
            file.write(line + '\n')

def write_in_fe(mof, mol_name=None, center_atom=None, temp=300, pressure=0.0):
    
    ## read the input file
    mof_styles = read_styles(mof)
    if mol_name is not None:
        mol_styles = read_styles(mol_name)
        merged_styles = merge_styles(mol_styles, mof_styles)
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

    job=mol_name+'_fe'
    modify_mof_data(f'emin_{mof}_{mol_name}', total_atm_num, total_b_num, total_a_num, total_d_num, total_im_num, 'fe')
    
    random_seed = random.randint(1, 10e5)
    os.makedirs(f'{mof}/', exist_ok=True)
    ## the bonded part
    os.makedirs(f'{mof}/bonded/', exist_ok=True)
    with open(f'in.BONDED_{mof}_{mol_name}', 'w') as f:
        f.write(bonded_part1_template.format(seed = random_seed, mof = mof,
                                             temp = temp, pressure = pressure, 
                                             pstyle = 'zero 12.5 nocoeff', 
                                             bstyle = merged_styles.get('bond_style'), 
                                             astyle = merged_styles.get('angle_style'), 
                                             dstyle = merged_styles.get('dihedral_style'), 
                                             istyle = merged_styles.get('improper_style'), 
                                             fram_atoms = mof_atm_num+1,
                                             job=job
                                  ))
        if mol_name is not None:
            f.write(f'group			cation type > {mof_atm_num}\n')
            f.write(f'group			center type {center_atom}\n')
            f.write(f'fix			SPRING2 center spring/self 1\n\n')
        
        f.write(bonded_part2_template.format(mof=mof))
    
    os.makedirs(f'{mof}/lj/', exist_ok=True)
    with open(f'in.LJ_{mof}_{mol_name}', 'w') as f:
        f.write(lj_part1_template.format(seed = random_seed, mof = mof, 
                                         temp = temp, pressure = pressure, 
                                         pstyle = 'lj/cut/soft 1 0.5 12.5', 
                                         bstyle = merged_styles.get('bond_style'), 
                                         astyle = merged_styles.get('angle_style'), 
                                         dstyle = merged_styles.get('dihedral_style'), 
                                         istyle = merged_styles.get('improper_style'), 
                                         fram_atoms = mof_atm_num+1,
                                         job=job
                                  ))
        if mol_name is not None:
            f.write(f'group			cation type > {mof_atm_num}\n')
            f.write(f'group			center type {center_atom}\n')
            f.write(f'fix			SPRING2 center spring/self 1\n\n')
        
        f.write(lj_part2_template.format(mof=mof))

    if 'coul' in merged_styles.get('pair_style'):
        os.makedirs(f'{mof}/q/', exist_ok=True)
        with open(f'in.Q_{mof}_{mol_name}', 'w') as f:
            f.write(q_part1_template.format(seed = random_seed, mof = mof, 
                                            temp = temp, pressure = pressure, 
                                            pstyle = merged_styles.get('pair_style'), 
                                            bstyle = merged_styles.get('bond_style'), 
                                            astyle = merged_styles.get('angle_style'), 
                                            dstyle = merged_styles.get('dihedral_style'), 
                                            istyle = merged_styles.get('improper_style'), 
                                            kstyle = merged_styles.get('kspace_style'),
                                            fram_atoms = mof_atm_num+1,
                                            job=job
                                            ))
            if mol_name is not None:
                f.write(f'group			cation type > {mof_atm_num}\n')
                f.write(f'group			center type {center_atom}\n')
                f.write(f'fix			SPRING2 center spring/self 1\n\n')
        
            f.write(q_part2_template.format(mof=mof))

    
    os.makedirs(f'{mof}/hr/', exist_ok=True)
    if mol_name is not None:
        with open(f'in.HR_{mof}_{mol_name}', 'w') as f:
            f.write(hr_with_mol_template.format(seed = random_seed, mof = mof, 
                                                temp = temp, pressure = pressure, 
                                                pstyle = merged_styles.get('pair_style'), 
                                                bstyle = merged_styles.get('bond_style'), 
                                                astyle = merged_styles.get('angle_style'), 
                                                dstyle = merged_styles.get('dihedral_style'), 
                                                istyle = merged_styles.get('improper_style'), 
                                                kstyle = merged_styles.get('kspace_style'),
                                                fram_atoms = mof_atm_num+1,
                                                cation_atoms=mof_atm_num,
                                                center_atom=center_atom,
                                                job=job
                                                ))
    else:
        with open(f'in.HR_{mof}_{mol_name}', 'w') as f:
            f.write(hr_no_mol_template.format(seed = random_seed, mof = mof, 
                                              temp = temp, pressure = pressure, 
                                              pstyle = merged_styles.get('pair_style'), 
                                              bstyle = merged_styles.get('bond_style'), 
                                              astyle = merged_styles.get('angle_style'), 
                                              dstyle = merged_styles.get('dihedral_style'), 
                                              istyle = merged_styles.get('improper_style'), 
                                              kstyle = merged_styles.get('kspace_style'),
                                              fram_atoms = mof_atm_num+1,
                                              job=job
                                              ))

    
    print('Inputs for Free energies has been generated.')

