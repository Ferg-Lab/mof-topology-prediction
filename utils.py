import os
from pathlib import Path
import lammpsio
import logging
import re
from ase.io import lammpsdata, cif, read

# read the styles from in files
def read_styles(in_mol):
    styles = {}
    try:
        with open(in_mol, 'r') as file:
            for line in file:
                if 'style' in line:
                    key, value = line.split(maxsplit=1)
                    styles[key] = value.strip()
    except FileNotFoundError:
        logging.error(f"{in_mol} not found.")
    return styles

## merge the styles from mof and mol
def merge_styles(styles1, styles2):
    merged_styles = {}
    hybrid_styles = {}
    
    for key in set(styles1.keys()).union(styles2.keys()):
        value1 = styles1.get(key)
        value2 = styles2.get(key)
        
        if value1 and value2:
            if value1 != value2:
                hybrid_entries = set(value1.split()) | set(value2.split())
                hybrid_entries.discard('hybrid')
                hybrid_styles[key] = f"hybrid {' '.join(hybrid_entries)}"
            else:
                merged_styles[key] = value1
        elif value1:
            merged_styles[key] = value1
        elif value2:
            merged_styles[key] = value2

    merged_styles.update(hybrid_styles)
    
    return merged_styles

## read data using lammpsio
def read_data(data_mol):
    if os.path.exists(data_mol):
        data_file = lammpsio.DataFile(data_mol, atom_style='full')
        mol_data = data_file.read()
        return mol_data
    else:
        raise ValueError("No data file found")
        return None


## read in the lammps parameters data
def read_param(data_mol):
    
    sections = {}
    current_section = None

    with open(data_mol, 'r') as file:
        lines = file.readlines()

        for i, line in enumerate(lines):
            if "Atoms" in line.strip():
                break
            line = line.split('#')[0].strip()
            
            if not line:
                continue

            if any(keyword in line for keyword in ["Masses", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs", "Pair Coeffs", "PairIJ Coeffs"]):
                current_section = line
                sections[current_section] = []
            elif current_section is not None:
                sections[current_section].append(line.split())

    return sections


def modify_mof_data(mof, total_atm_num, total_b_num, total_a_num, total_d_num, total_im_num, job):
    with open(mof, 'r') as f:
        lines = f.readlines()
    patterns = {
    'atom types': re.compile(r'^\s*\d+\s+atom types\s*$'),
    'bond types': re.compile(r'^\s*\d+\s+bond types\s*$'),
    'angle types': re.compile(r'^\s*\d+\s+angle types\s*$'),
    'dihedral types': re.compile(r'^\s*\d+\s+dihedral types\s*$'),
    'improper types': re.compile(r'^\s*\d+\s+improper types\s*$'),
    }

    new_lines = []
    skip_lines = False
    for line in lines:
        stripped_line = line.strip()
        if patterns['atom types'].match(stripped_line):
            new_lines.append(f"           {total_atm_num} atom types\n")
        elif patterns['bond types'].match(stripped_line):
            new_lines.append(f"           {total_b_num} bond types\n")
        elif patterns['angle types'].match(stripped_line):
            new_lines.append(f"           {total_a_num} angle types\n")
        elif patterns['dihedral types'].match(stripped_line):
            new_lines.append(f"           {total_d_num} dihedral types\n")
        elif patterns['improper types'].match(stripped_line):
            new_lines.append(f"           {total_im_num} improper types\n")
        elif "Masses" in stripped_line:
            skip_lines = True
        elif "Atoms" in stripped_line:
            skip_lines = False
            new_lines.append(line)
        elif not skip_lines:
            new_lines.append(line)

    with open(f'{mof}_{job}', 'w') as f:
        f.writelines(new_lines)

def lmp2cif(mof, mol_name):
    atomic_masses={
    1: 1.0, 2: 4.0, 3: 6.9, 4: 9.0, 5: 10.8, 6: 12.0, 7: 14.0, 8: 16.0, 9: 19.0, 10: 20.2,
    11: 23.0, 12: 24.3, 13: 27.0, 14: 28.1, 15: 31.0, 16: 32.1, 17: 35.5, 18: 39.9, 19: 39.1, 20: 40.1,
    21: 45.0, 22: 47.9, 23: 50.9, 24: 52.0, 25: 54.9, 26: 55.8, 27: 58.9, 28: 58.7, 29: 63.5, 30: 65.4,
    31: 69.7, 32: 72.6, 33: 74.9, 34: 79.0, 35: 79.9, 36: 83.8, 37: 85.5, 38: 87.6, 39: 88.9, 40: 91.2,
    41: 92.9, 42: 96.0, 43: 98.0, 44: 101.1, 45: 102.9, 46: 106.4, 47: 107.9, 48: 112.4, 49: 114.8, 50: 118.7,
    51: 121.8, 52: 127.6, 53: 126.9, 54: 131.3, 55: 132.9, 56: 137.3, 57: 138.9, 58: 140.1, 59: 140.9, 60: 144.2,
    61: 145.0, 62: 150.4, 63: 152.0, 64: 157.3, 65: 158.9, 66: 162.5, 67: 164.9, 68: 167.3, 69: 168.9, 70: 173.1,
    71: 175.0, 72: 178.5, 73: 180.9, 74: 183.8, 75: 186.2, 76: 190.2, 77: 192.2, 78: 195.1, 79: 197.0, 80: 200.6,
    81: 204.4, 82: 207.2, 83: 209.0, 84: 210.0, 85: 210.0, 86: 222.0, 87: 223.0, 88: 226.0, 89: 227.0, 90: 232.0,
    91: 231.0, 92: 238.0, 93: 237.0, 94: 244.0, 95: 243.0, 96: 247.0, 97: 247.0, 98: 251.0, 99: 252.0, 100: 257.0,
    101: 258.0, 102: 259.0, 103: 262.0, 104: 261.0, 105: 262.0, 106: 266.0, 107: 264.0, 108: 267.0, 109: 268.0, 110: 271.0,
    111: 272.0, 112: 285.0, 113: 284.0, 114: 289.0, 115: 288.0, 116: 292.0, 117: 295.0, 118: 294.0}
    if mol_name is not None:
        mof_data = f'data.emin_{mof}_{mol_name}'
        topo=f'emin_{mof}_{mol_name}'
    else:
        mof_data = f'data.emin_{mof}'
        topo=f'emin_{mof}'
    
    with open(mof_data, 'r') as file:
        full_data = file.read()
    lines = full_data.splitlines()
    start_index = full_data.find("Masses")
    end_index = full_data.find("Pair Coeffs # lj/cut", start_index) 
    end_index = next((full_data.find(line) for line in lines if "Pair" in line), len(full_data))
    masses_section = full_data[start_index:end_index].strip()
    masses_lines = masses_section.split('\n') 
   
    atom_type_dict = {}
    for line in masses_lines[1:]:
        parts = line.split()
        if len(parts) == 2:
            atom_type, mass = map(float, parts[:2])
            for atomic_number, mass_value in atomic_masses.items():
                if round(mass,1) == mass_value:
                    atom_type_dict[int(atom_type)] = atomic_number
                else:
                    continue
    
    print(f'Converting {topo} to CIF')
    cif_name = f"{topo}.cif"
    data = read(mof_data,format='lammps-data',Z_of_type=atom_type_dict)
    cif.write_cif(cif_name, data, cif_format=None, wrap=True, labels=None, loop_keys=None)
    

    
    
