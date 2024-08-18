import pormake as pm
import os
from contextlib import redirect_stdout, redirect_stderr

def build_mof(node_name, linker_name, topo_list, forcefield, replication):
    if node_name.endswith('.xyz'):
        node = pm.BuildingBlock(f'{node_name}')
        node_name = node_name.split('.')[0]
    else:
        node = pm.BuildingBlock(f'{node_name}.xyz')

    if linker_name.endswith('.xyz'):
        linker = pm.BuildingBlock(f'{linker_name}')
        linker_name = linker_name.split('.')[0]
    else:
        linker = pm.BuildingBlock(f'{linker_name}.xyz')
    
    node_c_number = node.atoms.symbols.count('X')
    linker_c_number = linker.atoms.symbols.count('X')
    
    database = pm.Database()
    builder = pm.Builder()
    topo_list = topo_list.split(',')
    for topo in topo_list:
        print("############### Working on "+topo+" ###############")
        topo_name = database.get_topo(topo)
        validaty = str(topo_name.check_validity).split('(')[1].split(')')[0].split(',')
        print("Description on the Topology")
        print(topo_name.check_validity)
        
        node_bbs = {}
        edge_bbs = {}
        
        if linker_c_number == 2: # it is an edge, else it is an organic node
            node_bbs[0] = node
            edge_bbs = {(0, 0): linker} # only a single node is supported 
        elif node_c_number == linker_c_number: # same connecting number for linker and node
            node_bbs[0] = node
            node_bbs[1] = linker
        else:
            for i, num in enumerate(validaty): # different connecting number for linker and node
                if int(num) == node_c_number:
                    node_bbs[i] = node
                else:
                    node_bbs[i] = linker
        try:
            with open(os.devnull, 'w') as fnull:
                with redirect_stdout(fnull), redirect_stderr(fnull):
                    mof = builder.build_by_type(topology=topo_name, node_bbs=node_bbs, edge_bbs=edge_bbs)
            
            mof.write_cif(f'./{topo}-{node_name}-{linker_name}.cif')
            if replication is not None:
                command = f'lammps-interface {topo}-{node_name}-{linker_name}.cif -ff {forcefield} --replication {replication} > tmp; rm tmp'
                os.system(command)
            else:
                command = f'lammps-interface {topo}-{node_name}-{linker_name}.cif -ff {forcefield} > tmp; rm tmp'
                os.system(command)
            print("Done! Check the cif and data files for further simulation")

        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            print("It is probably due to the incompatibility of building blocks and topology. Change to a correct topology.")
