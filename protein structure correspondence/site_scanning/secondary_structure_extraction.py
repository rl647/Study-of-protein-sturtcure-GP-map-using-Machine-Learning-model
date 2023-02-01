from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os
import numpy as np
import shutil
#%%
protein = '...'
paths = []
path = f'path_to_pdb'
for file in os.listdir(path):
    filename = os.fsdecode(file)
    if filename.endswith('.pdb'):
        paths.append(f'{path}/{filename}')
#%%

for i in paths:
    p = PDBParser()
    
    structure = p.get_structure("mutant2", f"{i}")
    model = structure[0]
    
    dssp = DSSP(model, f"{i}", dssp='mkdssp')
    
    sequence = ''
    sec_structure = ''
    for z in range(len(dssp)):
        a_key = list(dssp.keys())[z]
        sequence += dssp[a_key][1]
        if dssp[a_key][2] == '-' or dssp[a_key][2] == 'S' or dssp[a_key][2] == 'T'  or dssp[a_key][2] == 'I':
            sec_structure += 'C'
        elif dssp[a_key][2] == 'G' or dssp[a_key][2] == 'H':
            sec_structure += 'H'
        else:

            sec_structure += 'E'

    f=open(f'secondary_structure','a')
    f.write(i[47:-4]+'\t'+sec_structure+'\n')
    f.close()
