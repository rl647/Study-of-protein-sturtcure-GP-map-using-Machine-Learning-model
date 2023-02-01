

import subprocess as sp
import os
import random
import multiprocessing as mp
import re
import collections
from Bio.Emboss.Applications import NeedleCommandline
import shutil
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


#%%
protein = {}
# paths = []
path = '{path_to_pdb_file}'
paths = []
f=open('.../candidates')
for line in f:
    s=line.strip().split('\t')
    # paths.append(f'{path}/{s[0][:4]}')
    protein[s[0][:4]] = [s[0],f'{path}/{s[0][:4]}.pdb']
    paths.append(f'{path}/{s[0][:4]}.pdb')
f.close()    
#%%
asd = {}
paths=[]

random_paths = random.sample(paths, 1000)
#%%

from Bio import pairwise2
import multiprocessing as mp

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

from Bio import pairwise2
from Bio.pairwise2 import format_alignment            
secondary_structure = []
for val in random_paths[:]:


    p = PDBParser()
    
    structure = p.get_structure("mutant2", f"{val}")
    model = structure[0]
    
    dssp = DSSP(model, f"{val}", dssp='mkdssp')
    
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
    secondary_structure.append(sec_structure)
#%%
core = 76
a = ((len(secondary_structure)-1)*len(secondary_structure))/(2*core)
separate_dict = {}
n = 0
for i in range(1, core+1):
    # d[i]=[]
    if i == core:
        separate_dict[i] = [0, len(secondary_structure)-n]

    else:
        if i % 2 == 0:
            x = int((2*a+(n+1)**2-(n+1)+0.25)**0.5-0.5)+1
        else:
            x = int((2*a+(n+1)**2-(n+1)+0.25)**0.5-0.5)
        separate_dict[i] = [len(secondary_structure)-x, len(secondary_structure)-n]
        n = x


#%%
def ts_comparsions(k,paths=random_paths,sd=separate_dict):
    
    score = {}
    cp_paths = paths[sd[k][0]+1:]
    bench_list = paths[sd[k][0]:sd[k][1]]
    

    for bid, bench in enumerate(bench_list):
        t=0
        for cid, compare in enumerate(cp_paths):
            
            # print(compare)
            if compare in bench_list[:bid]  or compare==bench:
                pass
            else:
                aligned_residues = sp.getoutput(f'/home/runfeng/TMalign {bench} {compare}')
                    
                    
                a=aligned_residues.index("TM-score=")
                b=aligned_residues.index('(if normalized by length of Chain_2)')
                bench_score = aligned_residues[a+10:a+18]
                compare_score = aligned_residues[b-8:b]
                score[f'{bench[-10:-4]}_{compare[-10:-4]}']=[float(bench_score),float(compare_score)]
                    
                f=open(f'.../tertiary_score.txt','a')
                f.write(str(bench[-10:-4])+'\t'+str(compare[-10:-4])+'\t'+str(bench_score)+'\t'+str(compare_score)+'\n')
                f.close()
                t+=1
                if t%200 == 0:
                    print(str(cid)+'\n')
    return score
#%%
pool = mp.Pool(core)
pg = [pool.apply_async(ts_comparsions,args=(key,random_paths,separate_dict))for key in separate_dict.keys()]
pg = [p.get() for p in pg]

#%%



d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
seq = {}

path = 'Path_to_pdb'
# for line in f:
for file in os.listdir(path):
    filename = os.fsdecode(file)
    if filename.endswith('.pdb'):
        asd = {}
        seq[filename[:-4]] = []
        f=open(f'{path}/{filename}')
        for line in f:
            if int(line[22:26]) in asd:
                pass
            else:
                if line[17:20] in d3to1:
                    seq[filename[:-4]].append(d3to1[line[17:20]])
                else:
                    seq[filename[:-4]].append('X')
                asd[int(line[22:26])]=0
        f.close()
        seq[filename[:-4]] = ''.join(seq[filename[:-4]])

for k, v in seq.items():
    f=open(f'.../{k}.fasta','w')
    f.write('>'+k+'\n')
    f.write(v)
    f.close()
    
ps_paths = []
for i in random_paths:
    ps_paths.append(f'..../{i[-10:-4]}.fasta')
    




#%%

def ps_comparsions(k,paths=ps_paths,sd=separate_dict):
    path = f'...'
    score = {}
    cp_paths = paths[sd[k][0]+1:]
    bench_list = paths[sd[k][0]:sd[k][1]]
    

    for bid, bench in enumerate(bench_list):

        t=0
        for cid, compare in enumerate(cp_paths):

            if compare in bench_list[:bid]  or compare==bench:
                pass
            else:
                
                needle_cline = NeedleCommandline(asequence=bench, 
                                             bsequence=compare, 
                                             gapopen=10, gapextend=0.5, 
                                             endopen=10, endextend=0.5,
                                             outfile=f"{path}/needle_output/{bench[-12:-6]}_{compare[-12:-6]}.txt")
                                
                needle_cline()
                f=open(f"{path}/needle_output/{bench[-12:-6]}_{compare[-12:-6]}.txt")
                f = f.read()
                f = f[f.index('Identity:'):f.index('Identity:')+30]
                sim = float(f[f.index('(')+1:f.index('%')-1])
                os.remove(f"{path}/needle_output/{bench[-12:-6]}_{compare[-12:-6]}.txt")
                f=open(f'{path}/primary_score.txt','a')
                f.write(bench[-12:-6]+'\t'+compare[-12:-6]+'\t'+str(sim)+'\n')
                f.close()
                # shutil.rmtree(f"{path}/primary/{bench[-10:-6]}_{compare[-10:-6]}.txt")

#%%

pool = mp.Pool(core)
pg = [pool.apply_async(ps_comparsions,args=(key,ps_paths,separate_dict))for key in separate_dict.keys()]
pg = [p.get() for p in pg]



#%%


def ss_comparisons(k,paths, secondary_structure,sd):
    path = f'...'
    # print(sd[k])
    score = {}
    cp_paths = secondary_structure[sd[k][0]+1:]
    cp = paths[sd[k][0]+1:]
    bench_list = secondary_structure[sd[k][0]:sd[k][1]]
    # print(1)
    bp = paths[sd[k][0]:sd[k][1]]
    

    for bid, bench in enumerate(bench_list):
        t=0
        for cid, compare in enumerate(cp_paths):

            if compare in bench_list[:bid]  or compare==bench:
               pass
            else:
                
                ss1 = bench
                ss2 = compare
                # print(ss2)
                alignments = pairwise2.align.globalms(ss1, ss2,5,-4,-10,-0.5)
                if len(format_alignment(*alignments[0]).strip().split('\n'))>4:
                    alignments=list(filter(None,format_alignment(*alignments[0]).strip().split('\n')))[:3]
                else:
                    alignments=format_alignment(*alignments[0]).strip().split('\n')[:3]
                sim=alignments[1].count('|')/max([len(alignments[0]),len(alignments[2])])
                f=open(f'{path}/secondary_score.txt','a')
                f.write(bp[bid][-10:-4]+'\t'+cp[cid][-10:-4]+'\t'+str(sim)+'\n')
                f.close()
          

#%%

pool = mp.Pool(core)
pg = [pool.apply_async(ss_comparisons,args=(key,random_paths,secondary_structure,separate_dict))for key in separate_dict.keys()]
pg = [p.get() for p in pg]
