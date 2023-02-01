#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 15:41:10 2022

@author: runfeng
"""

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
protein='...'
primary = {}
secondary = {}
tertiary = {}

#%%
path = f'path_to_sequence'

for file in os.listdir(path):
    filename = os.fsdecode(file)
    if filename.endswith('.fasta'):
        primary[int(filename[:-6])] = f'{path}/{filename}'
        
path = f'path_to_pdb'

for file in os.listdir(path):
    filename = os.fsdecode(file)
    if filename.endswith('.pdb'):
        tertiary[int(filename[:-4])] = f'{path}/{filename}'
f=open(f'secondary_structure')
for line in f:
    s = line.strip().split('\t')
    secondary[int(s[0])] = s[1]    
        
#%%
path = f'.../needle_output'
for i in range(len(primary)-1):
    for e in range(i+1,len(primary)):
        
        needle_cline = NeedleCommandline(asequence=primary[i*100], 
                                     bsequence=primary[e*100], 
                                     gapopen=10, gapextend=0.5, 
                                     endopen=10, endextend=0.5,
                                     outfile=f"{path}/{i*100}_{e*100}.txt")
                        
        needle_cline()
        f=open(f"{path}/{i*100}_{e*100}.txt")
        f = f.read()
        f = f[f.index('Identity:'):f.index('Identity:')+30]
        sim = float(f[f.index('(')+1:f.index('%')-1])
        f=open(f'{path[:-14]}/primary_score.txt','a')
        f.write(str(i*100)+'\t'+str(e*100)+'\t'+str(sim)+'\n')
        f.close()
        os.remove(f"{path}/{i*100}_{e*100}.txt")
    

#%%
path = f'....'
for i in range(len(primary)-1):
    for e in range(i+1,len(primary)):
        ss1 = secondary[i*100]
        ss2 = secondary[e*100]
        alignments = pairwise2.align.globalms(ss1, ss2,5,-4,-10,-0.5)
        if len(format_alignment(*alignments[0]).strip().split('\n'))>4:
            alignments=list(filter(None,format_alignment(*alignments[0]).strip().split('\n')))[:3]
        else:
            alignments=format_alignment(*alignments[0]).strip().split('\n')[:3]
        sim=alignments[1].count('|')/max([len(alignments[0]),len(alignments[2])])
        f=open(f'{path}/secondary_score.txt','a')
        f.write(str(i*100)+'\t'+str(e*100)+'\t'+str(sim)+'\n')
        f.close()
        
#%%

path = f'...'
for i in range(len(primary)-1):
    for e in range(i+1,len(primary)):
        aligned_residues = sp.getoutput(f'/home/runfeng/TMalign {tertiary[i*100]} {tertiary[e*100]}')
            
            
        a=aligned_residues.index("TM-score=")
        b=aligned_residues.index('(if normalized by length of Chain_2)')
        bench_score = float(aligned_residues[a+10:a+18])
        compare_score = float(aligned_residues[b-8:b])
        sim=(bench_score+compare_score)/2
        f=open(f'{path}/tertiary_score.txt','a')
        f.write(str(i*100)+'\t'+str(e*100)+'\t'+str(sim)+'\n')
        f.close()
#%%






#%%
ps = {}
ss = {}
ts = {}



f=open(f'.../secondary_score.txt')
for line in f:
    s=line.strip().split('\t')
    if float(s[-1])>0.7:
        ss[s[0]+'_'+s[1]] = float(s[-1])
f.close()



f=open(f'.../primary_score.txt')
for line in f:
    s=line.strip().split('\t')
    if s[0]+'_'+s[1] in ss:
        # break
        
        ps[s[0]+'_'+s[1]] = float(s[-1])
f.close()

f=open(f'.../tertiary_score.txt')
for line in f:
    s=line.strip().split('\t')
    if s[0]+'_'+s[1] in ss:
            
        
        ts[s[0]+'_'+s[1]] = float(s[-1])
f.close()

#%%
x=[]
y=[]
z=[]
import matplotlib.pyplot as plt
for key,val in ps.items():
    x.append(val/100)
    y.append(ss[key])
    z.append(ts[key])
    
#%%
import numpy as np
width_height_1 = (10,5)
plt.figure(figsize=width_height_1)


plt.grid()

plt.xlabel('Similarity of Secondary Structure')
plt.ylabel('Similarity of Tertiary Structure')


plt.title('Secondary against Tertiary structure')
# plt.plot(y, z,'.')
cmap = plt.cm.rainbow

plt.scatter(y,z, c=x, cmap=cmap, edgecolor='none',s=50*np.ones(len(y)))
plt.xticks([0.1*i for i in range(7,11)])
plt.yticks([0.1*i for i in range(11)])
cbar = plt.colorbar()
cbar.set_label('Identity of seqeunces')














