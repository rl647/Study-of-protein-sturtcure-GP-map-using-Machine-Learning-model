#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 09:48:12 2022

@author: runfeng
"""

import os
import numpy as np
import random
import shutil

#%%
#protein sequences alphabets
alphabet=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' ]
#%%
protein='1b99'
kind = 'asymmetric'
#%%
# predict the original secondary structure
os.system(f'python3 .../Porter5.py -i .../{protein}.fasta --cpu 7 --fast')

#%%
# orginal secondary structures
ss = []
# orignal protein sequence
seq = []
f=open(f'.../{protein}.fasta.ss3')
f=f.readlines()[1:]
for line in f:
    s=line.strip().split('\t')
    ss.append(s[2])
    seq.append(s[1])
    
ss=''.join(ss)
seq=''.join(seq)

#%%
# mutating interface 
interface = {}
f=open(f'.../{kind}')

f=f.readlines()[1:]
for line in f:
    s=line.strip().split('\t')
    interface[int(s[0])]=float(s[2])
#%%
exclusion = []

target_position=[]

for i in interface.keys():
    if ss[i-1] =='E' or ss[i-1]=='H':
        target_position.append(i)

#%%

def mutations(protein=protein, kind=kind, seq=seq, interface=interface, target_loc=target_position,ss=ss, alphabet=alphabet):
    # for a certain score record all the possible phenotypes 
    score_ss = {}
    # for a certain ss record all the possible sequences
    ss_seuqence = {}

    # create path for the selected interface

    if not os.path.exists('...'):
        os.mkdir('...')

    mutated_seq = list(seq[:])
    all_score=[0]
    for i in interface.keys():

        # remove the original residues and shuffle the alphabet list
        ap=alphabet[:]
        ap.pop(ap.index(seq[i-1]))
        random.shuffle(ap)
        for index, elements in enumerate(ap):
            mutated_seq[i-1] = elements
            mutated_seq = ''.join(mutated_seq)
            os.mkdir(f'.../{i}')
            f=open(f'.../{protein}.fasta','w')
            f.write('>'+protein+'\n'+mutated_seq)
            f.close()
            os.system(f'python3 .../Porter5.py -i .../{i}/{protein}.fasta --cpu 7 --fast')
            current_ss = []
            
            f=open(f'.../{i}/{protein}.fasta.ss3')
            f=f.readlines()[1:]
            for line in f:
                s=line.strip().split('\t')
                current_ss.append(s[2])
                
                
            shutil.rmtree(f'.../{i}',ignore_errors=True)
            score = 0
            for (idx,mutant_sites), original_sites in zip(enumerate(current_ss), ss):
                if mutant_sites != original_sites  and idx+1 in target_loc and mutant_sites == 'C':
                    score+=interface[idx+1]
                elif mutant_sites != original_sites  and idx+1 not in interface:
                    score-=interface[i]
                
                    #score-=0.2*interface[i]
                elif mutant_sites != original_sites  and idx+1 in interface:
                    score+=0.1
            
                
            current_ss=''.join(current_ss)
            
            if current_ss == ss:

                score+=0.1
            print('score: ',score)
            
            if score in score_ss:
                score_ss[score].append(current_ss)
            else:
                score_ss[score]=[current_ss]
            if current_ss in ss_seuqence:
                ss_seuqence[current_ss].append(mutated_seq)
            else:
                ss_seuqence[current_ss]=[mutated_seq]
            
            sc = max(list(score_ss.keys()))
            if sc>0:
                
                mutated_seq = list(random.choice(ss_seuqence[random.choice(score_ss[sc])]))
            else:
                mutated_seq = list(seq[:])
                    
    for key, val in score_ss.items():
        score_ss[key] = list(set(val))
    return score_ss, ss_seuqence
#%%

score_ss, ss_seuqence=mutations()










                
                
                
            
            
        
        











    













































































