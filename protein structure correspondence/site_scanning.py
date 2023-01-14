#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 14:56:10 2022

@author: runfeng
"""

import numpy as np
import os
import random
import shutil
#%%
protein=''
ori_seq = []
f=open(f'.../{protein}.fasta')
ori_seq=f.readlines()[1][:-1]
f.close()
#%%
ori_ss = []
f=open(f'.../{protein}.fasta.ss3')
f=f.readlines()[1:]
for line in f:
    s=line.strip().split('\t')
    ori_ss.append(s[2])
ori_ss = ''.join(ori_ss)
#%%
alphabet=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' ]
def protein_site_scanning(ori_seq,ori_ss,alphabet,size,protein):
    # mutated_seq = list(ori_seq[:])
    mutated_seq = list(ori_seq[:])
    seq = {}
    seq[ori_seq] = 0 
    while len(seq)<=size+1:
        neutral_seq = list(seq.keys())[-1] 

        for indexes, residues in enumerate(neutral_seq):
            if len(seq)==size+1:
                break
            ap = alphabet[:]
            random.shuffle(ap)
            if residues != ori_seq[indexes]:
                ap.pop(ap.index(residues))
                ap.pop(ap.index(ori_seq[indexes]))
            else:
                ap.pop(ap.index(residues))
            for i, e in enumerate(ap):
                mutated_seq=list(neutral_seq)
                mutated_seq[indexes] = e
                mutated_seq=''.join(mutated_seq)
                os.mkdir(f'...')
                f=open(f'.../{protein}.fasta','w')
                f.write(f'> {protein}'+'\n'+mutated_seq)
                f.close()
                os.system(f'Porter5.py -i .../{protein}.fasta --cpu 7 --fast')
                mutated_ss = []
                f=open(f'.../{protein}.fasta.ss3')
                f=f.readlines()[1:]
                for line in f:
                    s=line.strip().split('\t')
                    mutated_ss.append(s[2])
                mutated_ss=''.join(mutated_ss)
                if mutated_ss == ori_ss:
                    print(len(seq))
                    seq[mutated_seq] = 1
                    shutil.rmtree(f'...')
                    break
                else:
                    shutil.rmtree(f'...')
    return seq

#%%
size=5000
seq = protein_site_scanning(ori_seq, ori_ss, alphabet, size, protein)
