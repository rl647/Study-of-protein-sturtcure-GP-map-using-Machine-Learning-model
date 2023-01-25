
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 21:43:11 2022

@author: runfeng
"""

import random
import numpy as np
import os
import multiprocessing as mp
import shutil
#%%

protein='...'
#%%
alphabet=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' ]

# the original sequence of the candidates

f=open(f'.../{protein}.fasta')
original_sequence=f.readlines()[1:][0][:-1]
f.close()


#%%

os.system(f'python3 .../Porter5.py -i .../{protein}.fasta --cpu 7 --fast')


#%%

mutated_sequence = {}
# paths to all the folder
paths = []

# all the candidates of mutated positions
loc = [i for i in range(len(original_sequence))]
# randomly select k starting points
choice =random.choices(loc,k=200000)
    
#%%
def mutation(c,choice):

    
    for i, e in enumerate(choice):
        if i%10000==0:
            print(i)
        
        
		# the positions that has not been mutated
        mpos = loc[:]
	
        mpos.pop(mpos.index(e))
		
		# the alphabet of amino acids
        a=alphabet[:]

	    
		# replace the selected sites by new amino acids 
        if original_sequence[e] in a:
            a.pop(a.index(original_sequence[e]))
        else:
            pass
        mutated_sequence[i+1]={}
        seq = list(original_sequence[:])
        seq[e] = random.choice(a)
	    
		# save the mutated sites to the dctionary
        mutated_sequence[i+1][1] = ''.join(seq)


        for t in range(len(original_sequence)-1):
            a=alphabet[:]
            # randomly selected new position from the positions that have not been mutated
            sp = random.choice(mpos)
            # update the list of positions that have not been mutated
            mpos.pop(mpos.index(sp))
            # the mutated sequences with 1 less positions changes
            lseq = list(mutated_sequence[i+1][t+1])   
            # remove the amino acid on the current selected postions
            if lseq[sp] in a:
                a.pop(a.index(lseq[sp]))
            else:
                pass
                # randomly select a aminod acid and replace the selected position by new amino acid
            lseq[sp] = random.choice(a)
            # save the mutant to the dictionary
            mutated_sequence[i+1][t+2] = ''.join(lseq)
            # save the mutant to the selected path
              

            f=open(f'.../mutant1.txt','a')
            f.write(''.join(lseq)+'\n')
            f.close()


            
#%%
asd = mutation(1,choice)

#%%

