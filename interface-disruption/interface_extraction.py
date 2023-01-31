#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 18:23:46 2022

@author: runfeng
"""

import os
import numpy as np
#%%
def interface_extraction(kind):
    interface = []
    BSA={}
    f=open(f'../pisa/{kind}')
    f=f.readlines()[1:]
    for line in f:
        if len(line)>1:
            s=line.strip().split('\t')
            
            s=list(filter(None,s))
            # print(s)
            if '##' not in s[0]:
                bsa=s[4].strip().split(' ')
                try:
                    float(s[4])
                except:
                    position = int(s[0])
                    
                    if position not in BSA:
                        
                        interface.append(position)
                        if len(interface)>0:
                            sp=interface[-1]
  
                            BSA[position]=float(bsa[0])
                    else:
                        
                        BSA[position]+=float(bsa[0])
        
    interface=sorted(list(set(interface)))
    
    
    interfaces=interface[:]
    for i,e in enumerate(interfaces):
        if i==0:
            pass
        else:
            if e-interfaces[i-1] <= 10:
                
                for x in range(e-interfaces[i-1]-1):
                    interface.append(interfaces[i-1]+x+1)
                interface=sorted(interface)
    interface=sorted(interface)
    # the complementary interface score is the average of that regions
    separated_interfaces={}
    z=1
    separated_interfaces[z]=[]
    for i,e in enumerate(interface):
        if len(separated_interfaces[z])==0:
            separated_interfaces[z].append(e)
        elif e-separated_interfaces[z][-1] <10:
            separated_interfaces[z].append(e)
        elif e-separated_interfaces[z][-1]>=10:
             z+=1
             separated_interfaces[z]=[e]

    score={}
    for keys,vals in separated_interfaces.items():
        s=0
        site=0
        for i in vals:
            if i in BSA:
                s+=BSA[i]
                site+=1
        score[keys]=s/site
    return BSA, interface, separated_interfaces, score

#%%

symmetric_BSA, symmetric_interface, symmetric_separated_interfaces, symmetric_score=interface_extraction('symmetric')

asymmetric_BSA, asymmetric_interface, asymmetric_separated_interfaces, asymmetric_score=interface_extraction('asymmetric')


si = symmetric_interface[:]
ai = asymmetric_interface[:]
for i in si:
    if i in ai:
        symmetric_interface.pop(symmetric_interface.index(i))

for i in ai:
    if i in si:
        asymmetric_interface.pop(asymmetric_interface.index(i))

for key,val in symmetric_separated_interfaces.items():
    vs = val[:]
    for i in val:
        if i not in symmetric_interface:
              vs.pop(vs.index(i))
    symmetric_separated_interfaces[key] = vs
    
for key,val in asymmetric_separated_interfaces.items():
    vs = val[:]
    for i in val:
        if i not in asymmetric_interface:
              vs.pop(vs.index(i))
    asymmetric_separated_interfaces[key] = vs

#%%


f=open(f'.../symmetric','w')
f.write('loc'+'\t'+'reg'+'\t'+'score'+'\n')
for key,val in symmetric_separated_interfaces.items():
    s=np.sqrt(symmetric_score[key])
    for i in val:
        f.write(str(i)+'\t'+str(key)+'\t'+f'{s:.2f}'+'\n')
f.close()

f=open(f'.../asymmetric','w')
f.write('loc'+'\t'+'reg'+'\t'+'score'+'\n')
for key,val in asymmetric_separated_interfaces.items():
    s=np.sqrt(asymmetric_score[key])
    for i in val:
        f.write(str(i)+'\t'+str(key)+'\t'+f'{s:.2f}'+'\n')
f.close()
