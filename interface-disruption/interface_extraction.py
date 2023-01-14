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
    f=open('....')
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


