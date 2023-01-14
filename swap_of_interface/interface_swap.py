#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 14:51:59 2022

@author: runfeng
"""
import subprocess as sp

import os

#%%
    
def interface_extraction(pisa):
    # extract the interface from pisa
    interface = []
    BSA={}
    f=open(f'{pisa}')
    f=f.readlines()[1:]
    for line in f:
        if len(line)>1:
            s=line.strip().split('\t')
            s=list(filter(None,s))
            # print(s)
            if '##' not in s[0]:
                bsa=s[4].strip().split(' ')
                
                # print(bsa)
                try:
                    float(s[4])
                except:
                    # if the current interface site deviate from the last interface site for less than 10 sites positions,
                    # the sites between them will be extracted as well
                    if len(interface)>0:
                        sep=interface[-1]
                        BSA[int(s[1][-8:])]=float(bsa[0])
                        if int(s[1][-8:])-sep<10:
                            for i in range(int(s[1][-8:])-sep):
                                interface.append(sep+1+i)
                        else:
                            interface.append(int(s[1][-8:]))
                            BSA[int(s[1][-8:])]=float(bsa[0])
                    else:
                        interface.append(int(s[1][-8:]))
                        BSA[int(s[1][-8:])]=float(bsa[0])
    
    
    # separate the interface into different regions
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
    return interface,separated_interfaces,BSA

#%%
def interface_switch(seq_path,apisa,bpisa,aname,bname,path,good):
 
    #apply TMalign
    
    aligned_residues = sp.getoutput(f'.../TMalign .../{aname}.pdb .../{bname}.pdb ')
    # print(aligned_residues)
    a=aligned_residues.index("TM-score=")
    b=aligned_residues.index('(if normalized by length of Chain_2)')
    bench_score = float(aligned_residues[a+10:a+18])
    compare_score = float(aligned_residues[b-8:b])
    # save the TMscore of alignment
    f=open(f'{path}/TM_score.txt','a')
    f.write(aname+'\t'+bname+'\t'+f'{bench_score:.3f}'+'\t'+f'{compare_score:.3f}'+'\n')
    f.close()
    
    
    s = int(aligned_residues.index('aligned residues)'))+len('aligned residues)')
    # extract the alignment residues from TMalign output
    aligned_residues=aligned_residues[s:]
    aligned_residues=aligned_residues.strip().split('\n')
    # print(aligned_residues)

    interface,separated_interfaces,BSA=interface_extraction(apisa)

    # extract the onterface of bsequence
    bseq_interface,bseq_separated_interfaces,bBSA=interface_extraction(bpisa)

    
    # record the interface in the TM alignment
    bsi = {}
    for key, val in bseq_separated_interfaces.items():
        for i,e in enumerate(aligned_residues[2]):
            if e=='-' and i+1<=val[0]:
                val=[x+1 for x in val]
            elif e=='-' and i+1<=val[-1] and i+1>val[0]:
                val.append(val[-1]+1)
            elif i+1 >val[-1]:
                break
        bsi[key] = val




    # sequence of heteromer
    aaligned=list(aligned_residues[0])
    # sequence of homomer
    baligned=list(aligned_residues[2])
    # indication of alignment between homomer and heteromer
    pair=list(aligned_residues[1])
    # regions that will be replaced and include interfaces
    sp_interface={}
    
    
    for keys,vals in separated_interfaces.items():
    
        
        for i, e in enumerate(aligned_residues[0]):
            if e=='-' and i+1<=vals[0]:
                vals=[x+1 for x in vals]
            elif e=='-' and i+1<=vals[-1] and i+1>vals[0]:
                vals.append(vals[-1]+1)

        vs=vals[:]
        s=0
        z=vs[0]-2-s
        v0=vs[0]
        good = good
        good_align=[':' for d in range(good)]
        # extend the regions until they are well aligned
        while pair[z:z+good] != good_align and z>=0:
            vs.insert(0,v0-1-s)
            s+=1
            z-=1
    
        s=0
        z=vs[-1]+s
        v1=vs[-1]
        while pair[z-good:z] != good_align:
            vs.append(v1+1+s)
            s+=1
            z+=1
    
            if z>=len(pair)-1:
                break

        sp_interface[keys] = vs   
        baligned[min(vs)-1:max(vs)] = aaligned[min(vs)-1:max(vs)] 
    aaligned=''.join(aaligned)
    baligned=''.join(baligned)
    


    # record the replaced region and its complementary region
    replaced_region={}
    replaced_region['a']={}
    replaced_region['b']={}
    domain_b = {}
    domain_region = {}
    domain_region['before'] = {}
    domain_region['after'] = {}
    
    other_interface={}
    other_interface['before'] = {}
    other_interface['after'] = {}
    
    # this for loop is to record the positions of interface of  b and  a
    si=list(sp_interface.values())
    for key, val in sp_interface.items():
        v = val[:]
        for i,e in enumerate(baligned):
            if e == '-' and i<min(val)-1:
                v=[s-1 for s in v]   
            elif e == '-' and i>=min(val)-1 and i<max(val):
                v.pop()
            elif i>=max(val):
                break
            
        replaced_region['b'][key]=v
        
        v = val[:]
        for i,e in enumerate(aaligned):
            if e == '-' and i<min(val)-1:
                v=[s-1 for s in v]
            elif e == '-' and i>=min(val)-1 and i<max(val):
                v.pop()
            elif i>=max(val):
                break
        replaced_region['a'][key]=v
        
        
        
        
    # this for loop is to record the complementary part of interface of b and a
    domain_b[0]=[i+1 for i in range(min(si[0])-1)]
    for i, e in enumerate(si):
        if i == len(si)-1:
            domain_b[i+1] = [s for s in range(max(e)+1,len(baligned)+1)]
        else:
            domain_b[i+1] = [s for s in range(max(e)+1,min(si[i+1]))]
    
    for key, val in domain_b.items():
        if len(val) == 0:
            pass
        else:
            v = val[:]
            for i,e in enumerate(baligned):
                if e == '-' and i<min(val)-1:
                    v=[s-1 for s in v]
                elif e == '-' and i>=min(val)-1 and i<max(val):
                    v.pop()
                elif i>=max(val):
                    break
    
            domain_region['after'][key]=v
            
            v = val[:]
            for i,e in enumerate(aligned_residues[2]):
                if e == '-' and i<min(val)-1:
                    v=[s-1 for s in v]
                elif e == '-' and i>=min(val)-1 and i<max(val):
                    v.pop()
                elif i>=max(val):
                    break
                
            domain_region['before'][key]=v
            
    
    other_interface={}
    other_interface['before'] = {}
    other_interface['after'] = {}
    # this loop is to record the interface postions in pdb files
    for key, val in bsi.items():
        v = val[:]
        for i,e in enumerate(baligned):
            if e == '-' and i<min(val)-1:
                v=[s-1 for s in v]
            elif e == '-' and i>=min(val)-1 and i<max(val):
                v.pop()
            elif i>=max(val):
                break
            
        other_interface['after'][key]=v
        v = val[:]
        for i,e in enumerate(aligned_residues[2]):
            if e == '-' and i<min(val)-1:
                v=[s-1 for s in v]
            elif e == '-' and i>=min(val)-1 and i<max(val):
                v.pop()
            elif i>=max(val):
                break
        other_interface['before'][key]=v 
    
    baligned = [x for x in baligned if x!='-']
    baligned=''.join(baligned)
        
    return domain_region,replaced_region,baligned





