#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: broszms
"""
import sys
import pandas as pd
import numpy as np
import MDAnalysis as mda
import math
#
fileName=str(sys.argv[1])
u=mda.Universe(fileName,fileName)
#
# Get Crosslinks for each Type
divalent=False
if divalent==True:
    cross1='L4Y'
    atom1='CE'
    cross2='L5Y'
    atom2='NZ' # closest connection between crosslinks
if divalent==False:
    cross1='LYX'
    atom1='C13'
    cross2='LY3'
    atom2='CG' # closest connection between crosslinks
#
data=pd.DataFrame(columns=['Idx','chainNr','chainId','resname','x','y','z'])
cnt=0
chainNr=0
chainBool=False
for a in u.atoms:
    if a.resname==cross1 and a.name==atom1 or a.resname==cross2 and a.name==atom2:
        chainNr=math.floor(a.segindex/3)
        data.loc[cnt] = [a.id,chainNr,a.chainID,a.resname,a.position[0],a.position[1],a.position[2]]
        cnt+=1
#
# Identify which Chains or Crosslinks are connected
chainIds=[]
pairs=[]
tmpPairs=[]
cnt=0
minDistAtom12=2.0
#
for k in range(len(data)):
    i=data.loc[k]
    vCA=np.array([i.x,i.y,i.z])
    chainIds.append(i.chainNr)
    cnt+=1
    for l in range(len(data[cnt:])):
        j=data.loc[l+cnt]
        if i.chainNr==j.chainNr and j.resname==i.resname:
            continue
        wCA=np.array([j.x,j.y,j.z])
        deltaCA=np.linalg.norm(wCA-vCA)
        if deltaCA < minDistAtom12:
            tmpPairs=np.array([i.Idx,i.chainNr,j.Idx,j.chainNr,deltaCA],dtype=object)
            break
    # Get Pairs of Triple Helices within a distance of deltaCA = 2 A
    pairs.append(tmpPairs) 
#
# Get List of Nearest Neighbors between Chains. 
line=[]
uniqueChains=[]
for i in pairs:
    tmp=[]
    for j in pairs:
        if j[1]==i[1]:
            tmp.append(j[3])
        if j[1]>i[1]:
            break
    if tmp!=line:
        tmp2=[[i[1]],list(np.unique(tmp))]
        # +1 to take iterator start from 0 and not from 1 into account
        uniqueChains.append([i+1 for k in tmp2 for i in k])
    line=tmp
#
# Get single chains without neighbours
chainIds=np.unique(chainIds)
for i in chainIds:
    out=True
    for j in pairs:
        if i in j: out=False
    if out==True:
        uniqueChains.append([i+1]) 
#
tmp=[]
mergedIds=[]
# Identify & Merge Chains with intersection
for ci,vi in enumerate(uniqueChains):
    for i in vi:
        for ck,vk in enumerate(uniqueChains[ci+1:]):
            if i in vk:
                tmp.append([vi,vk])
# 
mergedIds=[list(np.unique(i)) for i in tmp]
for k in mergedIds:
    for ci,vi in enumerate(uniqueChains):
        if k[0]==vi[0]:
            uniqueChains[ci]=k
        if k[1]==vi[0]:
            uniqueChains[ci]='Delete'
#
uniqueChains=[i for i in uniqueChains if i !='Delete']
with open ('AllTriplehelices.txt','w') as f:
    for l in uniqueChains:
        tmp=[str(i)+'.caps.pdb' for i in l]
        f.write(' '.join(tmp)+'\n')
#
with open ('triplehelices.txt','w') as f:
    for l in uniqueChains:
        tmp=[str(i)+'.caps.pdb' for i in l]
        if len(l)==1:
            continue
        f.write(' '.join(tmp)+'\n')
