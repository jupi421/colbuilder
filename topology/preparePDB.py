#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 11:59:04 2022

@author: broszms
"""
pdb=[]
with open('tmp.pdb','r') as f:
    pdb=[l for l in f]
#
iT=0
iTMap=0
tmp=0
start=int(pdb[1][22:26])-1
out=[]
outMap=[]
out.append(pdb[0])
outMap.append(pdb[0])
#
for i in range(1,len(pdb)):
    if pdb[i][0:3]=='END':
        break
    if pdb[i][0:3]=='TER':
        continue
    #
    tmp=(start-int(pdb[i][22:26]))
    if tmp<0:
        start=int(pdb[i][22:26])
        iT+=1
        iTMap+=1
    #
    if tmp>0 and pdb[i][21:22]=='B' or tmp>0 and pdb[i][21:22]=='C' or tmp>0 and pdb[i][21:22]=='A':
        start=int(pdb[i][22:26])
        iT=1
        iTMap+=1
    #
    if iT<10:
        out.append(pdb[i][:22]+'   '+str(int(iT))+pdb[i][26:])
    if iT>=10 and iT<100:
        out.append(pdb[i][:22]+'  '+str(int(iT))+pdb[i][26:])
    if iT>=100 and iT<1000:
        out.append(pdb[i][:22]+' '+str(int(iT))+pdb[i][26:])
    if iT>=1000 and iT<10000:
        out.append(pdb[i][:22]+''+str(int(iTMap))+pdb[i][26:])
    #
    if iTMap<10:
        outMap.append(pdb[i][:22]+'   '+str(int(iTMap))+pdb[i][26:])
    if iTMap>=10 and iTMap<100:
        outMap.append(pdb[i][:22]+'  '+str(int(iTMap))+pdb[i][26:])
    if iTMap>=100 and iTMap<1000:
        outMap.append(pdb[i][:22]+' '+str(int(iTMap))+pdb[i][26:])
    if iTMap>=1000 and iTMap<10000:
        outMap.append(pdb[i][:22]+''+str(int(iTMap))+pdb[i][26:])
#
out.append(pdb[-1])
with open('tmp.order.pdb','w') as f:
    for l in out:
        f.write(l)
#
with open('tmp.map.pdb','w') as f:
    for l in outMap:
        f.write(l)
        
