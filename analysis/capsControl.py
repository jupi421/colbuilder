#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 18:19:13 2022

@author: broszms
"""
import os
import numpy as np
#
os.chdir(os.getcwd())
mixed=[]
mixIdx=[]
with open('triplehelices.txt') as f:
    for l in f:
        mixed.append([i.replace('\n','') for i in l.split(' ') if i!='\n'])
        mixIdx.append(l.split('.')[0])
f.close()
#
ACE=[]
NME=[]
for i in mixIdx:
    with open('f'+str(int(i))+'.caps.pdb','r') as f:
        for l in f:
            if str(l[17:20])=='ACE':
                ACE.append(i)
            if str(l[17:20])=='NME':
                NME.append(i)
#
ACEfinal=np.unique(ACE)
NMEfinal=np.unique(NME)
#
notACE=[]
notNME=[]
#
for i in range(len(mixIdx)): 
    if mixIdx[i] not in ACEfinal:
        notACE.append(mixIdx[i])
    if mixIdx[i] not in NMEfinal:
        notNME.append(mixIdx[i])
#
print('Missing ACE: ')
print([str(i)+' ' for i in notACE])
print('Missing NME: ')
print([str(i)+' ' for i in notNME])
print('# ACE: '+str(len(ACEfinal))+' # NME: '+str(len(NMEfinal)))
