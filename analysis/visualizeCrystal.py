#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 17:05:44 2022

@author: broszms
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#
#
# P1 means the following rotation matrix: Triclinic unit cell
# Rot Matrix = [ [a, b*cos(gamma), c*cos(beta) ], 
#                [ 0, b*sin(gamma), c * (cos(alpha)-cos(beta)*cos(gamma)) / sin(gamma) ]
#                [0,0,c * sqrt(1-cos^2(beta)-(cos(alpha)-cos(beta)*cos(gamma)) /
#                 sin(gamma)*(cos(alpha)-cos(beta)*cos(gamma)) / sin(gamma)  ] ]
rotMatrix=np.matrix([[39.97,-7.238,-54.2489],[0,25.9598,-5.79142],[0,0,675.701]])
#
#
os.chdir(os.getcwd())
#
unpair=[]
with open('AllTriplehelices.txt','r') as f:
    for l in f:
        if len(l.split(' '))>1:
            continue
        unpair.append('Model '+str(float(l.split('.')[0])))
f.close()
#
# Get modified Chimera Symmetry Information
symContacts=[]
with open('CrystalContactsSym.txt') as f:
    symContacts=[l for l in f]
f.close()
#
symData=pd.DataFrame(columns=('modId','x','y','z'))
cnt=0
for i in range(0,int(len(symContacts)),4):
    symData.loc[cnt]=[symContacts[i].replace('\n',''),float(symContacts[i+1].split(' ')[-1].replace('\n','')),float(symContacts[i+2].split(' ')[-1].replace('\n','')),float(symContacts[i+3].split(' ')[-1].replace('\n',''))]
    cnt+=1
# Get original crystalcontact information
contacts=[]
with open('CrystalContacts.txt') as f:
    contacts=[l for l in f]
f.close()
#
data=pd.DataFrame(columns=('modId','x','y','z'))
cnt=0
for i in range(0,int(len(contacts)),4):
    data.loc[cnt]=[contacts[i].replace('\n',''),float(contacts[i+1].split(' ')[-1].replace('\n','')),float(contacts[i+2].split(' ')[-1].replace('\n','')),float(contacts[i+3].split(' ')[-1].replace('\n',''))]
    cnt+=1
#
# Transform Rot-Trans to Ids for modified Chimera
symChimera=pd.DataFrame(columns=('modId','ix','iy','iz'))
for i in range(0,len(symData)):
    b=[symData.loc[i].x,symData.loc[i].y,symData.loc[i].z]
    tmp=np.linalg.solve(rotMatrix,b)
    symChimera.loc[i]=[symData.loc[i].modId,np.round(tmp[0],0),np.round(tmp[1],0),np.round(tmp[2],0)]
#
# ... and for Chimera
chimera=pd.DataFrame(columns=('modId','ix','iy','iz'))
for i in range(0,len(data)):
    b=[data.loc[i].x,data.loc[i].y,data.loc[i].z]
    tmp=np.linalg.solve(rotMatrix,b)
    chimera.loc[i]=[data.loc[i].modId,np.round(tmp[0],0),np.round(tmp[1],0),np.round(tmp[2],0)]

# Get modified Post-Chimera model ids
cnt=0
pairModel=pd.DataFrame(columns=('modId','ix','iy','iz'))
for i in range(len(symChimera)):
    if str(symChimera.loc[i].modId) not in unpair:
        continue
    pairModel.loc[cnt]=[symChimera.loc[i].modId,symChimera.loc[i].ix,symChimera.loc[i].iy,symChimera.loc[i].iz]
    cnt+=1
#
cnt=0
subPairModel=pd.DataFrame(columns=('modId','ix','iy','iz'))
for i in range(len(pairModel)):
    for j in range(len(chimera)):
        if pairModel.loc[i].ix==chimera.loc[j].ix and pairModel.loc[i].iy==chimera.loc[j].iy and pairModel.loc[i].iz==chimera.loc[j].iz:
            subPairModel.loc[cnt]=[pairModel.loc[i].modId,pairModel.loc[i].ix,pairModel.loc[i].iy,pairModel.loc[i].iz]
            cnt+=1
            break
#
plt.rcParams["figure.figsize"] = (20,20)
plt.rcParams["font.size"] = (20)
fig=plt.figure(1)
ax=fig.add_subplot(projection='3d')
ax.scatter(subPairModel.ix,subPairModel.iy,subPairModel.iz,color='red',label='Post-Chimera',marker='o',s=200)
ax.scatter(chimera.ix,chimera.iy,chimera.iz,color='blue',label='Chimera',marker='o',s=200)
#ax.scatter(symChimera.ix,symChimera.iy,symChimera.iz,color='green',label='SymChimera',marker='o',s=200)
ax.set_xlabel('Dimension x')
ax.set_ylabel('Dimension y')
ax.set_zlabel('Dimension z')
ax.view_init(15,160)
plt.legend()
plt.savefig('FinalContacts.png',format='png')