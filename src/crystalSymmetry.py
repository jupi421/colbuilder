#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 14:36:31 2022

@author: broszms
"""
import os
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
#
path=os.getcwd()
os.chdir(path)
#
# P1 means the following rotation matrix: Triclinic unit cell
# Rot Matrix = [ [a, b*cos(gamma), c*cos(beta) ], 
#                [ 0, b*sin(gamma), c * (cos(alpha)-cos(beta)*cos(gamma)) / sin(gamma) ]
#                [0,0,c * sqrt(1-cos^2(beta)-(cos(alpha)-cos(beta)*cos(gamma)) /
#                 sin(gamma)*(cos(alpha)-cos(beta)*cos(gamma)) / sin(gamma)  ] ]
rotMatrix=np.matrix([[39.97,-7.238,-54.2489],[0,25.9598,-5.79142],[0,0,675.701]])
#
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
rotId=pd.DataFrame(columns=('modId','ix','iy','iz'))
for i in range(0,len(data)):
    b=[data.loc[i].x,data.loc[i].y,data.loc[i].z]
    tmp=np.linalg.solve(rotMatrix,b)
    rotId.loc[i]=[data.loc[i].modId,np.round(tmp[0],0),np.round(tmp[1],0),np.round(tmp[2],0)]
#
zMin=np.min(rotId.iz)
zMax=np.max(rotId.iz)
yMax=np.max(rotId.iy)
yMin=np.min(rotId.iy)
xMax=np.max(rotId.ix)
xMin=np.min(rotId.ix)
#
zSymmetry=[]
crystalRot=[[rotId.loc[i].ix,rotId.loc[i].iy,rotId.loc[i].iz] for i in range(len(rotId))]
#
for i in range(len(rotId)):
    #if rotId.loc[i].iy==yMin or rotId.loc[i].iy==yMax:
    #    continue
    # 
    if rotId.loc[i].iz==zMin or rotId.loc[i].iz==zMax or rotId.loc[i].iz==zMin+1 or rotId.loc[i].iz==zMax-1:
        tmp=[rotId.loc[i].ix,-rotId.loc[i].iy,rotId.loc[i].iz]
        if tmp not in crystalRot:
            zSymmetry.append(tmp)
#
# Add Extrema
if [xMin,yMax-1,zMin] not in zSymmetry and [xMin,yMax-1,zMin] not in crystalRot:
    zSymmetry.append([xMin,yMax-1,zMin])
if [xMax,yMax-1,zMax] not in zSymmetry and [xMax,yMax-1,zMax] not in crystalRot :
    zSymmetry.append([xMax,yMax-1,zMax])
if [xMin,yMin+1,zMin] not in zSymmetry and [xMin,yMin+1,zMin] not in crystalRot:
    zSymmetry.append([xMin,yMin+1,zMin])
if [xMax,yMin+1,zMax] not in zSymmetry and [xMax,yMin+1,zMax] not in crystalRot:
    zSymmetry.append([xMax,yMin+1,zMax])
#
# Manual addition of first and last layer
zSymmetry.append([-4,-3,zMin])
zSymmetry.append([-4,3,zMin])
zSymmetry.append([4,-3,zMax])
zSymmetry.append([4,3,zMax])
#
zSymmetry.append([-3,0,zMin])
zSymmetry.append([-8,0,zMin])
zSymmetry.append([3,0,zMax])
zSymmetry.append([8,0,zMax])
# Add Symmetry for y=0
zSymmetry.append([-1,0,zMin+1])
zSymmetry.append([-6,0,zMin+1])
zSymmetry.append([1,0,zMax-1])
zSymmetry.append([6,0,zMax-1])
# Add Symmetry for y=+/-1
zSymmetry.append([-8,-1,zMin])
zSymmetry.append([-8,1,zMin])
zSymmetry.append([-3,-1,zMin])
zSymmetry.append([-3,1,zMin])
zSymmetry.append([8,-1,zMax])
zSymmetry.append([8,1,zMax])
zSymmetry.append([3,-1,zMax])
zSymmetry.append([3,1,zMax])
#
zSymmetry.append([-1,-1,zMin+1])
zSymmetry.append([-1,1,zMin+1])
zSymmetry.append([1,-1,zMax-1])
zSymmetry.append([1,1,zMax-1])
zSymmetry.append([-6,-1,zMin+1])
zSymmetry.append([-6,1,zMin+1])
zSymmetry.append([6,-1,zMax-1])
zSymmetry.append([6,1,zMax-1])
#
#zSymmetry.append([1,-2,-1])
zSymmetry.append([1,-1,-1])
#zSymmetry.append([-3,2,-1])
#
#zSymmetry.append([-1,2,1])
zSymmetry.append([-1,1,1])
#zSymmetry.append([3,-2,1])
#
# Generate Rot-Trans-Matrices for symmetric models
rotSym=[[] for i in range(len(zSymmetry))]
for i in range(len(zSymmetry)):
    rotSym[i]=np.dot(rotMatrix,zSymmetry[i])
#
# Generate translation for symmetry candidates
cnt=0
lenRotId=0
crystalcontact=pd.DataFrame(columns=('modId','x','y','z'))
crystalcontactId=pd.DataFrame(columns=('modId','ix','iy','iz'))
for i in range(len(rotId)):
    #if rotId.loc[i].iy==yMin or rotId.loc[i].iy==yMax:
    #    continue
    #
    crystalcontact.loc[cnt]=['Model '+str(float(lenRotId)),data.loc[i]['x'],data.loc[i]['y'],data.loc[i]['z']]
    crystalcontactId.loc[cnt]=['Model '+str(float(lenRotId)),rotId.loc[i]['ix'],rotId.loc[i]['iy'],rotId.loc[i]['iz']]
    cnt+=1
    lenRotId+=1
#
pltCnt=cnt-1
for i in range(len(rotSym)):
    #if zSymmetry[i][1]==yMin or zSymmetry[i][1]==yMax:
    #    continue
    #
    crystalcontact.loc[cnt]=['Model '+str(float(lenRotId)),float(np.round(rotSym[i].item(0),2)),float(np.round(rotSym[i].item(1),2)),float(np.round(rotSym[i].item(2),2))]
    crystalcontactId.loc[cnt]=['Model '+str(float(lenRotId)),float(zSymmetry[i][0]),float(zSymmetry[i][1]),float(zSymmetry[i][2])]
    cnt+=1
    lenRotId+=1
#
with open('CrystalContactsSym.txt','w') as f:
    for i in range(len(crystalcontact)):
        f.write(crystalcontact.loc[i].modId+'\n')
        f.write('        1 0 0 '+str(crystalcontact.loc[i].x)+'\n')
        f.write('        0 1 0 '+str(crystalcontact.loc[i].y)+'\n')
        f.write('        0 0 1 '+str(crystalcontact.loc[i].z)+'\n')    
f.close()
#
with open('lenModelsSym.txt','w') as f:
    f.write(str(lenRotId))
#
# Check for overlap in output
listCrystal=[[crystalcontactId.loc[i].ix,crystalcontactId.loc[i].iy,crystalcontactId.loc[i].iz] for i in range(len(crystalcontactId.ix))]
tmp=[]
for i in range(len(listCrystal)):
    if listCrystal[i] in tmp:
        print('Duplicate at position: i = '+str(i))
    tmp.append(listCrystal[i])
#
plt.rcParams["figure.figsize"] = (20,20)
plt.rcParams["font.size"] = (20)
fig=plt.figure(1)
ax=fig.add_subplot(projection='3d')
#ax.scatter(crystalcontactId.loc[crystalcontactId.iz==0].ix,crystalcontactId.loc[crystalcontactId.iz==0].iy,crystalcontactId.loc[crystalcontactId.iz==0].iz,color='green',label='z=0',marker='o',s=200)
#ax.scatter(rotId.loc[rotId.iz==zMin].ix,rotId.loc[rotId.iz==zMin].iy,rotId.loc[rotId.iz==zMin].iz,color='red',label='z=0',marker='o',s=200)
#ax.scatter(rotId.loc[rotId.iz==zMax].ix,rotId.loc[rotId.iz==zMax].iy,rotId.loc[rotId.iz==zMax].iz,color='blue',label='z=0',marker='o',s=200)
#ax.scatter(rotId.loc[rotId.iz==zMin+1].ix,rotId.loc[rotId.iz==zMin+1].iy,rotId.loc[rotId.iz==zMin+1].iz,color='red',label='z=0',marker='o',s=200)
#ax.scatter(rotId.loc[rotId.iz==zMax-1].ix,rotId.loc[rotId.iz==zMax-1].iy,rotId.loc[rotId.iz==zMax-1].iz,color='blue',label='z=0',marker='o',s=200)
#ax.scatter(rotId.loc[rotId.iz==zMin+2].ix,rotId.loc[rotId.iz==zMin+2].iy,rotId.loc[rotId.iz==zMin+2].iz,color='red',label='z=0',marker='o',s=200)
#ax.scatter(rotId.loc[rotId.iz==zMax-2].ix,rotId.loc[rotId.iz==zMax-2].iy,rotId.loc[rotId.iz==zMax-2].iz,color='blue',label='z=0',marker='o',s=200)
#ax.scatter(rotId.loc[rotId.iz==zMin+3].ix,rotId.loc[rotId.iz==zMin+3].iy,rotId.loc[rotId.iz==zMin+3].iz,color='red',label='z=0',marker='o',s=200)
#ax.scatter(rotId.loc[rotId.iz==zMax-3].ix,rotId.loc[rotId.iz==zMax-3].iy,rotId.loc[rotId.iz==zMax-3].iz,color='blue',label='z=0',marker='o',s=200)
ax.scatter(crystalcontactId.loc[:pltCnt].ix,crystalcontactId.loc[:pltCnt].iy,crystalcontactId.loc[:pltCnt].iz,color='blue',label='Chimera',marker='o',s=200,alpha=0.2)
ax.scatter([i[0] for i in zSymmetry],[i[1] for i in zSymmetry],[i[2] for i in zSymmetry],color='red',label='synthetic',marker='o',s=200)
ax.set_xlabel('Dimension x')
ax.set_ylabel('Dimension y')
ax.set_zlabel('Dimension z')
ax.view_init(15,160)
plt.legend()
plt.savefig('CrystalContacts.png',format='png')
