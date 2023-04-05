#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:04:00 2023

@author: broszms
"""
import os
import numpy as np
import pandas as pd
#
# P1 means the following rotation matrix: Triclinic unit cell
# Rot Matrix = [ [a, b*cos(gamma), c*cos(beta) ], 
#                [ 0, b*sin(gamma), c * (cos(alpha)-cos(beta)*cos(gamma)) / sin(gamma) ]
#                [0,0,c * sqrt(1-cos^2(beta)-(cos(alpha)-cos(beta)*cos(gamma)) /
#                 sin(gamma)*(cos(alpha)-cos(beta)*cos(gamma)) / sin(gamma)  ] ]
#
def read_contact(file):
    # Reads crystal contacts from chimera output
    out=[]
    with open(file,'r') as f:
        out=[l for l in f]
    f.close()
    return out
#
def read_transform_matrix(contacts):
    # Every 4th element marks a new triple helical model
    out={'id':[],'id_x':[],'id_y':[],'id_z':[]}
    for i in range(0,int(len(contacts)),4):
        out['id'].append(contacts[i].replace('\n',''))
        out['id_x'].append(float(contacts[i+1].split(' ')[-1].replace('\n','')))
        out['id_y'].append(float(contacts[i+2].split(' ')[-1].replace('\n','')))
        out['id_z'].append(float(contacts[i+3].split(' ')[-1].replace('\n','')))
    return out
        
#
def get_shift_matrix(t_matrix):
    rotation_matrix=np.matrix([[39.97,-7.238,-54.2489],[0,25.9598,-5.79142],[0,0,675.701]])
    #
    out={'id':[],'id_x':[],'id_y':[],'id_z':[]}
    for i in range(0,len(t_matrix['id'])):
        b=[t_matrix['id_x'][i],t_matrix['id_y'][i],t_matrix['id_z'][i]]
        s_matrix=np.linalg.solve(rotation_matrix,b)
        out['id'].append(t_matrix['id'][i])
        out['id_x'].append(int(s_matrix[0]))
        out['id_y'].append(int(s_matrix[1]))
        out['id_z'].append(int(s_matrix[2]))
    return out
#
def get_transform_matrix(s_matrix):
    rotation_matrix=np.matrix([[39.97,-7.238,-54.2489],[0,25.9598,-5.79142],[0,0,675.701]])
    #
    out={'id':[],'id_x':[],'id_y':[],'id_z':[]}
    for i in range(0,len(s_matrix['id'])):
        b=[s_matrix['id_x'][i],s_matrix['id_x'][i],s_matrix['id_x'][i]]
        t_matrix=np.dot(rotation_matrix,b)
        out['id'].append(s_matrix['id'][i])
        out['id_x'].append(int(t_matrix[0]))
        out['id_y'].append(int(t_matrix[1]))
        out['id_z'].append(int(t_matrix[2]))
    return out
#
def symmetrize_shift_matrix(s_matrix):
    # Symmetrize Shift-Matrix to add new triple helices
    # Skip first line, since it containts the model Ids
    sl_matrix=[[s_matrix['id_x'][i],s_matrix['id_y'][i],s_matrix['id_z'][i]] for i in range(len(s_matrix['id']))]
    #
    out=[]
    xMax,xMin=np.max(s_matrix['id_x']),np.min(s_matrix['id_x'])
    yMax,yMin=np.max(s_matrix['id_y']),np.min(s_matrix['id_y'])
    zMax,zMin=np.max(s_matrix['id_z']),np.min(s_matrix['id_z'])
    for i in range(len(s_matrix['id'])):
        if s_matrix['id_x'][i]<xMax-2 or s_matrix['id_x'][i]>xMin+2:
            tmp=[s_matrix['id_x'][i],-s_matrix['id_y'][i],s_matrix['id_z'][i]]
            if tmp not in sl_matrix:
                out.append(tmp)
    #
    #     
    if [xMin,yMax-1,zMin] not in zSymmetry and [xMin,yMax-1,zMin] not in crystalRot:
        out.append([xMin,yMax-1,zMin])
    if [xMax,yMax-1,zMax] not in zSymmetry and [xMax,yMax-1,zMax] not in crystalRot :
        out.append([xMax,yMax-1,zMax])
    if [xMin,yMin+1,zMin] not in zSymmetry and [xMin,yMin+1,zMin] not in crystalRot:
        out.append([xMin,yMin+1,zMin])
    if [xMax,yMin+1,zMax] not in zSymmetry and [xMax,yMin+1,zMax] not in crystalRot:
        out.append([xMax,yMin+1,zMax])
                
    # fill_outter_layers
    out.append([yMin-1,yMin,zMin])
    out.append([yMin-1,yMax,zMin])
    out.append([yMax+1,yMin,zMax])
    out.append([yMax+1,yMax,zMax])
    #
    out.append([yMin,0,zMin])
    out.append([yMax,0,zMax])
    out.append([yMin-1,0,zMax])
    out.append([yMax+1,0,zMax])
    #
    out.append([-1,0,zMin+1])
    out.append([xMin,0,zMin+1])
    out.append([1,0,zMax-1])
    out.append([xMax,0,zMax-1])
    # Add Symmetry for y=+/-1
    out.append([xMin-2,-1,zMin])
    out.append([xMin-2,1,zMin])
    out.append([xMin+3,-1,zMin])
    out.append([xMin+3,1,zMin])
    out.append([xMax+2,-1,zMax])
    out.append([xMax+2,1,zMax])
    out.append([xMax-3,-1,zMax])
    out.append([xMax-3,1,zMax])
    #
    out.append([-1,-1,zMin+1])
    out.append([-1,1,zMin+1])
    out.append([1,-1,zMax-1])
    out.append([1,1,zMax-1])
    out.append([xMin,-1,zMin+1])
    out.append([xMin,1,zMin+1])
    out.append([xMax,-1,zMax-1])
    out.append([xMax,1,zMax-1])
    #
    out.append([1,-1,-1])
    out.append([-1,1,1])
    return out           

#def main():
    # Get all crystal contacts copies from chimera
    contacts=read_contact('CrystalContacts.txt')
    # Read transformation matrices
    transform_matrix=read_transform_matrix(contacts)
    # Get shift-matrices
    shift_matrix=get_shift_matrix(transform_matrix)
    # Symmetrize shift-matrices
    new_shift_matrix=symmetrize_shift_matrix(shift_matrix)
    
    
    import matplotlib.pyplot as plt
    plt.close('all')
    fig=plt.figure(1)
    ax=fig.add_subplot(projection='3d')
    ax.scatter(shift_matrix['id_x'],shift_matrix['id_y'],shift_matrix['id_z'],s=80)
    ax.scatter([i[0] for i in new_shift_matrix],[i[1] for i in new_shift_matrix],[i[2] for i in new_shift_matrix],s=80)

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