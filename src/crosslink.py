#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 15:42:38 2023

@author: broszms
"""
import numpy as np
from sklearn.metrics import pairwise_distances as pdist
#
def set_crosslink_topology(file,force_field):
    pdbs,coords=get_crosslink_coords(file)
    connect=get_crosslink_connect(coords,pdbs)
    bonds=set_crosslink_bonds(connect,force_field)
    return bonds
#
def get_crosslink_coords(file):
    # Identify Crosslinks from PDB and from Topology
    crosslinkCoords=[]
    crosslinkPDB=[]
    pdbIter=0
    with open(file,'r') as f:
        for l in f:
            if l[0:4]=='ATOM':
                pdbIter+=1
            # Trivalent Case
            if l[17:20]=='LYX' and l[12:15]=='SC4' or l[17:20]=='LY2' and l[12:15]=='SC1' or l[17:20]=='LY3' and l[12:15]=='SC1' or l[17:20]=='LYX' and l[12:15]=='SC5' :
                crosslinkPDB.append([str(pdbIter),l[17:20],l[12:15],l[21:26],float(l[29:38]),float(l[38:46]),float(l[46:56])])
                crosslinkCoords.append([float(l[29:38]),float(l[38:46]),float(l[46:56])])
            # Divalent Case
            if l[17:20]=='L4Y' and l[12:15]=='SC1' or l[17:20]=='L5Y' and l[12:15]=='SC1':
                crosslinkPDB.append([str(pdbIter),l[17:20],l[12:15],l[21:26],float(l[29:38]),float(l[38:46]),float(l[46:56])])
                crosslinkCoords.append([float(l[29:38]),float(l[38:46]),float(l[46:56])])
    f.close()
    return crosslinkPDB, crosslinkCoords

def get_crosslink_connect(crosslinkCoords,crosslinkPDB):
    # Get nearest crosslinks to determine connection 
    pairDist=pdist(crosslinkCoords)
    out=[]
    crosslinkConnect=[]
    #
    for p in pairDist:
        tmp=[]
        for k in np.argsort(p)[0:4]:
            if k not in out:
                tmp.append(crosslinkPDB[k])
                out.append(k)
        if tmp!=[]:
            crosslinkConnect.append(tmp)
    return crosslinkConnect

def set_crosslink_bonds(crosslinkConnect,force_field):
    #
    klyxly2=force_field[0]
    klyxly3=force_field[1]
    dlyxly2=format(force_field[2],'.3f')
    dlyxly3=format(force_field[3],'.3f')
    dly45=format(force_field[4],'.3f')
    kly45=force_field[5]
    # Setup Crosslink Topology
    crosslinkBonds=[]
    for c in crosslinkConnect:
        for clx in c:
            for cly in c:
                # Trivalent Case
                if clx[1]=='LYX' and cly[1]=='LY2' and clx[2]=='SC4':
                    crosslinkBonds.append([clx[0],cly[0],1,str(dlyxly2),str(klyxly2)])
                elif clx[1]=='LYX' and cly[1]=='LY3' and clx[2]=='SC5':
                    crosslinkBonds.append([clx[0],cly[0],str(1),str(dlyxly3),str(klyxly3)])
                # Divalent Case
                if clx[1]=='L4Y' and clx[2]=='SC1' and cly[1]=='L5Y' and cly[2]=='SC2':
                    crosslinkBonds.append([clx[0],cly[0],1,str(dly45),str(kly45)])
    return crosslinkBonds