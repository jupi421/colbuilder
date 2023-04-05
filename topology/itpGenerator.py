#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 14:08:00 2022

@author: broszms
"""
import os
import sys
import numpy as np
import crosslink
from analysis import triplehelix
#
def allocate(size):
    # allocate matrices to be filled by topology
    out=[]
    out=[[] for i in range(size)]
    return out
#
path=os.getcwd()
def run_itp(index,filename):
    j=int(sys.argv[1])
    filename=str(sys.argv[2])
    #
    # Define force field parameters for the crosslinks
    force_field=np.zeros(6)
    force_field[0]=9000   # klyxly2
    force_field[1]=9000   # klyxly3
    force_field[2]=0.290  # dlyxly2
    force_field[3]=0.240  # dlyxly3
    force_field[4]=0.350  # dly45
    force_field[5]=7500   # kly45
    #
    triplehelices=triplehelix.triplhelix_connect(filename)
    #
    if len(triplehelices[j-1])==1:
        sys.exit()
        #
    molecule=allocate(len(triplehelices[j-1]))
    molEnds=allocate(len(triplehelices[j-1]))
    atoms=allocate(len(triplehelices[j-1]))
    posres=allocate(len(triplehelices[j-1]))
    bonds=allocate(len(triplehelices[j-1]))
    table=allocate(len(triplehelices[j-1]))
    sites=allocate(len(triplehelices[j-1]))
    constraints=allocate(len(triplehelices[j-1]))
    angles=allocate(len(triplehelices[j-1]))
    exclusions=allocate(len(triplehelices[j-1]))
    dihedrals=allocate(len(triplehelices[j-1]))
    dihedralsOne=allocate(len(triplehelices[j-1]))
    dihedralsTwo=allocate(len(triplehelices[j-1]))
    virtualSites=allocate(len(triplehelices[j-1]))
    #
    for i in range(len(triplehelices[j-1])):
        #
        bondedType=''
        with open('itps/col_'+str(j)+'.'+str(i+1)+'.itp','r') as f:
            for l in f:
                if l[0]==';':
                    continue
                molecule[i].append(l)
                #
                if l=='[ atoms ]\n':
                    bondedType='atoms'
                elif l=='[ position_restraints ]\n':
                    molEnds[i]=int(molecule[i][len(molecule[i])-3].split(' ')[0])
                    bondedType='posres'
                elif l=='[ bonds ]\n':
                    bondedType='bonds'
                elif l=='[ constraints ]\n':
                    bondedType='constraints'
                elif l=='[ virtual_sitesn ]\n':
                        bondedType='virtualsites'
                elif l=='[ angles ]\n':
                    bondedType='angles'
                elif l=='[ dihedrals ]\n':
                    bondedType='dihedrals'
                #
                if bondedType=='atoms' and l.split(' ')[0]!='[':
                    atoms[i].append([k for k in l.split(' ') if k!=''])
                elif bondedType=='posres' and l.split(' ')[0]!='[' :
                    posres[i].append([k for k in l.split(' ') if k!=''])
                elif bondedType=='bonds' and l.split(' ')[0]!='[':
                    bonds[i].append([k for k in l.split(' ') if k!=''])
                elif bondedType=='constraints' and l.split(' ')[0]!='[' :
                    constraints[i].append([k for k in l.split(' ') if k!=''])
                elif bondedType=='virtualsites' and l.split(' ')[0]!='[':
                    virtualSites[i].append([k.replace('\n','') for k in l.split(' ') if k!=''])
                elif bondedType=='angles' and l.split(' ')[0]!='[':
                    angles[i].append([k for k in l.split(' ') if k!=''])
                elif bondedType=='exclusions' and l.split(' ')[0]!='[':
                    exclusions[i].append([k for k in l.split(' ') if k!=''])
                elif bondedType=='dihedrals' and l.split(' ')[0]!='[':
                    dihedrals[i].append([k for k in l.split(' ') if k!=''])
                if bondedType=='dihedrals' and l.split(' ')[0]!='[' and len([k for k in l.split(' ') if k!=' '])==8:
                    dihedralsOne[i].append([k for k in l.split(' ') if k!=''])
                elif bondedType=='dihedrals' and l.split(' ')[0]!='[' and len([k for k in l.split(' ') if k!=' '])==7:
                    dihedralsTwo[i].append([k for k in l.split(' ') if k!=''])
                #
        f.close()
        #
        with open('excl/col_'+str(j)+'.'+str(i+1)+'_go-excl.itp','r') as f:
            for l in f:
                exclusions[i].append([k for k in l.split(' ') if k!=''])
        f.close()
        #
        with open('sites/col_'+str(j)+'.'+str(i+1)+'_go-sites.itp','r') as f:
            for l in f:
                sites[i].append([k for k in l.split(' ') if k!=''])
        f.close()
        #
        with open('table/col_'+str(j)+'.'+str(i+1)+'_go-table.itp','r') as f:
            for l in f:
                table[i].append([k.replace('\n','') for k in l.split(' ') if k!=''])
        f.close()
    #
    #
    finalAtoms=[]
    finalPosRes=[]
    finalBonds=[]
    finalConstraints=[]
    finalAngles=[]
    finalExclusions=[]
    finalDihedrals=[]
    finalDihedralsOne=[]
    finalDihedralsTwo=[]
    # Manipulate Topology
    finalAtoms=atoms[0][:-1]
    finalPosRes=posres[0][:-2]
    finalBonds=bonds[0][:-2]
    finalConstraints=constraints[0][:-2]
    finalVirtualSites=virtualSites[0]
    finalAngles=angles[0]
    finalDihedrals=dihedrals[0]
    finalDihedralsOne=dihedralsOne[0]
    finalDihedralsTwo=dihedralsTwo[0]
    finalExclusions=exclusions[0]
    finalTable=table[0]
    finalSites=[i for i in sites[0]]
    deltaMerge=0
    #
    for i in range(1,len(triplehelices[j-1])):
        #
        deltaMerge+=molEnds[i-1]
        #
        tmp=[[int(k[0])+deltaMerge,k[1],k[2],k[3],k[4],int(k[5])+deltaMerge,k[6],k[7]] for k in atoms[i][:-1] if len(k)>3]
        for k in tmp:
            finalAtoms.append(k)
        #
        tmp=[[int(k[0])+deltaMerge,k[1],k[2],k[3],k[4]] for k in posres[i][:-1] if len(k)>3]
        for k in tmp:
            finalPosRes.append(k)
        #
        tmp=[[int(k[0])+deltaMerge,int(k[1])+deltaMerge,k[2],k[3],k[4]] for k in bonds[i][:-2] if len(k)>3]
        for k in tmp:
            finalBonds.append(k)
        #
        tmp=[[int(k[0])+deltaMerge,int(k[1])+deltaMerge,int(k[2])+deltaMerge,k[3],k[4],k[5]] for k in angles[i][:-1] if len(k)>3]
        for k in tmp:
            finalAngles.append(k)
        #
        tmp=[[int(k[0])+deltaMerge,int(k[1])+deltaMerge,int(k[2])+deltaMerge,int(k[3])+deltaMerge,k[4],k[5],k[6],k[7]] for k in dihedrals[i][:-1] if len(k)>4]
        for k in tmp:
            finalDihedrals.append(k)
        #
        tmp=[[int(k[0])+deltaMerge,int(k[1])+deltaMerge,int(k[2])+deltaMerge,int(k[3])+deltaMerge,k[4],k[5],k[6],k[7]] for k in dihedralsOne[i][:-1] if len(k)>4]
        for k in tmp:
            finalDihedralsOne.append(k)
        #
        tmp=[[int(k[0])+deltaMerge,int(k[1])+deltaMerge,int(k[2])+deltaMerge,int(k[3])+deltaMerge,k[4],k[5],k[6]] for k in dihedralsTwo[i][:-1] if len(k)>4]
        for k in tmp:
            finalDihedralsTwo.append(k)
        #
        tmp=[[int(k[0])+deltaMerge,int(k[1])+deltaMerge,k[2],k[3]] for k in constraints[i][1:-2] if len(k)>3]
        for k in tmp:
            finalConstraints.append(k)
        #
        tmp=[[str(int(k[0])+deltaMerge),str(int(1)),str(int(k[2])+deltaMerge)] for k in virtualSites[i][:-1]]
        for k in tmp:
            finalVirtualSites.append(k)
        #
        tmp=[[k[0],k[1],k[2],k[3],k[4],k[5],int(k[6])+deltaMerge,int(k[7])+deltaMerge,k[8]] for k in table[i][3:-2]]
        for k in tmp:
            finalTable.append(k)
        #
        tmp=[k for k in sites[i]]
        for k in tmp:
            finalSites.append(k)
        #
        tmp=[[int(k[0])+deltaMerge,int(k[1])+deltaMerge,k[3],int(k[4])+deltaMerge,int(k[5])+deltaMerge] for k in exclusions[i][2:-1]]
        for k in tmp:
            finalExclusions.append(k)
        #
    # 
    # Get Crosslink topology
    crosslinkBonds=crosslink.set_crosslink_topology('./pdbs/cross.'+str(j)+'.pdb',force_field)
    #
    # Write Topology     
    with open('itps/col_'+str(j)+'.itp','w') as f:
        f.write('; Merging of triple helices with regard to triplehelice.txt\n')
        f.write('[ moleculetype ]\n')
        f.write('col_'+str(j)+' 1\n')
        f.write('\n')
        #
        f.write('\n')
        f.write('[ atoms ]\n')
        for a in finalAtoms:
            f.write('{:>7}{:>18}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}'.format(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]))
        #
        f.write('\n')
        f.write('[ position_restraints ]\n')
        for a in finalPosRes:
            if len(a)<5:
                for i in a:
                    f.write(' '+i+' ')
                continue
            f.write('{:>7}{:>7}{:>7}{:>7}{:>7}'.format(a[0],a[1],a[2],a[3],a[4]))
        #
        f.write('#endif\n')
        f.write('\n')
        f.write('[ bonds ]\n')
        for b in finalBonds:
            if len(b)<5:
                f.write(';\n')
                continue            
            f.write('{:>8}{:>8}{:>8}{:>8}{:>10}'.format(b[0],b[1],b[2],b[3],b[4]))
        #
        f.write('; Crosslink bonds \n')
        for cb in crosslinkBonds:
            f.write('{:>8}{:>8}{:>8}{:>8}{:>10}\n'.format(cb[0],cb[1],cb[2],cb[3],cb[4]))
        #
        f.write('\n')
        f.write('[ constraints ]\n')
        for c in finalConstraints:
            if len(c)<4:
                for i in c:
                    f.write(' '+i+' ')
                continue  
            f.write('{:>7}{:>7}{:>7}{:>9}'.format(c[0],c[1],c[2],c[3]))
        #
        f.write('#endif\n')
        f.write('\n')
        f.write('[ virtual_sitesn ]\n')
        for v in finalVirtualSites:
            if len(v)<2:
                continue
            f.write('{:>8}{:>8}{:>8}\n'.format(v[0],v[1],v[2]))
        #
        f.write('\n')
        f.write('[ angles ]\n')
        for ag in finalAngles:
            if len(ag)<5:
                continue
            f.write('{:>7}{:>7}{:>7}{:>9}{:>9}{:>7}'.format(ag[0],ag[1],ag[2],ag[3],ag[4],ag[5]))
        #
        f.write('\n')
        f.write('[ dihedrals ]\n')
        for ag in finalDihedralsOne:
            if len(ag)<6:
                continue
            f.write('{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}'.format(ag[0],ag[1],ag[2],ag[3],ag[4],ag[5],ag[6],ag[7]))
        #
        f.write('\n')
        f.write('[ dihedrals ]\n')
        for ag in finalDihedralsTwo:
            if len(ag)<6:
                continue
            f.write('{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}'.format(ag[0],ag[1],ag[2],ag[3],ag[4],ag[5],ag[6]))
        #
    f.close()
    #
    with open('excl/col_'+str(j)+'_go-excl.itp','w') as f:
        f.write('[ exclusions ]\n')
        for ex in finalExclusions:
            if len(ex)<7:
                continue
            f.write('{:>7}{:>7};{:>7}{:>7}\n'.format(ex[0],ex[1],ex[4],ex[5]))
        #
    f.close()
    #
    # define virtual sites 
    with open('sites/col_'+str(j)+'_go-sites.itp','w') as f:
        f.write('[ atomtypes ]\n')
        for site in finalSites:
            if site[0]==';' or site[0]=='[':
                continue
            f.write('{:>12}{:>5}{:>7}{:>5}{:>5}{:>5}'.format(site[0],site[1],site[2],site[3],site[4],site[5]))
            f.write('\n')
        #
    f.close()
    #
    # write non-bonded interactions between virtual sites
    with open('table/col_'+str(j)+'_go-table.itp','w') as f:
        f.write('[ nonbond_params ]\n')
        for table in finalTable:
            if table[0]==';' or table[0]=='[':
                continue
            f.write('{:>15}{:>15}{:>5}{:>15}{:>15}{:>5}{:>7}{:>7}'.format(table[0],table[1],table[2],table[3],table[4],table[5],table[6],table[7]))
            f.write('\n')
        #
    f.close()

