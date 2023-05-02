#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 10:00:06 2023

@author: matthias
"""
import os
import sys
#
def read_go_itp(name):
    out=[]
    with open(name+'.itp','r') as f:
        for l in f:
            if l[0]=='[' or l[0]==';':
                continue
            out.append([k for k in l.split(' ') if k!=''  and k!='\t'])
    f.close()
    return out
#
def read_bb_pdb(bb_name):
    out=[]
    with open(bb_name,'r') as f:
        for l in f:
            if len(l)>3:
                tmp=[ i for i in l.split(' ') if i!='']
                out.append([tmp[1],tmp[3]])
    f.close()
    return out
#        
def generate_go_itp(go_name,bb_name):
    go_model={'sites':[],'table':[],'excl':[]}
    go_model['sites']=read_go_itp(go_name+'_go-sites')
    go_model['table']=read_go_itp(go_name+'_go-table')
    go_model['excl']=read_go_itp(go_name+'_go-excl')
    #
    bb_model=read_bb_pdb(bb_name)
    return go_model,bb_model


def run_map_go_itp(path_wd,filename,bbname):
    os.chdir(path_wd)
    #
    # Get Reference (Contact Map) + (BB File) for whole system
    go_ref_Model,bb_ref_model=generate_go_itp(filename,'BB_ref.pdb')
    #
    # Get Current BB File for current system
    bb_cur_model=read_bb_pdb(bbname)
    return go_ref_Model,bb_ref_model,bb_cur_model

if __name__=='main':
    filename='Rat'
    bbname='BB.pdb'
    path='/hits/fast/mbm/broszms/Collagen/colbuilder/tests/'
    goModel,bbrefModel,bbcurModel=run_map_go_itp(path,filename,bbname)