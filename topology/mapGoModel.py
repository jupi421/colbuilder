#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 10:00:06 2023

@author: matthias
"""
import os
import sys
#
def read_go_model(name):
    out=[]
    with open(name+'.itp','r') as f:
        for l in f:
            if l[0]=='[' or l[0]==';':
                continue
            out.append([k for k in l.split(' ') if k.replace(' ','')!='' ])
    f.close()
    return out
        
def generate_go_model(go_name):
    #
    goModel={'sites':[],'table':[],'excl':[]}
    goModel['sites']=read_go_model(go_name+'_go-sites')
    goModel['table']=read_go_model(go_name+'_go-table')
    goModel['excl']=read_go_model(go_name+'_go-excl')
    return goModel
#
def generate_bb_model(bb_name):
    
def run_go_model(filename):
    goModel=generate_go_model(filename)
    return goModel


if __name__ == '__main__':
    os.chdir('/store/go-martini/')
    filename='col'#str(sys.argv[1])
    run_go_model(filename)
    print('TODO')