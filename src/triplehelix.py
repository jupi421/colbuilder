#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 15:44:40 2023

@author: broszms
"""
import os
#
def triplhelix_connect(file):
    out=[]
    with open(file,'r') as f:
        for l in f:
            out.append([int(i.split('.')[0].replace('\n','')) for i in l.split(' ')])
    f.close()
    return out

