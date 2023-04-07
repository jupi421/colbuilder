#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 16:17:05 2023

@author: broszms
"""
import sys
from analysis import triplehelix
#
def write_top(file_size):
    # Write the final topology of the whole microfibril
    with open('system.top','w') as f:
        f.write('; This is the Topology for the Collagen Fibril\n')
        f.write('#define GO_VIRT\n')
        f.write('#include "martini_v3.0.0.itp"\n')
        f.write('#include "./sites/go-sites.itp"\n')
        f.write('#include "./table/go-table.itp"\n')
        f.write('\n')
        #
        for t in range(file_size):
            f.write('#include "./itps/col_'+str(t+1)+'.itp"\n')
            f.write('#include "./excl/col_'+str(t+1)+'_go-excl.itp"\n')
        #
        f.write('\n')
        f.write('#include "martini_v3.0.0_solvents_v1.itp"\n')
        f.write('#include "martini_v3.0.0_ions_v1.itp"\n')
        #
        f.write('\n')
        f.write('[ system ]\n')
        f.write('Collagen, Martini 3 and Go-Potentials \n')
        f.write('\n')
        #
        f.write('[ molecules ]\n')
        #
        for t in range(file_size):
            f.write('col_'+str(t+1)+'     1\n')
        #
    f.close()
#
def write_go_top(name,size,):
    # Write topology file for go-potentials
    # Here Go-sites contains all virtuals sites and 
    # Go-table contains the vs-vs interactions
    with open('go-'+name,'w') as f:
        for t in range(size):
            f.write('#include "col_'+str(t+1)+'_go-'+str(name)+'"\n')
    f.close()
#
def run_top(th_file):
    #
    triplehelices=triplehelix.triplhelix_connect(th_file)
    #
    write_top(len(triplehelices))
    #
    write_go_top('sites.itp',len(triplehelices))
    write_go_top('table.itp',len(triplehelices))
    
#
if __name__=='main':
    filename=str(sys.argv[1])
    triplehelices=triplehelix.triplhelix_connect(filename)
    #
    write_top(len(triplehelices))
    #
    write_go_top('sites.itp',len(triplehelices))
    write_go_top('table.itp',len(triplehelices))