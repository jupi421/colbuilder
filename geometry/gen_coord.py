#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:17:32 2023

@author: broszms
"""
import numpy as np
import geometry
import subprocess
#
class ContactSymmetry:
    """
    
    Converts Shift-Transformation matrix and fine-tunes shift-matrix
    
    ---
    
    input:  -f     crystal contacts file
   
    output: -o     symmetrized crystal contacts file      

    ---
        
    """
    def __init__(self,file_name):
        self.file_contact=file_name
        self.contacts=[]
        self.rotation_matrix=np.matrix([[39.97,-7.238,-54.2489],
                                        [0,25.9598,-5.79142],
                                        [0,0,675.701]]
                                       )
        self.t_matrix={
                        'model_id' : [],
                        't_x' : [],
                        't_y' : [],
                        't_z' : []
                        }
        
        self.s_matrix={
                        'model_id' : [],
                        's_x' : [],
                        's_y' : [],
                        's_z' : []
                        }
    #
    def read_contacts(self):
        with open(self.file_contact+'.txt','r') as f:
            self.contacts=f.readlines()
        f.close()
    #
    def set_t_matrix(self):
        for idx in range(0,len(self.contacts),4):
            self.t_matrix['model_id'].append(self.contacts[idx].replace('\n',''))
            self.t_matrix['t_x'].append(
                float(self.contacts[idx+1].split(' ')[-1])
                )
            self.t_matrix['t_y'].append(
                float(self.contacts[idx+2].split(' ')[-1])
                )
            self.t_matrix['t_z'].append(
                float(self.contacts[idx+3].split(' ')[-1])
                )
    #
    def get_s_matrix(self):
        for idx in range(len(self.t_matrix['model_id'])):
            print(idx)
            if self.t_matrix['model_id'][idx] in self.s_matrix['model_id']:
                continue
            #
            solve=np.linalg.solve(self.rotation_matrix,
                                [
                                    self.t_matrix['t_x'][idx],
                                    self.t_matrix['t_y'][idx],
                                    self.t_matrix['t_z'][idx],
                                    ]
                                )
            #
            self.s_matrix['model_id'].append(self.t_matrix['model_id'][idx])
            self.s_matrix['s_x'].append(int(solve[0]))
            self.s_matrix['s_y'].append(int(solve[1]))
            self.s_matrix['s_z'].append(int(solve[2]))
            
    #
    def get_t_matrix(self):
        for idx in range(len(self.s_matrix['model_id'])):
            if self.s_matrix['model_id'][idx] in self.t_matrix['model_id']:
                continue
            #
            solve=np.dot(self.rotation_matrix,
                         [
                             self.s_matrix['s_x'][idx],
                             self.s_matrix['s_y'][idx],
                             self.s_matrix['s_z'][idx],
                             ]
                         )
            #
            self.t_matrix['model_id'].append(self.s_matrix['model_id'][idx])
            self.t_matrix['t_x'].append(int(solve[0]))
            self.t_matrix['t_y'].append(int(solve[1]))
            self.t_matrix['t_z'].append(int(solve[2]))

    
    #
    def symmetrize(self):
        xMax,xMin=np.max(self.s_matrix['s_x']),np.min(self.s_matrix['s_x'])
        yMax,yMin=np.max(self.s_matrix['s_y']),np.min(self.s_matrix['s_y'])
        zMax,zMin=np.max(self.s_matrix['s_z']),np.min(self.s_matrix['s_z'])
        #
        print(xMax,' ',yMax,' ',zMax)
        #
       # for plane in range(zMin,zMax)
    
    #
    def get_plane(self,ip):
        self.plane={'p_x':[],'p_y':[],'p_z':[]}
        for idx in range(len(self.s_matrix['model_id'])):
            if self.s_matrix['z'][idx]==ip:
                self.plane['p_x'].append(self.s_matrix['s_x'][idx])
                self.plane['p_y'].append(self.s_matrix['s_y'][idx])
                self.plane['p_z'].append(self.s_matrix['s_z'][idx])
    #
    def run_system(self):
        # Read contacts and set transformation matrix 
        self.read_contacts()
        #
        self.set_t_matrix()
        # Set symmetrized matrix
        self.get_s_matrix()
        self.symmetrize()
        
        


class FibrilGenerator:
    """
    
    Generate collagen microfibril with the "crystal contacts" command from UCSF 
    Chimera. Chimera needs to be installed beforehand from:
    
    https://www.cgl.ucsf.edu/chimera/download.html
    
    Important: Install Chimera (python 2.7 based) and not ChimeraX (python 3.)
    
    ---
    
    input:  -f     pdb-file with single triple helix
            -nocc  number of crystal contacts
           
    output: -cc    crystal-contacts.txt      

    ---

    matrixget:  Calls Chimera via Python 2.7 script from terminal    
    
                -> Gets transformation matrices based on crystal contacts 
                   (e.g. no_cc = 60) to setup skeleton of microfibril
                
    matrixset:  Call Chimera via Python 2.7 script frpm terminal    
    
                -> set pdb-models based on updated, symmetriyzed transformation
                   matrices.
    
    
    """
    
    def __init__(self,path_geo,path,file_name,contact_distance,cut_off):
        self.path_geo=path_geo
        self.path=path
        self.pdb=file_name
        self.d_contact=contact_distance
        self.contacts='crystal_contacts'
        self.sym_contacts=''
        self.cut_off=cut_off
        self.nu_fibril=0   
    #
    def matrixget(self):
        return subprocess.run(
            'chimera --nogui --silent --script '+'"'+str(self.path_geo)+
            'matrixget.py '+str(self.path)+' '+str(self.pdb)+'.pdb '+
            str(self.d_contact)+' '+str(self.contacts)+'"',shell=True)
    #
    def symmetrize(self):
        self.sym_contacts='crystal_contacts_sym'
        return ContactSymmetry(self.contacts).run_system()
        
    #
    def matrixset(self):
        return subprocess.run(
            'chimera --nogui --silent --script '+'"'+str(self.path_geo)+
            'matrixset.py '+str(self.path)+' '+str(self.pdb)+'.pdb '+
            ' '+str(self.sym_contacts)+' '+str(self.nu_fibril)+' '+
            str(self.cut_off)+'"',shell=True)
    #
    def run_system(self):
        # Get Transformation Matrix
        #self.matrixget()
        # Symmetrize
        self.symmetrize()
        # Set symmetrized matrix
        #self.matrixset()
        

def run_gen_coord(path,file,contact_distance,cut_off):
    # TODO: Better solution with paths
    path_geo=geometry.__file__.replace('__init__.py','')
    # Get Crystal Contacts
    fibril = FibrilGenerator(path_geo,path,file,contact_distance,cut_off).run_system()
    return fibril