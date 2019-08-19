import chemistry as chem
import numpy as np
from math import *
from InternalMC import crystal_relax
import time
import os
import sys


def read_input(filename):
    '''
    Function for reading parameters for both MC into 2 dictionaries
    If expected keys are not in dictionary, defaults are set, and exceptions raised. 
    '''
    with open(filename, 'r') as f:
        WWW_param = {}
        relax_param = {}
        in_WWW = False
        in_internal = False
        for line in f:
            if(line[0:4].lower() == 'exit'):
                in_WWW = False
                in_internal = False
                continue
            if(line[0:3].upper() == 'WWW'):
                in_WWW = True
                continue 
            if(line[0:8].lower() == 'internal'):
                in_internal = True
                continue
            cols = line.split()
            if (len(cols) == 0): continue
            if (in_WWW):
                (key,val) = cols
                if key != 'out_dir':
                    WWW_param[key] = float(val)
                else:
                    WWW_param[key] = val
            elif (in_internal):
                (key,val) = cols
                if key != 'output':
                    relax_param[key] = float(val)
                else:
                    relax_param[key] = val
                
    #Set of defaults needed for program to run
    add_def = False            
    if 'total_moves' not in relax_param.keys():
        relax_param['total_moves'] = 10**5
        add_def = True
    else:
        relax_param['total_moves'] = int(10**(relax_param['total_moves']))          
    if 'kT' not in relax_param.keys():
        relax_param['kT'] = 0.5
        add_def = True
    if 'minT' not in relax_param.keys():
        relax_param['minT'] = 0.05
        add_def = True
    if 'max_move' not in relax_param.keys():
        relax_param['max_move'] = 1.0
        add_def = True
    if 'min_move' not in relax_param.keys():
        relax_param['min_move'] = 0.1
        add_def = True
    if 'k_12' not in relax_param.keys():
        relax_param['k_12'] = 27.0
        add_def = True
    if 'k3_O' not in relax_param.keys():
        relax_param['k3_O'] = 2.00
        add_def = True
    if 'k3_Si' not in relax_param.keys():
        relax_param['k3_Si'] = 4.32
        add_def = True                              
    if 'r_SiO' not in relax_param.keys():
        relax_param['r_SiO'] = 1.61
        add_def = True
    if 'theta_O' not in relax_param.keys():
        relax_param['theta_O'] = 180.
        add_def = True
    if 'theta_Si' not in relax_param.keys():
        relax_param['theta_Si'] = 109.471
        add_def = True
    if 'total_moves' not in WWW_param.keys():
        WWW_param['total_moves'] = 10**5
        add_def = True
    else:
        WWW_param['total_moves'] = int(10**(WWW_param['total_moves']))    
    if 'kT' not in WWW_param.keys():
        WWW_param['kT'] = 0.5
        add_def = True
    if 'minT' not in WWW_param.keys():
        WWW_param['minT'] = 0.05
        add_def = True
    if 'smallest_ring' not in WWW_param.keys():
        WWW_param['smallest_ring'] = 3
        add_def = True
    else:
        WWW_param['smallest_ring'] = int(WWW_param['smallest_ring'])   
    if 'out_dir' not in WWW_param.keys():
        WWW_param['out_dir'] = '../Output_Dump'
        add_def = True
    if add_def:
        print 'Not all required parameters present. Defaults incorperated '
    return (WWW_param,relax_param)
   

WWW_param, relax_param = read_input('../Input/test_internal.txt')

data_dir = WWW_param['out_dir']
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

SiO2 = chem.read_cif('../UniqueStart.cif')
chem.build_super(SiO2, 2,2,2)
SiO2.configure_bonding('Si','O',2.5)
SiO2.configure_bonding('O','Si',2.5)
sys.setrecursionlimit(10000)


cif_fname = '{}{}{}{}'.format('/Users/alggroup/Documents/Development/SquareCRN/test_internal/','progress','final','.cif')
chem.write_cif(cif_fname,SiO2) 
cif_fname = '{}{}{}{}'.format('/Users/alggroup/Documents/Development/SquareCRN/test_internal/','progress','final','.pdb')
chem.write_pdb(cif_fname, SiO2)