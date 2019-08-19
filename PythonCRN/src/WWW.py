'''
Created on 01 July 2016
Last Edited on 17 April 2017
@author: Phil Maffettone

TODO:
Remove test lines
'''
import chemistry as chem
import numpy as np
from math import *
from InternalMC import crystal_relax, start_E
import random
import copy
import os
import sys
from itertools import chain

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
            
if len(sys.argv) < 2: 
    print 'Please include input file as command argument...'
    sys.exit()
else:
    input_file = str(sys.argv[1])         
    WWW_param, relax_param = read_input(input_file)    
    
data_dir = WWW_param['out_dir']
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
  
data_fname = data_dir + 'output.dat'
with open(data_fname,'w') as f:
    f.write('%-12s%-12s%-13s%-12s\n' % ('Bond Swaps','InnerAccept','Temp', 'Energy'))
    f.write('------------------------------------------------------\n')

SiO2 = chem.read_cif('../UniqueStart.cif')
chem.build_super(SiO2, 2,2,2)
SiO2.configure_bonding('Si','O',2.5)
SiO2.configure_bonding('O','Si',2.5)
sys.setrecursionlimit(10000)

#Defining local parameters
total_moves = WWW_param['total_moves']
kT = WWW_param['kT']
minT = WWW_param['minT']
dT = -log(minT/kT)/total_moves
i = 0
N_accept = 0

E1 = start_E(SiO2,relax_param)
random.seed()

with open(data_fname,'a') as f:
    f.write('{:-8d}{:4s}{:-12.6f}{:-12.6f}{:-12.6f}\n'.format(0,'',0,kT,E1))
    
cif_fname = '{}{}'.format(data_dir,'start.cif')
chem.write_cif(cif_fname,SiO2)    
cif_fname = '{}{}'.format(data_dir,'start.pdb')
chem.write_pdb(cif_fname, SiO2)
while i < total_moves:
    reject = False
    MC_Results = {}
    a = random.choice([atom for atom in SiO2.atoms if atom.sym == 'Si'])
    bridge_ox = random.choice(a.bonded)
    b = random.choice([atom for atom in bridge_ox.bonded if atom is not a])
    a_ox = random.choice([atom for atom in a.bonded if atom is not bridge_ox])
    b_ox = random.choice([atom for atom in b.bonded if atom is not bridge_ox])
    a_next = random.choice([atom for atom in a_ox.bonded if atom is not a])
    b_next = random.choice([atom for atom in b_ox.bonded if atom is not b])    
        
    a_ox.bonded.append(b_next)
    b_next.bonded.append(a_ox)
    b_ox.bonded.append(a_next)
    a_next.bonded.append(b_ox)
    a_ox.bonded.remove(a_next)
    a_next.bonded.remove(a_ox)
    b_ox.bonded.remove(b_next)
    b_next.bonded.remove(b_ox)

    #Checking for rings at a
    k=1
    atom_tree = [] #Running list
    mem = [a] #Set of nodes from whence you came
    curr_branch =[] #Set of nodes you are going to, tuples which correspond to mem nodes
    curr_branch.append(tuple(a.bonded))
    atom_tree.extend(list(chain.from_iterable(curr_branch))) 
    while k<WWW_param['smallest_ring']:
        next_branch = []
        next_mem = []
        for j in xrange(len(curr_branch)):
            for node in curr_branch[j]:
                next_branch.append(tuple([atom for atom in node.bonded if atom is not mem[j]]))
                next_mem.append(node)
        curr_branch = next_branch
        mem = next_mem
        atom_tree.extend(list(chain.from_iterable(curr_branch)))

        k+=1
    if a in atom_tree:
        reject = True
    #Checking for rings at b
    k=1
    atom_tree = [] #Running list
    mem = [b] #Set of nodes from whence you came
    curr_branch =[] #Set of nodes you are going to, tuples which correspond to mem nodes
    curr_branch.append(tuple(b.bonded))
    atom_tree.extend(list(chain.from_iterable(curr_branch))) 
    while k<WWW_param['smallest_ring']:
        next_branch = []
        next_mem = []
        for j in xrange(len(curr_branch)):
            for node in curr_branch[j]:
                next_branch.append(tuple([atom for atom in node.bonded if atom is not mem[j]]))
                next_mem.append(node)
        curr_branch = next_branch
        mem = next_mem
        atom_tree.extend(list(chain.from_iterable(curr_branch)))
        k+=1
    if b in atom_tree:
        reject = True
    #Auto rejecting bond swaps which create bad rings
    
    if not reject:
        SiO2_cpy = copy.deepcopy(SiO2)
        MC_Results = crystal_relax(SiO2_cpy, relax_param)
        E2 = MC_Results['Total_E']
        dV = E2-E1
        dV = -1.0 #TEST ACCEPTS ALL MOVES
        try:
            w = exp(-1.*(dV)/(kT))
            r = np.random.rand()
            if w>=r:
                SiO2 = MC_Results['Crystal']
                E1=E2
                N_accept += 1
            else:
                reject = True
        except OverflowError:
            if dV/kT < 0:
                SiO2 = MC_Results['Crystal']
                E1=E2
                N_accept +=1
            else:
                reject = True
    
    
    #Move back if rejected for any reason
    if reject:
        a_ox.bonded.append(a_next)
        b_next.bonded.append(b_ox)
        b_ox.bonded.append(b_next)
        a_next.bonded.append(a_ox)
        a_ox.bonded.remove(b_next)
        a_next.bonded.remove(b_ox)
        b_ox.bonded.remove(a_next)
        b_next.bonded.remove(a_ox)
    
        
    if (i%100 == 0):
        cif_fname = '{}{}{:04d}{}'.format(data_dir,'progress',i/100,'.cif')
        chem.write_cif(cif_fname,SiO2)
        cif_fname = '{}{}{:04d}{}'.format(data_dir,'progress',i/100,'.pdb')
        chem.write_pdb(cif_fname, SiO2)
    i +=1
    
    try:
        with open(data_fname,'a') as f:
            f.write('{:-8d}{:4s}{:-12.6f}{:-12.6f}{:-12.6f}\n'.format(i,'',MC_Results['Acceptance'],kT,E1))
    except KeyError:
        #Case where no MC, so the dictionary is empty
        with open(data_fname, 'a') as f:
            f.write('{:-8d}{:4s}{:s}{:-12.6f}{:-12.6f}\n'.format(i,'','Ring Skip',kT,E1))
   
    kT = kT*(1-dT)

cif_fname = '{}{}'.format(data_dir,'final.cif')
chem.write_cif(cif_fname,SiO2)    
cif_fname = '{}{}'.format(data_dir,'final.pdb')
chem.write_pdb(cif_fname, SiO2)
