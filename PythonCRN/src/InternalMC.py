'''
Created on 2 Jul 2016
Modified on 3 Jul 2016
@author: PMM

TODO:
Generalize potentials to include attributes indicating usage.
Double check cart2frac

Speedups: 
Move + Change takes 5e-4 seconds, independent of size
Acceptance criteria takes 1e-5, independent of size
Total E takes 5e-2 for 2x2x2, 1e-1 for 3x3x3,  2e-1 for 4x4x4
The internal is mostly dependent on the number of moves 
'''
import chemistry as chem
import numpy as np
from math import *
import time
import os.path

def cart2frac(coords,cryst):
    '''
    A local conversion for cartesian to fractional coordinates
    to take the move from initial Angstrom range to fracitonal
    '''
    return  np.dot(cryst.ac2f,coords)

def frac2cart(coords,cryst):
    '''
    A local conversion for fractional to cartesian coordinates
    to take the move to initial Angstrom range
    '''
    return  np.dot(cryst.af2c,coords)
def start_E(cryst,param):
    '''
    A function similar to total_E which is accessible to 
    programs calling this module.
    Calculates total energy of SiO2, using bond harmonic and three potentials.
    Potentials are defined in the parameters
    '''
    k_12 = param['k_12']
    k3_O = param['k3_O']
    k3_Si = param['k3_Si']
    r_SiO = param['r_SiO']
    theta_O = param['theta_O']
    theta_Si = param['theta_Si']
    Si_O = chem.harmonic(k_12,r_SiO)
    Si_O_Si = chem.three_cos(k3_O,theta_O)
    O_Si_O = chem.three_cos(k3_Si,theta_Si)
    lat = cryst.lat
    
    V = 0.
    atoms = cryst.atoms
    gen = (atom for atom in atoms if atom.sym == 'Si')
    for a in gen:
        i=0
        j=0
        #Two Body Correlations
        for b in a.bonded:
            V += Si_O.pot(cryst,a,b)
        #Three Body Correlations
        while i< len(a.bonded):
            j=i+1
            while j < len(a.bonded):
                V += O_Si_O.pot(a,a.bonded[i],a.bonded[j])
                j += 1
            b = a.bonded[i]
            j=0
            while j < len(b.bonded):
                if b.bonded[j] is not a:
                    V += Si_O_Si.pot(b,a,b.bonded[j])
                j += 1
            i += 1
    return V
    
def crystal_relax(cryst,param):
    '''
    A complete crystal relaxation for the crystal after a bond swap.
    A second function should be written to read in a dicitonary of
    parameters for the MC relaxation.
    '''
    start = time.time()
    total_moves = param['total_moves']
    kT = param['kT']
    minT = param['minT']
    if (minT == 0): minT = 0.0001
    if (kT < minT) : kT = minT
    dT = -log(minT/kT)/total_moves
    max_move = param['max_move']
    min_move = param['min_move']
    k_12 = param['k_12']
    k3_O = param['k3_O']
    k3_Si = param['k3_Si']
    r_SiO = param['r_SiO']
    theta_O = param['theta_O']
    theta_Si = param['theta_Si']

    i=0
    N_atoms = len(cryst.atoms)
    N_accept = 0

    #Defining potentials
    Si_O = chem.harmonic(k_12,r_SiO)
    Si_O_Si = chem.three_cos(k3_O,theta_O)
    O_Si_O = chem.three_cos(k3_Si,theta_Si)
    lat = cryst.lat
   

    #Checking for output
    try:
        out_root = param['output']
        outputting = True
        i=0
        while True:
            fname = '{:s}{:04d}.dat'.format(out_root,i)
            if ( not os.path.isfile(fname)):
                with open(fname,'w') as f:
                    f.write('{:12s}{:12s}\n'.format('T','E'))
                    f.write('-----------------------------------------\n')
                break
            else:
                i+=1
    except KeyError:
        outputting = False

    def total_E(cryst):
        '''
        Calculates total energy of SiO2, using bond harmonic and three potentials.
        Potentials are defined in the parent function
        '''
        V = 0.
        atoms = cryst.atoms
        gen = (atom for atom in atoms if atom.sym == 'Si')
        for a in gen:
            i=0
            j=0
            #Two Body Correlations
            for b in a.bonded:
                V += Si_O.pot(cryst,a,b)
            #Three Body Correlations
            while i< len(a.bonded):
                j=i+1
                while j < len(a.bonded):
                    V += O_Si_O.pot(a,a.bonded[i],a.bonded[j])
                    j += 1
                b = a.bonded[i]
                j=0
                while j < len(b.bonded):
                    if b.bonded[j] is not a:
                        V += Si_O_Si.pot(b,a,b.bonded[j])
                    j += 1
                i += 1
        return V
    
    def local_change(move, a):
        '''
        Local change in energy associated with atom a being
        moved by length 3 cartesian vector, move.
        Returns the dV or change in energy.
        First the energy before the change including the
        bond harmonic, three, and neighboring harmonic is
        calculated. Then the energy after the change.
        '''
        V1 = 0.
        for b in a.bonded:
            V1 += Si_O.pot(cryst,a,b)
        i=0
        j=0
        if a.sym == 'Si':
            while i< len(a.bonded):
                j=i+1
                while j < len(a.bonded):
                    V1 += O_Si_O.pot(a,a.bonded[i],a.bonded[j])
                    j += 1
                b = a.bonded[i]
                j=0
                while j < len(b.bonded):
                    if b.bonded[j] is not a:
                        V1 += Si_O_Si.pot(b,a,b.bonded[j])
                    j += 1
                i += 1
        elif a.sym == 'O':
            while i< len(a.bonded):
                j=i+1
                while j < len(a.bonded):
                    V1 += Si_O_Si.pot(a,a.bonded[i],a.bonded[j])
                    j += 1
                b = a.bonded[i]
                j=0
                while j < len(b.bonded):
                    if b.bonded[j] is not a:
                        V1 += O_Si_O.pot(b,a,b.bonded[j])
                    j += 1
                i += 1

        V2 = 0.
        a.coords += move
        for b in a.bonded:
            V2 += Si_O.pot(cryst,a,b)
        i=0
        j=0
        if a.sym == 'Si':
            while i< len(a.bonded):
                j=i+1
                while j < len(a.bonded):
                    V2 += O_Si_O.pot(a,a.bonded[i],a.bonded[j])
                    j += 1
                b = a.bonded[i]
                j=0
                while j < len(b.bonded):
                    if b.bonded[j] is not a:
                        V2 += Si_O_Si.pot(b,a,b.bonded[j])
                    j += 1
                i += 1
        elif a.sym == 'O':
            while i< len(a.bonded):
                j=i+1
                while j < len(a.bonded):
                    V2 += Si_O_Si.pot(a,a.bonded[i],a.bonded[j])
                    j += 1
                b = a.bonded[i]
                j=0
                while j < len(b.bonded):
                    if b.bonded[j] is not a:
                        V2 += O_Si_O.pot(b,a,b.bonded[j])
                    j += 1
                i += 1
        return V2-V1

    #Defining initial parameters
    i=0
    while i < total_moves:
        #Producing move
        atom_idx = int(floor(np.random.rand()*N_atoms))
        a = cryst.atoms[atom_idx]
        r = np.random.rand()*(max_move-min_move)+min_move
        move = np.random.rand(3,)
        move = move/np.linalg.norm(move)*r
        move = cart2frac(move,cryst)
        

        #Acceptance probability
        dV = local_change(move, a)
        try:
            w = exp(-1.*(dV)/(kT))
            r = np.random.rand()
            if w>=r:
                N_accept += 1
            else:
                a.coords -= move
        except OverflowError:
            if dV/kT < 0:
                N_accept +=1
            else:
                a.coords -= move
        
        if(outputting and(i%100 == 0 or i==0)):
            with open(fname,'a') as f:
                f.write('{:12.3E}{:12.3E}\n'.format(kT,total_E(cryst)))
  
        kT = kT*(1-dT)
        i +=1
    result = {}
    result['Acceptance']  = float(N_accept)/float(total_moves)
    result['Crystal'] = cryst
    result['Total_E'] = total_E(cryst)
    result['Time'] = time.time() - start
    return result
