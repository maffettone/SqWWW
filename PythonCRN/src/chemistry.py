'''
Created on 16 Jun 2016
Last Edited on 28 Oct 2016
@author: Phil Maffettone

TODO:
PDB output that maintains connectivity
'''
from math import * 
import numpy as np

class Atom(object):
    '''
    Class for point atoms which contain relevant data
    coords...Coordinates
    label...Label or Name
    sym...Atom Symbol
    bonded...List of other atom instances which are in bonded to the atom
    '''
    def __init__(self, coords=None, label=None, sym=None, bonded=None):
        self.coords = coords
        self.label = label
        self.sym = sym
        self.bonded=bonded
        if self.coords is None:
            self.coords = np.zeros(3)
        if self.sym is None:
            if self.label is None:
                self.sym = 'H '
                self.label = 'H '
            else:
                self.sym = self.label[0]
                if self.label[1].isalpha():
                    self.sym+=self.label[1]
                else:
                    self.sym+=' '
        if self.label is None:
            self.label = self.sym
        if self.bonded is None:
            self.bonded = []

    def __repr__(self):
        return "%s:\t%s\t[%s]. Bonded to %i neighbors." % (self.__class__.__name__, self.label,
                                  ' '.join(map(lambda x: format(x,'.3f'),self.coords)),len(self.bonded))
        
    def __str__(self):
        return "%-14s%-10s %s" % (self.label,self.sym,
                                  '    '.join(map(lambda x: format(x,'.4f'),self.coords)))


class Crystal(object):
    '''
    Class for organization of atoms in a crystal. Includes atom list and 
    lattice parameter dictionary.
    atoms...Atoms are a separate class
    lat...Latice parameters
    isfract...logical dictating if the lattice parameters are in fractional coordinates
    ac2f... matrix for left multiplication to go from cartesian to fractional coords
    af2c... matrix for left multiplication to go from fractional to cartesian coords
    '''
    def __init__(self,atoms=[],lat={},isfract=True):
        self.atoms = atoms
        self.lat = lat
        self.isfract=isfract
        self.af2c= np.zeros((3,3))
        v = sqrt(1 - cos(lat['alpha'])**2 - cos(lat['beta'])**2 - cos(lat['gamma'])**2 +
             2*cos(lat['alpha'])*cos(lat['beta'])*cos(lat['gamma']))
        self.af2c[0,0] = self.lat['L_a']
        self.af2c[0,1] = self.lat['L_b']*cos(self.lat['gamma'])
        self.af2c[1,1] =self.lat['L_b']*sin(self.lat['gamma'])
        self.af2c[0,2] = self.lat['L_c']*cos(self.lat['beta'])
        self.af2c[1,2] = self.lat['L_c']*(cos(self.lat['alpha'])-cos(self.lat['beta'])*cos(self.lat['gamma']))/sin(self.lat['gamma'])
        self.af2c[2,2] = self.lat['L_c'] * v / sin(self.lat['gamma'])        
        self.ac2f = np.linalg.inv(self.af2c)
        
    def __str__(self):
        return ("Crystal Lattice containing %i atoms.\n"
                "a = %f\nb = %f\nc = %f\nalpha = %f\nbeta = %f\ngamma = %f\n" %
                (len(self.atoms),self.lat['L_a'],self.lat['L_b'],self.lat['L_c'],
                 degrees(self.lat['alpha']),degrees(self.lat['beta']),degrees(self.lat['gamma'])))
    
    def distance(self,atom1,atom2):
        '''
        First Calculates fractional distance, then converts to cartesian
        Works on fractional coordinates
        '''
        
        dist = np.subtract(atom2.coords,atom1.coords)
        if (dist[0] > 0.5):
            dist[0] -= 1.0
        elif (dist[0] < -0.5):
            dist[0] += 1.0
        if (dist[1] > 0.5):
            dist[1] -= 1.0
        elif (dist[1] < -0.5):
            dist[1] += 1.0
        if (dist[2] > 0.5):
            dist[2] -= 1.0
        elif (dist[2] < -0.5):
            dist[2] += 1.0 
        
        if self.isfract:
            dist = np.dot(self.af2c,dist)
        return np.linalg.norm(dist)        

    def configure_bonding(self,sym1,sym2,dmax):
        '''
        Builds pointers to nearest neighbors within dmax
        '''
        atoms = self.atoms
        for idx1,atom1 in enumerate(atoms):
            for idx2,atom2 in enumerate(atoms):
                if atom1.sym == sym1 and atom2.sym == sym2 and idx1 != idx2:
                    d = self.distance(atom1,atom2)
                    if d<dmax:
                        atom1.bonded.append(atom2)
                
        return 

    def frac2cart(self):
        if not self.isfract:
            return
        
        for atom in self.atoms:
            atom.coords = np.dot(self.af2c,atom.coords)            
        self.isfract = False
        
    def cart2frac(self):
        if self.isfract:
            return
        
        for atom in self.atoms:
            atom.coords = np.dot(self.ac2f,atom.coords)
        self.isfract  = True
    
 
class harmonic(object):
    '''
    Constructs the two body harmonic potential for a set of atom types
    from a force constant and equilibrium bond distance.
    '''
    def __init__(self,k,r0):
        self.k = k
        self.r0 = r0
    
    def __str__(self):
        return("Harmonic bond potential of the form 1/2 %.3f (r - %.3)^2" %
               (self.k,self.r0))
        
    def pot(self,cryst,atom1,atom2):
        ''' Computes harmonic potential for single atom pair in 
        fractional coordinates'''
        r = cryst.distance(atom1,atom2)
        V = 0.5*self.k*(r -self.r0)**2
        return V

class three(object):
    '''
    Constructs a three body harmonic potential for a set of atom types
    from a force constant and equilibrium bond angle (degrees).
    The angle is limited to (0,180)
    Works with fractional coordinates and degree inputs
    '''
    def __init__(self,k,theta0):
        self.k = k
        self.theta0 = radians(theta0)
        
    def __str__ (self):
        return("Harmonic bond angle potential of the form 1/2 %.3f (theta - %.3)^2" %
               (self.k,self.theta0))
    def pot(self,center,left,right):
        v1 = np.subtract(left.coords,center.coords)
        if (v1[0] > 0.5):
            v1[0] -= 1.0
        elif (v1[0] < -0.5):
            v1[0] += 1.0
        if (v1[1] > 0.5):
            v1[1] -= 1.0
        elif (v1[1] < -0.5):
            v1[1] += 1.0
        if (v1[2] > 0.5):
            v1[2] -= 1.0
        elif (v1[2] < -0.5):
            v1[2] += 1.0 
        v2 = np.subtract(right.coords,center.coords)
        if (v2[0] > 0.5):
            v2[0] -= 1.0
        elif (v2[0] < -0.5):
            v2[0] += 1.0
        if (v2[1] > 0.5):
            v2[1] -= 1.0
        elif (v2[1] < -0.5):
            v2[1] += 1.0
        if (v2[2] > 0.5):
            v2[2] -= 1.0
        elif (v2[2] < -0.5):
            v2[2] += 1.0 
        theta = np.dot(v1/np.linalg.norm(v1),v2/np.linalg.norm(v2))
        if theta<-1.0:
            theta = -1.0 + (-1.0-theta)
        theta = acos(theta)
        V = 0.5*self.k*(theta - self.theta0)**2
        return V   

class three_cos(object):
    '''
    Constructs a three body harmonic potential for a set of atom types
    from a force constant and equilibrium bond angle (degrees).
    The potential itself relates the cosine of the angles as opposed to
    the angles directly. 
    Works with fractional coordinates
    '''
    def __init__(self,k,theta0):
        self.k = k
        self.cos = cos(radians(theta0))
        
    def __str__ (self):
        return("Harmonic bond angle potential of the form 1/2 %.3f (cos(theta) - %.3)^2" %
               (self.k,self.cos))
    def pot(self,center,left,right):
        v1 = np.subtract(left.coords,center.coords)
        if (v1[0] > 0.5):
            v1[0] -= 1.0
        elif (v1[0] < -0.5):
            v1[0] += 1.0
        if (v1[1] > 0.5):
            v1[1] -= 1.0
        elif (v1[1] < -0.5):
            v1[1] += 1.0
        if (v1[2] > 0.5):
            v1[2] -= 1.0
        elif (v1[2] < -0.5):
            v1[2] += 1.0 
        v2 = np.subtract(right.coords,center.coords)
        if (v2[0] > 0.5):
            v2[0] -= 1.0
        elif (v2[0] < -0.5):
            v2[0] += 1.0
        if (v2[1] > 0.5):
            v2[1] -= 1.0
        elif (v2[1] < -0.5):
            v2[1] += 1.0
        if (v2[2] > 0.5):
            v2[2] -= 1.0
        elif (v2[2] < -0.5):
            v2[2] += 1.0 
        cos_thet = np.dot(v1/np.linalg.norm(v1),v2/np.linalg.norm(v2))
        V = 0.5*self.k*(cos_thet - self.cos)**2
        return V  
            
def read_cif(filename):
    '''
    Returns a dictionary of lattice parameters and a list of atoms_type atom
    Routine for reading .CIF files.
    Structured as a first attempt with i/o in python, using many comments
    Places cif file in new style dictionary using atom class
    '''
    data = {}
    cell_data = {}

    # Open the .CIF file
    with open(filename,'r') as f:
        reading_sym_ops = False
        reading_atom_sites = False
        line_num = 0

        # Read lines one by one
        for line in f:
            line_num += 1

            # Split into columns
            cols = line.split()
            if (len(cols) == 0): continue

            # ID the keywords
            # Keywords associated with cell parameters
            if (cols[0] == '_cell_length_a'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_length_a'] = float(cols[1])
            elif (cols[0] == '_cell_length_b'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_length_b'] = float(cols[1])
            elif (cols[0] == '_cell_length_c'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_length_c'] = float(cols[1])
            elif (cols[0] == '_cell_angle_alpha'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_angle_alpha'] = float(cols[1])
            elif (cols[0] == '_cell_angle_beta'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_angle_beta'] = float(cols[1])
            elif (cols[0] == '_cell_angle_gamma'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_angle_gamma'] = float(cols[1])
            elif (cols[0] == '_cell_volume'):
                cols[1:2] = cols[1].split('(',1)
                data['_cell_volume'] = float(cols[1])

            # Keywords associated with symmetry operations
            elif (cols[0] == '_symmetry_equiv_pos_as_xyz'):
                reading_sym_ops = True
                data['_symmetry_equiv_pos_as_xyz'] = []
            elif (reading_sym_ops):
                # Add the operation if the string is between single quotes.
                # Otherwise it's a sign we are done with the list.
                if (cols[0][0] == '\''  and  cols[0][-1] == '\''):
                    data['_symmetry_equiv_pos_as_xyz'].append(cols[0][1:-1])
                else:
                    reading_sym_ops = False

            # Keywords associated with atom sites
            elif (cols[0] == '_atom_site_label'):
                index_label = line_num
                data['_atom_site_label'] = []
                data['_atom_site_fract_x'] = []
                data['_atom_site_fract_y'] = []
                data['_atom_site_fract_z'] = []
                reading_atom_sites = True

            # Keep track of where the other labels are (order is important).
            elif (cols[0] == '_atom_site_fract_x'):
                index_x = line_num - index_label
            elif (cols[0] == '_atom_site_fract_y'):
                index_y = line_num - index_label
            elif (cols[0] == '_atom_site_fract_z'):
                index_z = line_num - index_label
            elif (cols[0] == '_atom_site_type_symbol'):
                index_sym = line_num - index_label
                data['_atom_site_type_symbol'] = []

            # If we are currently reading the atom sites...
            elif (reading_atom_sites):
                # Read the actual data if we have 4 columns or more of data.
                if (len(cols) >= 4):
                    data['_atom_site_label'].append(cols[0])
                    data['_atom_site_fract_x'].append(float(cols[index_x]))
                    data['_atom_site_fract_y'].append(float(cols[index_y]))
                    data['_atom_site_fract_z'].append(float(cols[index_z]))
                    if (data.has_key('_atom_site_type_symbol')): data['_atom_site_type_symbol'].append(cols[index_sym])

                # Stop reading atom sites if we found a line with fewer
                # columns, and which does not start with '_atom_site_'.
                elif (len(cols[0]) < 11  or  cols[0][:11] != '_atom_site_'):
                    reading_atom_sites = False

    #Creating a cleaner dictionary to return with a complete atom list
    cell_data['L_a'] =  float(data['_cell_length_a'])
    cell_data['L_b'] = float(data['_cell_length_b'])
    cell_data['L_c'] = float(data['_cell_length_c'])
    cell_data['alpha'] = radians(float(data['_cell_angle_alpha']))
    cell_data['beta'] = radians(float(data['_cell_angle_beta']))
    cell_data['gamma'] = radians(float(data['_cell_angle_gamma']))
    if (data.has_key('_cell_volume')) : cell_data['volume'] = float(data['_cell_volume'])

    #Creating atom list
    atoms = []
    sym = None
    label = data['_atom_site_label']
    coords = np.transpose(np.array([data['_atom_site_fract_x'],data['_atom_site_fract_y'],
                       data['_atom_site_fract_z']]))

    if (data.has_key('_atom_site_type_symbol')):
        sym = data['_atom_site_type_symbol']
        for i in xrange(len(data['_atom_site_label'])):
            atoms.append(Atom(coords[i,:],label[i],sym[i]))
    else:
        for i in xrange(len(data['_atom_site_label'])):
            atoms.append(Atom(coords[i,:],label[i]))

    # Fixing with float division and applying symmetry operations
    if (data.has_key('_symmetry_equiv_pos_as_xyz')):
        ops = data['_symmetry_equiv_pos_as_xyz']
        for i in xrange(len(ops)):
            ops[i] = ops[i].replace("/","./")

        # Two atoms are on top of each other if they are less than "eps" away.
        eps = 0.01  # in Angstrom
        imax = len(atoms)
        for i in xrange(imax):
            for op in ops:
                # Text evaluation of symmetry operation
                xn,yn,zn = eval(op)
                # Forcing into the unit cell
                xn = (xn + 10.0) % 1.0
                yn = (yn + 10.0) % 1.0
                zn = (zn + 10.0) % 1.0
                # Adding only unique new atoms
                new_atom = True
                for at in atoms:
                    if (abs(at.coords[0]-xn)<eps and abs(at.coords[1]-yn)<eps and abs(at.coords[2]-zn)<eps):
                        new_atom = False
                if (new_atom):
                    coords = np.array([xn,yn,zn])
                    atoms.append(Atom(coords,atoms[i].label,atoms[i].sym))

    #Sorting the atom list alphabetically
    atoms = sorted(atoms, key=lambda at: at.label)
    
    return Crystal(atoms,cell_data)

def write_cif(filename,crystal):
    '''
    Routine for writing cif file in P1 
    from lattice parameters and atom list.
    Works directly from the returned objects of read_cif()
    '''
    atoms = crystal.atoms
    lat = crystal.lat
    with open(filename,'w') as f:
        f.write('data_cif\n\n')
        f.write('%-32s%7.4f%s' % ('_cell_length_a',lat['L_a'],'(0)\n'))
        f.write('%-32s%7.4f%s' % ('_cell_length_b',lat['L_b'],'(0)\n'))
        f.write('%-32s%7.4f%s' % ('_cell_length_c',lat['L_c'],'(0)\n'))
        f.write('%-32s%7.4f%s' % ('_cell_angle_alpha',degrees(lat['alpha']),'(0)\n'))
        f.write('%-32s%7.4f%s' % ('_cell_angle_beta',degrees(lat['beta']),'(0)\n'))
        f.write('%-32s%7.4f%s' % ('_cell_angle_gamma',degrees(lat['gamma']),'(0)\n\n'))
        f.write('%-35s %s' % ('_symmetry_space_group_name_H-M',"'P1'\n"))
        f.write('%-35s %s' % ('_symmetry_Int_Tables_number','1\n'))
        f.write('%-35s %s' % ('_symmetry_cell_setting','triclinic\n\n'))
        f.write('loop_\n')
        f.write('_atom_site_label\n')
        f.write('_atom_site_type_symbol\n')
        f.write('_atom_site_fract_x\n')
        f.write('_atom_site_fract_y\n')
        f.write('_atom_site_fract_z\n')
        f.write('_atom_site_occupancy\n')
        
        for atom in atoms:
            f.write('%s     %s' % (atom, '1.000\n'))
            
def build_super(cryst,Nx=1,Ny=1,Nz=1):
    '''
    Builds an Nx x Ny x Nz super cell from a 
    list of atoms in fractional coordinates
    '''
    #Generates unique labels
    atomtypes = set([])
    for atom in cryst.atoms:
        atomtypes.add(atom.sym)
    atoms = cryst.atoms
    new_atoms = []
    lat = cryst.lat
    for i in xrange(Nx):
        for a in atoms:
            coords = np.add(a.coords,np.array([i,0.,0.]))
            new_atoms.append(Atom(coords,a.label,a.sym))
    atoms = new_atoms
    new_atoms=[]
    for i in xrange(Ny):
        for a in atoms:
            coords = np.add(a.coords,np.array([0.,i,0.]))
            new_atoms.append(Atom(coords,a.label,a.sym))  
    atoms=new_atoms
    new_atoms=[]   
    for i in xrange(Nz):
        for a in atoms:
            coords = np.add(a.coords,np.array([0.,0.,i]))
            new_atoms.append(Atom(coords,a.label,a.sym))
    for a in new_atoms:
        a.coords = np.divide(a.coords,np.array([Nx,Ny,Nz]))
    lat['L_a']= lat['L_a']*Nx
    lat['L_b'] = lat['L_b']*Ny
    lat['L_c'] = lat['L_c']*Nz
    cryst.atoms = new_atoms
    
    cryst.af2c= np.zeros((3,3))
    v = sqrt(1 - cos(lat['alpha'])**2 - cos(lat['beta'])**2 - cos(lat['gamma'])**2 +
         2*cos(lat['alpha'])*cos(lat['beta'])*cos(lat['gamma']))
    cryst.af2c[0,0] = cryst.lat['L_a']
    cryst.af2c[0,1] = cryst.lat['L_b']*cos(cryst.lat['gamma'])
    cryst.af2c[1,1] =cryst.lat['L_b']*sin(cryst.lat['gamma'])
    cryst.af2c[0,2] = cryst.lat['L_c']*cos(cryst.lat['beta'])
    cryst.af2c[1,2] = cryst.lat['L_c']*(cos(cryst.lat['alpha'])-cos(cryst.lat['beta'])*cos(cryst.lat['gamma']))/sin(cryst.lat['gamma'])
    cryst.af2c[2,2] = cryst.lat['L_c'] * v / sin(cryst.lat['gamma'])        
    cryst.ac2f = np.linalg.inv(cryst.af2c)
    
    #Generates Unique labels
    for sym in atomtypes:
        i=1
        for a in (atom for atom in cryst.atoms if atom.sym == sym):
            a.label = "{:s}{:d}".format(sym,i)
            i+=1
            
    return None

def write_pdb(filename,crystal):
    '''
    Routine for writing PDB file to include bonding
    '''
    crystal.frac2cart()
    atoms = crystal.atoms
    fmtstr= "{:6s}{:5d}{:1s}{:4s}{:1s}{:3s}{:1s}{:1s}{:4d}{:1s}{:3s}{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:10s}{:>2s}{:2s}\n"
    with open(filename,'w') as f:
        f.write('HEADER    Chemistry.py generation\n')
        for idx, atom in enumerate(atoms):
            f.write(fmtstr.format('ATOM',idx+1,'',atom.sym,'','MOL','','H',0,'','',atom.coords[0], 
                                  atom.coords[1],atom.coords[2],1.0,0.0,'',atom.sym,'0'))
        crystal.cart2frac()
        atoms = crystal.atoms
        for idx1,a in enumerate(atoms):
            conect = []
            for b in a.bonded:
                for idx2,c in enumerate(atoms):
                    d = crystal.distance(b,c)
                    if (d<0.001 and all(x < 0.5 for x in abs(np.subtract(a.coords,c.coords)))):
                        conect.append(idx2+1)
            f.write('{:6s}{:5d}'.format('CONECT',idx1+1))
            for idx2 in conect:
                f.write('{:5d}'.format(idx2))
            f.write('\n')               
        f.write('TER')    

def write_connectivity(filename,crystal):
    with open(filename,'w') as f:
        for a in crystal.atoms:
            f.write('{:s} {:s}\n'.format(a.label,[atom.label for atom in a.bonded]))
                
