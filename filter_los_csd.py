# version 0.04

################################### imports ###################################

import re, sys, os, io
import argparse

from itertools import product
from fractions import Fraction
from copy import deepcopy

import numpy as np
import pandas as pd

import CifFile # https://anaconda.org/conda-forge/pycifrw


############################## atomic properties ##############################

# main atomic symbols
_AS = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S',
       'Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga',
       'Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
       'Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm',
       'Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
       'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',
       'U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg',
       'Bh','Hs','Mt','Ds']

# atomic number - atomic symbol mapping
_AN2AS = {i+1: symbol for i, symbol in enumerate(_AS)}

# atomic symbol - atomic number mapping
_AS2AN = {symbol: i+1 for i, symbol in enumerate(_AS)}

# add isotopes and synonims to AS and AS2AN
_AS += ['D', 'T', 'Lw']
_AS2AN['D'] = 1
_AS2AN['T'] = 1
_AS2AN['Lw'] = 103

# atomic masses were extracted from CCDC "Elemental Data and Radii" 23.02.2019
# https://www.ccdc.cam.ac.uk/support-and-resources/ccdcresources/Elemental_Radii.xlsx
_AM = [1.008,4.003,6.941,9.012,10.811,12.011,14.007,15.999,18.998,20.18,22.991,24.305,
       26.982,28.086,30.974,32.066,35.453,39.948,39.098,40.078,44.956,47.867,50.942,
       51.996,54.938,55.845,58.933,58.693,63.546,65.39,69.723,72.61,74.922,78.96,
       79.904,83.8,85.468,87.62,88.906,91.224,92.906,95.94,98,101.07,102.906,106.42,
       107.868,112.411,114.818,118.71,121.76,127.6,126.904,131.29,132.905,137.327,
       138.906,140.116,140.908,144.24,145,150.36,151.964,157.25,158.925,162.5,164.93,
       167.26,168.934,173.04,174.967,178.49,180.948,183.84,186.207,190.23,192.217,
       195.078,196.967,200.59,204.383,207.2,208.98,210,210,222,223,226,227,232.038,
       231.036,238.029,237,244,243,247,247,251,252,257,258,259,262,261,262,266,264,
       269,268,271]
_AM = {i+1: m for i, m in enumerate(_AM)}

# covalent radii were extracted from CCDC "Elemental Data and Radii" 23.02.2019
# https://www.ccdc.cam.ac.uk/support-and-resources/ccdcresources/Elemental_Radii.xlsx
_RCOV = [0.23,1.5,1.28,0.96,0.83,0.68,0.68,0.68,0.64,1.5,1.66,1.41,1.21,1.2,1.05,1.02,
         0.99,1.51,2.03,1.76,1.7,1.6,1.53,1.39,1.61,1.52,1.26,1.24,1.32,1.22,1.22,1.17,
         1.21,1.22,1.21,1.5,2.2,1.95,1.9,1.75,1.64,1.54,1.47,1.46,1.42,1.39,1.45,1.54,
         1.42,1.39,1.39,1.47,1.4,1.5,2.44,2.15,2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,
         1.94,1.92,1.92,1.89,1.9,1.87,1.87,1.75,1.7,1.62,1.51,1.44,1.41,1.36,1.36,1.32,
         1.45,1.46,1.48,1.4,1.21,1.5,2.6,2.21,2.15,2.06,2,1.96,1.9,1.87,1.8,1.69,1.54,
         1.83,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5]
_RCOV = {i+1: rcov for i, rcov in enumerate(_RCOV)}

# CSD vdW radii (Bondi + 1.09 for H (Taylor) + 2.00 for unknown)
# values were extracted from CCDC "Elemental Data and Radii" 23.02.2019
# https://www.ccdc.cam.ac.uk/support-and-resources/ccdcresources/Elemental_Radii.xlsx
_RCSD = [1.09,1.4,1.82,2,2,1.7,1.55,1.52,1.47,1.54,2.27,1.73,2,2.1,1.8,1.8,1.75,
         1.88,2.75,2,2,2,2,2,2,2,2,1.63,1.4,1.39,1.87,2,1.85,1.9,1.85,2.02,2,2,
         2,2,2,2,2,2,2,1.63,1.72,1.58,1.93,2.17,2,2.06,1.98,2.16,2,2,2,2,2,2,2,
         2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1.72,1.66,1.55,1.96,2.02,2,2,2,2,2,2,2,
         2,2,1.86,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
_RCSD = {i+1: rcsd for i, rcsd in enumerate(_RCSD)}

# Bondi vdW
_RBND = deepcopy(_RCSD)
_RBND[1] = 1.2

# Alvarez vdW
_RALV = {1:1.2,2:1.43,3:2.12,4:1.98,5:1.91,6:1.77,7:1.66,8:1.5,9:1.46,10:1.58,
         11:2.5,12:2.51,13:2.25,14:2.19,15:1.9,16:1.89,17:1.82,18:1.83,19:2.73,
         20:2.62,21:2.58,22:2.46,23:2.42,24:2.45,25:2.45,26:2.44,27:2.4,28:2.4,
         29:2.38,30:2.39,31:2.32,32:2.29,33:1.88,34:1.82,35:1.86,36:2.25,37:3.21,
         38:2.84,39:2.75,40:2.52,41:2.56,42:2.45,43:2.44,44:2.46,45:2.44,46:2.15,
         47:2.53,48:2.49,49:2.43,50:2.42,51:2.47,52:1.99,53:2.04,54:2.06,55:3.48,
         56:3.03,57:2.98,58:2.88,59:2.92,60:2.95,62:2.9,63:2.87,64:2.83,65:2.79,
         66:2.87,67:2.81,68:2.83,69:2.79,70:2.8,71:2.74,72:2.63,73:2.53,74:2.57,
         75:2.49,76:2.48,77:2.41,78:2.29,79:2.32,80:2.45,81:2.47,82:2.6,83:2.54,
         89:2.8,90:2.93,91:2.88,92:2.71,93:2.82,94:2.81,95:2.83,96:3.05,97:3.4,
         98:3.05,99:2.7}
# set unknown radii as for CSD type
for an, rvdw in _RCSD.items():
    if an not in _RALV:
        _RALV[an] = rvdw

# Rowland & Taylor
_RR_T = {1:1.09,6:1.75,7:1.61,8:1.56,9:1.44,16:1.79,17:1.74,35:1.85,53:2.0}
# set unknown radii as for CSD type
for an, rvdw in _RCSD.items():
    if an not in _RR_T:
        _RR_T[an] = rvdw

# Chernyshov & Ananyev & Pidko
_RChAP = {1:1.21,5:1.91,6:1.87,7:1.76,8:1.74,9:1.55,15:2.17,16:1.95,17:1.91,
          33:2.09,34:2.03,35:2.00,53:2.17}
# set unknown radii as for CSD type
for an, rvdw in _RCSD.items():
    if an not in _RChAP:
        _RChAP[an] = rvdw

# list of vdW radii
_RVDW = {'csd'  : _RCSD,
         'bondi': _RBND,
         'rt'   : _RR_T,
         'alv'  : _RALV,
         'chap' : _RChAP}


def GetAtomicSymbol(atom):
    '''
    If input atom is atomic number, transforms it into atomic symbol
    Input:
        atom: atomic symbol or atomic number
    Output:
        atom: atomic symbol
    '''
    global _AN2AS, _AS2AN
    
    # check input and return atomic symbol or raise error
    if atom in _AS2AN:
        return atom
    elif atom in _AN2AS:
        return _AN2AS[atom]
    else:
        raise ValueError('Input value is not atomic symbol or atomic number: {0}'.format(atom))


def GetAtomicNumber(atom):
    '''
    If input atom is atomic symbol, transforms it into atomic number
    Input:
        atom: atomic symbol or atomic number
    Output:
        atom: atomic number
    '''
    global _AN2AS, _AS2AN
    
    # check input and return atomic symbol or raise error
    if atom in _AN2AS:
        return int(atom)
    elif atom in _AS2AN:
        return _AS2AN[atom]
    else:
        raise ValueError('Input value is not atomic symbol or atomic number: {0}'.format(atom))


def rCov(atom):
    '''
    Returns covalent radius of the input atom
    Input:
        atom: atomic symbol or atomic number
    Output:
        R: covalent radius in Angstroems
    '''
    global _RCOV
    
    return _RCOV[GetAtomicNumber(atom)]


def rvdW(atom, r_type = 'csd'):
    '''
    Returns covalent radius of the input atom
    Input:
        atom: atomic symbol or atomic number
        r_type: type of vdW radii:
            * 'csd' for CSD radii
            * 'bondi' for Bondi radii
            * 'alvarez' for Alvarez radii
    Output:
        R: covalent radius in Angstroems
    '''
    global _RVDW
    
    # check r_type
    if r_type not in _RVDW:
        raise ValueError('Unknown vdW radii type: {0}'.format(r_type))
    
    return _RVDW[r_type][GetAtomicNumber(atom)]



############################## Support functions ##############################

def _str2float(line):
    '''
    Removes precision given in parenthesis, e.g. '0.2532(12)' => '0.2532', and
    transforms obtained string to float
    '''
    return float(re.sub(r'\(.*\)', r'', line))


def _label2atom(line):
    '''
    Transforms atomic label to atomic symbol
    '''
    return re.sub(r'([a-zA-Z]+)\d+.*', r'\1', line)


def _pretify_label(line):
    '''
    Removes addends to atomic label, e.g. 'C12A' => 'C12'
    '''
    return re.sub(r'([a-zA-Z]+\d+).*', r'\1', line)


def _fraction2float(line):
    '''
    Wrapper for numbers like '3 1/2'
    '''
    return float(sum([Fraction(s) for s in line.split()]))


def _syms2matrix(line):
    '''
    Transforms symbolic symmetry operation to numpy matrix, e.g.:
                     | 1  1  0    0 |
    x+y,x-y,z+1/2 => | 1 -1  0    0 |
                     | 0  0  1  1/2 |
    '''
    sym = np.zeros( (3,4) )
    for idx, expr in enumerate(line.split(',')):
        if expr[0] not in ['+','-']:
            expr = '+' + expr
        # split expression to summons: x-y+1/2 => '+x', '-y', '+1/2'
        starts = []
        for i, letter in enumerate(expr):
            if letter in ['+','-']:
                starts.append(i)
        ends = starts[1:] + [len(expr)]
        parts = [expr[s:e] for s, e in zip(starts, ends)]
        for part in parts:
            if 'x' in part or 'y' in part or 'z' in part:
                # '+x' or '-x'
                if len(part) == 2:
                    if part[0] == '+':
                        val = 1
                    else:
                        val = -1
                # '+1/2x', etc
                else:
                    val = _fraction2float(part[:-1])
                # insert value to the matrix
                if part[-1] == 'x':
                    sym[idx,0] = val
                elif part[-1] == 'y':
                    sym[idx,1] = val
                elif part[-1] == 'z':
                    sym[idx,2] = val
                else:
                    raise ValueError('Wrong symmetry operation format')
            else:
                sym[idx,3] = _fraction2float(part)
    return sym



############################### Data functions ################################

def read_multiple_CIF(path):
    '''
    Reads multiple CIF and returns dictionary {Refcode: CIF lines}
    '''
    # use io because PyCifRW can not deal with full Windows paths like 'D:/wd/x.cif'
    with open(path, 'r') as inpf:
        text = [_.strip() for _ in inpf.readlines()]
    # find starts 
    starts = []
    ends = []
    refcodes = []
    for i, line in enumerate(text):
        if line[:13] == 'data_CSD_CIF_':
            starts.append(i)
            refcodes.append(line[13:])
    ends = starts[1:] + [len(text)]
    # split to cifs
    cifs = {r: text[s:e] for r, s, e in zip(refcodes, starts, ends)}
    
    return cifs


def extract_dXH(path):
    '''
    Extracts normalized lenghts of X-H bonds from the file
    '''
    dXH = {}
    # read the file
    with open(args.norm, 'r') as inpf:
        lines = [_.strip() for _ in inpf.readlines()]
        lines = [_ for _ in lines if _]
    # read X-H lengths
    for line in lines:
        parts = line.split()
        if len(parts) != 2:
            raise ValueError
        elem, dist = parts
        dXH[GetAtomicNumber(elem)] = float(dist)
    
    return dXH


def prepare_contacts(csv, LAB1, LAB2, DIST1):
    '''
    Prepares contacts for search inside of 'Crystal' object
    '''
    contacts = {}
    for idx in csv.index:
        refcode = csv.loc[idx,'Refcode']
        if refcode not in contacts:
            contacts[refcode] = []
        labs = sorted([csv.loc[idx,LAB1], csv.loc[idx,LAB2]])
        labs = [_pretify_label(_) for _ in labs]
        dist = csv.loc[idx,DIST1]
        cont = '{0}_{1}_{2:.3f}'.format(labs[0], labs[1], dist)
        contacts[refcode].append((cont, idx))
    
    return contacts


def get_max_dists(contacts):
    '''
    Returns maximal contact distance for cutting clusters: {plab: Rmax}
    '''
    max_dists = {}
    for contact in contacts:
        lab1, lab2, dist = contact.split('_')
        dist = float(dist)
        if lab1 in max_dists:
            if max_dists[lab1] < dist:
                max_dists[lab1] = dist
        else:
            max_dists[lab1] = dist
    
    return max_dists



################################ Math functions ###############################

def _sin(x):
    '''
    Sugar for x in degrees
    '''
    return np.sin(np.deg2rad(x))


def _cos(x):
    '''
    Sugar for x in degrees
    '''
    return np.cos(np.deg2rad(x))


def _coords_transform(data, M):
    '''
    Transforms crystallographic coordinates to Cartesian, or back
    '''
    df = deepcopy(data)
    df.loc[:,['x','y','z']] = np.dot(df.loc[:,['x','y','z']], M)
    return df


def _normalize_cell(cell):
    '''
    Normalizes cell to (0,0,0) unit cell
    '''
    for row, col in product(cell.index, ('x','y','z')):
        cell.loc[row,col] -= int(np.floor(cell.loc[row,col]))


def _find_duplicates(addend, cell, M, delta):
    '''
    Returns indexes of repeating atoms in addend in comparison to cell
    '''
    drop = []
    for i in addend.index:
        for j in cell.index:
            dr = [addend.loc[i,'x']-cell.loc[j,'x'],
                  addend.loc[i,'y']-cell.loc[j,'y'],
                  addend.loc[i,'z']-cell.loc[j,'z']]
            dr = np.array([[dr[0], dr[0]-1, dr[0]+1][np.argmin([abs(dr[0]), abs(dr[0]-1), abs(dr[0]+1)])],
                           [dr[1], dr[1]-1, dr[1]+1][np.argmin([abs(dr[1]), abs(dr[1]-1), abs(dr[1]+1)])],
                           [dr[2], dr[2]-1, dr[2]+1][np.argmin([abs(dr[2]), abs(dr[2]-1), abs(dr[2]+1)])]])
            d = sum(dr.dot(M)**2)**0.5
            if d < delta:
                drop.append(i)
                break
    return drop


def _dist_p2p(p1, p2, x):
    '''
    Returns distance between point x and plane defined by points (0,0,0), p1 and p2
    '''
    n = np.cross(p1, p2)
    n = n / np.linalg.norm(n)
    return abs(sum([a*b for a, b in zip(n, x)]))


def _rotation_matrix(axis, theta):
    '''
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    '''
    axis = np.asarray(axis)
    axis = axis / np.dot(axis, axis)**0.5
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]]).T



################################### Crystal ###################################

class Crystal():
    '''
    Represents molecular crystal obtained from CSD-formated CIF
    '''
    
    def _extract_cif_parameters(self):
        '''
        Extracts and prepares cell parameters, symops, coordinates and residue mapping
        'cif' is PyCifRW object corresponding to an analyzed crystal
        'mol2' is text lines of MOL2 file
        '''
        _, cif = CifFile.ReadCif(io.StringIO('\n'.join(self.cif_text)+'\n')).items()[0]
        self.cif = cif
        # check existance of needed fields
        for field in ['_cell_length_a', '_cell_length_b', '_cell_length_c',
                      '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma',
                      '_symmetry_equiv_pos_as_xyz',
                      '_atom_site_label', '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z',
                      '_geom_bond_atom_site_label_1', '_geom_bond_atom_site_label_2']:
            if field not in cif.block.keys():
                err = '{0} field is absent in {1}'.format(field, self.refcode)
                if field in ('_geom_bond_atom_site_label_1', '_geom_bond_atom_site_label_2'):
                    err += ' Please extract from ConQuest CIF containing bond info'
                raise KeyError(err)
        # cell parameters
        self.a = _str2float(cif.get('_cell_length_a'))
        self.b = _str2float(cif.get('_cell_length_b'))
        self.c = _str2float(cif.get('_cell_length_c'))
        self.alpha = _str2float(cif.get('_cell_angle_alpha'))
        self.beta = _str2float(cif.get('_cell_angle_beta'))
        self.gamma = _str2float(cif.get('_cell_angle_gamma'))
        # extract coordinates
        label = cif.get('_atom_site_label')
        self.N = len(label)
        x = [_str2float(_) for _ in cif.get('_atom_site_fract_x')]
        y = [_str2float(_) for _ in cif.get('_atom_site_fract_y')]
        z = [_str2float(_) for _ in cif.get('_atom_site_fract_z')]
        self.coords_raw = pd.DataFrame({'atom': [GetAtomicNumber(_label2atom(_)) for _ in label],
                                        'label': label,
                                        'plabel': [_pretify_label(_) for _ in label],
                                        'x': x, 'y': y, 'z': z})
        # extract symmetry operations
        self.syms_symbolic = cif.get('_symmetry_equiv_pos_as_xyz')
        self.syms = [_syms2matrix(_) for _ in self.syms_symbolic]
    
    def _generate_C2C_matrixes(self):
        '''
        Generates matrixes for Crystallographic <==> Cartesian coordinates
        transformation
        '''
        M = np.eye(3)
        M[0,0] = self.a
        M[1,0] = self.b*_cos(self.gamma)
        M[2,0] = self.c*_cos(self.beta)
        M[1,1] = self.b*_sin(self.gamma)
        M[2,1] = self.c*(_cos(self.alpha) - _cos(self.beta)*_cos(self.gamma))/_sin(self.gamma)
        M[2,2] = self.c*((1 - _cos(self.alpha)**2 - _cos(self.beta)**2 - _cos(self.gamma)**2 + 2*_cos(self.alpha)*_cos(self.beta)*_cos(self.gamma))**0.5)/_sin(self.gamma)
        M[abs(M)<10**-6] = 0
        self.cryst2cart = M
        self.cart2cryst = np.linalg.inv(M)
    
    def _find_XHs(self):
        '''
        Extracts X-H bonds from CIF bond block
        '''
        Hs = {}
        for lab1, lab2 in zip(self.cif.get('_geom_bond_atom_site_label_1'), self.cif.get('_geom_bond_atom_site_label_2')):
            if GetAtomicNumber(_label2atom(lab1)) == 1:
                if lab1 in Hs:
                    Hs[lab1] += [lab2]
                else:
                    Hs[lab1] = [lab2]
            elif GetAtomicNumber(_label2atom(lab2)) == 1:
                if lab2 in Hs:
                    Hs[lab2] += [lab1]
                else:
                    Hs[lab2] = [lab1]
        for lab, Xs in list(Hs.items()):
            if len(Xs) != 1:
                Hs.pop(lab)
            else:
                Hs[lab] = Xs[0]
        self.XHs = Hs
    
    def _find_dist(self, iX, iY):
        '''
        Returns distance between X and Y in coords (not cell)
        '''
        v = [self.xyz.loc[iY,axis]-self.xyz.loc[iX,axis] for axis in 'xyz']
        return sum([_**2 for _ in v])**0.5
    
    def _enlarge_dist(self, X, Y, D):
        '''
        Enlarges X-Y interatomic distance by moving Y if necessary
        '''
        iX = self.coords.loc[self.coords['label'] == X].index[0]
        iY = self.coords.loc[self.coords['label'] == Y].index[0]
        d = self._find_dist(iX, iY)
        x = self.coords.loc[iX,'x'] + D/d*(self.coords.loc[iY,'x']-self.coords.loc[iX,'x'])
        y = self.coords.loc[iX,'y'] + D/d*(self.coords.loc[iY,'y']-self.coords.loc[iX,'y'])
        z = self.coords.loc[iX,'z'] + D/d*(self.coords.loc[iY,'z']-self.coords.loc[iX,'z'])
        self.coords.loc[iY,'x'] = x
        self.coords.loc[iY,'y'] = y
        self.coords.loc[iY,'z'] = z
    
    def _normalize_XHs(self):
        '''
        Normalizes XH protons to crystallographic values
        '''
        self.coords = deepcopy(self.coords_raw)
        if not self.dXH:
            return None
        self.xyz = _coords_transform(self.coords, self.cryst2cart)
        for H, X in self.XHs.items():
            an = GetAtomicNumber(_label2atom(X))
            if an in self.dXH:
                self._enlarge_dist(X, H, self.dXH[an])
    
    def _generate_unit_cell(self):
        '''
        Applies symops to atomic coordinates, drops duplicated atoms, then
        normalizes to (0,0,0) cell
        '''
        cell = deepcopy(self.coords)
        cell.drop(cell.index, inplace = True)
        for idx, sym in enumerate(self.syms):
            addend = deepcopy(self.coords)
            for row, col in product(addend.index, ('x','y','z')):
                i = 'xyz'.index(col)
                addend.loc[row, col] = sym[i,0]*self.coords.loc[row,'x'] + sym[i,1]*self.coords.loc[row,'y'] + sym[i,2]*self.coords.loc[row,'z'] + sym[i,3]
            _normalize_cell(addend)
            drop = _find_duplicates(addend, cell, self.cryst2cart, self.dABmin)
            addend.drop(drop, inplace = True)
            cell = pd.concat([cell, addend], ignore_index = True)
        self.cell = cell
    
    def _find_axis_addends(self):
        '''
        Finds addends on a, b and c axes for search contacts
        '''
        axes = pd.DataFrame({'x': [1,0,0], 'y': [0,1,0], 'z': [0,0,1]})
        axes = _coords_transform(axes, self.cryst2cart)
        a = list(axes.loc[0]); b = list(axes.loc[1]); c = list(axes.loc[2])
        self.axis_addends = [(self.Rmax+self.Radd)/_dist_p2p(b,c,a), (self.Rmax+self.Radd)/_dist_p2p(a,c,b), (self.Rmax+self.Radd)/_dist_p2p(a,b,c)]
    
    def _enlarge_cell(self):
        '''
        Enlarges cell with axis addends
        '''
        self._find_axis_addends()
        cell = deepcopy(self.cell)
        # enlarge
        for i, axis in enumerate('xyz'):
            cell0 = deepcopy(cell)
            bot = int( 0 - np.ceil(self.axis_addends[i]) )
            top = int( 1 + np.ceil(self.axis_addends[i]) )
            for delta in range(bot, top):
                if delta == 0:
                    continue
                _cell = deepcopy(cell0)
                _cell.loc[:,axis] = [_ + delta for _ in _cell[axis]]
                cell = cell.append(_cell, ignore_index = True)
        # cut by precise cutoffs
        for i, axis in enumerate('xyz'):
            cell = cell.loc[(cell[axis] > -self.axis_addends[i]) & (cell[axis] < 1 + self.axis_addends[i])]
        cell.index = range(len(cell))
        # cartesian
        cell = _coords_transform(cell, self.cryst2cart)
        self.big_cell = cell
    
    def _cut_clusters(self):
        '''
        For each unique atom cuts cluster around
        '''
        self._enlarge_cell()
        self.clusters = {}
        for plabel, Rmax in self.Rmaxes.items():
            idx = self.big_cell.loc[self.big_cell['plabel'] == plabel].index[0]
            ori = [self.big_cell.loc[idx,axis] for axis in 'xyz']
            cluster = self.big_cell.loc[(self.big_cell['x'] > ori[0]-(Rmax+self.Radd)) & (self.big_cell['x'] < ori[0]+(Rmax+self.Radd)) & (self.big_cell['y'] > ori[1]-(Rmax+self.Radd)) & (self.big_cell['y'] < ori[1]+(Rmax+self.Radd)) & (self.big_cell['z'] > ori[2]-(Rmax+self.Radd)) & (self.big_cell['z'] < ori[2]+(Rmax+self.Radd))]
            cluster = cluster.reindex([idx] + [_ for _ in cluster.index if _ != idx])
            cluster.index = range(len(cluster))
            # normalize to ori
            for i, axis in enumerate('xyz'):
                cluster.loc[:,axis] = [_ - ori[i] for _ in cluster[axis]]
            # find distances
            ds = []
            for i in cluster.index:
                ds.append(sum([cluster.loc[i,axis]**2 for axis in 'xyz'])**0.5)
            cluster.loc[:,'D'] = ds
            cluster = cluster.sort_values(by = 'D')
            # drop atoms beyound 5A
            cluster = cluster.loc[cluster['D'] <= (Rmax+self.Radd)]
            cluster.index = range(len(cluster))
            # generate contact names
            contacts = []
            for i in cluster.index:
                labs = sorted([cluster.loc[cluster.index[0],'plabel'], cluster.loc[i,'plabel']])
                cont = '{0}_{1}_{2:.3f}'.format(labs[0], labs[1], cluster.loc[i,'D'])
                contacts.append(cont)
            contacts[0] = ''
            cluster.loc[:,'contact'] = contacts
            self.clusters[plabel] = cluster
    
    def __init__(self, refcode, cif_text, Rmaxes, Radd = 2,
                 dXH = {6: 1.089, 7: 1.015, 8: 0.993},
                 dABmin = 0.005, r_type = 'csd'):
        '''
        Extracts crystal parameters from CIF and label->residue mapping from MOL2;
        then prepares data for searching line-of-sight intermolecular contacts
        Input:
            - refcode: CSD entry refcode
            - cif_text: list of CIF lines
            - normXH: if True, X-H bonds are normalized
            - dABmin: minimal possible distance between similar atoms
            - r_type: type of vdW radii used
        Output:
            - Crystal object
        '''
        # initial data
        self.error = None
        self.refcode = refcode
        self.cif_text = cif_text
        self.Rmaxes = Rmaxes
        self.Rmax = max(Rmaxes.values())
        self.Radd = Radd
        self.dXH = dXH
        self.dABmin = dABmin
        self.r_type = r_type
        # extract and prepare info
        self._extract_cif_parameters()
        self._generate_C2C_matrixes()
        # normalize Hs if needed
        self._find_XHs()
        self._normalize_XHs()
        # cut clusters around unique atoms
        self._generate_unit_cell()
        self._cut_clusters()
    
    def get_shielding(self, contact):
        '''
        Returns contact shielding
        '''
        lab1, lab2, D = contact.split('_')
        D = float(D)
        X = deepcopy(self.clusters[lab1])
        # find contacting atoms
        idx1 = X.index[0]
        idx2s = list(X.loc[X['contact'] == contact].index)
        if len(idx2s) == 0:
            return None, None
        idx2 = idx2s[0]
        # find rotation matrix
        a0 = np.array([X.loc[idx2,axis] for axis in 'xyz'])
        a0 = a0 / np.dot(a0, a0)**0.5
        a1 = np.array([1.0, 0, 0])
        axis = np.cross(a0, a1)
        if np.dot(axis, axis) != 0:
            axis = axis / np.dot(axis, axis)**0.5
            angle = np.arccos(sum(a0*a1))
            M = _rotation_matrix(axis, angle)
            # reorient Cartesian coordinates so that A-B || Ox
            X = _coords_transform(X, M)
        elif a0[0] < 0:
            X.loc[:,'x'] = [-x for x in X['x']]
        # find points on the atomic surfaces
        ori1 = 0.0
        ori2 = X.loc[idx2,'x']
        # drop atoms that are too far from the A-B line (Ox axis)
        drop = []
        rABmax = max(rvdW(X.loc[idx1,'atom'], self.r_type), rvdW(X.loc[idx2,'atom'], self.r_type))
        for idx in [_ for _ in X.index if _ not in (idx1, idx2)]:
            x = X.loc[idx,'x']
            y = X.loc[idx,'y']
            z = X.loc[idx,'z']
            r = rvdW(X.loc[idx,'atom'], self.r_type) + rABmax + 0.5
            if x < ori1 - r or x > ori2 + r:
                drop.append(idx)
                continue
            d = (y**2 + z**2)**0.5
            if d > r:
                drop.append(idx)
        X = X.drop(drop)
        # collect shielding info
        info = {'val': float('inf'), 'C': None}
        As = (rvdW(X.loc[idx1,'atom'], self.r_type), 0, -rvdW(X.loc[idx1,'atom'], self.r_type))
        Bs = (D - rvdW(X.loc[idx2,'atom'], self.r_type), D, D + rvdW(X.loc[idx2,'atom'], self.r_type))
        for idx in [_ for _ in X.index if _ not in (idx1, idx2)]:
            x = X.loc[idx,'x']
            d = (X.loc[idx,'y']**2 + X.loc[idx,'z']**2)**0.5
            # case of long contact
            if As[0] < Bs[0]:
                if x < As[0]:
                    R0 = (d**2 + (As[0]-x)**2)**0.5
                elif x > Bs[0]:
                    R0 = (d**2 + (x-Bs[0])**2)**0.5
                else:
                    R0 = d
            # case of short contact
            else:
                if x < Bs[0]:
                    R0 = (d**2 + (As[0]-x)**2)**0.5
                elif x > As[0]:
                    R0 = (d**2 + (x-Bs[0])**2)**0.5
                else:
                    delta = max(abs(x-Bs[0]), abs(As[0]-x))
                    R0 = (d**2 + delta**2)**0.5
            # shielding
            dr = R0 - rvdW(X.loc[idx,'atom'], self.r_type)
            if dr < info['val']:
                info['val'] = dr
                info['C'] = X.loc[idx,'plabel']
        
        return info['val'], info['C']



################################## Main code ##################################

if __name__ == '__main__':
    
    # cmd arguments
    parser = argparse.ArgumentParser(description = 'This script classifies contacts found with CCDC ConQuest program as line-of-sight (LoS) or non-LoS. If you use it in your research, please cite generously: DOI ???')
    parser.add_argument('path_csv',
                        help = 'path to .csv containing info about contacts')
    parser.add_argument('path_cif',
                        help = 'path to .cif containing entries from .csv file')
    parser.add_argument('-r', '--radii', default = 'csd',
                        choices = ['csd', 'bondi', 'rt', 'alv', 'chap'],
                        help = 'type of used vdW radii')
    parser.add_argument('--lab1', default = 'LAB1',
                        help = 'name of column containing atomic label of the first atom')
    parser.add_argument('--lab2', default = 'LAB2',
                        help = 'name of column containing atomic label of the second atom')
    parser.add_argument('--dist', default = 'DIST1',
                        help = 'name of column containing contact distance')
    parser.add_argument('--norm', default = 'csd',
                        help = 'type of X-H bonds normalization: "csd"/"no"/path to the file containing normalized X-H bond distances')
    parser.add_argument('--tol', default = 0.005, type = float,
                        help = 'minimal possible distance between different atoms')
    
    # handle cmd arguments
    args = parser.parse_args()
    # paths
    if not os.path.exists(args.path_cif):
        print('Error: bad .CIF filepath: {0}'.format(args.path_cif))
    if not os.path.exists(args.path_csv):
        print('Error: bad .CSV filepath: {0}'.format(args.path_csv))
    # normalization
    if args.norm == 'csd':
        dXH = {6: 1.089, 7: 1.015, 8: 0.993}
    elif args.norm == 'no':
        dXH = {}
    else:
        if not os.path.exists(args.norm):
            print('Error: bad normalization filepath: {0}'.format(args.path_cif))
        try:
            dXH = extract_dXH(args.norm)
        except:
            print('Error: bad normalization file format, see manual for the details')
            sys.exit()
    
    # tolerance
    if args.tol <= 0 or args.tol > 0.5:
        print('Error: tolerance (--tol) must be in (0,0.5] interval')
    
    # read input files
    cifs = read_multiple_CIF(args.path_cif)
    csv = pd.read_csv(args.path_csv)
    csv.columns = [_.strip() for _ in csv.columns]
    
    # check csv
    for col in ('Refcode', args.lab1, args.lab2, args.dist):
        if col not in csv.columns:
            print('Error: there is no {0} column in .csv file'.format(col))
            sys.exit()
    for val in csv.loc[:,args.lab1]:
        if not re.match('([a-zA-Z]+)\d+', val):
            print('Error: column {0} contains bad atomic label: {1}'.format(args.lab1, val))
            sys.exit()
    for val in csv.loc[:,args.lab2]:
        if not re.match('([a-zA-Z]+)\d+', val):
            print('Error: column {0} contains bad atomic label: {1}'.format(args.lab2, val))
            sys.exit()
    try:
        val = csv.loc[0,args.dist]
        val = float(csv.loc[0,args.dist])
    except:
        print('Error: column {0} contains non-numeric value: {1}'.format(args.dist, val))
        sys.exit()
    
    # prepare contacts
    contacts = prepare_contacts(csv, args.lab1, args.lab2, args.dist)
    csv.loc[:,'LOS'] = ['?']*csv.shape[0]
    csv.loc[:,'SHIELDING'] = [np.nan]*csv.shape[0]
    csv.loc[:,'SHIELD_ATOM'] = ['']*csv.shape[0]
    
    # print absent refcodes in file
    no_cif = [_ for _ in contacts if _ not in cifs]
    if no_cif:
        print('These entries are absent in CIF file:')
        for refcode in no_cif:
            print(refcode)
        print('')
    
    # test contacts
    print('Start contact searching:')
    refcodes = [_ for _ in contacts if _ in cifs]
    for refcode in list(cifs.keys()):
        if refcode not in refcodes:
            cifs.pop(refcode)
    N = len(refcodes)
    for ii, refcode in enumerate(refcodes):
        print('{0}: {1} of {2}'.format(refcode, ii+1, N))
        # find maximal dist
        Rmaxes = get_max_dists([contact for contact, idx in contacts[refcode]])
        # generate cell
        try:
            X = Crystal(refcode, cifs[refcode], Rmaxes, 2,
                        dXH, args.tol, args.radii)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            for contact, idx in contacts[refcode]:
                csv.loc[idx,'LOS'] = '!'
            print('{0}: error with cell generation'.format(refcode))
            continue
        # find contacts
        for contact, idx in contacts[refcode]:
            try:
                shielding, C = X.get_shielding(contact)
                if shielding is not None:
                    los = '+' if shielding > 0 else '-'
                    csv.loc[idx,'LOS'] = los
                    csv.loc[idx,'SHIELDING'] = shielding
                    csv.loc[idx,'SHIELD_ATOM'] = C
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                csv.loc[idx,'LOS'] = '!'
                print('{0}: error with contact search: {1}'.format(refcode, contact))
    
    # save csv with LOS info
    path_out = os.path.splitext(args.path_csv)[0] + '_los.csv'
    csv.to_csv(path_out, index = False, float_format = '%.3f')


