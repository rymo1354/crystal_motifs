# -*- coding: utf-8 -*-
import os, sys
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
        
class CiftoStructure:
    # Convert .cif file to pymatgen Structure object
    pass

class AseAtomstoStructure:
    # Convert an ASE atoms object to pymatgen Structure object
    pass

class SymmetrizeStructure(object):
    # Useful class for symmetrizing and reducing site-symmetry broken, large structures prior to graph analysis
    def __init__(self, sym_prec, angle_tolerance):
        self.sym_prec = sym_prec
        self.angle_tolerance= angle_tolerance
        return
    
    def symmetrize_structure(self, structure):
        sga = SpacegroupAnalyzer(structure, self.sym_prec, self.angle_tolerance)
        symmetrized_structure = sga.get_symmetrized_structure()
        primitive_structure = symmetrized_structure.get_primitive_structure()
        return primitive_structure

class CoordinateConverter(object):
    def __init__(self):
        # Used to convert fractional, oblique coordinates to Cartesian, orthogonal coordinates 
        # and vice-versa for 3D plotting, using the a, b, and c lattice vectors.
        # The origins are fixed between the coordinate systems using these conversion matrices
        return
    
    def __repr__(self):
        return "Coordinate Converter"
    
    def __str__(self):
        return self.__repr__()
    
    def angle_between_two_vectors(self, vec1, vec2):
        inner = np.inner(vec1, vec2)
        norms = np.multiply(np.linalg.norm(vec1), np.linalg.norm(vec2))

        cos = np.divide(inner, norms)
        rad = np.arccos(np.clip(cos, -1.0, 1.0))
        return rad
        
    def cartesian_to_fractional_cm(self, cart_coords, a_vec, b_vec, c_vec):
        
        vol = np.dot(a_vec, np.cross(b_vec, c_vec))
        a = np.linalg.norm(a_vec)
        b = np.linalg.norm(b_vec)
        c = np.linalg.norm(c_vec)
        alpha = self.angle_between_two_vectors(b_vec, c_vec) # between b and c vectors
        beta = self.angle_between_two_vectors(a_vec, c_vec) # between a and c vectors
        gamma = self.angle_between_two_vectors(a_vec, b_vec) # between a and b vectors
        
        cm11 = np.divide(1, a)
        cm12 = np.negative(np.divide(np.cos(gamma), np.multiply(a, np.sin(gamma))))
        cm13_top = np.subtract(np.multiply(np.cos(alpha), np.cos(gamma)), np.cos(beta))
        cm13_bottom = np.multiply(vol, np.sin(gamma))
        cm13 = np.multiply(np.multiply(b, c), np.divide(cm13_top, cm13_bottom))
        cm1 = np.array([cm11, cm12, cm13])
        
        cm21 = 0
        cm22 = np.divide(1, np.multiply(b, np.sin(gamma)))
        cm23_top = np.subtract(np.multiply(np.cos(beta), np.cos(gamma)), np.cos(alpha))
        cm23_bottom = np.multiply(vol, np.sin(gamma))
        cm23 = np.multiply(np.multiply(a, c), np.divide(cm23_top, cm23_bottom))
        cm2 = np.array([cm21, cm22, cm23])
        
        cm31 = 0
        cm32 = 0
        cm33 = np.divide(np.multiply(np.multiply(a, b), np.sin(gamma)), vol)
        cm3 = np.array([cm31, cm32, cm33])
        
        conversion_matrix = np.array([cm1, cm2, cm3])
        
        return conversion_matrix
    
    def cartesian_to_fractional(self, cart_coords, a_vec, b_vec, c_vec):
        conversion_matrix = self.cartesian_to_fractional_cm(cart_coords, a_vec, b_vec, c_vec)
        frac_coords = np.dot(conversion_matrix, cart_coords)
        return frac_coords
    
    def fractional_to_cartesian_cm(self, frac_coords, a_vec, b_vec, c_vec):
        
        vol = np.dot(a_vec, np.cross(b_vec, c_vec))
        a = np.linalg.norm(a_vec)
        b = np.linalg.norm(b_vec)
        c = np.linalg.norm(c_vec)
        alpha = self.angle_between_two_vectors(b_vec, c_vec)
        beta = self.angle_between_two_vectors(a_vec, c_vec)
        gamma = self.angle_between_two_vectors(a_vec, b_vec)
        
        cm11 = a
        cm12 = np.multiply(b, np.cos(gamma))
        cm13 = np.multiply(c, np.cos(beta))
        cm1 = np.array([cm11, cm12, cm13])
        
        cm21 = 0
        cm22 = np.multiply(b, np.sin(gamma))
        cm23_top = np.subtract(np.cos(alpha), np.multiply(np.cos(beta), np.cos(gamma)))
        cm23_bottom = np.sin(gamma)
        cm23 = np.multiply(c, np.divide(cm23_top, cm23_bottom))
        cm2 = np.array([cm21, cm22, cm23])
        
        cm31 = 0
        cm32 = 0
        cm33 = np.divide(vol, np.multiply(np.multiply(a, b), np.sin(gamma)))
        cm3 = np.array([cm31, cm32, cm33])
        
        conversion_matrix = np.array([cm1, cm2, cm3])
                                          
        return conversion_matrix
    
    def fractional_to_cartesian(self, frac_coords, a_vec, b_vec, c_vec):
        conversion_matrix = self.fractional_to_cartesian_cm(frac_coords, a_vec, b_vec, c_vec)
        cart_coords = np.dot(conversion_matrix, frac_coords)
        return cart_coords

class PolyhedronTypes:
    # Rules for polyhedron assignments and checking their identities
    # To be loaded into Polyhedron3D geometry class
    # Number of faces, shape(s) of faces, connectivities(?) of faces
    pass
