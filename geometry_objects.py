# -*- coding: utf-8 -*-
import sys
import numpy as np
from helpers import CoordinateConverter
from pymatgen.core.structure import Structure
from helpers import HiddenPrints

class Point3D(object):    
    def __init__(self, coordinates, parent_structure=None, index=None, image_from_index=np.array([0, 0, 0])):
        # Creates a Point in 3-D space from Cartesian Coordinates in 3D space
        # Can associate with indexed site of a pymatgen.core.structure.Structure obj and/or periodic image
        if len(coordinates) != 3: # different dimensions
            print('coordinates must be in three dimensions')
            sys.exit(1)
            
        self.coordinates = coordinates
        
        if type(parent_structure) is Structure:
            self.parent_structure = parent_structure
            possible_indices = [i for i in range(len(parent_structure))]
            if index in possible_indices:
                self.index = index
                self.specie = self.parent_structure[self.index].specie
            else:
                self.index = None
                self.specie = None
            if len(image_from_index) == 3:
                if type(image_from_index) is np.ndarray:
                    if np.issubdtype(image_from_index.dtype, np.integer) == True:
                       self.image_from_index = image_from_index
                else:
                    self.image_from_index = np.array([0, 0, 0])
            else:
                self.image_from_index = np.array([0, 0, 0])
        else:
            self.parent_structure = None
            self.index = None
            self.specie = None
    
    def __repr__(self):
        return "Point3D %s %s" % (self.specie, self.coordinates)
    
    def __str__(self):
        return self.__repr__()
    
    def is_point(self, other_point, periodic=False):
        same_point = False
        
        if self.parent_structure == other_point.parent_structure:
            if self.index == other_point.index:
                if periodic == False: # Check if same Cartesian point
                    if np.array_equal(self.coordinates, other_point.coordinates):
                        same_point = True
                else: # Check if same Periodic point
                    if self.parent_structure != None and self.index != None: # Must be specified
                        s = self.parent_structure
                        ops = other_point.parent_structure
                        if np.array_equal(s[self.index].frac_coords, ops[other_point.index].frac_coords):
                            same_point = True
                    else:
                        print('Must specify structure, index and image to check periodic equivalence')
                        sys.exit(1)
        return same_point
        
    def point_shift(self, other_point):
        # shift that makes current point equal to another periodic point in Cartesian space
        cc = CoordinateConverter()
        if self.is_point(other_point, periodic=True) == True: # Check that points are same, with periodicity
            matrix = self.parent_structure.lattice.matrix
            frac_coords_1 = cc.cartesian_to_fractional(self.coordinates, matrix[0], matrix[1], matrix[2])
            frac_coords_2 = cc.cartesian_to_fractional(other_point.coordinates, matrix[0], matrix[1], matrix[2])
            return np.rint(np.subtract(frac_coords_2, frac_coords_1))
        else:
            print('Not the same point considering periodic structure')
            return np.array([None, None, None])
        
    def shift_point(self, shift=np.array([0, 0, 0])):
        # Returns point shifted by image provided
        cc = CoordinateConverter()
        matrix = self.parent_structure.lattice.matrix
        frac_coords_old = cc.cartesian_to_fractional(self.coordinates, matrix[0], matrix[1], matrix[2])
        frac_coords_new = np.add(frac_coords_old, shift)
        cartesian_coords_new = cc.fractional_to_cartesian(frac_coords_new, matrix[0], matrix[1], matrix[2])
        shifted_point = Point3D(cartesian_coords_new, self.parent_structure, self.index, 
                                np.add(self.image_from_index, shift))
        return shifted_point
    
class Edge3D(object):
    def __init__(self, point1, point2):
        # Creates an Edge defined by two 3-D Points
        self.point1 = point1
        self.point2 = point2
        self.points = np.array([self.point1, self.point2])
        
    def __repr__(self):
        return "Edge3D %s" % (self.points)
    
    def __str__(self):
        return self.__repr__()
    
    def __distance__(self):
        # Length of the edge in Cartesian coordinates supplied
        v = np.subtract(np.array(self.point1.coordinates), np.array(self.point2.coordinates))
        return np.linalg.norm(v)
    
    def order_edge(self, other_edge, periodic=False):
        arranged_edges = [Edge3D(other_edge.point1, other_edge.point2), Edge3D(other_edge.point2, other_edge.point1)]
        for arranged_edge in arranged_edges:
            checked_points = [False for i in range(len(self.points))]
            checked_shifts = [None for i in range(len(self.points))]
            for i, point in enumerate(self.points):
                if self.points[i].is_point(arranged_edge.points[i], periodic) == True:
                    checked_points[i] = True
                    checked_shifts[i] = self.points[i].point_shift(arranged_edge.points[i])
            if np.array_equal(checked_points, [True for i in range(len(self.points))]): # Same points
                if np.all(checked_shifts == checked_shifts[0]): # Same periodic shifts for all points
                    return arranged_edge
        print('Not the same edge')
        return None
        
    def edge_shift(self, other_edge):
        ordered_edge = self.order_edge(other_edge, periodic=True)
        if ordered_edge != None:
            pt1_shift = self.point1.point_shift(ordered_edge.point1)
            pt2_shift = self.point2.point_shift(ordered_edge.point2)
            null_shift = np.array([None, None, None])
            if np.array_equal(pt1_shift, pt2_shift) and not np.array_equal(pt1_shift, null_shift):
                return pt1_shift
        else:
            print('Not the same edge with periodic=True')
            return np.array([None, None, None])
        
    def is_edge(self, other_edge, periodic=False):
        same_edge = False
        ordered_edge = self.order_edge(other_edge, periodic)
        if ordered_edge != None:
            same_edge = True
        return same_edge
    
    def shift_edge(self, shift=np.array([0, 0, 0])):
        shifted_points = []
        for point in self.points:
            shifted_point = point.shift_point(shift)
            shifted_points.append(shifted_point)
        shifted_edge = Edge3D(shifted_points[0], shifted_points[1])
        return shifted_edge

class Face3D(object):
    def __init__(self, ordered_points):
        # Creates a Face defined by an array of 3-D Points
        # Points assumed to be connected from left to right, with beginning and end of array connecting
        # to add: area calculator
        self.vertices = ordered_points
        self.edges = [Edge3D(self.vertices[i], self.vertices[i+1]) for i in range(len(self.vertices)-1)] + [Edge3D(self.vertices[-1], self.vertices[0])]
        self.shape = self.__shape__()
        
    def __repr__(self):
        return "Face3D %s" % (self.vertices)
    
    def __str__(self):
        return self.__repr__()
    
    def __shape__(self):
        if len(self.vertices) == 3:
            return 'triangle'
        elif len(self.vertices) == 4:
            return 'quadrilateral'
        elif len(self.vertices) == 5:
            return 'pentagon'
        elif len(self.vertices) == 6:
            return 'hexagon'
        elif len(self.vertices) == 7:
            return 'heptagon'
        elif len(self.vertices) == 8:
            return 'octagon'
        else:
            return 'N-polygon'
    
    def __vector_normal__(self):
        # Vector normal to the Face's plane in Cartesian coordinates
        vector1 = np.subtract(self.vertices[1].coordinates, self.vertices[0].coordinates)
        vector2 = np.subtract(self.vertices[2].coordinates, self.vertices[1].coordinates)
        normal_vector = np.cross(vector1, vector2)
        return normal_vector
    
    def face_arrangements(self, face):
        arranged_faces = [] # All possible arrangements of vertices preserving connectivity
        arr1 = list(face.vertices) # starting point, i.e., [1, 2, 3, 4]
        arr2 = arr1[::-1] # reverse order, i.e., [4, 3, 2, 1]
        for i in range(len(arr1)):
            arr1 = [arr1[-1]] + arr1[0:-1] # now [4, 1, 2, 3]
            arranged_faces.append(Face3D(np.array(arr1)))
            arr2 = [arr2[-1]] + arr2[0:-1] # now [1, 4, 3, 2]
            arranged_faces.append(Face3D(np.array(arr2)))
        return arranged_faces
    
    def order_face(self, other_face, periodic=False):
        arranged_faces = self.face_arrangements(other_face)
        for arranged_face in arranged_faces:
            checked_points = [False for i in range(len(self.vertices))]
            checked_shifts = [None for i in range(len(self.vertices))]
            for i, vertice in enumerate(self.vertices):
                if self.vertices[i].is_point(arranged_face.vertices[i], periodic) == True:
                    checked_points[i] = True
                    checked_shifts[i] = self.vertices[i].point_shift(arranged_face.vertices[i])
            if np.array_equal(checked_points, [True for i in range(len(self.vertices))]): # Same points
                if np.all(checked_shifts == checked_shifts[0]): # Same periodic shifts for all points
                    return arranged_face
        print('Not the same face')
        return None
    
    def face_shift(self, other_face):
        ordered_face = self.order_face(other_face, periodic=True)
        if ordered_face != None:
            pts_shift = np.array([self.vertices[i].point_shift(ordered_face.vertices[i]) for i in range(len(self.vertices))])
            if np.all(pts_shift[0] == pts_shift):
                return pts_shift[0]
        else:
            print('Not the same face with periodic=True')
            return np.array([None, None, None]) 
    
    def is_face(self, other_face, periodic=False):
        same_face = False
        if len(self.vertices) != len(other_face.vertices): # Check that different sized faces aren't compared
            return same_face
        
        ordered_face = self.order_face(other_face, periodic)
        if ordered_face != None:
            same_face = True
        return same_face
    
    def shift_face(self, shift=np.array([0, 0, 0])):
        shifted_points = []
        for point in self.vertices:
            shifted_point = point.shift_point(shift)
            shifted_points.append(shifted_point)
        shifted_face = Face3D(shifted_points)
        return shifted_face
    
class Polyhedron3D(object):
    def __init__(self, points, edges, faces, center):
        self.points = points
        self.edges = edges
        self.faces = faces
        self.center = center
        self.poly_type = self.get_poly_type()
        self.hp = HiddenPrints()
        
        # Get the connectors Polyhedron3D shares with itself within the lattice
        ssp, sse, ssf = set(), set(), set()
        
        self.self_shared_points = [p for p in self.shared_points(self, periodic=True) if p in ssp or (ssp.add(p) or False)]
        self.self_shared_edges = [e for e in self.shared_edges(self, periodic=True) if e in sse or (sse.add(e) or False)]
        self.self_shared_faces = [f for f in self.shared_faces(self, periodic=True) if f in ssf or (ssf.add(f) or False)]
        
    def __repr__(self):
        return "Polyhedron3D %s %s" % (self.center, self.points)
    
    def __str__(self):
        return self.__repr__()
    
    def get_poly_type(self):
        # Can add in other classifications here, or not
        if len(self.faces) == 4:
            return 'tetrahedron'
        elif len(self.faces) == 5:
            return 'pentahedron'
        elif len(self.faces) == 6:
            return 'hexahedron'
        elif len(self.faces) == 7:
            return 'heptahedron'
        elif len(self.faces) == 8:
            return 'octahedron'
        elif len(self.faces) == 9:
            return 'enneahedron'
        elif len(self.faces) == 10:
            return 'decahedron'
        elif len(self.faces) == 11:
            return 'hendecahedron'
        elif len(self.faces) == 12:
            return 'dodecahedron'
        else:
            return 'N-sided polyhedron'
            
    def shared_points(self, other_polyhedron, periodic=False):
        with self.hp:
            shared_points = []
            for point in self.points:
                for other_point in other_polyhedron.points:
                    if point.is_point(other_point, periodic) == True:
                        shared_points.append(point)   
        return shared_points
        
    def shared_edges(self, other_polyhedron, periodic=False, duplicate=False):
        with self.hp:
            shared_edges = []
            for edge in self.edges:
                for other_edge in other_polyhedron.edges:
                    if edge.is_edge(other_edge, periodic) == True:
                        if edge.is_edge(other_edge) == True or duplicate == False: # Real space edges are shared
                            shared_edges.append(edge)
                        else: # Real space edges NOT shared; only periodic
                            shared_edges.append(edge)
                            shared_edges.append(other_edge)       
        return shared_edges
    
    def shared_faces(self, other_polyhedron, periodic=False, duplicate=False):
        with self.hp:
            shared_faces = []
            for face in self.faces:
                for other_face in other_polyhedron.faces:
                    if face.is_face(other_face, periodic) == True:
                        if face.is_face(other_face) == True or duplicate == False: # Real space faces are shared
                            shared_faces.append(face)
                        else: # Real space faces NOT shared; only periodic
                            shared_faces.append(face)
                            shared_faces.append(other_face)
        return shared_faces 
    
    def shared_center(self, other_polyhedron, periodic=False):
        with self.hp:
            if self.center.is_point(other_polyhedron.center, periodic):
                return self.center
            else:
                return None
    
    def is_poly(self, other_polyhedron, periodic=False):
        same_poly = False
        shared_faces = self.shared_faces(other_polyhedron, periodic, duplicate=False)
        shared_center = self.shared_center(other_polyhedron, periodic)
        if len(shared_faces) == len(self.faces) and shared_center != None:
            same_poly = True
        return same_poly
    
    def shift_poly(self, shift=np.array([0, 0, 0])):
        shifted_center = self.center.shift_point(shift)
        shifted_points = [point.shift_point(shift) for point in self.points]
        shifted_edges = [edge.shift_edge(shift) for edge in self.edges]
        shifted_faces = [face.shift_face(shift) for face in self.faces]
        shifted_poly = Polyhedron3D(shifted_points, shifted_edges, shifted_faces, shifted_center)
        return shifted_poly