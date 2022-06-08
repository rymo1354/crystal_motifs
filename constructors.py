# -*- coding: utf-8 -*-
import sys
import numpy as np
from copy import deepcopy
from geometry_objects import Point3D, Edge3D, Face3D, Polyhedron3D
from helpers import CoordinateConverter
from scipy.spatial import ConvexHull

class Face3DfromSimplicesConstructor(object):
    def __init__(self, sig_figs=1):
        # Atom positions of structures should not be within sig_figs of each other
        self.sig_figs = sig_figs
        
    def check_forms_plane(self, three_coordinates):
        form_plane = True
        units = np.round(np.array([np.divide(coord, self.vector_normalization(coord)) for coord in three_coordinates]), 
                         self.sig_figs)
        # Check if the three coordinates are on same line / are the same point, within sig_figs
        # Check to make sure no points are equal to zero
        if len(np.unique(units, axis=0)) == 1 or len(np.unique(three_coordinates, axis=0)) == 1: 
            form_plane = False
        return form_plane
    
    def vector_normalization(self, vector):
        if np.array_equal(vector, np.array([0, 0, 0])):
            return 1 # Avoid np errors with dividing by 0
        else:
            return np.linalg.norm(vector)
        
    def coplanar_coord(self, plane_coords, coord):
        coplanar = False
        if self.check_forms_plane(plane_coords) == False:
            print("Points %s cannot form plane" % plane_coords)
            sys.exit(1)
        else:
            distances = [self.get_distance(plane_coord, coord) for plane_coord in plane_coords]
            planes_coord_ind = distances.index(max(distances))
            plane_vector = self.get_plane_vector(plane_coords, planes_coord_ind)
            normal_component = np.dot(plane_vector, np.subtract(plane_coords[planes_coord_ind], coord))
            if np.round(np.abs(normal_component), self.sig_figs) == 0:
                coplanar = True
        return coplanar
    
    def coplanar_coords(self, plane_coords, coords):
        in_plane = [False for i in range(len(coords))]
        for coord_ind, coord in enumerate(coords):
            is_in_plane = self.coplanar_coord(plane_coords, coord)
            in_plane[coord_ind] = is_in_plane
        if np.array_equal(in_plane, np.array([True for i in range(len(coords))])):
            return True
        else:
            return False
    
    def get_plane_vector(self, coords_in_plane, coords_point_ind):
        all_inds = [i for i in range(len(coords_in_plane)) if i != coords_point_ind]
        
        vector1 = np.subtract(coords_in_plane[coords_point_ind], coords_in_plane[all_inds[0]])
        vector2 = np.subtract(coords_in_plane[coords_point_ind], coords_in_plane[all_inds[1]])
        vector_normal = np.cross(vector1, vector2)
        return vector_normal 
    
    def get_distance(self, coord1, coord2):
        vector = np.subtract(coord1, coord2)
        return np.linalg.norm(vector)
    
    def check_merged(self, simplice, merged_list):
        if merged_list != []:
            if np.any(np.all(simplice == merged_list, axis=1)):
                return True
            else:
                return False
            return False
        
    def merge_simplices(self, parent_simplice, simplice_coordinates, merge_simplice):
        new_simplice = deepcopy(parent_simplice)
        for c_ind in merge_simplice:
            distances = [self.get_distance(simplice_coordinates[c_ind], 
                                           simplice_coordinates[pc_ind]) for pc_ind in parent_simplice]
            min_inds_dist = np.argpartition(distances, 2) # indices of two smallest distances
            min_ind_ps = parent_simplice.index(parent_simplice[min_inds_dist[0]])
            min2_ind_ps = parent_simplice.index(parent_simplice[min_inds_dist[1]])
            bigger_ind_ps = np.max([min_ind_ps, min2_ind_ps])
            if np.abs(np.subtract(min_ind_ps, min2_ind_ps)) == 1: # Adjacent simplices connected in parent_simplice
                new_simplice = new_simplice[:bigger_ind_ps] + [c_ind] + new_simplice[bigger_ind_ps:]
            else: # indices wrap around to preserve connectivity
                new_simplice = new_simplice + [c_ind]
        return new_simplice
    
    def get_true_faces(self, original_simplices, simplice_coordinates):
        # Simplices and their associated coordinates from ConvexHull
        new_simplices = []
        merged = [] # Do not consider these anymore, have already been merged
        
        for s_i, simplice in enumerate(list(original_simplices)):
            if self.check_merged(simplice, merged) == True:
                continue # If simplice is already merged, pass it over in the for loop
            new_simplice = list(deepcopy(simplice)) # Create the new simplice to be formed
            for s_i2, simplice2 in enumerate(list(original_simplices)[s_i+1:]):
                if self.check_merged(simplice2, merged) == True:
                    continue # If simplice to check against is already merged, pass it over        
                coords1 = [simplice_coordinates[s] for s in simplice] # To check plane, from 3 points
                coords2 = [simplice_coordinates[s2] for s2 in simplice2]    
                if self.coplanar_coords(coords1, coords2) == True:
                    # Use s_i2 indices to remove from copy_original_simplices
                    different_lst = [s for s in simplice2 if s not in new_simplice] # compare to new simplice
                    new_simplice = self.merge_simplices(new_simplice, simplice_coordinates,
                                                       different_lst)
                    merged.append(simplice2)
            new_simplices.append(new_simplice)
        return new_simplices
    
class PMGPeriodicStructurePolyhedron3DConstructor:
    def __init__(self, nn_finder, sig_figs=1, **kwargs):
        self.nn_finder = nn_finder(**kwargs)
        self.sig_figs = sig_figs
        self.cc = CoordinateConverter()
        # Could be error(s) assuming matrix a, b, and c vectors are 0, 1, 2, 
        # elements from structure.lattice.matrix; should check this
        
    def __repr__(self):
        return "PMGPeriodicStructurePolyhedron3DConstructor with %s" % self.nn_finder.__class__
    
    def __str__(self):
        return self.__repr__()
        
    def get_points(self, structure, center_site_index, center_image=np.array([0., 0., 0.])):
        nn = self.nn_finder.get_nn_info(structure, center_site_index)
        neighbor_points = []
        for i in range(len(nn)):
            cartesian_coordinates = self.cc.fractional_to_cartesian(np.add(nn[i]['site'].frac_coords, center_image), 
                                                               structure.lattice.matrix[0], 
                                                               structure.lattice.matrix[1], 
                                                               structure.lattice.matrix[2])
            neighbor_point = Point3D(cartesian_coordinates, structure, nn[i]['site_index'], # PMG Version specific 
                                     center_image + nn[i]['image'])
            neighbor_points.append(neighbor_point)
        return neighbor_points
        
    def get_edges(self, simplices, neighbor_points):
        edge_indices_sets = []
        for simplice in simplices: # list of simplices
            for i in range(len(simplice)):
                if i == len(simplice)-1:
                    edge_indices = set([simplice[i], simplice[0]])
                else:
                    edge_indices = set([simplice[i], simplice[i+1]])    
                if edge_indices not in edge_indices_sets:
                    edge_indices_sets.append(edge_indices)
        edge_indices_lsts = [list(edge_indices_set) for edge_indices_set in edge_indices_sets]
        
        edges = []
        for edge_indices in edge_indices_lsts:
            edge = Edge3D(neighbor_points[edge_indices[0]], neighbor_points[edge_indices[1]])
            edges.append(edge)
        
        return edges 
    
    def get_faces(self, simplices, neighbor_points):
        faces = []
        for simplice in simplices: # list of simplices
            face = Face3D(np.array([neighbor_points[i] for i in simplice]))
            faces.append(face)
        return faces
    
    def get_center_point(self, structure, center_site_index, center_image=np.array([0., 0., 0.])):
        center_site = structure[center_site_index]
        cartesian_coordinates = self.cc.fractional_to_cartesian(np.add(center_site.frac_coords, center_image), 
                                                           structure.lattice.matrix[0], 
                                                           structure.lattice.matrix[1], 
                                                           structure.lattice.matrix[2])
        center_point = Point3D(cartesian_coordinates, structure, center_site_index, center_image)
        return center_point
    
    def polyhedron_constructor(self, structure, center_site_index, center_image=np.array([0., 0., 0.])):
        if center_site_index > len(structure):
            print('Site index for center must be < %s' % str(len(structure)))
            sys.exit(1)
        
        center_point = self.get_center_point(structure, center_site_index, center_image)
        neighbor_points = self.get_points(structure, center_site_index, center_image)
        neighbor_coordinates = [n.coordinates for n in neighbor_points]
        
        ch = ConvexHull(neighbor_coordinates) # still use this for now
        ffsc = Face3DfromSimplicesConstructor(self.sig_figs) # for numerical stability
        refaced_simplices = ffsc.get_true_faces(ch.simplices, ch.points)
        
        edges = self.get_edges(refaced_simplices, neighbor_points)
        faces = self.get_faces(refaced_simplices, neighbor_points)
        polyhedron = Polyhedron3D(neighbor_points, edges, faces, center_point)
        
        return polyhedron
    
class TranslatedPolyhedron3DConstructor:
    # WARNING: This is currently untested code
    def __init__(self):
        return
    
    def __repr__(self):
        return "TranslatedPolyhedron3DConstructor"
    
    def __str__(self):
        return self.__repr__()
    
    def get_shared_polyhedrons(self, polyhedrons, shared_object, periodic=False):
        # Need to compute the center index and shift for the polyhedron along the shared edge
        shared_polyhedrons = []
        shared_connectors = [[] for i in range(len(polyhedrons))]
        for polyhedron_ind, polyhedron in enumerate(polyhedrons):
            shared = False
            if isinstance(shared_object, Point3D):
                for point in polyhedron.points:
                    if shared_object.is_point(point, periodic=periodic):
                        shared = True
                        shared_connectors[polyhedron_ind].append(point)
            elif isinstance(shared_object, Edge3D):
                for edge in polyhedron.edges:
                    if shared_object.is_edge(edge, periodic=periodic):
                        shared = True
                        shared_connectors[polyhedron_ind].append(edge)
            elif isinstance(shared_object, Face3D):
                for face in polyhedron.faces:
                    if shared_object.is_face(face, periodic=periodic):
                        shared = True
                        shared_connectors[polyhedron_ind].append(face)
            else:
                continue
            
            if shared == True:
                shared_polyhedrons.append(polyhedron)
            else:
                shared_polyhedrons.append(None)
        return shared_polyhedrons, shared_connectors
    
    def shared_polyhedrons_constructor(self, polyhedrons, shared_object, periodic=False):
        
        shared_polys, connectors = self.get_shared_polyhedrons(polyhedrons, shared_object, periodic=True)
        shifted_polys = []
        for poly_ind, poly in enumerate(shared_polys):
            shifts = []
            if poly != None:
                conns = connectors[poly_ind]
                for conn in conns:
                    if isinstance(conn, Point3D):
                        shift = shared_object.shift_point(conn)
                        if not np.array_equal(shift, np.array([0, 0, 0])): # to ignore same polyhedrons generated
                            shifts.append(shift)
                    elif isinstance(conn, Edge3D):
                        shift = shared_object.shift_edge(conn)
                        if not np.array_equal(shift, np.array([0, 0, 0])):
                            shifts.append(shift)
                    elif isinstance(conn, Face3D):
                        shift = shared_object.shift_face(conn)
                        if not np.array_equal(shift, np.array([0, 0, 0])):
                            shifts.append(shift)
            for shift in shifts:
                shifted_poly = poly.shift_poly(shift)
                shifted_polys.append(shifted_poly)
                
        return shifted_polys