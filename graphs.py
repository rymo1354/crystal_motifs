# -*- coding: utf-8 -*-

import sys
from crystal_motifs.geometry_objects import Point3D, Edge3D, Face3D
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Specie, Element
from crystal_motifs.constructors import PMGPeriodicStructurePolyhedron3DConstructor
import numpy as np

class GraphNode(object):
    def __init__(self, polyhedron):
        self.polyhedron = polyhedron
    
    def __repr__(self):
        return "GraphNode %s(%s) %s" % (self.polyhedron.center.specie, self.polyhedron.center.index, self.polyhedron.poly_type)
    
    def __str__(self):
        return self.__repr__()

class GraphEdge(object):
    def __init__(self, connector):
        self.connector = connector
    
    def __repr__(self):
        if isinstance(self.connector, Point3D):
            return 'GraphEdge Point3D'
        elif isinstance(self.connector, Edge3D):
            return 'GraphEdge Edge3D'
        elif isinstance(self.connector, Face3D):
            return 'GraphEdge Face3D'
        else:
            return 'GraphEdge %s' % type(self.connector)
    
    def __str__(self):
        return self.__repr__()

class GraphGenerator(object):
    def __init__(self, graph_type, nn_finder, sig_figs=1, **kwargs):
        # High default sig_figs tolerance (only compared to single sig fig)
        # This results in generally good performance for rattled structures
        self.graph_type = graph_type
        self.poly_generator = PMGPeriodicStructurePolyhedron3DConstructor(nn_finder, sig_figs, **kwargs)
        
    def subset_check(self, check_connectors, check_against_connectors, periodic=False):
        # return the check_connectors that are not in check_against_connectors
        # ordered by lower degree connector, higher degree connector i.e., (points, faces)
        # check real-space similarities, not periodic space
        reduced_connectors = []
        for check_connector in check_connectors:
            connector_present = []
            for check_against_connector in check_against_connectors:
                if type(check_connector) != type(check_against_connector):
                    if isinstance(check_connector, Point3D) and isinstance(check_against_connector, Edge3D):
                        truth_array = [check_connector.is_point(cacp, periodic) for cacp in check_against_connector.points]
                    elif isinstance(check_connector, Point3D) and isinstance(check_against_connector, Face3D):
                        truth_array = [check_connector.is_point(cacp, periodic) for cacp in check_against_connector.vertices]
                    elif isinstance(check_connector, Edge3D) and isinstance(check_against_connector, Face3D):
                        truth_array = [check_connector.is_edge(cacp, periodic) for cacp in check_against_connector.edges]
                    else:
                        print('Invalid connector types and/or connector order; should be (lower degree, higher degree)')
                        sys.exit(1)
                else:
                    if isinstance(check_connector, Point3D) and isinstance(check_against_connector, Point3D):
                        truth_array = [check_connector.is_point(check_against_connector, periodic)]
                    elif isinstance(check_connector, Edge3D) and isinstance(check_against_connector, Edge3D):
                        truth_array = [check_connector.is_edge(check_against_connector, periodic)]
                    elif isinstance(check_connector, Face3D) and isinstance(check_against_connector, Face3D):
                        truth_array = [check_connector.is_face(check_against_connector, periodic)]
                    else:
                        print('Invalid connector types; should be lists of Point3D, Edge3D and Face3D')
                        sys.exit(1)
                connector_present += truth_array
                
            if np.all(connector_present == False): # Check that connector isn't present in any other same or higher-level connectors
                reduced_connectors.append(check_connector)
                        
        return reduced_connectors
    
    def get_polyhedrons(self, structure, species=None, indices=None, **kwargs):
        polyhedrons = []
        if type(structure) == Structure: # check that a PMG structure object is passed
            if species != None:
                if all(isinstance(x, Specie) for x in species):
                    indices = [i for i in range(len(structure)) if structure[i].specie in species]
                    for indice in indices:
                        polyhedron = self.poly_generator.polyhedron_constructor(structure, indice)
                        polyhedrons.append(polyhedron)
                elif all(isinstance(x, Element) for x in species):
                    indices = [i for i in range(len(structure)) if structure[i].element in species]
                    for indice in indices:
                        polyhedron = self.poly_generator.polyhedron_constructor(structure, indice)
                        polyhedrons.append(polyhedron)
                else:
                    print('Incorrect species; must be list of pmg.core.periodic_table Species or Element types')
                    sys.exit(1)
            elif indices != None:
                if all(isinstance(x, int) and x < len(structure) for x in indices):
                    for indice in indices:
                        polyhedron = self.poly_generator.polyhedron_constructor(structure, indice)
                        polyhedrons.append(polyhedron)
            else:
                print('Incorrect indices; must be list of integers less than the length of the Structure object')
                sys.exit(1)
        return polyhedrons
        
    def get_nodes_edges_dct(self, polyhedrons):
        ### fmt: {(node_pair): {'point': [], 'edge': [], 'face': []}}
        all_nodes_edges_dct = {}
        graph_nodes = [GraphNode(p) for p in polyhedrons]
        for gn1_idx, gn1 in enumerate(graph_nodes):
            for gn2 in graph_nodes[gn1_idx:]: # Don't repeat checked connections
                node_key = (gn1, gn2)
                all_nodes_edges_dct[node_key] = {}
                if gn1.polyhedron.is_poly(gn2.polyhedron, periodic=False) == False: # Handles the periodicity problem
                    shared_points = gn1.polyhedron.shared_points(gn2.polyhedron, periodic=True)
                    shared_edges = gn1.polyhedron.shared_edges(gn2.polyhedron, periodic=True)
                    shared_faces = gn1.polyhedron.shared_faces(gn2.polyhedron, periodic=True)
                else:
                    shared_points = gn1.polyhedron.self_shared_points
                    shared_edges = gn1.polyhedron.self_shared_edges
                    shared_faces = gn1.polyhedron.self_shared_faces
                # Get the reducible connectors, checking to make sure the comparison is ok
                if shared_points != [] and shared_edges != []:
                    shared_points = self.subset_check(shared_points, shared_edges, periodic=True)
                if shared_edges != [] and shared_faces != []:
                    shared_edges = self.subset_check(shared_edges, shared_faces, periodic=True)
                
                all_nodes_edges_dct[node_key]['point'] = [GraphEdge(sp) for sp in shared_points]
                all_nodes_edges_dct[node_key]['edge'] = [GraphEdge(se) for se in shared_edges]
                all_nodes_edges_dct[node_key]['face'] = [GraphEdge(sf) for sf in shared_faces]
        
        return all_nodes_edges_dct
    
    def generate_complete_graph(self, structure, species=None, indices=None, **kwargs):
        ### specify **kwargs for self.graph_type
        ### specify either a list of pymatgen Species objects to build polyhedrons for, or indices of these sites
        ### get graph of shared Point3Ds, Edge3Ds and Face3Ds between polyhedrons of a structure and species/indices
        polyhedrons = self.get_polyhedrons(structure, species, indices, **kwargs)        
        all_nodes_edges_dct = self.get_nodes_edges_dct(polyhedrons)
        complete_graph = self.graph_type(**kwargs)
        for node_pair in list(all_nodes_edges_dct.keys()):
            for connector_type in list(all_nodes_edges_dct[node_pair]):
                for graph_edge in all_nodes_edges_dct[node_pair][connector_type]:
                    complete_graph.add_edge(node_pair[0], node_pair[1], graph_edge)
        return complete_graph
    
class GraphComparator(object):
    # Check to see if one graph is isomorphic to another, based on comparison rules
    pass
    
        

