# -*- coding: utf-8 -*-

import json
from scipy.spatial.qhull import QhullError
from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import DummySpecies, Species, Specie, Element
from pymatgen.analysis.local_env import CrystalNN, BrunnerNN_relative
from networkx import MultiGraph
from crystal_motifs.graphs import GraphGenerator
from crystal_motifs.geometry_objects import Point3D
from tqdm import tqdm
from time import process_time, time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import functools
import argparse
from crystal_motifs.helpers import HiddenPrints, MolecularStructure
from pathlib import Path
import os
import xml
import sys
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.structure import Structure

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', '--directory', help='directory to check', required=True)
    parser.add_argument(
        '-w', '--write_file_name', help='.json that is written', required=True)
    parser.add_argument(
        '-c', '--consider_unparsable', help='Consider unparsable files; default = False', action='store_true')
    parser.add_argument(
        '-a', '--anions_AB', help='Anions assigned to A/B sites of anti-perovskites; default = False', action='store_true')
    args = parser.parse_args()

    return args

def get_paths(directory, consider_unparsable):
    paths = []
    for root, dirs, files in os.walk(directory):
        if consider_unparsable is True:
            contcar_path = os.path.join(root, 'CONTCAR')
            if os.path.exists(contcar_path):
                paths.append(contcar_path)
        else:
            vasprun_path = os.path.join(root, 'vasprun.xml')
            if os.path.exists(vasprun_path):
                try:
                    v = Vasprun(vasprun_path)
                    if v.converged is True:
                        contcar_path = os.path.join(root, 'CONTCAR')
                        if os.path.exists(contcar_path):
                            paths.append(contcar_path)
                    else:
                        continue
                except xml.etree.ElementTree.ParseError:
                    continue
    return paths

def get_data(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
    return data

def write_data(write_file, write_dict):
    with open(write_file, 'w') as f:
        json.dump(write_dict, f)

def species_mapping(unique_species, anions_ab):
    anions = [s for s in unique_species if s.oxi_state < 0]
    cations = [s for s in unique_species if s.oxi_state > 0]

    def get_mappings(lst):
        lst1 = []
        lst2 = []
        split = len(lst) - 1
        if split == 1:
            lst1.append([lst[0]])
            lst2.append([lst[1]])
        else:
            s_indices = ''.join([str(i) for i in range(len(lst))])
            while split >= 0.5*len(lst):
                s_combos = list(combinations(s_indices, split))
                spec_combos = []
                o_spec_combos = []
                for tup in s_combos:
                    i_lst = []
                    for i in tup:
                        i_lst.append(int(i))
                    spec_combos.append([lst[i] for i in i_lst])
                    o_spec_combos.append([lst[i] for i in range(len(lst)) if i not in i_lst])
                lst1 += spec_combos
                lst2 += o_spec_combos
                split += -1
        return lst1, lst2

    if anions_ab is False and len(anions) >= 1 and len(cations) > 1: # true for perovskites
        X_sites = anions
        A_sites, B_sites = get_mappings(cations)
    elif anions_ab is True and len(cations) >= 1 and len(anions) > 1: # true for anti-perovskites
        X_sites = cations
        A_sites, B_sites = get_mappings(anions)
    else: # not enough unique species to form perovskite
        print('%s does not have enough unique species to form a perovskite or anti-perovksite' % str(unique_species))
        sys.exit(1)
    return A_sites, B_sites, X_sites

def check_3D_perovskite(cmpd, graph, a_species, b_species):
    is_3d_perovskite = True
    nodes = list(graph.nodes)
    for node in nodes:
        if node.polyhedron.center.specie in b_species:
            # 1. Check that all B-site nodes are octahedrons with only triangular faces
            if node.polyhedron.poly_type == 'Nothing coordinating':
                print('%s: B-site(s) with no neighbors; check oxidation states and/or structure' % cmpd)
                is_3d_perovskite = False
                break
            if node.polyhedron.poly_type != '8-sided polyhedron':
                print('%s: Not all B-site nodes are octahedrons' % cmpd)
                is_3d_perovskite = False
                break
            face_shapes = np.array([face.shape for face in node.polyhedron.faces])
            if np.all(face_shapes == '3-sided polygon') is False:
                print('%s: Not all B-site faces are triangles' % cmpd)
                is_3d_perovskite = False
                break
            # 2. Check that all B-site nodes share six Point3D edges
            all_edge_point3Ds = []
            for node2 in nodes:
                if node2.polyhedron.center.specie in b_species:
                    edge_data = graph.get_edge_data(node, node2)
                    if edge_data is not None:
                        shared_edges = list(edge_data.keys())
                        edges_point3Ds = [isinstance(shared_edge.connector, Point3D) for shared_edge in shared_edges]
                        all_edge_point3Ds += edges_point3Ds
            all_point3Ds = np.all(np.array(all_edge_point3Ds) == True)
            if all_point3Ds is False:  # Not all edges are Point3Ds
                print('%s: Not all B-site graph edges are Point3Ds' % cmpd)
                is_3d_perovskite = False
                break
            if len(all_edge_point3Ds) != 6:  # Six graph edges shared
                print('%s: Not all B-site octahedrons share exactly six Point3Ds' % cmpd)
                is_3d_perovskite = False
                break

        elif node.polyhedron.center.specie in a_species: # if in a_species
            # 3. Check that all A-site nodes have between 8 and 12 coordinating neighbors
            if node.polyhedron.poly_type == 'Nothing coordinating':
                print('%s: A-site(s) with no neighbors; check oxidation states and/or structure' % cmpd)
                is_3d_perovskite = False
                break
            if len(node.polyhedron.points) < 8 or len(node.polyhedron.points) > 12:
                print('%s: A-site(s) with %s neighbors' % (cmpd, str(len(node.polyhedron.points))))
                is_3d_perovskite = False
                break
        else:
            print('Inappropriate nodes, improper species')
            sys.exit(1)

    return is_3d_perovskite

def map_organics(specie, dct):
    if isinstance(specie, DummySpecies):
        organic = '(' + list(dct.keys())[list(dct.values()).index(specie._symbol)] + ')'
        if specie.oxi_state > 0:
            return organic + str(specie.oxi_state) + '+'
        else:
            return organic + str(np.abs(specie.oxi_state)) + '-'
    else:
        return str(specie) 

def classifier(anions_ab, path):
    hp = HiddenPrints()

    # Get structure and add oxidation states
    struct = Structure.from_file(path) 
    mso = MolecularStructure(struct) # check if molecules can be identified in the structure & add oxidation states
    
    # Designate NNfinder to be used
    if len(mso.structure) == len(mso.molecular_structure): # no molecules identified
        gg = GraphGenerator(CrystalNN, 8, cation_anion=True)
        sp_dct = {} 
    else:
        gg = GraphGenerator(BrunnerNN_relative) # works for DummySpecies when molecules were identified
        sp_dct = mso.formulas_dct

    # Get the species assignments and generate the graph
    unique_species = list(np.unique(mso.molecular_structure.species))
    A_sites, B_sites, X_sites = species_mapping(unique_species, anions_ab)
    graph_species = [s for s in unique_species if s not in X_sites] # X_sites are vertices, edges and faces, to be nodes in graph
    cmpd = mso.structure.composition.reduced_formula # return this, without DummySpecies  
    with hp:
        graph = gg.generate_complete_graph(mso.molecular_structure, MultiGraph, species=graph_species)
    
    # Iterate over possible site assignments and designate perovskite/non-perovskite
    assignment_dct = {}
    assignment_dct[cmpd] = {}
    count = 1
    for i, assgn in enumerate(A_sites):
        assignment = 'assignment_%s' % str(count)
        is_3d_perovskite = check_3D_perovskite(cmpd, graph, A_sites[i], B_sites[i])
        to_assign = {'is_3d_perovskite': is_3d_perovskite, 'A': [map_organics(a, sp_dct) for a in A_sites[i]], 'B': [map_organics(b, sp_dct) for b in B_sites[i]], 'X': [map_organics(x, sp_dct) for x in X_sites]}
        assignment_dct[cmpd][assignment] = to_assign
        count += 1

        assignment = 'assignment_%s' % str(count)
        is_3d_perovskite = check_3D_perovskite(cmpd, graph, B_sites[i], A_sites[i])
        to_assign = {'is_3d_perovskite': is_3d_perovskite, 'A': [map_organics(b, sp_dct) for b in B_sites[i]], 'B': [map_organics(a, sp_dct) for a in A_sites[i]], 'X': [map_organics(x, sp_dct) for x in X_sites]}
        assignment_dct[cmpd][assignment] = to_assign
        count += 1

    return cmpd, assignment_dct 

def pool_classifier_map(paths, anions_ab, nprocs):
    partial_classifier = functools.partial(classifier, anions_ab)
    with ProcessPoolExecutor(max_workers=nprocs) as executor:
        data_dict = {cmpd: a_dct[cmpd] for cmpd, a_dct in executor.map(partial_classifier, paths)}
    return data_dict


if __name__ == '__main__':
    t1_p_start = process_time()
    t1_start = time()
    args = argument_parser()

    print('Running analysis...')
    paths = get_paths(args.directory, args.consider_unparsable)
    classifier_dict = pool_classifier_map(paths, args.anions_AB, mp.cpu_count())
    print(classifier_dict)
    print('Writing data...')
    write_data(args.write_file_name, classifier_dict)

    t1_p_finish = process_time()
    t1_finish = time()
    print('Elapsed time:', t1_finish-t1_start)
    print('Process time:', t1_p_finish-t1_p_start)
