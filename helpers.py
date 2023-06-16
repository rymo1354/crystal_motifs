# -*- coding: utf-8 -*-
import os, sys
import numpy as np
import re
from copy import deepcopy
from pymatgen.core.structure import Structure
from pymatgen.core.structure import Molecule
from pymatgen.core.periodic_table import Element, Species, DummySpecies
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.composition import Composition
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.analysis.local_env import CutOffDictNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase.data.pubchem import pubchem_conformer_search
import openbabel.pybel as pb
import pubchempy as pcp

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

class MolecularStructure(object):
    
    def __init__(self, structure,  organics=[Element('C'), Element('N'), Element('H')]):
        self.structure = structure
        self.organics = organics
        organic_neighbors = self.get_organic_neighbors(self.structure, self.organics)
        if len(organic_neighbors) == 0:
            self.molecular_structure = deepcopy(self.structure)
            self.molecular_structure.add_oxidation_state_by_guess()
        else:
            molecule_dcts = self.molecular_indices(organic_neighbors)
            estimated_molecules = self.make_molecules(self.structure, organic_neighbors, molecule_dcts)
            self.molecules, self.smiles = self.identify_molecules(estimated_molecules)
            self.molecular_structure, self.formulas_dct = self.dummy_structure(self.structure, self.molecules, self.smiles)

    def neighbor_dct(self, tolerance=0.2):
        cutoff_dct = {}
        for el1_ind, el1 in enumerate(self.organics):
            for el2 in self.organics[el1_ind:]:
                bond_length = el1.atomic_radius + el2.atomic_radius
                cutoff_dct[(str(el1), str(el2))] = bond_length*(1 + tolerance)
        return cutoff_dct

    def estimate_formal_charge(self, site, nn_dict):
        ''' Incomplete formal charge estimator for C, N and H organics'''
        neighbor_elements = [s.specie for s in nn_dict[site]['sites']]
        if site.specie == Element.H:
            valence = site.specie.valence[1]
            charge = valence - len(neighbor_elements)
        else:
            valence = site.specie.valence[1] + 2 # add in s-orbital electrons
            if (8 - valence) == len(neighbor_elements):
                charge = 0
            elif (8 - valence) < len(neighbor_elements):
                charge = len(neighbor_elements) - (8 - valence)
            else:
                print('Negative charge estimate for %s not supported at this time' % str(site))
                sys.exit(1)
        return charge

    def get_organic_neighbors(self, structure, organics):
        organic_neighbors = {}
        conn = CutOffDictNN(cut_off_dict=self.neighbor_dct())
        for site_ind, site in enumerate(structure):
            if site.specie in organics:
                organic_neighbors[site] = {'index': None, 'sites': [], 'indices': [], 'images': [], 'charge': None}
                nn_info = conn.get_nn_info(structure, site_ind)
                for nn in nn_info:
                    if nn['site'].specie in organics or nn['site'].specie is Element('H'):
                        organic_neighbors[site]['index'] = site_ind
                        organic_neighbors[site]['sites'].append(nn['site'])
                        organic_neighbors[site]['indices'].append(nn['site_index'])
                        organic_neighbors[site]['images'].append(nn['image'])
                organic_neighbors[site]['charge'] = self.estimate_formal_charge(site, organic_neighbors)
        return organic_neighbors

    def molecular_indices(self, organic_neighbors):
        site_indices = [organic_neighbors[k]['index'] for k in list(organic_neighbors.keys())]
        molecule_dcts = []

        def match_indices(organic_neighbors, index_dct):
            new_dct = deepcopy(index_dct)
            for s in list(organic_neighbors.keys()):
                if organic_neighbors[s]['index'] in list(new_dct.keys()):
                    for ii, i in enumerate(organic_neighbors[s]['indices']):
                        if i not in list(new_dct.keys()):
                            new_dct[i] = tuple(np.add(new_dct[organic_neighbors[s]['index']], organic_neighbors[s]['images'][ii])) 
            if set(new_dct.keys()) == set(index_dct.keys()):
                return new_dct
            else:
                return match_indices(organic_neighbors, new_dct)

        while site_indices:
            molecule_dct = match_indices(organic_neighbors, {site_indices[0]: (0, 0, 0)})
            molecule_dcts.append(molecule_dct)
            site_indices = [i for i in site_indices if i not in list(molecule_dct.keys())]
        return molecule_dcts

    def make_molecules(self, structure, organic_neighbors, molecule_dcts):
        molecules = []
        for molecule_dct in molecule_dcts:
            indices = [i for i in list(molecule_dct.keys())]
            charges = [organic_neighbors[structure[i]]['charge'] for i in indices]
            images = [molecule_dct[i] for i in indices]
            species = [structure[i].species for i in indices]
            coords = []
            for ind in indices:
                site = PeriodicSite(species=structure[ind].species, lattice=structure[ind].lattice, 
                                    coords=np.add(structure[ind].frac_coords, molecule_dct[ind]), coords_are_cartesian=False)
                coords.append(site.coords)
            molecule = Molecule(species=species, coords=coords, charge=np.sum(charges), site_properties={'indices': indices, 'images': images})
            molecules.append(molecule)
        return molecules
                
    def smiles_string(self, molecule):
        bma_obj = BabelMolAdaptor(molecule)
        pbm_obj = pb.Molecule(bma_obj.openbabel_mol)
        smi_obj = pb.readstring('can', pbm_obj.write("can"))
        smiles = "{}".format(smi_obj.write("can"))
        smiles_stripped = re.sub('\s+', '', smiles)
        return smiles_stripped 
 
    def identify_molecules(self, molecules):
        identified_molecules = []
        smiles_strings = []
        for id, molecule in enumerate(molecules):
            smiles_string = self.smiles_string(molecule)
            smiles_strings.append(smiles_string)
            try:
                db_query = pubchem_conformer_search(smiles=smiles_string)
                pubchem_id = db_query[0].data['PUBCHEM_COMPOUND_CID'] # if multiple, only take first entry
                charge = pcp.Compound.from_cid(pubchem_id).charge
                molecule.set_charge_and_spin(charge=charge)
            except ValueError:
                print('No pubchem entry for %s' % smiles_string)
            identified_molecules.append(molecule)
        return identified_molecules, smiles_strings

    def dummy_dct(self, formulas):
        formulas_dct = {}
        count = 0
        for formula in formulas:
            if formula in formulas_dct:
                pass
            else:
                formulas_dct[formula] = 'Z%s' % str(count)
                count += 1
        return formulas_dct

    def dummy_structure(self, s, molecule_objects, formulas):
        formulas_dct = self.dummy_dct(formulas)
        s_new = deepcopy(s)
        s_new.remove_sites(np.array([mo.site_properties['indices'] for mo in molecule_objects]).flatten())
        
        for molecule_id, molecule in enumerate(molecule_objects):
            if molecule.charge == 1: 
                s_new.append(species=Species('Cs', 1), coords=molecule.center_of_mass, coords_are_cartesian=True)
            elif molecule.charge == 2:
                s_new.append(species=Species('Ba', 2), coords=molecule.center_of_mass, coords_are_cartesian=True)
            elif molecule.charge == 3:
                s_new.append(species=Species('La', 3), coords=molecule.center_of_mass, coords_are_cartesian=True)
            elif molecule.charge == 4:
                s_new.append(species=Species('Si', 4), coords=molecule.center_of_mass, coords_are_cartesian=True)
            elif molecule.charge == 5:
                s_new.append(species=Species('Bi', 5), coords=molecule.center_of_mass, coords_are_cartesian=True)
            else:
                print('Molecule charge of %s not supported at this time;' % str(molecule.charge))
        s_new.add_oxidation_state_by_guess()

        start = len(s_new)-len(molecule_objects)
        for proxy_site_ii, proxy_site_ind in enumerate(range(start, len(s_new), 1)):
            s_new.replace(i=proxy_site_ind, species=DummySpecies(symbol=formulas_dct[formulas[proxy_site_ii]], 
                          oxidation_state=molecule.charge), coords=molecule_objects[proxy_site_ii].center_of_mass, coords_are_cartesian=True)
        return s_new, formulas_dct

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
