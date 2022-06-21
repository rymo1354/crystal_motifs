# -*- coding: utf-8 -*-
import unittest
from crystal_motifs.helpers import CoordinateConverter, SymmetrizeStructure
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from crystal_motifs.geometry_objects import Point3D, Edge3D, Face3D
from crystal_motifs.constructors import PMGPeriodicStructurePolyhedron3DConstructor
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import CrystalNN

class PointBuilder(object):
    def __init__(self):
        self.structure = Structure.from_file('test_poscars/mp-5924_MgSiO3_POSCAR.txt')
        self.matrix = self.structure.lattice.matrix
        self.shift = np.array([0, 0, -1])
        self.cc = CoordinateConverter()
        
        index1 = 12 # index 12 of structure, oxygen anion
        index2 = 13 # index 13 of structure, oxygen anion
        point1_coords = self.cc.fractional_to_cartesian(self.structure[index1].frac_coords, self.matrix[0], 
                                                        self.matrix[1], self.matrix[2])
        point2_coords = self.cc.fractional_to_cartesian(self.structure[index2].frac_coords, self.matrix[0], 
                                                        self.matrix[1], self.matrix[2])
        point3_coords = self.cc.fractional_to_cartesian(np.add(self.structure[index1].frac_coords, self.shift), 
                                                        self.matrix[0], self.matrix[1], self.matrix[2])
        self.point1 = Point3D(point1_coords, self.structure, index1)
        self.point2 = Point3D(point2_coords, self.structure, index2)
        self.point3 = Point3D(point3_coords, self.structure, index1)

class TestCoordinateConverter(unittest.TestCase):
    def setUp(self):
        self.pb = PointBuilder()
        
    def test_coordinate_conversion(self):
        # Conversion from fractional to cartesian and back again is the same, within floating point
        cc = CoordinateConverter()
        fractional_coordinates = cc.cartesian_to_fractional(self.pb.point1.coordinates, self.pb.matrix[0], 
                                                           self.pb.matrix[1], self.pb.matrix[2])
        cartesian_coordinates = cc.fractional_to_cartesian(fractional_coordinates, self.pb.matrix[0], 
                                                           self.pb.matrix[1], self.pb.matrix[2])
        self.assertIsNone(assert_array_almost_equal(self.pb.point1.coordinates, cartesian_coordinates))
        with self.assertRaises(AssertionError):
            assert_array_almost_equal(self.pb.point2.coordinates, cartesian_coordinates)

class TestPoint3D(unittest.TestCase):
    def setUp(self):
        self.pb = PointBuilder()
    
    def test_is_point(self):
        self.assertTrue(self.pb.point1.is_point(self.pb.point1, periodic=False)) # Same is same
        self.assertFalse(self.pb.point1.is_point(self.pb.point2, periodic=False)) # Different is not same
        self.assertTrue(self.pb.point1.is_point(self.pb.point3, periodic=True)) # Same in Periodic are same in Periodic
        self.assertFalse(self.pb.point1.is_point(self.pb.point3, periodic=False)) # Same in Periodic are not same in Cartesian
        
    def test_point_shift(self):
        self.assertIsNone(assert_array_equal(self.pb.shift, self.pb.point1.point_shift(self.pb.point3))) # Recover shift
        self.assertIsNone(assert_array_equal(np.array([None, None, None]), self.pb.point1.point_shift(self.pb.point2)))
    
    def test_shift_point(self):
        shifted_point3 = self.pb.point3.shift_point(np.negative(self.pb.shift))
        self.assertIsNone(assert_array_almost_equal(shifted_point3.coordinates, self.pb.point1.coordinates))
    
class TestEdge3D(unittest.TestCase):
    def setUp(self):
        self.pb = PointBuilder()
        self.edge1 = Edge3D(self.pb.point1, self.pb.point2)
        self.edge2 = Edge3D(self.pb.point2, self.pb.point1)
        self.edge3 = Edge3D(self.pb.point1, self.pb.point3)
        self.edge4 = Edge3D(self.pb.point3, self.pb.point2.shift_point(self.pb.shift))
        
    def test_order_edge(self):
        for point_ind, point in enumerate(self.edge1.points):
            self.assertFalse(point.is_point(self.edge2.points[point_ind]))
            self.assertTrue(point.is_point(self.edge1.order_edge(self.edge2).points[point_ind]))
    
    def test_edge_shift(self):
        self.assertIsNone(assert_array_equal(self.pb.shift, self.edge1.edge_shift(self.edge4))) # Shifted edge
        self.assertIsNone(assert_array_equal(np.array([None, None, None]), self.edge1.edge_shift(self.edge3))) # Different edge
        
    def test_is_edge(self):
        self.assertTrue(self.edge1.is_edge(self.edge4, periodic=True))
        self.assertFalse(self.edge1.is_edge(self.edge4, periodic=False))
        self.assertFalse(self.edge1.is_edge(self.edge3, periodic=True))
    
    def test_shift_edge(self):
        shifted_edge1 = self.edge1.shift_edge(self.pb.shift)
        self.assertTrue(shifted_edge1.is_edge(self.edge4, periodic=False))
        
    def test_distance(self):
        # These do not always agree with the pymatgen distances between points - determine why this is
        pass
    
class TestFace3D(unittest.TestCase):
    def setUp(self):
        self.pb = PointBuilder()
        self.face1 = Face3D([self.pb.point1, self.pb.point2, self.pb.point3])
        self.face2 = Face3D([self.pb.point1.shift_point(self.pb.shift), self.pb.point2.shift_point(self.pb.shift), 
                             self.pb.point3.shift_point(self.pb.shift)])
        self.face3 = Face3D([self.pb.point1.shift_point(self.pb.shift), self.pb.point2, self.pb.point3])
        self.face4 = Face3D([self.pb.point2, self.pb.point1, self.pb.point3])
    
    def test_get_shape(self):
        self.assertEqual(['3-sided polygon', '3-sided polygon', '3-sided polygon'], [self.face1.shape, self.face2.shape, self.face3.shape])
    
    def test_order_face(self):
        for point_ind, point in enumerate(self.face1.vertices):
            self.assertFalse(point.is_point(self.face2.vertices[point_ind])) # Not the same in real space
            self.assertTrue(point.is_point(self.face1.order_face(self.face4).vertices[point_ind], periodic=True))
            self.assertTrue(point.is_point(self.face1.order_face(self.face4).vertices[point_ind]))
    
    def test_face_shift(self):
        self.assertIsNone(assert_array_equal(self.pb.shift, self.face1.face_shift(self.face2))) # Shifted face
        self.assertIsNone(assert_array_equal(np.array([None, None, None]), self.face1.face_shift(self.face3))) # Different edge
    
    def test_is_face(self):
        self.assertTrue(self.face1.is_face(self.face2, periodic=True))
        self.assertFalse(self.face1.is_face(self.face2, periodic=False))
        self.assertFalse(self.face1.is_face(self.face3, periodic=True)) # Different shifts, not the same face
    
    def test_shift_face(self):
        shifted_face4 = self.face4.shift_face(self.pb.shift)
        self.assertTrue(shifted_face4.is_face(self.face2, periodic=False))
        
class PolyhedronBuilder(object):
    def __init__(self):
        self.perovskite = Structure.from_file('test_poscars/mp-22013_LaTiO3_POSCAR.txt')
        self.perovskite.add_oxidation_state_by_element({'La': 3, 'Ti': 3, 'O': -2})
        self.cubic_perovskite = Structure.from_file('test_poscars/mp-5827_CaTiO3_POSCAR.txt')
        self.cubic_perovskite.add_oxidation_state_by_element({'Ca': 2, 'Ti': 4, 'O': -2})
        self.fluorite = Structure.from_file('test_poscars/mp-20194_CeO2_POSCAR.txt')
        self.fluorite.add_oxidation_state_by_element({'Ce': 4, 'O': -2})
        
        self.rattled_fluorite = Structure.from_file('test_poscars/rattled_mp-20194_CeO2_POSCAR.txt')
        self.rattled_fluorite.add_oxidation_state_by_element({'Ce': 4, 'O': -2})
        ss = SymmetrizeStructure(sym_prec=0.1, angle_tolerance=5)
        self.sym_rattled_fluorite = ss.symmetrize_structure(self.rattled_fluorite)
        
        self.constructor = PMGPeriodicStructurePolyhedron3DConstructor(CrystalNN, weighted_cn=True, cation_anion=True)
        
        self.ti_octahedron = self.constructor.polyhedron_constructor(self.perovskite, 4)
        self.ce_cube = self.constructor.polyhedron_constructor(self.fluorite, 0)
        self.ti_cubic_octahedron = self.constructor.polyhedron_constructor(self.cubic_perovskite, 1)
        self.rattled_ce_non_cube = self.constructor.polyhedron_constructor(self.rattled_fluorite, 0)
        self.sym_rattled_ce_cube = self.constructor.polyhedron_constructor(self.sym_rattled_fluorite, 0)
        
class TestPolyhedron3DConstructor(unittest.TestCase):
    def setUp(self):
        self.polyb = PolyhedronBuilder()
        
    def test_shapes(self):
        self.assertEqual(self.polyb.ti_octahedron.poly_type, '8-sided polyhedron')
        self.assertEqual(self.polyb.ce_cube.poly_type, '6-sided polyhedron')
        self.assertEqual(self.polyb.sym_rattled_ce_cube.poly_type, '6-sided polyhedron') # Expected behavior when symmetrized
        self.assertNotEqual(self.polyb.rattled_ce_non_cube.poly_type, '6-sided polyhedron') # Unexpected behavior when not symmetrized
        
    ### Could add more tests to this, but is probably sufficient for now ###
    
class TestPolyhedron3D(unittest.TestCase):
    def setUp(self):
        self.polyb = PolyhedronBuilder()
        self.la_octahedron_image = self.polyb.constructor.polyhedron_constructor(self.polyb.perovskite, 4, np.array([0, 0, -1]))
        self.ce_cube_image = self.polyb.constructor.polyhedron_constructor(self.polyb.fluorite, 0, np.array([1, 0, 0]))
        
    def get_all_shared(self, polyhedron, structure, indices, type_shared):
        shared = []
        for index in indices:
            other_polyhedron = self.polyb.constructor.polyhedron_constructor(structure, index)
            if type_shared == Point3D:
                if polyhedron.is_poly(other_polyhedron) == False:
                    shared += polyhedron.shared_points(other_polyhedron, periodic=True)
                else:
                    shared += polyhedron.self_shared_points
            elif type_shared == Edge3D:
                if polyhedron.is_poly(other_polyhedron) == False:
                    shared += polyhedron.shared_edges(other_polyhedron, periodic=True)
                else:
                    shared += polyhedron.self_shared_edges
            elif type_shared == Face3D:
                if polyhedron.is_poly(other_polyhedron) == False:
                    shared += polyhedron.shared_faces(other_polyhedron, periodic=True)
                else:
                    shared += polyhedron.self_shared_faces
            else:
                pass
        return list(set(shared)) # To not consider repeated connectors for units that are repeated
        
    def test_shared_points(self):
        shared_perovskite = self.get_all_shared(self.polyb.ti_octahedron, self.polyb.perovskite, [4, 5, 6, 7], Point3D)
        shared_cubic_perovskite = self.get_all_shared(self.polyb.ti_cubic_octahedron, self.polyb.cubic_perovskite, 
                                                 [1], Point3D)
        shared_fluorite = self.get_all_shared(self.polyb.ce_cube, self.polyb.fluorite, [0, 1, 2, 3], Point3D)
        shared_sr_fluorite = self.get_all_shared(self.polyb.sym_rattled_ce_cube, self.polyb.sym_rattled_fluorite, [0], Point3D)
        
        self.assertEqual(len(shared_perovskite), 6) # six vertices shared in periodic space by the octahedron
        self.assertEqual(len(shared_cubic_perovskite), 6) # six vertices shared in periodic space by the octahedron
        self.assertEqual(len(shared_fluorite), 8) # eight vertices shared in periodic space by the high-sym cube
        self.assertEqual(len(shared_sr_fluorite), 8) # eight vertices sared in periodic space by the lower sym cube
        
    def test_shared_edges(self):
        shared_perovskite = self.get_all_shared(self.polyb.ti_octahedron, self.polyb.perovskite, [4, 5, 6, 7], Edge3D)
        shared_cubic_perovskite = self.get_all_shared(self.polyb.ti_cubic_octahedron, self.polyb.cubic_perovskite, 
                                                 [1], Edge3D)
        shared_fluorite = self.get_all_shared(self.polyb.ce_cube, self.polyb.fluorite, [0, 1, 2, 3], Edge3D)
        shared_sr_fluorite = self.get_all_shared(self.polyb.sym_rattled_ce_cube, self.polyb.sym_rattled_fluorite, [0], Edge3D)
        
        self.assertEqual(len(shared_perovskite), 0) # no edges shared in periodic space by the octahedron
        self.assertEqual(len(shared_cubic_perovskite), 0) # no edges shared in periodic space by the octahedron
        self.assertEqual(len(shared_fluorite), 12) # twelve edges shared in periodic space by the high-sym cube
        self.assertEqual(len(shared_sr_fluorite), 12) # twelve edges shared in periodic space by the lower sym cube
        
    def test_shared_faces(self):
        shared_perovskite = self.get_all_shared(self.polyb.ti_octahedron, self.polyb.perovskite, [4, 5, 6, 7], Face3D)
        shared_cubic_perovskite = self.get_all_shared(self.polyb.ti_cubic_octahedron, self.polyb.cubic_perovskite, 
                                                 [1], Face3D)
        shared_fluorite = self.get_all_shared(self.polyb.ce_cube, self.polyb.fluorite, [0, 1, 2, 3], Face3D)
        shared_sr_fluorite = self.get_all_shared(self.polyb.sym_rattled_ce_cube, self.polyb.sym_rattled_fluorite, [0], Face3D)
        
        self.assertEqual(len(shared_perovskite), 0) # no faces shared in periodic space by the octahedron
        self.assertEqual(len(shared_cubic_perovskite), 0) # no faces shared in periodic space by the octahedron
        self.assertEqual(len(shared_fluorite), 0) # no faces shared in periodic space by the high-sym cube
        self.assertEqual(len(shared_sr_fluorite), 0) # no faces shared in periodic space by the lower sym cube
    
if __name__ == '__main__':
    unittest.main()
