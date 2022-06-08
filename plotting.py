# -*- coding: utf-8 -*-
import numpy as np
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from geometry_objects import Point3D, Edge3D, Face3D
from sympy import Polygon
import matplotlib.pyplot as plt
from matplotlib import cm
import networkx as nx
import sys

class PMGPeriodicStructurePolyhedron3DPlotting:
    def __init__(self):
        return
    
    def get_unit_cell_shift_centroid(self, matrix, centroid):
        xyz_prime = np.add(matrix[0], np.add(matrix[1], matrix[2]))
        xyz_prime_half = np.divide(xyz_prime, 2)
        shift = np.subtract(centroid, xyz_prime_half)
        return shift

    def get_unit_cell_from_matrix(self, matrix, shift=np.array([0, 0, 0])):
        origin = np.array([0, 0, 0])
        x_prime = matrix[0]
        y_prime = matrix[1]
        z_prime = matrix[2]
    
        xy_prime = np.add(matrix[0], matrix[1])
        xz_prime = np.add(matrix[0], matrix[2])
        yz_prime = np.add(matrix[1], matrix[2])
        xyz_prime = np.add(xy_prime, z_prime)

        linex_prime = np.array([origin, x_prime])
        liney_prime = np.array([origin, y_prime])
        linez_prime = np.array([origin, z_prime])

        linexxy_prime = np.array([xy_prime, x_prime])
        lineyyz_prime = np.array([yz_prime, y_prime])
        linexxz_prime = np.array([xz_prime, x_prime])

        linexyy_prime = np.array([xy_prime, y_prime])
        lineyzz_prime = np.array([yz_prime, z_prime])
        linexzz_prime = np.array([xz_prime, z_prime])

        linexzxyz_prime = np.array([xz_prime, xyz_prime])
        lineyzxyz_prime = np.array([yz_prime, xyz_prime])
        linexyxyz_prime = np.array([xy_prime, xyz_prime])

        all_lines = np.array([linex_prime, liney_prime, linez_prime, linexxy_prime, 
                              lineyyz_prime, linexxz_prime, linexyy_prime, lineyzz_prime, 
                              linexzz_prime, linexzxyz_prime, lineyzxyz_prime, linexyxyz_prime])

        return np.add(shift, all_lines)
    
    def centroid(self, polyhedra_ib):
        centroid_points = []
        for poly in polyhedra_ib:
            center = poly.center.coordinates
            centroid_points.append(center)
        kmeans = KMeans(n_clusters=1).fit(centroid_points)
        centroid = kmeans.cluster_centers_[0]
        return centroid
    
    def plot_unit_cell(self, ax, matrix, centroid, color='grey', alpha=0.5, linestyle='--'):
        
        shift = self.get_unit_cell_shift_centroid(matrix, centroid)
        lines = self.get_unit_cell_from_matrix(matrix, shift)
        
        for line in lines:
            X = line.T[0]
            Y = line.T[1]
            Z = line.T[2]
            ax.plot3D(X, Y, Z, linestyle=linestyle, color=color, alpha=alpha)
        
        return 
    
    def plot_label(self, ax, polyhedra_ib, size=20):
        for poly in polyhedra_ib:
            center = poly.center.coordinates
            specie = poly.center.specie
            ax.text(center[0], center[1], center[2],  
                    '$%s$' % str(specie.element), 
                    size=size, zorder=2, color='k', ha='center', va='center') 
        return 
    
    def plot_faces(self, ax, faces, color, alpha):
        for face in faces:
            coordinates = []
            for all_vertice in face.vertices:
                coords = all_vertice.coordinates
                coordinates.append(coords)
            coordinates = np.array(coordinates)
            plot_faces = [list(zip(coordinates.T[0], coordinates.T[1], coordinates.T[2]))]
            ax.add_collection3d(Poly3DCollection(plot_faces, facecolors=color, alpha=alpha))
        return 
    
    def plot_edges(self, ax, edges, color, alpha, linestyle='-'):
        for edge in edges:
            for i in range(len(edges)):
                poly_line = np.array([edges[i].point1.coordinates, edges[i].point2.coordinates])
                ax.plot3D(poly_line.T[0], poly_line.T[1], poly_line.T[2], color=color, alpha=alpha,
                         linestyle=linestyle)   
        return 
    
    def plot_polyhedrons(self, matrix, polyhedra_ib, polyhedra_oob, ax_range=10, colors=None, show_cell=True):
        if colors == None:
            colors = ['grey' for poly in polyhedra_ib]
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111, projection="3d")
        centroid = self.centroid(polyhedra_ib)
        if show_cell == True:
            self.plot_unit_cell(ax, matrix, centroid)
        else:
            pass
        
        all_shared_faces = []
        all_shared_edges = []
        
        for poly_ind, poly in enumerate(polyhedra_ib):
            #print(colors[poly_ind])
            self.plot_faces(ax, poly.faces, colors[poly_ind], 0.08)
            self.plot_edges(ax, poly.edges, 'grey', 0.1, '-')
            for poly2_ind, poly2 in enumerate(polyhedra_ib[poly_ind+1:]):
                shared_faces = poly.shared_faces(poly2, periodic=True, duplicate=True)
                shared_edges = poly.shared_edges(poly2, periodic=True, duplicate=True)
                if shared_faces != []:
                    all_shared_faces += list(shared_faces)
                if shared_edges != []:
                    all_shared_edges += list(shared_edges)
        if all_shared_faces == []:
            self.plot_edges(ax, all_shared_edges, 'red', 0.2, '--')
        else:
            self.plot_faces(ax, all_shared_faces, 'red', 0.5)
        self.plot_label(ax, polyhedra_ib, 15)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        #ax.set_box_aspect([1,1,1])
        ax.autoscale(enable=False, axis = 'both')  #you will need this line to change the Z-axis
        ax.set_xbound(centroid[0]-(ax_range/2), centroid[0]+(ax_range/2))
        ax.set_ybound(centroid[1]-(ax_range/2), centroid[1]+(ax_range/2))
        ax.set_zbound(centroid[2]-(ax_range/2), centroid[2]+(ax_range/2))

        # make the panes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.set_axis_off()
        
    def plot_structure(self, constructor_3d, structure, colors, species=None, inds=None, oob=False):
        all_pss = []
        all_colors = []
        oob_pss = []
        if species != None:
            for i, specie in enumerate(species):
                indices = [i for i in range(len(structure)) if structure[i].specie == specie]
                ps = []
                ps_colors = []
                for ind in indices:
                    p = constructor_3d.polyhedron_constructor(structure, ind)
                    ps.append(p)
                ps_colors.append(colors[i])

            all_pss += ps
            all_colors += ps_colors
        
        elif inds != None:
            ps = []
            ps_colors = []
            for ind_ind, ind in enumerate(inds):
                p = constructor_3d.polyhedron_constructor(structure, ind)
                ps.append(p)
                ps_colors.append(colors[ind_ind])
    
            all_pss += ps
            all_colors += ps_colors
        
        else:
            print('Can only specify one')
            sys.exit(1)
        
        self.plot_polyhedrons(structure.lattice.matrix, all_pss, oob_pss, ax_range=20, colors=all_colors)
        return 
    
class GraphNetworkPlotting():
    def __init__(self):
        return
    
    def divide_edges(self, graph):
        ### Divides edges into node-to-node and node-to-self, and returns their counts (for self loops)
        graph_edges = list(graph.edges)
        
        node_to_self_edge_sets = []
        node_to_self_edges = []
        
        node_to_node_edge_sets = []
        node_to_node_edges = []
        
        for edge in graph_edges:
            check_set = set([edge[0], edge[1]])
            if len(check_set) == 1: # node to self
                if check_set not in node_to_self_edge_sets: # Not currently counted
                    node_to_self_edge_sets.append(check_set)
                    node_to_self_edges.append([edge])
                else:
                    index = node_to_self_edge_sets.index(check_set)
                    node_to_self_edges[index] += [edge]
            else:
                if check_set not in node_to_node_edge_sets: # Not currently counted
                    node_to_node_edge_sets.append(check_set)
                    node_to_node_edges.append([edge])
                else:
                    index = node_to_node_edge_sets.index(check_set)
                    node_to_node_edges[index] += [edge]
        return node_to_node_edge_sets, node_to_node_edges, node_to_self_edge_sets, node_to_self_edges 
    
    def draw_nodes(self, graph, pos, colors, default=cm.tab20.colors):
        ### Draw the nodes first - will rework this to change shapes and colors 
        species = []
    
        for node in list(graph.nodes):
            specie = node.polyhedron.center.specie
            poly_type = node.polyhedron.poly_type
            if poly_type == 'octahedron':
                shape = 'D'
            else:
                shape = 'o'
            if specie not in species:
                species.append(specie)
            if colors is None:
                color = default[species.index(specie)]
                color = np.reshape(np.array(color), (-1, 3))
            else:
                color = colors[species.index(specie)]
            poly_type = node.polyhedron.poly_type
            nx.draw_networkx_nodes(graph, pos, node_size=500, node_shape=shape, nodelist=[node],
                                   node_color=color, zorder=1)
            
        nx.draw_networkx_labels(graph, pos, labels={node: node.polyhedron.center.index for node in graph.nodes()}, 
                                                    zorder=2)
        
    def edge_style(self, edge):
        if isinstance(edge[2].connector, Point3D):
            linestyle = 'dotted'
        elif isinstance(edge[2].connector, Edge3D):
            linestyle = 'dashed'
        elif isinstance(edge[2].connector, Face3D):
            linestyle = 'solid'
        else:
            linestyle = 'dashdot' # Should not be used here
        return linestyle
    
    def get_loop_centers(self, poly_center, n, scaling=0.015):
        if n == 1:
            numpy_vertice = np.add(poly_center, np.array([0, 1]))
            numpy_vertices = np.reshape(numpy_vertice, (-1, 2))
        elif n == 2:
            numpy_vertices = np.add(poly_center, np.array([[0, 1], [0, -1]]))
        else:
            polygon = Polygon(poly_center, 1, n=n)
            sympy_vertices = polygon.vertices
            numpy_vertices = np.array(sympy_vertices).astype(np.float64)
        return np.multiply(scaling, numpy_vertices)
    
    def draw_self_loop(self, center, radius, linestyle):
        ax = plt.gca()
        ring = mpatches.Wedge(center, radius, -180, 180)
        p = PatchCollection([ring], edgecolor='black', linestyle=linestyle, facecolor='none')
        ax.add_collection(p)
    
    def plot_graph(self, graph, colors=None):
        pos = nx.spring_layout(graph)
        plt.figure()
        ax = plt.gca()
        self.draw_nodes(graph, pos, colors) ### Draw the nodes; will rework this
        
        ### For the connections
        ntn_edge_sets, ntn_edges, nts_edge_sets, nts_edges = self.divide_edges(graph)
        
        ### For the node-to-self connections
        for es_ind, es in enumerate(nts_edges):
            node_pos = pos[list(nts_edge_sets[es_ind])[0]] ### get the position of the single node
            vertices = self.get_loop_centers(node_pos, len(es)) ### get the centers of the loops
            for e_ind, e in enumerate(es):
                linestyle = self.edge_style(e)
                radii = np.linalg.norm(np.array(node_pos)-np.array(vertices[e_ind]))
                self.draw_self_loop(vertices[e_ind], radii, linestyle)
                
        ### For the node-to-node connections
        ntn_counts = [1 if (len(ntn_edges[i]) % 2) == 0 else 0.5 for i in range(len(ntn_edge_sets))]
        ntn_signs = [1 for i in range(len(ntn_counts))]
        for es_ind, es in enumerate(ntn_edges):
            for e_ind, e in enumerate(es):
                linestyle = self.edge_style(e)
                edge_set = ntn_edge_sets[es_ind]
                set_index = ntn_edge_sets.index(edge_set) # works for repeats
                set_count = ntn_counts[set_index]
                
                if set_count == 0.5: # Draw direct line for odd number of total connections
                    rad = 0
                    ntn_counts[set_index] = 1
                else:
                    rad = str(0.1*ntn_counts[set_index]*ntn_signs[set_index])
                    if ntn_signs[set_index] == -1:
                        ntn_counts[set_index] += 1
                    ntn_signs[set_index] = ntn_signs[set_index] * -1
                        
                ax.annotate("", xy=pos[e[0]], xycoords='data', xytext=pos[e[1]], textcoords='data', zorder=0,
                            arrowprops=dict(arrowstyle="-", color="black", shrinkA=5, shrinkB=5, 
                            patchA=None, patchB=None, linestyle=linestyle, connectionstyle="arc3, rad=rrr".replace('rrr', rad)))
        
        ### Add legend patch based on identities of the nodes
        plt.axis('off')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.tight_layout()
        plt.show()
