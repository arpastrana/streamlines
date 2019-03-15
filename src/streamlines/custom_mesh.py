'''
A mesh wrapper that contains rs and compas information.
'''

__name__ = "Structural Mesh"
__author__ = "Rafael Pastrana"
__version__ = "0.0.4"
__creation__ = "2018.11.12"
__date__ = "2018.11.12"


import compas
import compas_rhino
import rhinoscriptsyntax as rs

from compas.datastructures import Mesh
from compas import geometry as cg
from compas_rhino import helpers
from math import degrees
from math import acos
from utilities import Utilities
from compas.geometry import KDTree
from compas.topology import dijkstra_distances


ut = Utilities()


class StructuralMesh():

    def __init__(self, gh_mesh, unify=False):
        self.c_mesh = helpers.mesh.mesh_from_guid(Mesh, gh_mesh)
        self.gh_mesh = gh_mesh

        self.adj = None
        self.e_weights = None
        self.avg_edge = None

        if unify:
            compas.topology.mesh_unify_cycles(self.c_mesh, root=None)

        self.set_face_planes()
        self.set_edge_planes()
        self.set_adjacency()
        self.set_avg_edge_length()
        # self.set_dijkstra_distances()

    def set_avg_edge_length(self):
        lengths = [self.c_mesh.edge_length(u, v) for u, v in self.c_mesh.edges()]
        self.avg_edge = sum(lengths) / len(lengths)

    def set_adjacency(self):
        adjacency = {key: self.c_mesh.vertex_neighbors(key) for key in self.c_mesh.vertices()}
        weight = {(u, v): self.c_mesh.edge_length(u, v) for u, v in self.c_mesh.edges()}
        weight.update({(v, u): weight[(u, v)] for u, v in self.c_mesh.edges()})

        self.adj = adjacency
        self.e_weights = weight

    def set_dijkstra_distances(self):
        for v in self.c_mesh.vertices():
            dijkstra_d = dijkstra_distances(self.adj, self.e_weights, v)
            self.c_mesh.set_vertex_attribute(v, 'dijkstra', dijkstra_d)

    def set_face_vectors(self, vectorField, name, normalize=True):
        for f_key in self.c_mesh.faces():

            vector_a = vectorField[int(f_key)]
            if normalize is True:
                vector_a = cg.normalize_vector(vector_a)  # normalize or not?
            vector_b = cg.scale_vector(vector_a, -1.)
            self.c_mesh.set_face_attribute(f_key, str(name) + '_a', vector_a)
            self.c_mesh.set_face_attribute(f_key, str(name) + '_b', vector_b)

    def set_face_planes(self):
        for f_key in self.c_mesh.faces():
            ct = self.c_mesh.face_centroid(f_key)
            normal = self.c_mesh.face_normal(f_key)
            self.c_mesh.set_face_attribute(f_key, 'plane', [ct, normal])

    def set_edge_planes(self):
        for u, v in self.c_mesh.edges(data=False):

            f_keys = self.c_mesh.edge_faces(u, v)
            f_keys = [f_key for f_key in f_keys if f_key is not None]

            normals = [self.c_mesh.face_normal(f_key) for f_key in f_keys]
            normals = list(map(lambda x: ut.align_vector(x, normals[0]),
                               normals
                               )
                           )

            vec = reduce(lambda x, y: cg.add_vectors(x, y), normals)

            vec = cg.cross_vectors(self.c_mesh.edge_vector(u, v), vec)
            plane = [self.c_mesh.edge_midpoint(u, v), vec]
            self.c_mesh.set_edge_attribute((u, v), 'plane', plane)

    def set_vertex_vectors(self, name):
        for v in self.c_mesh.vertices():
            # 1. find keys of vertex' neighbouring faces
            face_idxs = self.c_mesh.vertex_faces(v, ordered=True)

            # 2. get attributes of the found faces and calculate weights
            vec_a = self.c_mesh.get_faces_attribute(name=str(name) + '_a',
                                                   keys=face_idxs
                                                   )
            vec_b = self.c_mesh.get_faces_attribute(name=str(name) + '_b',
                                                   keys=face_idxs
                                                   )

            # 3. calculate weights
            pt_cloud = [self.c_mesh.face_centroid(k) for k in face_idxs]
            weights = ut.get_dist_weights(self.c_mesh.vertex_coordinates(v),
                                          pt_cloud
                                          )

            # 4. multiply vectors by weights
            nd_vec_a = ut.vectors_weight_reduce(vec_a, weights)
            nd_vec_b = ut.vectors_weight_reduce(vec_b, weights)

            # 6. assign vector attributes to vertices
            self.c_mesh.set_vertex_attribute(v, str(name) + '_a', nd_vec_a)
            self.c_mesh.set_vertex_attribute(v, str(name) + '_b', nd_vec_b)

    def set_vertex_vectors_angles(self, name):
        for v in self.c_mesh.vertices():
            # 1. find keys of vertex' neighbouring faces
            face_idxs = self.c_mesh.vertex_faces(v, ordered=True)

            # 2. get attributes of the found faces and calculate weights
            vec_a = self.c_mesh.get_faces_attribute(name=str(name) + '_a',
                                                   keys=face_idxs
                                                   )
            vec_b = self.c_mesh.get_faces_attribute(name=str(name) + '_b',
                                                   keys=face_idxs
                                                   )

            # 3. get connected edges
            angles = []
            for face_idx in face_idxs:
                edges = []
                for edge in self.c_mesh.face_halfedges(face_idx):
                    if v in edge:
                        edges.append(edge)
                # if len(edges) > 1:
                e_1 = self.c_mesh.edge_vector(edges[0][0], edges[0][1])
                e_2 = self.c_mesh.edge_vector(edges[1][0], edges[1][1])
                angle = cg.angle_vectors(e_1, e_2, deg=True)
                angles.append(angle)

            mass_angle = reduce(lambda x, y: x+y, angles)
            weights = list(map(lambda x: x/mass_angle, angles))

            # 4. multiply vectors by weights
            nd_vec_a = ut.vectors_weight_reduce(vec_a, weights)
            nd_vec_b = ut.vectors_weight_reduce(vec_b, weights)

            # 6. assign vector attributes to vertices
            self.c_mesh.set_vertex_attribute(v, str(name) + '_a', nd_vec_a)
            self.c_mesh.set_vertex_attribute(v, str(name) + '_b', nd_vec_b)

    def get_vector_on_face(self, point, f_key, name, vec=[0, 0, 0]):
        '''
        Inputs:
        - point: compas point
        - f_key: mesh face key
        - name: vector field name, as stored in custom mesh

        Optional:
        - vec: reference vector for alignment
        '''

        # getters
        v_keys = self.c_mesh.face_vertices(f_key)
        pt_cloud = self.c_mesh.face_coordinates(f_key)
        v_vectors_a = self.c_mesh.get_vertices_attribute(str(name) + '_a',
                                                        keys=v_keys
                                                        )
        # v_vectors_b = self.c_mesh.get_vertices_attribute(str(name) + '_b',
        #                                                 keys=v_keys
        #                                                 )

        # v_vectors = ut.filter_aligned_vectors(vec, v_vectors_a, v_vectors_b)

        v_vectors = list(map(lambda x: ut.align_vector(x, vec), v_vectors_a))
        weights = ut.get_dist_weights(point, pt_cloud)
        new_vectors = list(map(lambda x, y: x * y, v_vectors, weights))

        # new_vectors = []
        # for idx, vec in enumerate(v_vectors):
        #     new_vector = cg.scale_vector(vec, weights[idx])
        #     new_vectors.append(new_vector)

        return reduce(lambda x, y: cg.add_vectors(x, y), new_vectors)

    def get_vector_on_face_ext(self, point, f_key, name, vec=[0, 0, 0]):
        f_keys = [f_key]
        f_keys.extend(self.c_mesh.face_neighbours(f_key))
        v_keys = []
        pt_cloud = []

        for f_key in f_keys:
            v_keys.extend(self.c_mesh.face_vertices(f_key))

        v_keys = list(set(v_keys))
        v_vectors_a = self.c_mesh.get_vertices_attribute(str(name) + '_a',
                                                        keys=v_keys
                                                        )
        v_vectors_b = self.c_mesh.get_vertices_attribute(str(name) + '_b',
                                                        keys=v_keys
                                                        )

        v_vectors = ut.filter_aligned_vectors(vec, v_vectors_a, v_vectors_b)

        for v_key in v_keys:
            pt_cloud.append(self.c_mesh.vertex_coordinates(v_key))

        weights = ut.get_dist_weights(point, pt_cloud)
        new_vectors = []
        for idx, vec in enumerate(v_vectors):
            new_vector = cg.scale_vector(vec, weights[idx])
            new_vectors.append(new_vector)

        return reduce(lambda x, y: cg.add_vectors(x, y), new_vectors)

    def get_vector_on_face_faces(self, point, f_key, name, vec=[0, 0, 0]):
        pt_cloud = []

        f_keys = [f_key]
        f_keys.extend(self.c_mesh.face_neighbors(f_key))
        f_keys = list(set(f_keys))

        f_vectors_a = self.c_mesh.get_faces_attribute(f_keys,
                                                      str(name) + '_a'
                                                      )
        f_vectors_b = self.c_mesh.get_faces_attribute(f_keys,
                                                      str(name) + '_b',
                                                      )
        f_vectors = ut.filter_aligned_vectors(vec, f_vectors_a, f_vectors_b)

        for f_key in f_keys:
            pt_cloud.append(self.c_mesh.face_centroid(f_key))

        weights = ut.get_dist_weights(point, pt_cloud)
        new_vectors = []
        for idx, vec in enumerate(f_vectors):
            new_vector = cg.scale_vector(vec, weights[idx])
            new_vectors.append(new_vector)

        return reduce(lambda x, y: cg.add_vectors(x, y), new_vectors)

    def get_data_on_nodes(self, data_tag, mode='max'):
        n_data = []

        for v in self.c_mesh.vertices():
            faces = self.c_mesh.vertex_faces(v)
            f_data = self.c_mesh.get_faces_attribute(name=data_tag, keys=faces)

            if mode == 'max':
                val = max(f_data)
            elif mode == 'min':
                val = min(f_data)
            elif mode == 'avg':
                val = sum(f_data) / len(f_data)
            elif mode == 'angle_avg':
                weights = ut.get_angle_weights(self.c_mesh, v, faces)
                val = ut.values_weight_reduce(f_data, weights)

            n_data.append(val)
        return n_data

    def get_edge_labels(self, name, tol):
        for index, (u, v, attr) in enumerate(self.c_mesh.edges(True)):
            vec_u_a = self.c_mesh.get_vertex_attribute(u, str(name) + '_a')
            vec_v_a = self.c_mesh.get_vertex_attribute(v, str(name) + '_a')

            dot = cg.dot_vectors(vec_u_a, vec_v_a)

            if dot < 0 - tol:
                label = -1
            elif dot > 0 + tol:
                label = 1
            else:
                label = 0

            self.c_mesh.set_edge_attribute((u, v), 'label', label)

    def get_face_labels(self, name, tol=0.6):
        self.get_edge_labels(name, tol)

        for f_key in self.c_mesh.faces():
            valency = 1
            for edge in self.c_mesh.face_halfedges(f_key):
                valency *= self.c_mesh.get_edge_attribute(edge, 'label')
            self.c_mesh.set_face_attribute(f_key, 'label', valency)

    def find_vector_on_face(self, point, f_key, name):
        v_keys = self.c_mesh.face_vertices(f_key)
        vectors_a = self.c_mesh.get_vertices_attribute(str(name) + '_a',
                                                      keys=v_keys
                                                      )
        vectors_b = self.c_mesh.get_vertices_attribute(str(name) + '_b',
                                                      keys=v_keys
                                                      )

        weights = ut.get_dist_weights(point,
                                      self.c_mesh.face_coordinates(f_key)
                                      )

        vec_a = ut.vectors_weight_reduce(vectors_a, weights)
        vec_b = ut.vectors_weight_reduce(vectors_b, weights)
        return vec_a, vec_b

    def find_vector_face_centers(self, name):
        vectors = []
        for f_key in self.c_mesh.faces():
            ct = self.c_mesh.face_centroid(f_key)
            vec_a, vec_b = self.find_vector_on_face(ct, f_key, name)
            vectors.append((vec_a, vec_b))
        return vectors

    def stream_face_labels(self):
        l_neg = []
        l_null = []
        l_pos = []

        for f_key in self.c_mesh.faces():
            label = self.c_mesh.get_face_attribute(f_key, 'label')

            if label == -1:
                l_neg.append(self.c_mesh.face_centroid(f_key))
            elif label == 0:
                l_null.append(self.c_mesh.face_centroid(f_key))
            elif label == 1:
                l_pos.append(self.c_mesh.face_centroid(f_key))

        l_neg = list(map(lambda x: rs.AddPoint(x), l_neg))
        l_null = list(map(lambda x: rs.AddPoint(x), l_null))
        l_pos = list(map(lambda x: rs.AddPoint(x), l_pos))

        return l_neg, l_null, l_pos

    def closest_point(self, point, maxdist=None):
        maxdist = maxdist
        point, face = rs.MeshClosestPoint(self.gh_mesh,
                                          rs.AddPoint(*point),
                                          maxdist
                                          )
        point = [point.X, point.Y, point.Z]
        return point, face
