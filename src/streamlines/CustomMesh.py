'''
A mesh wrapper that contains rs and compas information.
'''


__name__ = "Structural Mesh"
__author__ = "Rafael Pastrana"
__version__ = "0.0.3"
__creation__ = "2018.08.03"
__date__ = "2018.08.03"


import compas
import compas_rhino
import rhinoscriptsyntax as rs

from compas.datastructures import Mesh
from compas import geometry as cg
from math import degrees
from math import acos
from Utilities import Utilities
from compas.geometry import KDTree
from compas.topology import dijkstra_distances


ut = Utilities()


class StructuralMesh():

    def __init__(self, ghMesh, unify=False):
        self.cMesh = compas_rhino.helpers.mesh_from_guid(Mesh, ghMesh)
        self.ghMesh = ghMesh

        self.adj = None
        self.e_weights = None
        self.avg_edge = None

        if unify:
            compas.topology.mesh_unify_cycles(self.cMesh, root=None)

        self.set_face_planes()
        self.set_edge_planes()
        self.set_adjacency()
        self.set_avg_edge_length()
        # self.set_dijkstra_distances()

    def set_avg_edge_length(self):
        lengths = [self.cMesh.edge_length(u, v) for u, v in self.cMesh.edges()]
        self.avg_edge = sum(lengths) / len(lengths)

    def set_adjacency(self):
        adjacency = {key: self.cMesh.vertex_neighbours(key) for key in self.cMesh.vertices()}
        weight = {(u, v): self.cMesh.edge_length(u, v) for u, v in self.cMesh.edges()}
        weight.update({(v, u): weight[(u, v)] for u, v in self.cMesh.edges()})

        self.adj = adjacency
        self.e_weights = weight

    def set_dijkstra_distances(self):
        for v in self.cMesh.vertices():
            dijkstra_d = dijkstra_distances(self.adj, self.e_weights, v)
            self.cMesh.set_vertex_attribute(v, 'dijkstra', dijkstra_d)

    def set_face_vectors(self, vectorField, name, normalize=True):
        for f_key in self.cMesh.faces():

            vector_a = vectorField[int(f_key)]
            if normalize is True:
                vector_a = cg.normalize_vector(vector_a)  # normalize or not?
            vector_b = cg.scale_vector(vector_a, -1.)
            self.cMesh.set_face_attribute(f_key, str(name) + '_a', vector_a)
            self.cMesh.set_face_attribute(f_key, str(name) + '_b', vector_b)

    def set_face_planes(self):
        for f_key in self.cMesh.faces():
            ct = self.cMesh.face_centroid(f_key)
            normal = self.cMesh.face_normal(f_key)
            self.cMesh.set_face_attribute(f_key, 'plane', [ct, normal])

    def set_edge_planes(self):
        for u, v in self.cMesh.edges(data=False):

            f_keys = self.cMesh.edge_faces(u, v)
            f_keys = [f_key for f_key in f_keys if f_key is not None]

            normals = [self.cMesh.face_normal(f_key) for f_key in f_keys]
            normals = list(map(lambda x: ut.align_vector(x, normals[0]),
                               normals
                               )
                           )

            vec = reduce(lambda x, y: cg.add_vectors(x, y), normals)

            vec = cg.cross_vectors(self.cMesh.edge_vector(u, v), vec)
            plane = [self.cMesh.edge_midpoint(u, v), vec]
            self.cMesh.set_edge_attribute((u, v), 'plane', plane)

    def set_vertex_vectors(self, name):
        for v in self.cMesh.vertices():
            # 1. find fkeys of vertex' neighbouring faces
            face_idxs = self.cMesh.vertex_faces(v, ordered=True)

            # 2. get attributes of the found faces and calculate weights
            vec_a = self.cMesh.get_faces_attribute(str(name) + '_a',
                                                   fkeys=face_idxs
                                                   )
            vec_b = self.cMesh.get_faces_attribute(str(name) + '_b',
                                                   fkeys=face_idxs
                                                   )

            # 3. calculate weights
            pt_cloud = [self.cMesh.face_centroid(k) for k in face_idxs]
            weights = ut.get_dist_weights(self.cMesh.vertex_coordinates(v),
                                          pt_cloud
                                          )

            # 4. multiply vectors by weights
            nd_vec_a = ut.vectors_weight_reduce(vec_a, weights)
            nd_vec_b = ut.vectors_weight_reduce(vec_b, weights)

            # 6. assign vector attributes to vertices
            self.cMesh.set_vertex_attribute(v, str(name) + '_a', nd_vec_a)
            self.cMesh.set_vertex_attribute(v, str(name) + '_b', nd_vec_b)

    def set_vertex_vectors_angles(self, name):
        for v in self.cMesh.vertices():
            # 1. find fkeys of vertex' neighbouring faces
            face_idxs = self.cMesh.vertex_faces(v, ordered=True)

            # 2. get attributes of the found faces and calculate weights
            vec_a = self.cMesh.get_faces_attribute(str(name) + '_a',
                                                   fkeys=face_idxs
                                                   )
            vec_b = self.cMesh.get_faces_attribute(str(name) + '_b',
                                                   fkeys=face_idxs
                                                   )

            # 3. get connected edges
            angles = []
            for face_idx in face_idxs:
                edges = []
                for edge in self.cMesh.face_halfedges(face_idx):
                    if v in edge:
                        edges.append(edge)
                # if len(edges) > 1:
                e_1 = self.cMesh.edge_vector(edges[0][0], edges[0][1])
                e_2 = self.cMesh.edge_vector(edges[1][0], edges[1][1])
                angle = cg.angle_vectors(e_1, e_2, deg=True)
                angles.append(angle)

            mass_angle = reduce(lambda x, y: x+y, angles)
            weights = list(map(lambda x: x/mass_angle, angles))

            # 4. multiply vectors by weights
            nd_vec_a = ut.vectors_weight_reduce(vec_a, weights)
            nd_vec_b = ut.vectors_weight_reduce(vec_b, weights)

            # 6. assign vector attributes to vertices
            self.cMesh.set_vertex_attribute(v, str(name) + '_a', nd_vec_a)
            self.cMesh.set_vertex_attribute(v, str(name) + '_b', nd_vec_b)

    def get_vector_on_face(self, point, f_key, name, vec=[0, 0, 0]):
        v_keys = self.cMesh.face_vertices(f_key)
        v_vectors_a = self.cMesh.get_vertices_attribute(str(name) + '_a',
                                                        keys=v_keys
                                                        )
        v_vectors_b = self.cMesh.get_vertices_attribute(str(name) + '_b',
                                                        keys=v_keys
                                                        )

        v_vectors = ut.filter_aligned_vectors(vec, v_vectors_a, v_vectors_b)
        pt_cloud = self.cMesh.face_coordinates(f_key)
        weights = ut.get_dist_weights(point, pt_cloud)

        new_vectors = []
        for idx, vec in enumerate(v_vectors):
            new_vector = cg.scale_vector(vec, weights[idx])
            new_vectors.append(new_vector)
        return reduce(lambda x, y: cg.add_vectors(x, y), new_vectors)

    def get_vector_on_face_ext(self, point, f_key, name, vec=[0, 0, 0]):
        f_keys = [f_key]
        f_keys.extend(self.cMesh.face_neighbours(f_key))
        v_keys = []
        pt_cloud = []

        for f_key in f_keys:
            v_keys.extend(self.cMesh.face_vertices(f_key))

        v_keys = list(set(v_keys))
        v_vectors_a = self.cMesh.get_vertices_attribute(str(name) + '_a',
                                                        keys=v_keys
                                                        )
        v_vectors_b = self.cMesh.get_vertices_attribute(str(name) + '_b',
                                                        keys=v_keys
                                                        )

        v_vectors = ut.filter_aligned_vectors(vec, v_vectors_a, v_vectors_b)

        for v_key in v_keys:
            pt_cloud.append(self.cMesh.vertex_coordinates(v_key))

        weights = ut.get_dist_weights(point, pt_cloud)
        new_vectors = []
        for idx, vec in enumerate(v_vectors):
            new_vector = cg.scale_vector(vec, weights[idx])
            new_vectors.append(new_vector)

        return reduce(lambda x, y: cg.add_vectors(x, y), new_vectors)

    def get_data_on_nodes(self, data_tag, mode='max'):
        n_data = []

        for v in self.cMesh.vertices():
            faces = self.cMesh.vertex_faces(v)
            f_data = self.cMesh.get_faces_attribute(data_tag, fkeys=faces)

            if mode == 'max':
                val = max(f_data)
            elif mode == 'min':
                val = min(f_data)
            elif mode == 'avg':
                val = sum(f_data) / len(f_data)
            elif mode == 'angle_avg':
                weights = ut.get_angle_weights(self.cMesh, v, faces)
                val = ut.values_weight_reduce(f_data, weights)

            n_data.append(val)
        return n_data

    def get_edge_labels(self, name, tol):
        for index, (u, v, attr) in enumerate(self.cMesh.edges(True)):
            vec_u_a = self.cMesh.get_vertex_attribute(u, str(name) + '_a')
            vec_v_a = self.cMesh.get_vertex_attribute(v, str(name) + '_a')

            dot = cg.dot_vectors(vec_u_a, vec_v_a)

            if dot < 0 - tol:
                label = -1
            elif dot > 0 + tol:
                label = 1
            else:
                label = 0

            self.cMesh.set_edge_attribute((u, v), 'label', label)

    def get_face_labels(self, name, tol=0.6):
        self.get_edge_labels(name, tol)

        for f_key in self.cMesh.faces():
            valency = 1
            for edge in self.cMesh.face_halfedges(f_key):
                valency *= self.cMesh.get_edge_attribute(edge, 'label')
            self.cMesh.set_face_attribute(f_key, 'label', valency)

    def find_vector_on_face(self, point, f_key, name):
        v_keys = self.cMesh.face_vertices(f_key)
        vectors_a = self.cMesh.get_vertices_attribute(str(name) + '_a',
                                                      keys=v_keys
                                                      )
        vectors_b = self.cMesh.get_vertices_attribute(str(name) + '_b',
                                                      keys=v_keys
                                                      )

        weights = ut.get_dist_weights(point,
                                      self.cMesh.face_coordinates(f_key)
                                      )

        vec_a = ut.vectors_weight_reduce(vectors_a, weights)
        vec_b = ut.vectors_weight_reduce(vectors_b, weights)
        return vec_a, vec_b

    def find_vector_face_centers(self, name):
        vectors = []
        for f_key in self.cMesh.faces():
            ct = self.cMesh.face_centroid(f_key)
            vec_a, vec_b = self.find_vector_on_face(ct, f_key, name)
            vectors.append((vec_a, vec_b))
        return vectors

    def stream_face_labels(self):
        l_neg = []
        l_null = []
        l_pos = []

        for f_key in self.cMesh.faces():
            label = self.cMesh.get_face_attribute(f_key, 'label')

            if label == -1:
                l_neg.append(self.cMesh.face_centroid(f_key))
            elif label == 0:
                l_null.append(self.cMesh.face_centroid(f_key))
            elif label == 1:
                l_pos.append(self.cMesh.face_centroid(f_key))

        l_neg = list(map(lambda x: rs.AddPoint(x), l_neg))
        l_null = list(map(lambda x: rs.AddPoint(x), l_null))
        l_pos = list(map(lambda x: rs.AddPoint(x), l_pos))

        return l_neg, l_null, l_pos

    def closest_point(self, point, maxdist=None):
        maxdist = maxdist
        point, face = rs.MeshClosestPoint(self.ghMesh,
                                          rs.AddPoint(*point),
                                          maxdist
                                          )
        point = [point.X, point.Y, point.Z]
        return point, face
