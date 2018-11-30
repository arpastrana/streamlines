'''
A set of functions for putting everyday tasks at ease.
'''


__name__ = "Utilities"
__author__ = "Rafael Pastrana"
__version__ = "0.0.2"
__creation__ = "2018.07.24"
__date__ = "2018.07.27"


import compas
import compas_rhino
import compas.geometry as cg
import math
import rhinoscriptsyntax as rs

from compas.datastructures import Mesh
from compas.geometry import Line
from compas.geometry import KDTree
from compas.topology import delaunay_from_points


class Utilities():

    def __init__(self):
        pass

    def is_point_close(self, point, objs, distance, exclude=None, d_out=False):
        try:
            tree = KDTree(objs)
            nnbr, label, dst = tree.nearest_neighbour(point, exclude)
            if dst is not None:
                if dst <= distance:
                    # print('distance is {}'.format(dst))
                    if d_out is True:
                        return (True, dst)
                    else:
                        return True
                else:
                    if d_out is True:
                        return (False, dst)
                    else:
                        return False
        except Exception:
            return None
        return None

    def align_vector(self, vector_a, vector_b):
        if cg.dot_vectors(vector_a, vector_b) < 0:
            return cg.scale_vector(vector_a, -1.)
        else:
            return vector_a

    def get_dist_weights(self, point, point_cloud):
        dst = []
        weigths = []

        for c_point in point_cloud:
            dst.append(cg.distance_point_point(point, c_point))

        if len(dst) > 1:
            massDistance = reduce(lambda x, y: x+y, dst)
            return list(map(lambda x: x/massDistance, dst))

        else:
            return [1.0]

    def get_angle_weights(self, c_mesh, v_key, f_keys):
        angles = []
        weigths = []

        angles = []
        for f_key in f_keys:
            edges = []
            for edge in c_mesh.face_halfedges(f_key):
                if v_key in edge:
                    edges.append(edge)

            e_1 = c_mesh.edge_vector(edges[0][0], edges[0][1])
            e_2 = c_mesh.edge_vector(edges[1][0], edges[1][1])

            try:
                angle = cg.angle_vectors(e_1, e_2, deg=True)
            except Exception:
                angle = 1e-6
            angles.append(angle)

        if len(angles) > 1:
            mass_angle = reduce(lambda x, y: x+y, angles)
            if mass_angle > 0:
                return list(map(lambda x: x/mass_angle, angles))
        return [1.0]

    def vectors_weight_reduce(self, vectors, weights):
        newVectors = []
        try:
            guideVec = vectors[1]
        except Exception:
            guideVec = vectors[0]

        for idx, vec in enumerate(vectors):
            vec = self.align_vector(vec, guideVec)
            newVector = cg.scale_vector(vec,
                                        weights[idx]
                                        )
            newVectors.append(newVector)
        return reduce(lambda x, y: cg.add_vectors(x, y), newVectors)

    def values_weight_reduce(self, values, weights):
        new_values = [val * weights[idx] for idx, val in enumerate(values)]
        return reduce(lambda x, y: x + y, new_values)

    def angle_vectors(self, u, v, deg=True):
        a = cg.dot_vectors(u, v) / (cg.length_vector(u) * cg.length_vector(v))
        if deg:
            return math.degrees(math.acos(a))
        return math.acos(a)

    def filter_aligned_vectors(self, vec, vecs_a, vecs_b):
        vec_tuples = zip(vecs_a, vecs_b)
        aligned_vecs = []

        for vec_t in vec_tuples:
            if cg.dot_vectors(vec, vec_t[0]) < 0:
                aligned_vecs.append(vec_t[1])

            else:
                aligned_vecs.append(vec_t[0])
        return aligned_vecs

    def delaunay_is_point_close(self, pt, points, distance):
        if pt in points:
            idx = points.index(pt)
            try:
                faces = delaunay_from_points(points)
            except Exception:
                return None

            if len(faces) > 0:
                triplets = [face for face in faces if idx in face]

                if len(triplets) > 0:
                    for triplet in triplets:
                        try:
                            circle = cg.circle_from_points(points[triplet[0]],
                                                           points[triplet[1]],
                                                           points[triplet[2]]
                                                           )
                            if circle:
                                if circle[1] * 2 < distance:
                                    print('Distance Delaunay is {}'.format(circle[1]))
                                    return True
                        except Exception:
                            return None
        return None

    def limit_vector(self, vector, limit):
        if cg.length_vector(vector) > limit:
            return cg.scale_vector(cg.normalize_vector(vector), limit)
        else:
            return vector

    def get_points_in_distance(self, point, points, distance):
        pass

    def normalise_data(self, data, cmin=None, cmax=None):
        fmax = cmax if cmax is not None else max([abs(x) for x in data])
        fmin = cmin if cmin is not None else min([abs(x) for x in data])
        fabs = max([abs(fmin), abs(fmax)])

        fscaled = []
        for value in data:
            new_val = value / fabs
            if new_val > 1:
                new_val = 1
            elif new_val < -1:
                new_val = -1
            fscaled.append(new_val)

        return fscaled, fabs

    def colorbar(self, data, type=255):

        colors = []
        for fsc in data:
            r = +abs(fsc + 0.25) * 2 - 0.5
            g = -abs(fsc - 0.25) * 2 + 1.5
            b = -(fsc - 0.25) * 2

            r = max([0, min([1, r])])
            g = max([0, min([1, g])])
            b = max([0, min([1, b])])

            colors.append(rs.CreateColor([i * type for i in [r, g, b]]))
        return colors

    def reconstruct_mesh_from_faces(self, mesh, fkeys):
        active_keys = {}
        new_faces = {}
        count = 0

        # get new faces
        for fkey in fkeys:
            new_face = []
            for vkey in mesh.face_vertices(fkey):

                if vkey not in active_keys:
                    active_keys[vkey] = count
                    count += 1

                new_face.append(active_keys.get(vkey))
            new_faces[fkey] = new_face

        # get new vertices
        new_vertices = []
        for a_key, idx in active_keys.items():
            new_vertices.append((mesh.vertex_coordinates(a_key), idx))
        new_vertices = sorted(new_vertices, key=lambda x: x[1])
        new_vertices = [x[0] for x in new_vertices]

        # retrieve face list
        new_faces = [face for fkey, face in new_faces.items()]

        # create new mesh
        new_mesh = Mesh.from_vertices_and_faces(new_vertices, new_faces)
        return new_mesh

    def make_gh_mesh(self, c_mesh):
        rs_points = []
        for v in c_mesh.vertices():
            rs_points.append(rg.Point3d(*c_mesh.vertex_coordinates(v)))

        faces = []
        for idx, f in enumerate(c_mesh.faces()):
            f_vts = c_mesh.face_vertices(f)
            if len(f_vts) < 4:
                f_vts.append(f_vts[-1])
            faces.append(f_vts)

        new_gh_mesh = rs.AddMesh(rs_points, faces, vertex_colors=None)
        return new_gh_mesh

    def slice_list_flags(self, my_list, flags):
        if all(flags) is True:
            print('all true')
            return [my_list]

        indices = [idx for idx, fl in enumerate(flags) if fl is False]
        indices.append(list(range(len(my_list)))[-1])
        indices = set(indices)
        indices = sorted(list(indices))
        print('slicing indices are {}'.format(indices))

        start = 0
        new_list = []
        for count, index in enumerate(indices):

            # if index == list(range(len(my_list)))[-2]:
            #     index += 1

            slice = my_list[start:index]
            new_list.append(slice)
            start = index + 1

        new_list = [slice for slice in new_list if len(slice) > 0]
        return new_list

    def intersect_(self):
        pass
