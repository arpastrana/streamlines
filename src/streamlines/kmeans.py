#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
A routine that implements a K-means routine to cluster similar vectors.

Based on the work of Variational Shape Approximation by Cohen-Steiner and on
the python script developed by Jesus Galvez.

The clustering routing non-overlapping regions according to a specified metric.

I. Algorithm:

1. Make a Mesh datastructure with the following properties embedded:
     A. Face Centers
     B. Face Normals
     C. Vertex Normals
     D. Face Count
     E. Face Indexes
     F. Vertices
     G. Area per Face
     H. Edge Adjacency
     I. Face Adjacency

2. Compute Weighed Normals (Face normal * face area)
3. Create Initial Seeds (if not specified) with *create seeds*
4. Run K-Means Method
'''

__author__ = 'Rafael Pastrana'
__name__ = 'K-Means Clustering'
__version__ = '0.0.5'
__date__ = '17.04.20'


import heapq
import itertools
import math
import random
import time

import numpy as np
import compas.geometry as cg

from functools import total_ordering
from functools import reduce

from streamlines.utilities import Utilities


ut = Utilities()


def k_means(clusters, faces, iters, mergesplit=False, callback=None):
    '''
    2. Run for N given iterations only.
    1. Create Cluster Colors. Make global *counter*.
    2B. Create Proxys with *get_proxy*
    3. Test whether it is the 1st iteration or not with global *counter*.
    4. If 1st, get proxies from clusters through from *get_proxy_seed*
    4B. If not, proxies are the regions. Or the other way around.
    5. Build a queue with the seeds' adjacent faces through *build_queue*
    6. Grow up a cluster with *assign_to_region* method.
    7. Create new proxies from created regions with *grow_seeds* method.
    8. New proxies become the proxies.
    9. Repeat
    '''
    all_clusters = []
    print('seeds are: {}'.format([cl.seed for cl_key, cl in clusters.items()]))

    for it in range(iters):
        print('#####')
        s1 = time.time()

        new_clusters, new_faces = get_new_clusters(clusters, faces)

        q = Queue(new_clusters, new_faces)
        q.init_queue()
        q.assign()
        clusters = q.get_clusters()
        all_clusters.append(output_cls(clusters))

        if mergesplit is True:
            if it < iters-1:
                clusters = merge_split(clusters)

        if callback:
            callback(k=it, clusters=clusters)

        e1 = time.time()
        print('Iteration {} execution time was {}'.format(it, (e1-s1)))
    print('End Clusters are: {}'.format(clusters))
    return all_clusters


def furthest_init(num, faces, callback=None):  # no dependency
    print('#####')
    print('Furthest init started')
    s0 = time.time()
    clusters = {0: Cluster(0, 0)}
    all_clusters = []

    for i in range(num):
        new_clusters, new_faces = get_new_clusters(clusters, faces)

        q = Queue(new_clusters, new_faces)
        q.init_queue()
        q.assign()
        clusters = q.get_clusters()
        all_clusters.append(output_cls(clusters))

        if i < num-1:
            t_s = get_cluster_to_split(clusters)
            clusters = split_cluster(t_s, clusters)

        if callback:
            callback(k=i, clusters=clusters)

        e0 = time.time()

    print('Furthest init execution time was {}'.format(e0-s0))
    return all_clusters


def output_cls(clusters):
    new_clusters = {}
    for c_key, cluster in clusters.items():
        new_cl = Cluster(cluster.id, cluster.seed)
        new_cl.copy_cluster(cluster)
        new_clusters[c_key] = new_cl
    return new_clusters


def get_new_clusters(clusters, faces):
    n_clusters = {}

    for key, cluster in clusters.items():
        cluster.harvest_faces(faces)
        cluster.set_vector_proxy()

        n_cluster = Cluster(cluster.id, cluster.get_new_seed())
        n_clusters[n_cluster.id] = n_cluster
        cluster.clear_faces()  # clears out cluster and face relationship

    return n_clusters, faces


def clear_clusters(faces):
    for f_key, face in faces.items():
        face.clear_cluster()


def get_random_seeds(maximum, num):
    return random.sample(range(0, maximum), num)


def get_random_color(num):
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)
    return (r, g, b)


def get_random_colors(values):
    return list(map(get_random_color, values))


def make_faces(s_mesh, tag, weight=False):  # no dep
    faces = {}

    for f_key in s_mesh.cMesh.faces():
        face = Face(f_key)
        halfedges = s_mesh.cMesh.face_halfedges(f_key)
        vector = s_mesh.cMesh.get_face_attribute(f_key, tag)

        face.set_halfedges([tuple(sorted(h)) for h in halfedges])
        face.set_vertices(s_mesh.cMesh.face_vertices(f_key))
        face.set_area(s_mesh.cMesh.face_area(f_key))
        face.set_neighbours(s_mesh.cMesh.face_neighbors(f_key))

        face.set_vector(vector)
        face.set_weighed_vector(weight)

        faces[f_key] = face

    return faces


def cluster_adjacency(clusters):
    comb = itertools.combinations(range(len(clusters)), 2)
    cl_comb = [(clusters.get(x[0]), clusters.get(x[1])) for x in comb]
    return list(filter(is_adjacent, cl_comb))


def is_adjacent(cluster_pair):
    vert_1 = cluster_pair[0].get_faces_halfedges()
    vert_2 = cluster_pair[1].get_faces_halfedges()
    return len(vert_1.intersection(vert_2)) > 0


def get_clusters_to_merge(adj_clusters):
    return min(adj_clusters, key=lambda x: simulate_merge(x[0], x[1]))


def get_cluster_to_split(clusters):
    return max(clusters.items(), key=lambda x: x[1].get_distortion())[1]


def simulate_merge(cluster_1, cluster_2):
    t_faces = set(cluster_1.faces + cluster_2.faces)
    errors = get_errors(t_faces, get_proxy(t_faces))
    return get_distortion(errors)


def merge_clusters(t_m, clusters):
    resilient = t_m[0]
    resilient.absorb_cluster(t_m[1])
    new_clusters = {v.id: v for k, v in clusters.items() if v.id != t_m[1].id}
    new_clusters[resilient.id] = resilient
    return new_clusters, t_m[1].id


def split_cluster(s_cluster, clusters, new_id=None):
    new_fkey = s_cluster.get_worst_seed()
    s_cluster.remove_face(new_fkey)

    if new_id is None:
        new_id = max(clusters.keys()) + 1

    clusters[new_id] = Cluster(new_id, new_fkey)
    clusters[s_cluster.id] = s_cluster
    return clusters


def merge_split(clusters):
    adj_cl = cluster_adjacency(clusters)  # returns objects list
    to_merge = get_clusters_to_merge(adj_cl)  # returns objects tuple
    to_split = get_cluster_to_split(clusters)  # returns object single

    if execute_merge_split(to_merge, to_split) is True:
        clusters, new_id = merge_clusters(to_merge, clusters)
        clusters = split_cluster(to_split, clusters, new_id)
    return clusters


def execute_merge_split(t_m, t_s):
    to_merge_err = reduce(lambda x, y: x+y, [x.get_distortion() for x in t_m])
    merged_err = simulate_merge(t_m[0], t_m[1])
    dif = merged_err - to_merge_err
    worst_err = t_s.get_distortion()

    if math.fabs(dif) < 0.5 * worst_err:  # 0.5
        print('merge-split is True')
        return True

    else:
        print('merge-split is False')
        return False


def get_proxy(faces):
    w_ve = [face.w_vector for face in faces]
    w_ve = list(map(lambda x: ut.align_vector(x, w_ve[0]), w_ve))
    r_ve = reduce(lambda x, y: cg.add_vectors(x, y), w_ve)
    return cg.normalize_vector(r_ve)


def get_errors(faces, proxy):
    return [face.get_error(proxy) for face in faces]


def get_distortion(errors):
    return reduce(lambda x, y: x+y, errors)


class Queue():
    def __init__(self, clusters, faces):
        self.clusters = clusters
        self.faces = faces
        self.queue = []
        heapq.heapify(self.queue)

    def init_queue(self):
        for c_key, cluster in self.clusters.items():
            cluster.add_seed(self.faces)
            cluster.set_vector_proxy()
            n_faces = self.get_neighbour_faces(cluster.seed)
            self.update_queue(n_faces, c_key)

    def update_queue(self, n_faces, c_key):
        for f in n_faces:
            if f.cluster is None:
                error = f.get_error(self.clusters.get(c_key).proxy)
                entry = {'fkey': f.fkey, 'cluster': c_key}
                heapq.heappush(self.queue, KeyDict(error, entry))

    def assign(self):
        while len(self.queue) > 0:
            entry = heapq.heappop(self.queue)
            cu_f = entry.dct.get('fkey')
            face = self.faces.get(cu_f)

            if face.cluster is None:
                c_key = entry.dct.get('cluster')
                cluster = self.clusters.get(c_key)
                cluster.add_face(face)
                self.update_queue(self.get_neighbour_faces(cu_f), c_key)

    def get_neighbour_faces(self, f_key):
        for nf in self.faces.get(f_key).neighbours:
            yield self.faces.get(nf)

    def get_clusters(self):
        for ckey, cluster in self.clusters.items():
            cluster.set_vector_proxy()
            cluster.set_distortion()
        return self.clusters


class Cluster():
    def __init__(self, id, f):
        self.id = id
        self.seed = f

        self.faces = []
        self.faces_keys = []

        self.proxy = None
        self.distortion = None
        self.add_face_key(self.seed)

    def remove_face(self, fkey):
        self.faces_keys = [k for k in self.faces_keys if k != int(fkey)]
        self.faces = [f for f in self.faces if f.fkey != fkey]
        self.set_vector_proxy()

    def absorb_cluster(self, other_cluster):
        for o_face in other_cluster.faces:
            self.add_face(o_face)
        self.set_vector_proxy()
        self.set_faces_in_cluster()

    def copy_cluster(self, other_cluster):
        self.faces_keys = list(set(other_cluster.faces_keys))
        self.proxy = other_cluster.proxy
        self.distortion = other_cluster.distortion

    def add_face_key(self, f_key):
        if f_key not in self.faces_keys:
            self.faces_keys.append(f_key)

    def add_seed(self, faces):
        seed_face = faces.get(self.seed)
        seed_face.set_cluster(self.id)
        self.faces.append(seed_face)

    def add_face(self, face):
        if face.fkey not in self.faces_keys:
            face.set_cluster(self.id)
            self.faces_keys.append(face.fkey)
            self.faces.append(face)

    def harvest_faces(self, faces):
        for key in self.faces_keys:
            self.faces.append(faces.get(key))

    def get_weighed_vectors(self):
        return [face.w_vector for face in self.faces]

    def get_vectors(self):
        return [face.vector for face in self.faces]

    def set_faces_in_cluster(self):
        for face in self.faces:
            face.cluster = self.id

    def set_vector_proxy_normal(self):  # NORMAL
        w_ve = self.get_weighed_vectors()
        r_ve = reduce(lambda x, y: cg.add_vectors(x, y), w_ve)
        self.proxy = cg.normalize_vector(r_ve)

    def set_vector_proxy(self):  # NEW
        w_ve = self.get_weighed_vectors()
        w_ve = list(map(lambda x: ut.align_vector(x, w_ve[0]), w_ve))
        r_ve = reduce(lambda x, y: cg.add_vectors(x, y), w_ve)
        self.proxy = cg.normalize_vector(r_ve)

    def get_errors(self):
        return [face.get_error(self.proxy) for face in self.faces]

    def get_new_seed(self):  # to improve timing
        return min(self.faces, key=lambda x: x.get_error(self.proxy)).fkey

    def get_worst_seed(self):
        return max(self.faces, key=lambda x: x.get_error(self.proxy)).fkey

    def get_face_keys(self):
        return [f.fkey for f in self.faces]

    def get_faces_halfedges(self):
        face_halfedges = set()
        for f in self.faces:
            face_halfedges.update(f.halfedges)
        return face_halfedges

    def get_faces_vertices(self):
        face_vertices = set()
        for f in self.faces:
            face_vertices.update(f.vertices)
        return face_vertices

    def clear_faces(self):
        for face in self.faces:
            face.clear_cluster()
        self.faces[:] = []

    def get_distortion(self):
        return reduce(lambda x, y: x+y, self.get_errors())

    def set_distortion(self):
        self.distortion = self.get_distortion()

    def __repr__(self):
        f = len(self.faces)
        fk = len(self.faces_keys)
        s = self.seed
        return 'id:{0} seed:{1} faces:{2}, keys:{3}'.format(self.id, s, f, fk)


class Face():
    def __init__(self, fkey):
        self.fkey = fkey
        self.vertices = None
        self.halfedges = None

        self.vector = None
        self.vector_length = None
        self.w_vector = None
        self.area = None
        self.neighbours = None
        self.error = None

        self.cluster = None

    def set_vertices(self, vertices):
        self.vertices = vertices

    def set_halfedges(self, halfedges):
        self.halfedges = halfedges

    def set_cluster(self, cluster):
        self.cluster = cluster

    def set_vector(self, vector):
        self.vector = cg.normalize_vector(vector)
        self.vector_length = cg.length_vector(vector)

    def set_weighed_vector(self, area_weight=False):
        if area_weight is True:
            self.w_vector = cg.scale_vector(self.vector, self.area)
        else:
            self.w_vector = self.vector

    def set_area(self, area):
        self.area = area

    def set_neighbours(self, neighbours):
        self.neighbours = [n for n in neighbours if n is not None]

    def get_error(self, proxy, area_weight=False):
        ali_vec = ut.align_vector(self.vector, proxy)
        difference = cg.subtract_vectors(ali_vec, proxy)

        # error = cg.length_vector_sqrd(difference)  # original
        error = cg.length_vector(difference)

        # error = cg.angle_vectors(ali_vec, proxy)
        # error = error ** 2
        
        # if area_weight is True:
        #     error = self.area * cg.length_vector_sqrd(difference)

        # do something about weights
        # w_1 = 0.3
        # w_2 = 0.7
        # w_error = w_1 * error + w_2 * self.vector_length

        return error

    def set_error(self, proxy):  # NEW
        self.error = self.get_error(proxy)

    def clear_cluster(self):
        self.cluster = None


@total_ordering
class KeyDict(object):
    def __init__(self, key, dct):
        self.key = key
        self.dct = dct

    def __lt__(self, other):
        return self.key < other.key

    def __eq__(self, other):
        return self.key == other.key

    def __repr__(self):
        return '{0.__class__.__name__}(key={0.key}, dct={0.dct})'.format(self)
