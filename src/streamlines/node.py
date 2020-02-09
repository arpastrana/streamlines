#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
A node class with boid temperament.
'''


__name__ = "Node"
__author__ = "Rafael Pastrana"
__version__ = "0.0.4"
__creation__ = "2018.11.12"
__date__ = "2018.11.12"

import sys
import imp

import compas
import compas.geometry as cg
import compas_rhino
import math

from streamlines.boid import Boid
from compas.datastructures import Mesh
from compas.geometry import Point
from compas.geometry import Line
from compas.topology import dijkstra_distances


TOL = 1e-6


class Node(Boid):

    def __init__(self):
        self.x = 0
        self.y = 0
        self.z = 0

        self.pos = [self.x, self.y, self.z]
        self.vel = [0, 0, 0]

        self.f_id = None
        self.edge = None

        self.is_intersection = False
        self.is_start = False

        self.point_right = None
        self.point_left = None

    def pull_to_mesh(self, s_mesh, name):
        self.pos, self.f_id = s_mesh.closest_point(self.pos)
        self.update_pos(self.pos)
        edge = self.pull_to_edge(s_mesh, self.f_id)

        if edge is not None:
            self.edge = edge
            self.is_intersection = True

        self.vel = s_mesh.get_vector_on_face(self.pos, self.f_id, name)

    def update_pos(self, point):
        try:
            self.x = point[0]
            self.y = point[1]
            self.z = point[2]
            self.pos = point
        except Exception:
            rh_pos = [point.X, point.Y, point.Z]
            self.pos = rh_pos
            self.update_pos(rh_pos)

    def update_pos_rh(self, rh_point):
        rh_pos = [rh_point.X, rh_point.Y, rh_point.Z]
        self.update_pos(rh_pos)

    def update_vel(self, vector):
        self.vel = vector

    def clone_from(self, node):
        self.x = node.x
        self.y = node.y
        self.z = node.z

        self.pos = node.pos
        self.vel = node.vel

        self.f_id = node.f_id
        self.edge = node.edge

        self.is_intersection = node.is_intersection

    def pull_to_edge(self, s_mesh, f_id):
        for u, v in s_mesh.c_mesh.face_halfedges(f_id):
            line = s_mesh.c_mesh.edge_coordinates(u, v)

            if cg.distance_point_point(line[0], line[1]) > TOL:
                if cg.distance_point_line(self.pos, line) < TOL:
                    return (u, v)
        return None

    def find_neighbors(self, c_mesh, objs, rings=1, pos=True):
        count = 0
        f_keys = set()
        f_keys.add(self.f_id)

        while count < rings:
            temp_f_keys = []
            for f_key in f_keys:
                temp_f_keys.extend(c_mesh.face_neighbors(f_key))

            f_keys.update(temp_f_keys)
            count += 1

        n_objs = [obj.pos for obj in objs if obj.f_id in f_keys]
        if len(n_objs) > 0:
            if pos is True:
                return n_objs
            else:
                return [obj for obj in objs if obj.f_id in f_keys]
        else:
            return None

    def find_dijkstra_neighbors(self, s_mesh, objs, length, pos=True):
        f_keys = set()
        v_keys = set()
        vertices = s_mesh.c_mesh.face_vertices(self.f_id)

        for v in s_mesh.c_mesh.face_vertices(self.f_id):
            dst = s_mesh.c_mesh.get_vertex_attribute(v, 'dijkstra')
            v_keys.update([v_key for v_key, d in dst.items() if d <= length])

        for v_key in v_keys:
            f_keys.update(s_mesh.c_mesh.vertex_faces(v_key))

        n_objs = [obj.pos for obj in objs if obj.f_id in f_keys]

        if len(n_objs) > 0:
            if pos is True:
                return n_objs
            else:
                return [obj for obj in objs if obj.f_id in f_keys]
        else:
            return None

    def set_links(self, pt_right, pt_left):
        self.pt_right = pt_right
        self.pt_left = pt.left
