#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
A wrapper for a bunch of nodes with boid temperament.
'''


__name__ = "Streamline"
__author__ = "Rafael Pastrana"
__version__ = "0.0.4"
__creation__ = "2018.11.12"
__date__ = "2018.11.12"


import compas
import compas_rhino
import math

from compas.datastructures import Mesh
from compas.geometry import Line
from compas.geometry import Polyline
from compas.geometry import KDTree

from streamlines.node import Node

try:
    import rhinoscriptsyntax as rs
except:
    if compas.IPY:
        raise

class Streamline():

    def __init__(self, start):
        self.start = start
        self.grow = True
        self.id = None
        self.parent = None
        self.child = None

        self.tree = None
        self.polyline = None

        self.nodes = []
        self.nodes_forth = []
        self.nodes_back = []
        self.nd_search = []
        self.face_trace = []

    def __repr__(self):
        s = self.start.f_id
        num = len(self.nodes)
        try:
            le = self.polyline.length
        except Exception:
            le = 0
        return 'f_id={0}, nodes={1}, length={2})'.format(s, num, le)

    def create_polyline(self):
        if len(self.nodes) > 2:
            points = [node.pos for node in self.nodes]
            self.polyline = Polyline(points)

    def add_node(self, node):
        self.nodes.append(node)

    def arrange_nodes(self):
        nodes_back = self.nodes_back
        self.nodes[:] = []
        nodes_back.reverse()
        self.nodes.extend(nodes_back)
        self.nodes.append(self.start)
        self.nodes.extend(self.nodes_forth)

    def activate_growth(self):
        self.grow = True

    def update_search_tree(self):
        # self.tree = KDTree([nd.pos for nd in self.nd_search if nd is not None])
        pass

    def is_point_close(self, point, distance, n):
        # exclude = set()
        # for i in range(n):
        #     nnbr = self.tree.nearest_neighbor(point, exclude)
        #     exclude.add(nnbr[1])

        # nnbr, label, dst = self.tree.nearest_neighbor(point, exclude)
        # if dst is not None:
        #     if dst < distance:
        #         return True
        # return False
        pass

    def resample_polyline(self, length):
        try:
            rs_poly = rs.AddPolyline(self.polyline.points)
            if length <= self.polyline.length:
                a = rs.DivideCurveLength(rs_poly, length, False)
                num_div = len(rs.DivideCurveLength(rs_poly, length, False, False))
                new_pts = rs.DivideCurve(rs_poly, num_div, False, True)

                new_rs_poly = rs.AddPolyline(new_pts)
                segs = rs.ExplodeCurves(new_rs_poly)

            else:
                print('it is a line!')
                segs = [rs.AddLine(rs.CurveStartPoint(rs_poly), rs.CurveEndPoint(rs_poly))]
                print(segs)

            out_pts = []
            out_vec = []

            for seg in segs:
                new_pt = rs.CurveMidPoint(seg)
                v = rs.VectorCreate(rs.CurveEndPoint(seg), rs.CurveStartPoint(seg))
                new_pt = [new_pt.X, new_pt.Y, new_pt.Z]
                new_vec = [v.X, v.Y, v.Z]
                new_vec = compas.geometry.normalize_vector(new_vec)
                out_pts.append(new_pt)
                out_vec.append(new_vec)
            # print('succesfully resampled')
            return out_pts, out_vec

        except Exception:
            print('Polyline could not be resampled.')
            return None, None

    def resample_curve(self, le):
        try:
            # rs_poly = rs.AddPolyline(self.polyline.points)
            crv = rs.AddInterpCurve(self.polyline.points)

            if le <= self.polyline.length:
                a = rs.DivideCurveLength(crv, le, False)
                div = rs.DivideCurveLength(crv, le, False, False)
                new_pts = rs.DivideCurve(crv, len(div), False, True)
                new_par = rs.DivideCurve(crv, len(div), False, False)

            else:
                print('it is a line!')
                new_pts = [rs.CurveStartPoint(crv), rs.CurveEndPoint(crv)]
                new_par = [0.0, rs.CurveLength(crv)]
                print(segs)

            out_pts = []
            out_vec = []

            for idx, pt in enumerate(new_pts):

                pt = [pt.X, pt.Y, pt.Z]

                v = rs.CurveTangent(crv, new_par[idx])
                vec = [v.X, v.Y, v.Z]
                # vec = compas.geometry.normalize_vector(vec)

                out_pts.append(pt)
                out_vec.append(vec)
            # print('succesfully resampled')

            return out_pts, out_vec

        except Exception:
            print('Interpolated Curve could not be resampled.')
            return None, None
