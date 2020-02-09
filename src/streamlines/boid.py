#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
A boid object that seeks, aligns, coheres, and separates.
'''


__name__ = "Boid"
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

from compas.datastructures import Mesh
from compas.geometry import Point
from compas.geometry import Line
from compas.geometry import subtract_vectors

from streamlines.utilities import Utilities

try:
    import rhinoscriptsyntax as rs
except:
    if compas.IPY:
        raise

TOL = 1e-6
ut = Utilities()


class Boid():

    def __init__(self):
        self.x = 0
        self.y = 0
        self.z = 0

        self.pos = [self.x, self.y, self.z]
        self.vel = [0, 0, 0]
        self.acc = [0, 0, 0]

        self.max_speed = 1.0
        self.max_force = 1.0

    def seek(self, target):
        desired = cg.subtract_vectors(target, self.pos)
        desired = cg.normalize_vector(desired)
        desired = cg.scale_vector(desired, self.max_speed)

        # Reynold's steering formula
        steer = cg.substract_vectors(desired, self.vel)

        ut.limit_vector(steer, self.max_force)
        return steer

    def separate(self, dst, boids, mxf=0.02, mxs=0.02):
        count = 0
        sum_vec = [0, 0, 0]
        tol = 1e-6

        for boid in boids:
            d = cg.distance_point_point(self.pos, boid.pos)
            print('distance is {}'.format(d))

            if d > tol and d < dst:
                dif = subtract_vectors(self.pos, boid.pos)
                dif = cg.normalize_vector(cg.scale_vector(dif, 1/d))
                sum_vec = cg.add_vectors(dif, sum_vec)
                count += 1

            if count > 0:
                sum_vec = cg.normalize_vector(cg.scale_vector(sum_vec,
                                                              1/count)
                                              )
                sum_vec = cg.scale_vector(sum_vec, mxs)

                # Reynold's steering formula
                steer = subtract_vectors(sum_vec, self.vel)
                steer = ut.limit_vector(steer, mxf)
                return steer
        return sum_vec

    def align(self, dst, boids, mxf=0.02, mxs=0.02):
        count = 0
        sum_vec = [0, 0, 0]

        for boid in boids:
            d = cg.distance_point_point(self.pos, boid.pos)
            if d > 0 and d < dst:
                sum_vec = cg.add_vectors(boid.vel, sum_vec)
                count += 1

        if count > 0:
            sum_vec = cg.normalize_vector(cg.scale_vector(sum_vec, 1/count))
            sum_vec = cg.scale_vector(sum_vec, mxs)
            # Reynold's steering formula
            steer = cg.subtract_vectors(sum_vec, self.vel)
            steer = ut.limit_vector(steer, mxf)
            return steer
        return sum_vec

    def flock(self, boids, dst, f_s, f_al):
        sep = cg.scale_vector(self.separate(dst, boids), f_s)
        ali = cg.scale_vector(self.align(dst, boids), f_al)

        self.apply_force(sep)
        self.apply_force(ali)

    def apply_force(self, vector):
        self.acc = cg.add_vectors(self.acc, vector)

    def update(self):
        self.vel = cg.add_vectors(self.vel, self.acc)
        self.vel = ut.limit_vector(self.vel, self.max_speed)
        self.pos = cg.add_vectors(self.pos, self.vel)
        self.acc = cg.scale_vector(self.acc, 0.)

    def run(self, boids, dst, f_s, f_ali):
        self.flock(boids, dst, f_s, f_ali)
        self.update()
