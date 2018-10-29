'''
A compas implementation of a streamline system.
'''


__name__ = "Streamsystem Compas"
__author__ = "Rafael Pastrana"
__version__ = "0.0.6"
__creation__ = "2018.08.19"
__date__ = "2018.08.19"

import sys
import imp
import compas
import compas_rhino
import heapq
import math
import compas.geometry as cg
import rhinoscriptsyntax as rs
import Rhino.Geometry as rg
import time

import Node
import Streamline
import Utilities

imp.reload(Node)
imp.reload(Streamline)
imp.reload(Utilities)

from Node import Node
from Streamline import Streamline
from compas.datastructures import Mesh
from System import DateTime
from Utilities import Utilities

ut = Utilities()


class StreamSystem():

    def __init__(self,
                 sMesh,
                 as_tag,
                 dL=0.02,
                 min_sp=0.05,
                 uni_sp=False
                 ):

        self.sMesh = sMesh
        self.dL = dL
        self.min_sp = min_sp / 2
        self.uni_sp = uni_sp
        self.s_tag = as_tag

        self.v_tag = None
        self.w = [0, 0, 0]
        self.min_length = 1e-6

        self.streamlines = []
        self.lines = []
        self.int_points = []
        self.objs = []

        self.seed_pts = []
        self.offsets = []

    def makeStreamlines3D(self, seeds):
        for rh_seed in seeds:
            seed = Node()
            seed.update_pos_rh(rh_seed)
            seed.pull_to_mesh(self.sMesh, self.v_tag)

            poly_pts = self.make_streamline(seed)

    def set_tracing_data(self, vector_tag, w, min_length):
        self.v_tag = vector_tag
        self.w = w
        self.min_length = min_length

    def make_streamlines_mebarki(self, seeds):
        seeds = seeds
        while len(seeds) > 0:
            pt_seed = seeds.pop(0)
            seed = Node()
            valid = True

            try:
                seed.update_pos_rh(pt_seed)
            except Exception:
                seed.update_pos(pt_seed)

            seed.pull_to_mesh(self.sMesh, self.v_tag)

            if self.uni_sp is False:
                sp = self.sMesh.cMesh.get_face_attribute(seed.f_id, self.s_tag)
            else:
                sp = self.min_sp

            rings = math.ceil(self.sMesh.avg_edge / sp)
            nbrs = seed.find_neighbours(self.sMesh.cMesh, self.objs, rings)
            if nbrs:
                if ut.is_point_close(seed.pos, nbrs, sp) is True:
                    valid = False

            if valid is True:
                poly_pts = self.make_streamline(seed)

    def make_streamlines_jobard(self):
        '''
        """"""" PART A """"""
        1. create queue list with streamlines. store them as they appear.
        2. find the most stressed point on the design space
        3. make streamline
        4. append to the queue

        """"""" PART B """"""
        5. pop streamline-i from queue
        6. sample points on the streamline-i at min spacing
        7. along streamline-i, find point at most stressed face
        8. select point and offset to the right and left of streamline-i
        9. make point on the right. make point on the left
        10. take point on the right. make streamline-j. append to queue
        11. take point on the left. make streamline-k. append to queue
        12. select next sampling point on (second most stressed point)

        """"""" PART C """"""
        13. repeat points 8 until 11 until no more samplint points exist

        """"""" PART D """"""
        14. repeat parts B and C until queue is empty

        15. finish
        '''

        # 1. create queue list with streamlines. store them as they appear
        queue = []
        heapq.heapify(queue)

        # 2. find the most stressed point on the design space faces
        f_keys = [f for f in self.sMesh.cMesh.faces()]
        as_ = self.sMesh.cMesh.get_faces_attribute(self.s_tag, fkeys=f_keys)
        max_ = min(zip(as_, f_keys), key=lambda x: x[0])[1]
        start_pt = self.sMesh.cMesh.face_centroid(max_)

        start_nd = Node()
        start_nd.update_pos(start_pt)
        start_nd.pull_to_mesh(self.sMesh, self.v_tag)

        # 3. make streamline
        strm = self.make_streamline(start_nd, rs_pt=False)
        self.objs[:] = []  # clear objs search
        self.streamlines[:] = []  # clear first streamline

        # 4. append to the queue
        queue.append(strm)
        # heapq.heappush(queue, KeyDict(error, entry))
        count = 0

        # 5. pop streamline-i from queue
        #  while len(queue) > 0:
        for i in range(40):
            if len(queue) > 0:
                cu_strm = queue.pop(0)
                count += 1

            # 6. sample points on the streamline-i at min spacing
                if cu_strm is not None:
                    s_pts, s_vels = cu_strm.resample_polyline(0.2)
                    if s_pts is not None:

                        s_ = []
                        heapq.heapify(s_)

                        for idx, pt in enumerate(s_pts):
                            nd = Node()
                            nd.update_pos(s_pts[idx])
                            nd.pull_to_mesh(self.sMesh, self.v_tag)
                            nd.update_vel(s_vels[idx])
                            sep_ = self.sMesh.cMesh.get_face_attribute(nd.f_id,
                                                                       self.s_tag
                                                                       )

                            heapq.heappush(s_, (sep_, nd))

                    # 8. select point and offset to the right and left of streamline-i
                        # while len(s_) > 0:
                        for i in range(4):
                            if len(s_) > 0:
                                cu_tup = heapq.heappop(s_)
                                cu_sep, cu_nd = cu_tup
                                cu_sep = self.min_sp

                                if count != 1:
                                    cu_sep = cu_sep * 2

                                self.seed_pts.append(rs.AddPoint(*cu_nd.pos))
                                nrml = self.sMesh.cMesh.face_normal(cu_nd.f_id)
                                of_vec = cg.cross_vectors(cu_nd.vel, nrml)
                                of_vec = cg.normalize_vector(of_vec)

                        # # # 9B. make point on the right with particle tracing
                                t_nd_r = Node()
                                t_nd_r.update_pos(cu_nd.pos)
                                t_nd_r.pull_to_mesh(self.sMesh, self.v_tag)
                                t_nd_r.update_vel(of_vec)
                                strm_right = Streamline(t_nd_r)
                                pt_right = self.offset(strm_right, cu_sep)

                        # # 9C. make point on the left with particle tracing
                                t_nd_l = Node()
                                t_nd_l.update_pos(cu_nd.pos)
                                t_nd_l.pull_to_mesh(self.sMesh, self.v_tag)
                                t_nd_l.update_vel(cg.scale_vector(of_vec, -1))
                                strm_left = Streamline(t_nd_l)
                                pt_left = self.offset(strm_left, cu_sep)

                        # 10. take point on the right. make streamline-j. append to queue


                                if pt_right is not None:
                                    nd_ri = Node()
                                    nd_ri.update_pos(pt_right)
                                    nd_ri.pull_to_mesh(self.sMesh, self.v_tag)
                                    valid = True

                                    if self.uni_sp is False:
                                        sp = self.sMesh.cMesh.get_face_attribute(nd_ri.f_id, self.s_tag)
                                    else:
                                        sp = self.min_sp

                                    rings = math.ceil(self.sMesh.avg_edge / sp)
                                    nbrs = nd_ri.find_neighbours(self.sMesh.cMesh, self.objs, rings)
                                    if nbrs:
                                        if ut.is_point_close(nd_ri.pos, nbrs, sp) is True:
                                            valid = False

                                    if valid is True:
                                        strm_right = self.make_streamline(nd_ri, rs_pt=False)
                                        if strm_right is not None:
                                            queue.append(strm_right)

                        # 11. take point on the left. make streamline-k. append to queue
                                if pt_left is not None:
                                    nd_left = Node()
                                    nd_left.update_pos(pt_left)
                                    nd_left.pull_to_mesh(self.sMesh, self.v_tag)
                                    valid = True

                                    if self.uni_sp is False:
                                        sp = self.sMesh.cMesh.get_face_attribute(nd_left.f_id, self.s_tag)
                                    else:
                                        sp = self.min_sp

                                    rings = math.ceil(self.sMesh.avg_edge / sp)
                                    nbrs = nd_left.find_neighbours(self.sMesh.cMesh, self.objs, rings)
                                    if nbrs:
                                        if ut.is_point_close(nd_left.pos, nbrs, sp) is True:
                                            valid = False

                                    if valid is True:
                                        strm_left = self.make_streamline(nd_left, rs_pt=False)
                                        if strm_left is not None:
                                            queue.append(strm_left)

    def make_streamline(self, seed, rs_pt=True):
        strm = Streamline(seed)
        strm.nd_search.append(seed)  # append to search

        self.grow(strm, self.v_tag, self.w, 'forth')
        strm.activate_growth()
        self.grow(strm, self.v_tag, self.w, 'back')

        strm.arrange_nodes()
        strm.create_polyline()

        if strm.polyline is not None:
            if strm.polyline.length > self.min_length:

                self.streamlines.append(strm.polyline)
                self.objs.extend(strm.nodes)

                if rs_pt is True:
                    pl = [rs.AddPoint(*pt) for pt in strm.polyline.points]
                else:
                    pl = strm
                return pl
        return None

    def grow(self,
             streamline,
             vtag,
             weights,
             direction='forth',
             exit=True,
             max_iter=1000,
             prox_factor=0.2,
             self_prox_factor=1.0,
             agent_prox_factor=2.0,
             other_prox_factor=1.0
             ):

        # 1. set major parameters
        proj_dist = self.dL  # scale for vector to be projected on plane
        prox = prox_factor * self.dL  # for checking if it's inside the gh mesh

        dL = self.dL  # step size from growth
        ds_own = self.dL * self_prox_factor  # distance to stop agains itself

        count = 0
        nodes = [streamline.start]

        # do forward integration
        while streamline.grow is True:

            # 0. get seed point to grow from and initial vector
            nd = nodes[-1]
            if count == 0:
                if direction == 'back':
                    nd.vel = cg.scale_vector(nd.vel, -1.)

            # 1. get threshold distances
            d = self.min_sp
            if self.uni_sp is False:
                d = self.sMesh.cMesh.get_face_attribute(nd.f_id, self.s_tag)
            r = math.ceil(d / self.sMesh.avg_edge)
            od = d * other_prox_factor

            # 2. check if too close to other curves before growing. if so, stop
            nbrs = nd.find_neighbours(self.sMesh.cMesh, self.objs, r, False)
            if nbrs:
                nbrs_pos = [nbr.pos for nbr in nbrs]
                fl, sep_d = ut.is_point_close(nd.pos, nbrs_pos, od, d_out=True)

                if fl is True:
                    continue
                    streamline.grow = False

                sep_vec = nd.separate(d * agent_prox_factor, nbrs)
                ali_vec = nd.align(d * agent_prox_factor, nbrs)

            # 3 check if too close to itself with delaunay triangulation
            o_nbrs = nd.find_neighbours(self.sMesh.cMesh, streamline.nd_search)
            if o_nbrs:
                flag = ut.delaunay_is_point_close(nd.pos, o_nbrs, ds_own)
                if flag is True:
                    print('delaunay ON')
                    streamline.grow = False

            # 3. find vector to follow
            c_vec = self.sMesh.get_vector_on_face(nd.pos, nd.f_id,
                                                  vtag,
                                                  nd.vel
                                                  )
            fol_vec = cg.scale_vector(cg.normalize_vector(c_vec), weights[0])

            # 4. get alignment and separation vectors
            try:
                sep_vec = cg.scale_vector(sep_vec, sep_d * weights[1] / d)
                ali_vec = cg.scale_vector(sep_vec, weights[2])
            except Exception:
                sep_vec = cg.scale_vector([0, 0, 0], weights[1])
                ali_vec = cg.scale_vector([0, 0, 0], weights[2])

            # 5. add up follow, separation and alignment vectors
            new_cvec = cg.add_vectors([0, 0, 0], fol_vec)
            new_cvec = cg.add_vectors(new_cvec, sep_vec)
            new_cvec = cg.add_vectors(new_cvec, ali_vec)
            new_cvec = cg.normalize_vector(new_cvec)
            new_cvec = ut.align_vector(new_cvec, nd.vel)
            c_vec = new_cvec  # applying the agent based vector

            # 6. project resultant vector to face plane
            plane = self.sMesh.cMesh.get_face_attribute(nd.f_id, 'plane')
            temp_pt = cg.translate_points([nd.pos],
                                          cg.scale_vector(c_vec,
                                                          proj_dist)
                                          )[0]
            proj_pt = cg.project_point_plane(temp_pt, plane)
            c_vec = cg.normalize_vector(cg.vector_from_points(nd.pos, proj_pt))
            c_vec = ut.align_vector(cg.scale_vector(c_vec, dL), nd.vel)

            # 7. create new point by displacing along projected vector
            new_pos = cg.translate_points([nd.pos], c_vec)[0]
            new_nd = Node()
            new_nd.clone_from(nd)
            new_nd.update_pos(new_pos)
            new_nd.update_vel(c_vec)

            # 8. check for edge intersections
            face_edges = self.sMesh.cMesh.face_halfedges(nd.f_id)
            line = [nd.pos, new_nd.pos]
            ints = self.find_edge_intersections(line, face_edges, new_nd)

            if ints is not None:
                int_pos, u, v, nextNd = ints
                new_nd.edge = (u, v)
                new_nd.is_intersection = True  # it is an intersection
                self.int_points.append(rs.AddPoint(*int_pos))

                # find new face index
                new_faces = list(filter(lambda x: x != nd.f_id,
                                        self.sMesh.cMesh.edge_faces(u, v)
                                        )
                                 )
                # test against found faces
                if new_faces[0] is None:
                    new_nd.update_pos(int_pos)
                    new_nd.f_id = new_faces[0]
                    streamline.nd_search.append(new_nd)  # append to search
                    nodes.append(new_nd)
                    print('It hit the boundary')
                    streamline.grow = False
                    break

                else:
                    new_nd.update_pos(int_pos)
                    new_nd.f_id = new_faces[0]

            else:
                tpt, tid = self.sMesh.closest_point(new_nd.pos)
                dst = cg.distance_point_point(tpt, new_nd.pos)

                if dst < prox:
                    new_nd = Node()
                    new_nd.clone_from(nd)
                    new_nd.update_pos(new_pos)
                    new_nd.f_id = tid
                    new_nd.update_vel(c_vec)
                else:
                    print('dst is: {}'.format(dst))
                    print('OUT. distance is larger than tol')
                    streamline.grow = False
                    break

            # 9. store new node
            streamline.nd_search.append(new_nd)  # append to search
            nodes.append(new_nd)

            # 10. increase the count
            count += 1

            # 11. additional break flag for security
            if exit is True:
                if count >= max_iter:
                    print('count was exceeded')
                    streamline.grow = False

        # 12. the while loop is done, extend correspoinding nodes list
        nodes.pop(0)
        if direction == 'forth':
            streamline.nodes_forth.extend(nodes)
        elif direction == 'back':
            streamline.nodes_back.extend(nodes)

    def find_edge_intersections(self, line, edges, node):
        tol = 1e-6
        start = line[0]
        intersections = []

        for u, v in edges:
            plane = self.sMesh.cMesh.get_edge_attribute((u, v), 'plane')
            compas_int = cg.intersection_segment_plane(line, plane)

            if compas_int is not None:
                cpt = cg.closest_point_on_segment(compas_int, line)
                dst = cg.distance_point_point(compas_int, cpt)

                if dst < tol:
                    if node.edge is None:
                        return compas_int, u, v, node

                    elif node.edge[0] in (u, v) or node.edge[1] in (u, v):
                        if cg.distance_point_point(start, compas_int) > tol:
                            return compas_int, u, v, node

                    else:
                        intersections.append((compas_int, u, v, node))
        if len(intersections) > 0:
            return intersections[0]  # take the first entry only

        return None

    def offset(self,
               streamline,
               target_length,
               exit=True,
               growth_rate=0.01,
               max_iter=2000,
               prox_factor=0.5,
               ):

        # 1. set major parameters
        proj_dist = growth_rate  # scale for vector to be projected on plane
        prox = prox_factor * growth_rate  # for checking if it's inside
        dL = growth_rate  # step size from growth

        count = 0
        nodes = [streamline.start]
        out = None

        while streamline.grow is True:

            # 1. create polyline to check length
            streamline.create_polyline()
            if streamline.polyline is not None:
                if streamline.polyline.length > target_length:
                    streamline.grow = False

            # 2. get seed point to grow from and initial vector
            nd = nodes[-1]

            # 3. find vector to follow
            c_vec = nd.vel

            # 4. project resultant vector to face plane
            plane = self.sMesh.cMesh.get_face_attribute(nd.f_id, 'plane')
            temp_pt = cg.translate_points([nd.pos],
                                          cg.scale_vector(c_vec,
                                                          proj_dist)
                                          )[0]
            proj_pt = cg.project_point_plane(temp_pt, plane)
            c_vec = cg.normalize_vector(cg.vector_from_points(nd.pos, proj_pt))
            c_vec = ut.align_vector(cg.scale_vector(c_vec, dL), nd.vel)

            # 5. create new point by displacing along projected vector
            new_pos = cg.translate_points([nd.pos], c_vec)[0]
            new_nd = Node()
            new_nd.clone_from(nd)
            new_nd.update_pos(new_pos)
            new_nd.update_vel(c_vec)

            # 8. check for edge intersections
            face_edges = self.sMesh.cMesh.face_halfedges(nd.f_id)
            line = [nd.pos, new_nd.pos]
            ints = self.find_edge_intersections(line, face_edges, new_nd)

            if ints is not None:
                int_pos, u, v, nextNd = ints
                new_nd.edge = (u, v)
                new_nd.is_intersection = True  # it is an intersection
                self.int_points.append(rs.AddPoint(*int_pos))

                # find new face index
                new_faces = list(filter(lambda x: x != nd.f_id,
                                        self.sMesh.cMesh.edge_faces(u, v)
                                        )
                                 )
                # test against found faces
                if new_faces[0] is None:
                    new_nd.update_pos(int_pos)
                    new_nd.f_id = new_faces[0]
                    streamline.nd_search.append(new_nd)  # append to search
                    nodes.append(new_nd)
                    print('It hit the boundary')
                    streamline.grow = False
                    break

                else:
                    new_nd.update_pos(int_pos)
                    new_nd.f_id = new_faces[0]

            else:
                tpt, tid = self.sMesh.closest_point(new_nd.pos)
                dst = cg.distance_point_point(tpt, new_nd.pos)

                if dst < prox:
                    new_nd = Node()
                    new_nd.clone_from(nd)
                    new_nd.update_pos(new_pos)
                    new_nd.f_id = tid
                    new_nd.update_vel(c_vec)
                else:
                    print('dst is: {}'.format(dst))
                    print('OUT. distance is larger than tol')
                    streamline.grow = False
                    break

            # 9. store new node
            nodes.append(new_nd)
            streamline.nodes.append(new_nd)

            # 10. increase the count
            count += 1

            # 11. additional break flag for security
            if exit is True:
                if count >= max_iter:
                    print('count was exceeded')
                    streamline.grow = False

            try:
                out = streamline.nodes[-1].pos
            except Exception:
                continue

        return out
