'''
A compas implementation of a streamline system.
'''


__name__ = "Streamsystem"
__author__ = "Rafael Pastrana"
__version__ = "0.0.5"
__creation__ = "2018.11.12"
__date__ = "2019.03.12"


import imp
import heapq
import compas.geometry as cg
import rhinoscriptsyntax as rs

import node
import streamline
import utilities

from compas.geometry import KDTree
from compas.geometry import Polyline

imp.reload(node)
imp.reload(streamline)
imp.reload(utilities)

from utilities import Utilities
from node import Node
from streamline import Streamline
from timeit import default_timer as timer

ut = Utilities()


class Streamsystem:

    def __init__(self,
                 s_mesh,
                 as_tag,
                 dL=0.02,
                 min_sp=0.05,
                 uni_sp=False
                 ):

        self.s_mesh = s_mesh
        self.dL = dL
        self.min_sp = min_sp
        self.uni_sp = uni_sp
        self.s_tag = as_tag

        self.v_tag = None
        self.w = [0, 0, 0]
        self.min_length = 1e-6

        self.streamlines = []
        self.short_streamlines = []
        self.lines = []
        self.int_points = []
        self.objs = []
        self.tree = KDTree()

        self.seed_pts = []
        self.offsets = []
        self.first_pt = None
        self.extra_streamlines = []
        self.links_right = []
        self.other_points = []
        self.offset_vectors = []
        self.intersection_points = []

        self.grow_1_timing = 0.0
        self.grow_2_timing = 0.0
        self.grow_3_timing = 0.0
        self.grow_4_timing = 0.0
        self.grow_5_timing = 0.0

        self.timing_6 = 0.0

    def clone_streamsystem(self, streamsystem):
        self.streamlines = streamsystem.streamlines
        self.objs = streamsystem.objs
        self.tree = streamsystem.tree

    def set_tracing_data(self, vector_tag, w, min_length):
        self.v_tag = vector_tag
        self.w = w
        self.min_length = min_length

    # @ut.timing
    def get_new_streamline(self, point, o_prox=1.0, st_o_prox=0.8, s_prox=4.0):
        seed_node = self.make_new_node(point)
        sep = self.get_threshold_distance(seed_node) * st_o_prox
        start = timer()
        invalid = self.get_neighbor_proximity(seed_node.pos, sep)
        end = timer()
        thing = end - start
        self.timing_6 += thing

        if invalid is False:

            new_streamline = self.make_streamline(seed_node, o_prox, s_prox)
            return new_streamline

        # print('seed was at start too close to others')
        return None

    # @ut.timing
    def make_streamline(self, seed, o_prox=1.0, s_prox=1.0):
    	start = timer()
    	strm = Streamline(seed)
    	strm.nd_search.append(seed)  # append to search
    	end = timer()
    	thing = end - start
    	self.grow_3_timing += thing
    	
    	start = timer()
    	self.grow(strm, vtag=self.v_tag, direction='forth', dL=self.dL, o_prox=o_prox, s_prox=s_prox)
    	end = timer()
    	thing = end - start
    	self.grow_1_timing += thing
    	strm.activate_growth()
    	
    	start = timer()
    	self.grow(strm, vtag=self.v_tag, direction='back', dL=self.dL, o_prox=o_prox, s_prox=s_prox)
    	end = timer()
    	thing = end - start
    	self.grow_2_timing += thing
    	strm.activate_growth()
    	
    	start = timer()
    	strm.arrange_nodes()
    	strm.create_polyline()
    	end = timer()
    	thing = end - start
    	self.grow_4_timing += thing

    	start = timer()
    	if strm.polyline is not None:
    		poly_length = strm.polyline.length
    		if poly_length > self.min_length:
    			self.streamlines.append(strm)
    			self.objs.extend(strm.nodes)
    		else:
    			self.short_streamlines.append(strm)
    			print('length is: {} but wont be added'.format(poly_length))
    	else:
    		print('creation of polyline failed')
    	end = timer()
    	thing = end - start
    	self.grow_5_timing += thing
    	return strm

    def make_streamlines_mebarki(self, seeds, o_prox, st_o_prox):
        while len(seeds) > 0:
            print('******* new streamline started *********')
            self.update_search_tree()
            pt_seed = seeds.pop(0)
            streamline = self.get_new_streamline(pt_seed, o_prox=o_prox, st_o_prox=st_o_prox)
            print('streamline was {}'.format(streamline))

    @ut.timing
    def make_streamlines_jobard(self, strat_f, o_prox, st_o_prox, s_prox, num_samples, ite):
        # 1. create queue list with streamlines. store them as they appear
        queue = []
        heapq.heapify(queue)

        # 2. find the most stressed point on the design space faces
        min_sp = self.get_max_face_spacing() * strat_f
        start_nd = self.make_new_node(self.get_max_face_centroid())
        self.first_pt = rs.AddPoint(*start_nd.pos)

        # 3. make first streamline
        strm = self.make_streamline(start_nd)
        self.objs[:] = []  # clear objs search
        self.streamlines[:] = []  # clear first streamline

        # 4. append first streamline to the queue
        # print('******* initial streamline created *********')
        heapq.heappush(queue, (0, strm))

        count = 0
        node_count = 0
        timing_1 = 0.
        timing_2 = 0.
        timing_3 = 0.
        # 5. pop streamline-i from queue
        for i in range(ite):
            if len(queue) > 0:
                # print('******* new streamline popped *********')
                as_, cu_strm = heapq.heappop(queue)
                count += 1

            # 6. sample points on the streamline-i at min spacing OK
                s_ = self.make_sampling_nodes(cu_strm, min_sp)

                # 8. select node and offset to the right and left
                # while len(s_) > 0:
                
                if s_ is not None:
                    # while len(s_) > 0:
                    for i in range(num_samples):
                        if len(s_) > 0:
                            cu_sep, cu_nd = heapq.heappop(s_)
                            node_count += 1

                            if self.uni_sp is True:
                                cu_sep = self.min_sp

                            if count == 1:
                                cu_sep = cu_sep / 2

                        # 9. make points with particle tracing
                            # pt_right = self.get_offset_point(cu_nd, cu_sep)
                            # pt_left = self.get_offset_point(cu_nd, cu_sep, True)
                            _start = timer()
                            pt_right, pt_left = self.get_side_pts(cu_nd, cu_sep)
                            _end = timer()
                            _thing = _end - _start
                            timing_1 += _thing 

                        # 10. make streamline-j. append to queue
                            for pt_ in [pt_right, pt_left]:
                                new_streamline = None
                                if pt_ is not None:
                                	_start = timer()
                                	new_streamline = self.get_new_streamline(pt_, o_prox, st_o_prox, s_prox)
                                	_end = timer()
                                	_thing = _end - _start
                                	timing_2 += _thing

                                if new_streamline is not None:
                                    # heapq.heappush(queue, (count, new_streamline))
                                    heapq.heappush(queue, (cu_sep, new_streamline))

        # print('Number of Processed Nodes was: {}'.format(node_count))
        print('get side pts took: {} ms'.format(timing_1 * 1000))
        print('make streamlines took: {} ms'.format(timing_2 * 1000))
        print('grow 1 took: {} ms'.format(self.grow_1_timing * 1000))
        print('grow 2 took: {} ms'.format(self.grow_2_timing * 1000))
        print('instantiating streamlines took: {} ms'.format(self.grow_3_timing * 1000))
        print('polyline making took: {} ms'.format(self.grow_4_timing * 1000))
        print('polyline misc took: {} ms'.format(self.grow_5_timing * 1000))
        print('neighbor start took: {} ms'.format(self.timing_6 * 1000))

    # @ut.timing
    def get_side_pts(self, my_node, my_sep):
    	pt_r = self.get_offset_point(my_node, my_sep)
    	pt_l = self.get_offset_point(my_node, my_sep, True)
    	return pt_r, pt_l

    # @ut.timing
    def grow(self,
             strm,
             vtag=None,
             direction='forth',
             dL=0.02,
             check_proximity=True,
             exit=True,
             max_iter=2000,
             pull_factor=0.2,
             o_prox=1.0,
             s_prox=1.0,
             output=False,
             target_length=None
             ):

        # 1. set major parameters
        proj_dist = dL  # scale for vector to be projected on plane
        pr = pull_factor * dL  # for checking if it's inside the gh mesh
        ds_self = dL * s_prox  # distance to stop agains itself

        count = 0
        nodes = [strm.start]

        # do forward integration
        while True:
            # 0. check length
            if target_length is not None:
                if self.check_length(nodes, target_length) is False:
                    strm.grow = False
                    break

            # 1. get seed point to grow from and initial vector
            nd = nodes[-1]
            strm.update_search_tree()
            self.get_growth_direction(nd, count, direction)

            # 2. get threshold distances
            max_dist = self.get_threshold_distance(nd)

            # 3. check for proximity to other streamlines
            if check_proximity is True:
                if self.check_proximity(strm, nd, max_dist, o_prox, ds_self) is False:
                    strm.grow = False
                    break

            # 3b. temporary for simmetry:  # temporary
            alpha = 0.00
            if nd.x <= alpha or nd.y >= alpha:
                # print('larger than {}'.format(alpha))
                break

            # 4. find vector to follow
            if vtag is None:
                vec = nd.vel
            else:
                vec = self.s_mesh.get_vector_on_face(nd.pos, nd.f_id, vtag, nd.vel)
                vec = ut.align_vector(cg.normalize_vector(vec), nd.vel)

            # 5. project vector to face plane and create new point
            vec = self.project_vector_on_mesh(nd, vec, proj_dist)
            new_pos = cg.translate_points([nd.pos], vec)[0]

            # 6. create new displace node
            new_nd = self.clone_node(nd, new_pos, vec)

            # 7. check for intersections with mesh edges
            n_nd, ap, flag = self.intersect(nd, new_nd, strm, new_pos, vec, pr)

            # 8. store new node
            if ap is True:
                strm.nd_search.append(n_nd)  # append to search
                nodes.append(n_nd)

            # 9, test for growth after intersection
            if flag is False:
                strm.grow = False
                break

            # 10. increase the count
            count += 1

            # 11. additional break flag for security
            if exit is True:
                if count >= max_iter:
                    print('count was exceeded')
                    strm.grow = False
                    break

        # 12. the while loop is done, extend correspoinding nodes list
        nodes.pop(0)
        if direction == 'forth':
            strm.nodes_forth.extend(nodes)
        elif direction == 'back':
            strm.nodes_back.extend(nodes)

        if output is True:
            if len(nodes) > 0:
                return nodes[-1].pos
            return None

    def get_growth_direction(self, node, count, direction):
        if count == 0:
            if direction == 'back':
                node.vel = cg.scale_vector(node.vel, -1.)

    def get_threshold_distance(self, node):
        d = self.min_sp
        if self.uni_sp is False:
            d = self.s_mesh.c_mesh.get_face_attribute(node.f_id, self.s_tag)
        return d

    # @ut.timing
    def check_proximity(self, streamline, node, max_dist, other_prox, ds_self):
        if self.is_point_close(node.pos, max_dist * other_prox) is True:
            # print('too close to others')
            return False

        elif streamline.is_point_close(node.pos, ds_self, 5) is True:
            # print('too close to itself')
            return False

        # o_nbrs = node.find_neighbors(self.s_mesh.c_mesh, streamline.nd_search)
        # if o_nbrs:
        #     if ut.delaunay_is_point_close(node.pos, o_nbrs, ds_self) is True:
        #         print('threshold to itself: {}'.format(max_dist * other_prox))
        #         print('too close to itself. delaunay on')
        #         return False
        return True

    def check_length(self, nodes, target_length):
        # print('target length is {}'.format(target_length))
        points = [node.pos for node in nodes]
        polyline = Polyline(points)
        poly_length = polyline.length

        # print('polyline length is: {}'.format(poly_length))
        if poly_length > target_length:
            # print('it is too long')
            return False
        return True

    def is_point_close(self, point, distance, exclude=None):
        nnbr, label, dst = self.tree.nearest_neighbor(point, exclude)
        if dst is not None:
            if dst < distance:
                # print('Distance to KDTree is {}'.format(dst))
                return True
        return False

    def update_search_tree(self):
        self.tree = KDTree([nd.pos for nd in self.objs])

    def project_vector_on_mesh(self, node, vector, length):
        plane = self.s_mesh.c_mesh.get_face_attribute(node.f_id, 'plane')
        temp = cg.translate_points([node.pos], cg.scale_vector(vector, 1))
        proj = cg.project_points_plane(temp, plane)[0]
        proj = cg.normalize_vector(cg.vector_from_points(node.pos, proj))
        return cg.scale_vector(proj, length)

    def clone_node(self, node, new_position, new_velocity):
        new_node = Node()
        new_node.clone_from(node)
        new_node.update_pos(new_position)
        new_node.update_vel(new_velocity)
        return new_node

    # @ut.timing
    def intersect(self, node, new_node, streamline, new_position, vector, pr):
        ints = self.get_edge_intersections(node, new_node)

        if ints is not None:
            int_pos, u, v, a_node = ints
            new_node.edge = (u, v)
            new_node.is_intersection = True  # it is an intersection

            # find new face index
            new_faces = list(filter(lambda x: x != node.f_id,
                                    self.s_mesh.c_mesh.edge_faces(u, v)
                                    )
                             )
            # test against found faces
            if new_faces[0] is None:  # it hit the boundary
                # print('It hit the boundary. It stopped.')
                new_node.update_pos(int_pos)
                new_node.f_id = new_faces[0]
                streamline.nd_search.append(new_node)  # append to search

                # nodes.append(new_node)   # to return?
                # streamline.grow = False  # to return?
                # break

                return new_node, True, False

            else:  # it is an edge intersection
                new_node.update_pos(int_pos)
                new_node.f_id = new_faces[0]  # to return?
                # streamline.grow = True

                return new_node, True, True

        else:
            tpt, tid = self.s_mesh.closest_point(new_node.pos)
            dst = cg.distance_point_point(tpt, new_node.pos)

            if dst < pr:  # it is on a mesh face, but no edge intersection
                new_node = Node()
                new_node.clone_from(node)
                new_node.update_pos(new_position)  # input
                new_node.f_id = tid
                new_node.update_vel(vector)  # input

                # streamline.grow = True

                return new_node, True, True

            else:  # it is outside of the mesh
                # print('dst is: {}'.format(dst))
                # print('OUT. distance is larger than tol')

                # streamline.grow = False
                # break

                return None, False, False

    def get_edge_intersections(self, node, other_node):
        face_edges = self.s_mesh.c_mesh.face_halfedges(node.f_id)
        line = [node.pos, other_node.pos]
        ints = self.find_edge_intersections(line, face_edges, other_node)
        return ints

    def find_edge_intersections(self, line, edges, node):
        tol = 1e-6
        start = line[0]
        intersections = []

        for u, v in edges:
            plane = self.s_mesh.c_mesh.get_edge_attribute((u, v), 'plane')
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

    def make_new_node(self, point, velocity=None):
        new_node = Node()
        new_node.update_pos(point)
        new_node.pull_to_mesh(self.s_mesh, self.v_tag)

        if velocity is not None:
            new_node.update_vel(velocity)

        return new_node

    # @ut.timing
    def make_sampling_nodes(self, strm, threshold_sp, length=0.2, heap=True):
        # samp_pts, samp_vels = strm.resample_polyline(length)
        samp_pts, samp_vels = strm.resample_curve(length)

        if samp_pts is not None:
            sampling = []
            if heap is True:
                heapq.heapify(sampling)

            for idx, point in enumerate(samp_pts):
                s_nd = self.make_new_node(point, samp_vels[idx])
                s = self.s_mesh.c_mesh.get_face_attribute(s_nd.f_id, self.s_tag)

                if s <= threshold_sp:
                    if heap is True:
                        heapq.heappush(sampling, (s, s_nd))
                    else:
                        sampling.append((s, s_nd))
                else:
                    # print('spacing is {}'.format(s))
                    # print('threshold is {}'.format(threshold_sp))
                    print('spacing larger than given threshold')
            return sampling
        return None

    def get_offset_vector(self, node):
        self.seed_pts.append(rs.AddPoint(*node.pos))
        normal = self.s_mesh.c_mesh.face_normal(node.f_id)
        return cg.normalize_vector(cg.cross_vectors(node.vel, normal))

    # @ut.timing
    def get_offset_point(self, node, offset, reverse=False):
        # print('offset distance is {}'.format(offset))
        offset_vector = self.get_offset_vector(node)
        if reverse is True:
            offset_vector = cg.scale_vector(offset_vector, -1.0)
        offset_node = self.make_new_node(node.pos, offset_vector)
        offset_streamline = Streamline(offset_node)
        offset_point = self.grow(offset_streamline,
                                 dL=0.01,
                                 check_proximity=False,
                                 output=True,
                                 target_length=offset
                                 )

        offset_streamline.arrange_nodes()
        offset_streamline.create_polyline()

        if offset_streamline.polyline is not None:
            self.offsets.append(offset_streamline)
        return offset_point

    def get_neighbor_proximity(self, point, distance, exclude=None):
        self.update_search_tree()
        nnbr, label, dst = self.tree.nearest_neighbor(point, exclude)
        if dst is not None:
            if dst < distance:
                # print('Distance at start to KDTree is {}'.format(dst))
                return True
        return False

    def get_max_face_centroid(self):
        f_keys = [f for f in self.s_mesh.c_mesh.faces()]
        sep = self.s_mesh.c_mesh.get_faces_attribute(name=self.s_tag, keys=f_keys)
        min_1 = min(zip(sep, f_keys), key=lambda x: x[0])[1]

        min_nbrs = self.s_mesh.c_mesh.face_neighbors(min_1)
        sep_nbrs = self.s_mesh.c_mesh.get_faces_attribute(name=self.s_tag, keys=min_nbrs)
        min_2 = min(zip(sep_nbrs, min_nbrs), key=lambda x: x[0])[1]

        min_1 = self.s_mesh.c_mesh.face_centroid(min_1)
        min_2 = self.s_mesh.c_mesh.face_centroid(min_2)

        return cg.Line(min_1, min_2).midpoint

    def get_max_face_spacing(self):
        f_keys = [f for f in self.s_mesh.c_mesh.faces()]
        sep = self.s_mesh.c_mesh.get_faces_attribute(name=self.s_tag, keys=f_keys)
        return max(sep)
