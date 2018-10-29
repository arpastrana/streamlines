# ints = self.get_edge_intersections(nd, new_nd)
#
# if ints is not None:
#     int_pos, u, v, nextNd = ints
#     new_nd.edge = (u, v)
#     new_nd.is_intersection = True  # it is an intersection
#     self.int_points.append(rs.AddPoint(*int_pos))
#
#     # find new face index
#     new_faces = list(filter(lambda x: x != nd.f_id,
#                             self.sMesh.cMesh.edge_faces(u, v)
#                             )
#                      )
#     # test against found faces
#     if new_faces[0] is None:
#         new_nd.update_pos(int_pos)
#         new_nd.f_id = new_faces[0]
#         streamline.nd_search.append(new_nd)  # append to search
#         nodes.append(new_nd)
#         print('It hit the boundary')
#         streamline.grow = False
#         break
#
#     else:
#         new_nd.update_pos(int_pos)
#         new_nd.f_id = new_faces[0]
#
# else:
#     tpt, tid = self.sMesh.closest_point(new_nd.pos)
#     dst = cg.distance_point_point(tpt, new_nd.pos)
#
#     if dst < prox:
#         new_nd = Node()
#         new_nd.clone_from(nd)
#         new_nd.update_pos(new_pos)
#         new_nd.f_id = tid
#         new_nd.update_vel(vec)
#     else:
#         print('dst is: {}'.format(dst))
#         print('OUT. distance is larger than tol')
#         streamline.grow = False
#         break
#


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
    start_nd = self.make_new_node(self.get_max_face_centroid())

    # 3. make first streamline
    strm = self.make_streamline(start_nd)
    self.objs[:] = []  # clear objs search
    self.streamlines[:] = []  # clear first streamline

    # 4. append first streamline to the queue
    heapq.heappush(queue, (max_[0], strm))
    count = 0

    # 5. pop streamline-i from queue
    # while len(queue) > 0:
    for i in range(40):
        if len(queue) > 0:
            as_, cu_strm = heapq.heappop(self.queue)
            count += 1

        # 6. sample points on the streamline-i at min spacing OK
            s_ = self.make_sampling_points(cu_strm)

            # 8. select node and offset to the right and left
            # while len(s_) > 0:
            for i in range(4):
                if len(s_) > 0:
                    cu_sep, cu_nd = heapq.heappop(s_)
                    cu_sep = self.min_sp

                    if count != 1:
                        cu_sep = cu_sep * 2

                # 9. make points with particle tracing
                    pt_right = self.get_offset_point(cu_nd, cu_sep)
                    pt_left = self.get_offset_point(cu_nd, cu_sep, True)

                # 10. make streamline-j. append to queue
                    for pt_ in [pt_right, pt_left]:
                        if pt is not None:
                            new_streamline = self.get_new_streamline(pt_)
                        if new_streamline is not None:
                            heapq.heappush(queue, (error, new_streamline))
