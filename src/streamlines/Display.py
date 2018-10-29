'''
A bunch of functions tailored to display a variety of Stream objects.
'''


__name__ = "Display"
__author__ = "Rafael Pastrana"
__version__ = "0.0.1"
__creation__ = "2018.08.01"
__date__ = "2018.08.01"


import compas
import compas_rhino
import imp
import Rhino.Geometry as rg
import rhinoscriptsyntax as rs
import Utilities

imp.reload(Utilities)

from compas.datastructures import Mesh
from compas import geometry as cg
from Utilities import Utilities


ut = Utilities()


class Display():

    def __init__(self):
        pass

    def plot_mesh(self, s_mesh, data_tag, mode='max'):
        # modes: max, min, avg, angle_avg
        data = s_mesh.get_data_on_nodes(data_tag, mode)

        rs_points = []
        for v in s_mesh.cMesh.vertices():
            rs_points.append(s_mesh.cMesh.vertex_coordinates(v))

        faces = []
        for idx, f in enumerate(s_mesh.cMesh.faces()):
            f_vts = s_mesh.cMesh.face_vertices(f)
            if idx == 0:
                print('initial indices are: {}'.format(f_vts))

            if len(f_vts) < 4:
                f_vts.append(f_vts[-1])

            faces.append(f_vts)

        n_data = ut.normalise_data(data)
        color_data = ut.colorbar(n_data[0])

        new_gh_mesh = rs.AddMesh(rs_points,
                                 faces,
                                 vertex_normals=None,
                                 texture_coordinates=None,
                                 vertex_colors=None
                                 )

        return new_gh_mesh, data, color_data

    def draw_vectors(self, points, vectors, scale):
        rs_lines = []
        for idx, point in enumerate(points):
            vector_a, vector_b = vectors[idx]

            vector_a = cg.scale_vector(vector_a, scale)
            vector_b = cg.scale_vector(vector_b, scale)

            pt_a = cg.add_vectors(point, vector_a)
            pt_b = cg.add_vectors(point, vector_b)

            line = rs.AddLine(rs.AddPoint(pt_a),    # rs
                              rs.AddPoint(pt_b)
                              )
            rs_lines.append(line)  # rs
        return rs_lines

    def draw_face_vectors(self, s_mesh, name, scale=0.02):
        cts = [s_mesh.cMesh.face_centroid(f) for f in s_mesh.cMesh.faces()]
        vectors_a = s_mesh.cMesh.get_faces_attribute(str(name) + '_a')
        vectors_b = s_mesh.cMesh.get_faces_attribute(str(name) + '_b')
        return self.draw_vectors(cts, zip(vectors_a, vectors_b), scale)

    def draw_node_vectors(self, s_mesh, name, scale=0.02):
        vts = [s_mesh.cMesh.vertex_coordinates(v) for v in s_mesh.cMesh.vertices()]
        vectors_a = s_mesh.cMesh.get_vertices_attribute(str(name) + '_a')
        vectors_b = s_mesh.cMesh.get_vertices_attribute(str(name) + '_b')
        return self.draw_vectors(vts, zip(vectors_a, vectors_b), scale)
