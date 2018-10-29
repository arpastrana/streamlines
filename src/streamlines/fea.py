'''
A bunch of utilities that help to communicate with abaqus through compas FEA.
'''

__author__ = 'Rafael Pastrana, Andrew Liew'
__name__ = 'FEA Communication'
__version__ = '0.0.2'
__date__ = '16.08.27'


import Rhino
import compas.geometry as cg
import rhinoscriptsyntax as rs
from compas_fea.cad import rhino
from compas_fea.structure import Structure
from compas_fea.structure import ShellSection
from compas_fea.structure import Concrete
from compas_fea.structure import Steel
from compas_fea.structure import ElementProperties as Properties


def clean_vectors(vectors):
    cleaned_vectors = []
    for vector in vectors:
        try:
            vector = [vector.X, vector.Y, vector.Z]
        except Exception:
            vector = None
        cleaned_vectors.append(vector)
    return cleaned_vectors


def steel_area_to_spacing(s_area, diameter):
    TOL = 1e-3
    bar_area = get_bar_area(diameter)

    if s_area > TOL:
        spacing = s_area / bar_area
        spacing = 1 / spacing
    else:
        return None

    def get_bar_area(diameter):
        return math.pi * math.pow((diameter / 10), 2) / 4

    return spacing


def get_rebar_angles(x_dir, ori):
    angle = cg.angle_vectors(ori, x_dir, deg=True)
    angle_2 = angle + 90.0

    # if angle > 90.0:
    #     angle = angle-90
    #     angle_2 = angle - 90.0

    return angle, angle_2


def get_rebar_spacing(mesh, fkey):
    TAG = 'spacing_'
    AS_TAGS = ['asx_', 'asy_']
    SF = 1.0

    diam_x = mesh.get_face_attribute(fkey, 'bar_diameter_x') / 1000.
    diam_y = mesh.get_face_attribute(fkey, 'bar_diameter_y') / 1000.


    sp_x = mesh.get_face_attribute(fkey, TAG + AS_TAGS[0]) * SF
    sp_y = mesh.get_face_attribute(fkey, TAG + AS_TAGS[1]) * SF
    return [diam_x, diam_y], [sp_x, sp_y]


def get_rc_property(fkey, s_sec, m_conc, m_reb, diam, spac, he, ang):
    pr_n = 'ep_' + str(fkey)
    rebar = get_rebar(m_reb, diam, spac, he, ang)
    ep = Properties(name=pr_n, material=m_conc, section=s_sec, elements=[fkey], reinforcement=rebar)

    if fkey == 0:
        # print(angle, t, cover, sp_x_b, sp_x_t, diam_x)
        print(rebar)
    return ep


def get_rebar(mat_rebar, diameter, spacing, height, angle):
    rebar = {
            'r_x': {'pos': height[0],
                    'spacing': spacing[0],
                    'material': mat_rebar,
                    'dia': diameter[0],
                    'angle': angle[0]
                    },
            'r_y': {'pos': height[1],
                    'spacing': spacing[1],
                    'material': mat_rebar,
                    'dia': diameter[1],
                    'angle': angle[1]
                    }
            }
    return rebar


def add_nodes_elements_from_mesh(stru, mesh, s_sec, m_conc, m_reb, angles, heights, elset):

    for key in sorted(list(mesh.vertices()), key=int):
        stru.add_node(mesh.vertex_coordinates(key))

    ek = []
    for fkey in mesh.faces():

        f = []
        for i in mesh.face[fkey]:
            f.append(stru.check_node_exists(mesh.vertex_coordinates(i)))
        ek.append(stru.add_element(nodes=f, type='ShellElement'))

        # add rc properties
        diam, spac = get_rebar_spacing(mesh, fkey)
        ang = angles[fkey]
        he = heights[fkey]
        ep = get_rc_property(fkey, s_sec, m_conc, m_reb, diam, spac, he, ang)
        stru.add_element_properties(ep)

    stru.add_set(name=elset, type='element', selection=ek)

    return ek


def add_nodes_elements_from_mesh_simple(stru, mesh, elset):
    for key in sorted(list(mesh.vertices()), key=int):
        stru.add_node(mesh.vertex_coordinates(key))
    ek = []

    for fkey in mesh.faces():

        f = []
        for i in mesh.face[fkey]:
            f.append(stru.check_node_exists(mesh.vertex_coordinates(i)))

        ek.append(stru.add_element(nodes=f, type='ShellElement'))

    if elset:
        stru.add_set(name=elset, type='element', selection=ek)

    return ek


def add_element_set(structure, guids, name):
    added_ele = set()

    for guid in guids:
        if rs.IsMesh(guid):
            vertices = rs.MeshVertices(guid)
            faces = rs.MeshFaceVertices(guid)
            nodes = [structure.add_node(vertex) for vertex in vertices]

            for f in rs.MeshFaceVertices(guid):
                nodes = [structure.check_node_exists(vertices[i]) for i in f]

                if nodes[-1] == nodes[-2]:
                    del nodes[-1]

                ekey = structure.add_element(nodes=nodes, type='ShellElement')
                if ekey is not None:
                    added_ele.add(ekey)

    structure.add_set(name=name, type='element', selection=list(added_ele))


def add_points_sets(structure, points, names):
    for idx, point_list in enumerate(points.Branches):
        name = names[idx]
        check_points = [rs.IsPoint(pt) for pt in point_list]
        if all(check_points):
            add_node_set(structure=structure, pt_guids=point_list, name=name)
        else:
            print('*****.set not created *****'.format(name))


def add_node_set(structure, pt_guids, name):
    nodes = []
    for pt_guid in pt_guids:
        if rs.IsPoint(pt_guid):
            node = structure.check_node_exists(rs.PointCoordinates(pt_guid))
            if node is not None:
                nodes.append(node)
    structure.add_set(name=name, type='node', selection=nodes)
