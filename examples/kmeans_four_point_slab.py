import math

from random import choice
from time import time
from functools import partial

from compas.geometry import scale_vector
from compas.geometry import add_vectors
from compas.geometry import normalize_vector
from compas.geometry import length_vector
from compas.geometry import cross_vectors


def vector_lines_on_faces(mesh, vector_tag, uniform=True, factor=0.02):
    '''
    '''
    def line_sdl(start, direction, length):
        direction = normalize_vector(direction[:])
        a = add_vectors(start, scale_vector(direction, -length))
        b = add_vectors(start, scale_vector(direction, +length))

        return a, b

    lines = []
    for fkey, attr in mesh.faces(data=True):
        vector = attr.get(vector_tag)

        if not vector:
            raise ValueError('Vector {} not defined on face {}'.format(vector_tag, fkey))

        if uniform:
            vec_length = factor
        else:
            vec_length = length_vector(vector) * factor

        pt = mesh.face_centroid(fkey)
        lines.append(line_sdl(pt, vector, vec_length))
    return lines


def line_tuple_to_dict(line):
    '''
    '''
    a, b = line
    return {'start': a, 'end': b}


if __name__ == '__main__':

    from compas.datastructures import Mesh
    from compas.datastructures import mesh_unify_cycles

    from compas.utilities import i_to_rgb
    from compas.utilities import i_to_black

    from compas_plotters import MeshPlotter

    from streamlines.custom_mesh import StructuralMesh

    from streamlines.kmeans import make_faces
    from streamlines.kmeans import furthest_init
    from streamlines.kmeans import k_means

    # ==========================================================================
    # Constants
    # ==========================================================================

    HERE = '/Users/arpj/code/libraries/streamlines/examples/four_point_slab.json'

    tags = [
        'n_1',
        'n_2',
        'm_1',
        'm_2',
        'ps_1_top',
        'ps_1_bot',
        'ps_1_mid',
        'ps_2_top',
        'ps_2_bot',
        'ps_2_mid',
        'custom_1',
        'custom_2'
        ]

    # (odd numbers only!) (after 11, starts to get confused!) but at 19, kind of works again
    NUM = 13  # number of clusters 
    ITERS = 30  # number of iterations
    MERGESPLIT = True  # merge split in k means. True is good for this example, but not for knitcandela!
    EXPORT_PNG = False

    THERE = '/Users/arpj/code/libraries/streamlines/examples/gif_{0}_{1}/kmeans_{0}_{1}_'
    THERE = THERE.format(NUM, ITERS)

    vector_tag_1 = 'ps_1_top'  # ps_1_top
    vector_tag_2 = 'ps_2_top'  # ps_1_top
    vector_tag = 'ps_12_top'  # ps_1_top

    # ==========================================================================
    # Import mesh
    # ==========================================================================

    mesh = Mesh()
    mesh.load(HERE)
    mesh_unify_cycles(mesh)

    # ==========================================================================
    # rebuild mesh
    # ==========================================================================

    # new_mesh = Mesh()

    # all_vertices = set()
    # for idx, tup in enumerate(mesh.faces(True)):
    #     fkey, attr = tup

    #     # 4.5 x 6.0 m rectancle
    #     if mesh.face_centroid(fkey)[0] < -0.05:  # x - mesh deleter by symmetry
    #         continue
    #     if mesh.face_centroid(fkey)[1] < 0.:  # y - mesh deleter by symmetry
    #         continue

    #     attr_dict = {k:v for k, v in attr.items()}
    #     face = mesh.face_vertices(fkey)
    #     new_mesh.add_face(key=idx, vertices=face, attr_dict=attr_dict)
    #     all_vertices.update(face)

    # for vkey, attr in mesh.vertices(True):
    #     if vkey not in all_vertices:
    #         continue
    #     attr_dict = {k:v for k, v in attr.items()}
    #     new_mesh.add_vertex(vkey, attr_dict=attr_dict)

    # mesh = new_mesh

    # ==========================================================================
    # 45 degrees field
    # ==========================================================================

    for fkey, attr in mesh.faces(True):
        vec_1 = attr[vector_tag_1]
        y = 1.0 / math.tan(math.radians(45.0))
        x_vec = vec_1
        y_vec = cross_vectors(x_vec, [0.0, 0.0, 1.0])  # global Z
        y_vec = scale_vector(y_vec, y)
        vec_3 = normalize_vector(add_vectors(x_vec, y_vec))
        
        mesh.set_face_attribute(fkey, name=vector_tag, value=vec_3)

    # ==========================================================================
    # Create PS vector lines
    # ==========================================================================

    lines = vector_lines_on_faces(mesh, vector_tag, True, factor=0.05)

    lines = [line for line in map(line_tuple_to_dict, lines)]
    for line in lines:
        line['width'] = 0.60

    # ==========================================================================
    # Instantiate StructuralMesh()
    # ==========================================================================

    str_mesh = StructuralMesh(mesh)

    for tag in tags:
        vector_field = mesh.get_faces_attribute(keys=list(mesh.faces()), name=tag)
        str_mesh.set_face_vectors(vector_field, tag, normalize=True)
        str_mesh.set_vertex_vectors_angles(tag)

    # ==========================================================================
    # Define Callback
    # ==========================================================================

    def callback(k, plotter, clusters, filepath, export):
        num = len(list(clusters.keys()))

        facedict = {}

        for idx, cluster in clusters.items():
            color = [i / 255 for i in i_to_rgb(idx / num)]
            for fkey in cluster.faces_keys:
                facedict[fkey] = color

        facecolors = sorted(facedict.items(),  key=lambda x: x[0])
        facecolors = [x[1] for x in facecolors]
        plotter.facecollection.set_facecolors(facecolors)

        if export:
            plotter.save(THERE + '{}_{}.png'.format(time(), k))
        plotter.update(pause=0.50)

    # ==========================================================================
    # Set up Plotter
    # ==========================================================================

    plotter = MeshPlotter(mesh, figsize=(12, 9))
    plotter.draw_lines(lines)
    plotter.draw_faces()
    # plotter.update(pause=0.5)

    callback = partial(callback, plotter=plotter, filepath=THERE, export=EXPORT_PNG)

    # ==========================================================================
    # Set up K-Means algorithm
    # ==========================================================================

    faces = make_faces(str_mesh, vector_tag, weight=False)
    clusters = furthest_init(NUM, faces, callback=None)
    
    sel_clusters = clusters[-1]
    all_clusters = k_means(sel_clusters, faces, ITERS, MERGESPLIT, callback=None)

    final_clusters = all_clusters[-1]

    callback(1, clusters=final_clusters)
    # ==========================================================================
    # Visualization
    # ==========================================================================
    plotter.update(pause=0.5)
    plotter.show()
