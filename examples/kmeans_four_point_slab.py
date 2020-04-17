from compas.geometry import scale_vector
from compas.geometry import add_vectors
from compas.geometry import normalize_vector
from compas.geometry import length_vector


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

    NUM = 7  # number of clusters
    ITERS = 50  # number of iterations

    # ==========================================================================
    # Import mesh
    # ==========================================================================

    mesh = Mesh()
    mesh.load(HERE)
    mesh_unify_cycles(mesh)

    # ==========================================================================
    # Create PS vector lines
    # ==========================================================================

    vector_tag = 'ps_1_bot'
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
    # Set up K-Means algorithm
    # ==========================================================================

    faces = make_faces(str_mesh, vector_tag, weight=False)

    print('faces made')

    clusters = furthest_init(NUM, faces)

    print('clusters made')
    sel_clusters = clusters[-1]
    all_clusters = k_means(sel_clusters, faces, ITERS)

    print('k means was run')

    final_clusters = all_clusters[-1]

    # ==========================================================================
    # Visualization
    # ==========================================================================

    plotter = MeshPlotter(mesh, figsize=(12, 9))

    num = len(final_clusters)
    for idx, cluster in final_clusters.items():
        plotter.draw_faces(keys=cluster.faces_keys, facecolor=i_to_rgb(idx/num))

    plotter.draw_lines(lines)
    plotter.show()
