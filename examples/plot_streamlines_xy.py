import numpy as np

from random import seed
from random import choice

from compas.geometry import scale_vector
from compas.geometry import add_vectors
from compas.geometry import normalize_vector
from compas.geometry import length_vector
from compas.geometry import closest_point_on_plane
from compas.geometry import distance_point_point
from compas.geometry import cross_vectors

from compas.geometry import is_point_in_triangle
from compas.geometry import is_point_in_triangle_xy

from compas.utilities import flatten


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


def barycentric_coordinates(point, triangle, clamp=False):
    '''
    '''
    a, b, c = triangle
    pt = point

    def clamper(value):
        if not clamp:
            return value
        if value < 0.0:
            return 0.0
        elif value > 1.0:
            return 1.0
        return value
    
    def barycentric_1():
        numerator = (b[1] - c[1]) * (pt[0] - c[0]) + (c[0] - b[0]) * (pt[1] - c[1])
        denominator = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
        return clamper(numerator / denominator)

    def barycentric_2():
        numerator = (c[1] - a[1]) * (pt[0] - c[0]) + (a[0] - c[0]) * (pt[1] - c[1])
        denominator = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
        return clamper(numerator / denominator)

    bar_1 = barycentric_1()
    bar_2 = barycentric_2()
    bar_3 = 1 - bar_1 - bar_2
    return [bar_1, bar_2, bar_3]


def barycentric_to_cartesian(barycentric, triangle):
    '''
    '''
    bar_1, bar_2, bar_3 = barycentric
    a, b, c = triangle

    x = a[0] * bar_1 + b[0] * bar_2 + c[0] * bar_3
    y = a[1] * bar_1 + b[1] * bar_2 + c[1] * bar_3

    return [x, y, 0.0]


def is_point_in_triangle_barycentric_xy(point, triangle):
    '''
    '''
    if length_vector(barycentric_coordinates(point, triangle)) <= 1.0:
        return True
    return False


def trimesh_closest_point_xy(mesh, point):
    '''
    '''
    closest_pts = []

    for fkey in mesh.faces():
        triangle = mesh.face_coordinates(fkey)
        plane = (mesh.face_centroid(fkey), mesh.face_normal(fkey))
        closest_pt = closest_point_on_plane(point, plane)
        
        bars = barycentric_coordinates(closest_pt, triangle, clamp=True)
        closest_pt = barycentric_to_cartesian(bars, triangle)

        closest_pts.append((closest_pt, fkey, distance_point_point(point, closest_pt)))

    if not closest_pts:
        return None, None, None

    return sorted(closest_pts, key=lambda x: x[2])[0]


if __name__ == '__main__':

    import time

    from compas.datastructures import Mesh
    from compas.datastructures import mesh_unify_cycles

    from compas_plotters import MeshPlotter

    from streamlines.streamsystem import Streamsystem
    from streamlines.custom_mesh import StructuralMesh


    # HERE = '/Users/arpj/code/princeton/directional_clustering/data/json_files/four_point_slab.json'
    
    # HERE = '/Users/arpj/code/princeton/directional_clustering/data/json_files/two_point_wall_res_005.json'

    # HERE =
    # '/Users/arpj/code/princeton/directional_clustering/data/json_files/two_point_wall_res_005_k_5.json'
    
    # HERE = "/Users/arpj/code/princeton/directional_clustering/data/json_files/square_wall_cantilever_res_005_k_3.json"
    HERE = "/Users/arpj/code/libraries/directional_clustering/data/json/four_point_slab.json"


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


    # vector_tag = 'n_2_k_3' # for walls
    vector_tag = 'm_1'  # for slabs

    # ==========================================================================
    # Import mesh
    # ==========================================================================

    mesh = Mesh()
    mesh = Mesh.from_json(HERE)
    mesh_unify_cycles(mesh)

    # ==========================================================================
    # Instantiate StructuralMesh()
    # ==========================================================================

    str_mesh = StructuralMesh(mesh)

    # ==========================================================================
    # Create closest-point seeds
    # ==========================================================================

    test_pt = [3.0, 2.5, 0.5]
    test_pt = [0.5, 0.5, 0.0]
    output = trimesh_closest_point_xy(mesh, test_pt)
    closest_pt, fkey, dist = output

    seed(0)

    seeds = [closest_pt]

    # ==========================================================================
    # Create linear space grid
    # ==========================================================================

    # num_x = 3  # 8
    # num_y = 10  # 8

    # P = np.array([mesh.vertex_coordinates(v) for v in mesh.vertices()])

    # max_x = np.amax(P[:, 0])
    # min_x = np.amin(P[:, 0])
    # max_y = np.amax(P[:, 1])
    # min_y = np.amin(P[:, 1])

    # P = np.zeros((num_x, num_y, 3))

    # x_space = np.linspace(min_x, max_x, num_x).reshape((-1, 1))
    # y_space = np.linspace(min_y, max_y, num_y).reshape((1, -1))

    # P[:, :, 0] = x_space
    # P[:, :, 1] = y_space

    # seeds = list(flatten(P.tolist()))
    # seed_points = seeds[:]

    # ==========================================================================
    # Other seeds
    # ==========================================================================

    # seeds = [mesh.vertex_coordinates(500)]
    # seeds = [choice(str_mesh.face_centroids) for i in range(10)]
    # seeds = str_mesh.boundary_polygon[:]

    # ==========================================================================
    # Second field seeds
    # ==========================================================================

    seeds_2 = seeds[:]
    
    # ==========================================================================
    # Set up Streamsystem()
    # ==========================================================================
    
    streamsystem = Streamsystem(str_mesh, dL=0.10, min_sp=0.20, uni_sp=True)
    streamsystem.set_tracing_data(vector_tag, [0, 0, 0], min_length=0.10)

    # ==========================================================================
    # Execute tracing routine
    # ==========================================================================
    
    start_time = time.time()
    streamsystem.make_streamlines_mebarki(seeds, o_prox=0.2, st_o_prox=0.2)
    end_time = time.time()
    print('elapsed time: {} seconds'.format(time.time() - start_time))

    # ==========================================================================
    # Visualization
    # ==========================================================================
    
    polylines = []
    control_points = []
    for streamline in streamsystem.streamlines:
        polylines.append({'points': streamline.polyline.points,
                         'color': (0, 0, 255)})

        for idx, xyz in enumerate(streamline.polyline.points):
            control_points.append({'pos': xyz, 'facecolor': (255, 255, 255), 'radius': 0.0025, 'text': str(idx)})

    # ==========================================================================
    # Rotate field
    # ==========================================================================

    for fkey in str_mesh.c_mesh.faces():
        vector = str_mesh.c_mesh.face_attribute(key=fkey, name=vector_tag)
        global_z = [0.0, 0.0, 1.0]
        other_vector = cross_vectors(vector, global_z)
        o_vector_tag = "x_" + vector_tag
        str_mesh.c_mesh.face_attribute(key=fkey, name=o_vector_tag, value=other_vector)

    # ==========================================================================
    # Set up Streamsystem() 2
    # ==========================================================================
    
    streamsystem_2 = Streamsystem(str_mesh, dL=0.10, min_sp=0.20, uni_sp=True)
    streamsystem_2.set_tracing_data(o_vector_tag, [0, 0, 0], min_length=0.10)

    # ==========================================================================
    # Execute tracing routine 2
    # ==========================================================================
    
    start_time = time.time()
    streamsystem_2.make_streamlines_mebarki(seeds_2, o_prox=0.2, st_o_prox=0.2)
    end_time = time.time()
    print('elapsed time: {} seconds'.format(time.time() - start_time))

    # ==========================================================================
    # Visualization 2
    # ==========================================================================
    
    polylines_2 = []
    for streamline in streamsystem_2.streamlines:
        polylines_2.append({'points': streamline.polyline.points, 'color': (255, 0, 0)})

    # ==========================================================================
    # Create PS vector lines
    # ==========================================================================
    
    lines = vector_lines_on_faces(mesh, vector_tag, True, factor=0.02)
    lines = [line for line in map(line_tuple_to_dict, lines)]
    for line in lines:
        line['width'] = 0.50

    # ==========================================================================
    # Visualization
    # ==========================================================================
    
    points = []

    for pt in seeds:
        points.append({"pos": pt, "radius": 0.05, "facecolor": (255, 255, 255)})

    # ==========================================================================
    # Visualization
    # ==========================================================================

    plotter = MeshPlotter(mesh, figsize=(12,9))

    # plotter.draw_edges(color=(10, 10, 10))
    plotter.draw_faces()
    # plotter.draw_lines(lines)

    plotter.draw_points(points)

    plotter.draw_polylines(polylines)
    plotter.draw_polylines(polylines_2)

    plotter.show()
