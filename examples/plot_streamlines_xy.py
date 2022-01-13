import numpy as np

from random import choice

from compas.geometry import scale_vector
from compas.geometry import add_vectors
from compas.geometry import normalize_vector
from compas.geometry import length_vector
from compas.geometry import closest_point_on_plane
from compas.geometry import distance_point_point

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


    HERE = "/Users/arpj/code/libraries/directional_clustering/data/json/four_point_slab.json"


    TAGS = [
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

    INTEGRATOR = "mebarki"  # mebarki, jobard
    VF_NAMES = ["m_1", "m_2"]
    COLORS = iter([(255, 0, 0), (0, 0, 255)])

    STEP_SIZE = 0.1
    LENGTH_MIN = 0.1
    SPACING_MIN = 0.2
    SPACING_UNIFORM = True

    # spacing
    MIN_DIST_OTHERS = SPACING_MIN
    MIN_DIST_OTHERS_START = SPACING_MIN
    MIN_DIST_SELF = SPACING_MIN

    # jobard integrator settings
    NUM_SAMPLES = 5
    ITERS = 10
    FIRST_SEED = [0.0, 0.0, 0.0]
    HEAP = False

    # plotter settings
    DRAW_EDGES = True
    DRAW_FACES = False
    DRAW_VF = False

    # ==========================================================================
    # Import mesh
    # ==========================================================================

    mesh = Mesh()
    mesh = Mesh.from_json(HERE)
    mesh_unify_cycles(mesh)

    # ==========================================================================
    # Visualization
    # ==========================================================================

    plotter = MeshPlotter(mesh, figsize=(12,9))

    if DRAW_EDGES:
        plotter.draw_edges(color=(100, 100, 100), width=0.1)

    if DRAW_FACES:
        plotter.draw_faces()

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

    seeds = [closest_pt]

    # ==========================================================================
    # Create linear space grid
    # ==========================================================================

    num_x = 3  # 8
    num_y = 10  # 8

    P = np.array([mesh.vertex_coordinates(v) for v in mesh.vertices()])

    max_x = np.amax(P[:, 0])
    min_x = np.amin(P[:, 0])
    max_y = np.amax(P[:, 1])
    min_y = np.amin(P[:, 1])

    P = np.zeros((num_x, num_y, 3))

    x_space = np.linspace(min_x, max_x, num_x).reshape((-1, 1))
    y_space = np.linspace(min_y, max_y, num_y).reshape((1, -1))

    P[:, :, 0] = x_space
    P[:, :, 1] = y_space

    seeds = list(flatten(P.tolist()))

    # ==========================================================================
    # Other seeds
    # ==========================================================================

    seeds = [choice(str_mesh.face_centroids) for i in range(10)]
    seeds = str_mesh.boundary_polygon[:]

    # ==========================================================================
    # Trace streamlines
    # ==========================================================================

    for vf_name in VF_NAMES:

        print(f"Integrating vector field: {vf_name}")

    # ==========================================================================
    # Set up Streamsystem()
    # ==========================================================================

        streamsystem = Streamsystem(str_mesh, dL=STEP_SIZE, min_sp=SPACING_MIN, uni_sp=SPACING_UNIFORM)
        streamsystem.set_tracing_data(vf_name, [0, 0, 0], min_length=LENGTH_MIN)

    # ==========================================================================
    # Execute tracing routine
    # ==========================================================================

        start_time = time.time()
        if INTEGRATOR == "mebarki":
            streamsystem.make_streamlines_mebarki(seeds, o_prox=MIN_DIST_OTHERS, st_o_prox=MIN_DIST_OTHERS_START)
        elif INTEGRATOR == "jobard":
            streamsystem.make_streamlines_jobard(min_sp=SPACING_MIN,
                                                 o_prox=MIN_DIST_OTHERS,
                                                 st_o_prox=MIN_DIST_OTHERS_START,
                                                 s_prox=MIN_DIST_SELF,
                                                 num_samples=NUM_SAMPLES,
                                                 iters=ITERS,
                                                 start_pt=FIRST_SEED,
                                                 heap=HEAP)
        else:
            raise ValueError
        end_time = time.time()
        print('elapsed time: {} seconds'.format(time.time() - start_time))

    # ==========================================================================
    # Visualization
    # ==========================================================================

        polylines = []
        color = next(COLORS)
        for streamline in streamsystem.streamlines:
            polylines.append({'points': streamline.polyline.points, 'color': color})
        plotter.draw_polylines(polylines)

    # ==========================================================================
    # Create PS vector lines
    # ==========================================================================

        if DRAW_VF:
            lines = vector_lines_on_faces(mesh, vf_name, True, factor=0.02)
            lines = [line for line in map(line_tuple_to_dict, lines)]
            for line in lines:
                line['width'] = 0.50

            plotter.draw_lines(lines)

    # ==========================================================================
    # Show
    # ==========================================================================

    plotter.show()
