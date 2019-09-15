import heapq


from compas.geometry import constrained_delaunay_triangle
from compas.geometry import scale_vector
from compas.geometry import add_vectors
from compas.geometry import normalize_vector
from compas.geometry import length_vector
from compas.geometry import closest_point_on_plane
from compas.geometry import distance_point_point

from compas.geometry import is_point_in_triangle
from compas.geometry import is_point_in_triangle_xy


HERE = '/Users/arpj/code/libraries/streamlines/examples/four_point_slab.json'

'''
tags = [
    n_1
    n_2
    m_1
    m_2
    ps_1_top
    ps_1_bot
    ps_1_mid
    ps_2_top
    ps_2_bot
    ps_2_mid
    custom_1
    custom_2
]
'''


def vector_lines_on_faces(mesh, vector_tag, length=0.02):
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

        pt = mesh.face_centroid(fkey)
        lines.append(line_sdl(pt, vector, length))

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
    # heapq.heapify(closest_points)

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

    from compas.datastructures import Mesh
    from compas.datastructures import mesh_unify_cycles

    from compas_plotters import MeshPlotter


    mesh = Mesh()
    mesh.load(HERE)
    mesh_unify_cycles(mesh)

    lines = vector_lines_on_faces(mesh, 'ps_1_bot', length=0.05)
    lines = [line for line in map(line_tuple_to_dict, lines)]

    test_pt = [3.0, 2.5, 0.5]
    # test_pt = [2.11956516901652, 2.9333333174387612, 0.0]
    output = trimesh_closest_point_xy(mesh, test_pt)
    closest_pt, fkey, dist = output
    print(output)

    plotter = MeshPlotter(mesh, figsize=(12,9))
    plotter.draw_faces()
    plotter.draw_lines(lines)
    plotter.draw_points([{'pos': closest_pt, 'facecolor': (255, 0, 0), 'radius': 0.02},
                        {'pos': test_pt, 'facecolor': (0, 0, 255), 'radius': 0.02}
                        ]
                        )
    plotter.show()
