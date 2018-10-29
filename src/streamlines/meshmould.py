'''
An text interface to communicate with the Mesh Mould technology.
'''


__name__ = "Mesh Mould Interface"
__author__ = "Rafael Pastrana"
__version__ = "0.0.1"
__creation__ = "2018.09.16"
__date__ = "2018.09.16"

import rhinoscriptsyntax as rs
import compas
import compas.geometry as cg
import os.path


def import_text_file(path, name):
    filepath = path + name + '.txt'

    with open(filepath, 'r') as f:
        raw_data = f.readlines()

    raw_data = [x[1:-2] for x in raw_data]
    raw_data = [x.split(', ') for x in raw_data]
    data = [x for x in raw_data if len(x) > 2]

    booleans = [x for x in raw_data if len(x) < 3]
    booleans = [x for x in booleans if len(x) > 1]

    points = []
    for line in data:
        polygon = []
        index = 0
        lim = 3

        for num in range(4):
            polygon.append(list(map(lambda x: float(x), line[index:index + lim])))
            index += lim

        points.append(polygon)
    return points, booleans


def export_cells(columns, points, filepath, tol=0.001):
    search_tree = make_search_tree(points)
    with open(filepath, "w") as f:

        for idx, cells in enumerate(columns.Branches):
            for j, cell in enumerate(cells):

                vts = get_vertices_coordinates(cell)
                f.write(str(vts))
                f.write("\n")

                bool = get_booleans(cell, search_tree)
                f.write(str(bool))
                f.write("\n")

            f.write("\n")
        f.close()
    print('exported!')

# write text
def export_text_file(self,my_list,filepath):
    # write points
    # write booleans
    # if branch is over, insert an empty row
    # save
    with open(filepath, "w", encoding='utf-8') as f:
        f.write(str(my_list))


# make compas point from rhino point
def make_c_pt(rh_pt):
    if rh_pt is not None:
        return [rh_pt.X, rh_pt.Y, rh_pt.Z]


# make kd-tree
def make_search_tree(rh_pts):
    c_pts = list(map(lambda x: make_c_pt(x), rh_pts))
    tree = cg.KDTree(c_pts)
    return tree


# make booleans
def get_booleans(curve, tree, tol=0.001):
    # explode curve and get segments
    segs = rs.ExplodeCurves(curve, delete_input=False)
    mid = list(map(lambda x: rs.CurveMidPoint(x), segs))
    target = list(map(lambda x: make_c_pt(x), [mid[0], mid[2]]))

    booleans = [None, None]
    for idx, t in enumerate(target):
        nbr, key, d = tree.nearest_neighbour(t)
        bool = 0

        if d <= tol:
            bool = 1
        booleans[idx] = bool

    return booleans


# transform cells to points
def get_vertices_coordinates(curve):
    # explode curve and get vertices
    coordinates = []
    vertices = rs.PolylineVertices(curve)
    vertices = vertices[:-1]

    for vertex in vertices:
        c_pt = make_c_pt(vertex)
        coordinates.extend(c_pt)

    print('number of points is {}, number of coordinates is {}'.format(len(vertices), len(coordinates)))
    return coordinates
