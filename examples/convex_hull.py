import random

from compas.datastructures import Mesh

from compas.geometry import distance_point_point
from compas.geometry import convex_hull_numpy
from compas.topology import unify_cycles

from compas_viewers.meshviewer import MeshViewer

radius = 5
origin = (0., 0., 0.)
count = 0
points = []

while count < 10:
    x = (random.random() - 0.5) * radius * 2
    y = (random.random() - 0.5) * radius * 2
    z = (random.random() - 0.5) * radius * 2
    pt = x, y, z

    if distance_point_point(origin, pt) <= radius:
        points.append(pt)
        count += 1

vertices, faces = convex_hull_numpy(points)

i_index = {i: index for index, i in enumerate(vertices)}

vertices = [points[index] for index in vertices]
faces = [[i_index[i] for i in face] for face in faces]
faces = unify_cycles(vertices, faces)

mesh = Mesh.from_vertices_and_faces(vertices, faces)

viewer = MeshViewer()
# viewer.setup()
viewer.mesh = mesh
viewer.show()