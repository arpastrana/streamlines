import math
import compas
import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial import distance
from scipy.sparse.linalg import eigs
from scipy.sparse.linalg import eigsh
from scipy.sparse.csgraph import laplacian
from scipy.sparse import csr_matrix

from sklearn.metrics.pairwise import euclidean_distances
from sklearn.neighbors import kneighbors_graph
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler

from time import time
from functools import partial

from compas.geometry import scale_vector
from compas.geometry import add_vectors
from compas.geometry import normalize_vector
from compas.geometry import cross_vectors
from compas.geometry import length_vector
from compas.geometry import angle_vectors
from compas.geometry import distance_point_point

from compas.datastructures import Mesh
from compas.datastructures import Network
from compas.datastructures import mesh_unify_cycles
from compas.datastructures import mesh_dual

from compas.utilities import i_to_rgb

from compas_plotters import MeshPlotter
from compas_plotters import NetworkPlotter

from streamlines.utilities import Utilities


ut = Utilities()

# ==========================================================================
# Globals
# ==========================================================================

HERE = '/Users/arpj/code/libraries/streamlines/examples/three_point_circle_slab.json'

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

vector_tag_1 = 'ps_1_top'  # ps_1_top
vector_tag_2 = 'ps_2_top'  # ps_1_top
vector_tag = 'ps_12_top'  # ps_1_top
smooth_iters = 0


# ==========================================================================
# Import mesh
# ==========================================================================

mesh = Mesh()
mesh.load(HERE)

# ==========================================================================
# rebuild mesh
# ==========================================================================

new_mesh = Mesh()

all_vertices = set()
for idx, tup in enumerate(mesh.faces(True)):
	fkey, attr = tup

	attr_dict = {k:v for k, v in attr.items()}
	face = mesh.face_vertices(fkey)
	new_mesh.add_face(key=idx, vertices=face, attr_dict=attr_dict)
	all_vertices.update(face)

for vkey, attr in mesh.vertices(True):
	if vkey not in all_vertices:
		continue
	attr_dict = {k:v for k, v in attr.items()}
	new_mesh.add_vertex(vkey, attr_dict=attr_dict)

mesh = new_mesh

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
# Average smooth vector field
# ==========================================================================

for _ in range(smooth_iters):
	averaged_vectors = {}
	for fkey in mesh.faces():
		nbrs = mesh.face_neighbors(fkey)
		vectors = mesh.get_faces_attribute(keys=nbrs, name=vector_tag)
		vectors.append(mesh.get_face_attribute(fkey, name=vector_tag))

		vectors = list(map(lambda x: ut.align_vector(x, vectors[0]), vectors))

		vectors = np.array(vectors)
		avg_vector = np.mean(vectors, axis=0).tolist()
		averaged_vectors[fkey] = avg_vector

	for fkey in mesh.faces():
		mesh.set_face_attribute(fkey, name=vector_tag, value=averaged_vectors[fkey])

# ==========================================================================
# Process PS vectors
# ==========================================================================

angles = {}
centroids = {}
for fkey, attr in mesh.faces(data=True):
	vector = attr.get(vector_tag)
	angle = angle_vectors([1.0, 1.0, 0.0], vector, deg=True)
	angles[fkey] = angle
	centroids[fkey] = np.array(mesh.face_centroid(fkey))

print('max angle', min(angles.values()))
print('min angle', max(angles.values()))

for idx, angle in angles.items():
	if angle >= 90.0:
		angle = 180.0 - angle
	angles[idx] = angle

for idx, angle in angles.items():
	if angle >= 45:
		angle = 90.0 - angle
	angles[idx] = angle

print('max angle', min(angles.values()))
print('min angle', max(angles.values()))


# ==========================================================================
# Create Face Adjacency Network - keys from 0 to N!
# ==========================================================================

n = mesh.number_of_faces()

network = Network()

for fkey in mesh.faces():
	xyz = centroids[fkey]
	attr_dict = {k: v for k, v in zip('xyz', xyz)}
	network.add_vertex(key=fkey, attr_dict=attr_dict)


for fkey in mesh.faces():	
	nbrs = mesh.face_neighbors(fkey)
	for nbr in nbrs:
		if fkey == nbr:
			continue

		angle_diff = math.fabs(angles[fkey] - angles[nbr])

		try:
			ad = network.get_edge_attribute((fkey, nbr), 'angle_diff')
			if ad:
				continue
		except:
			network.add_edge(fkey, nbr, attr_dict={'angle_diff': angle_diff})

# # ==========================================================================
# # color up
# # ==========================================================================

anglemax = max(network.get_edges_attribute('angle_diff'))
print('angle diff max', anglemax)

colors = {}
for u, v, attr in network.edges(True):
	angle_diff = attr['angle_diff']
	color = i_to_rgb(angle_diff / anglemax)
	colors[(u, v)] = color

# # ==========================================================================
# # Set up Plotter
# # ==========================================================================

plotter = NetworkPlotter(network, figsize=(12, 9))
# plotter.draw_faces(facecolor=colors)
plotter.draw_vertices(radius=0.01)
plotter.draw_edges(color=colors)
plotter.show()
