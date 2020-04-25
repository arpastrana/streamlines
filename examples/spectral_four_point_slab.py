import math
import compas
import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial import distance
from scipy.sparse.linalg import eigs
from scipy.sparse.linalg import eigsh

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
from compas.geometry import length_vector
from compas.geometry import angle_vectors
from compas.geometry import distance_point_point

from compas.datastructures import Mesh
from compas.datastructures import Network
from compas.datastructures import mesh_unify_cycles
from compas.datastructures import mesh_dual

from compas.utilities import i_to_rgb

from compas_plotters import MeshPlotter


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

# ==========================================================================
# Import mesh
# ==========================================================================

mesh = Mesh()
mesh.load(HERE)
mesh_unify_cycles(mesh)

# ==========================================================================
# Create PS vector lines
# ==========================================================================

vector_tag = 'ps_1_top'  # ps_1_top

angles = {}
centroids = {}
for fkey, attr in mesh.faces(data=True):
	vector = attr.get(vector_tag)
	angle = angle_vectors([1.0, 0.0, 0.0], vector, deg=True)
	angles[fkey] = angle
	centroids[fkey] = np.array(mesh.face_centroid(fkey))


print('max angle', min(angles.values()))
print('min angle', max(angles.values()))

for idx, angle in angles.items():
	if angle <= 90.0:
		continue
	angles[idx] = 180.0 - angle

print('max angle', min(angles.values()))
print('min angle', max(angles.values()))

anglemax = max(angles.values())
colors = {}
for idx, angle in angles.items():
	color = i_to_rgb(angle / anglemax)
	colors[idx] = color

vectors = {}
for fkey, angle in angles.items():
	y = 1.0 / math.tan(math.radians(angle))
	x_vec = [1.0, 0.0, 0.0]
	y_vec = [0.0, y, 0.0]
	vec = normalize_vector(add_vectors(x_vec, y_vec))
	vectors[fkey] = vec

# ==========================================================================
# Create Face Adjacency Network - keys from 0 to N!
# ==========================================================================

n = mesh.number_of_faces()

# A = np.zeros((n, n))
# print('Shape of Adjacency Matrix', A.shape)

# for fkey in mesh.faces():
# 	nbrs = mesh.face_neighbors(fkey)

# 	for nbr in nbrs:
# 		# A[fkey, nbr] = 1.0
# 		angle_diff = math.fabs(angles[fkey] - angles[nbr])
# 		A[fkey, nbr] = angle_diff

# # all options
# A_all = np.zeros((n, n))

# ANGLES = sorted(list(angles.items()), key=lambda x: x[0])
# ANGLES = np.array([x[1] for x in ANGLES])
# B = ANGLES.reshape((1, -1)) - ANGLES.reshape((-1, 1))


# CENT = sorted(list(centroids.items()), key=lambda x: x[0])
# CENT = np.array([x[1] for x in CENT])

# C = euclidean_distances(CENT, CENT)
# C /= np.amax(C)

# B /= np.amax(B)

# A_all = B - np.exp(1 / C)

# for fkey in mesh.faces():
# 	for okey in mesh.faces():
# 		if fkey == okey:
# 			A_all[fkey, okey] = 0.0
# 		else:
# 			angle_diff = math.fabs(angles[fkey] - angles[okey])

# 			pta = centroids[fkey]
# 			ptb = centroids[okey]

# 			# dist = distance_point_point(pta, ptb)

# 			# dist = np.linalg.norm(pta - ptb)
# 			dist = distance.euclidean(pta, ptb)
			
# 			A_all[fkey, okey] = angle_diff + dist


# X = np.zeros((n, 3))
# for fkey, vector in vectors.items():
# 	for i in range(len(vector)):
# 		if i != 2:
# 			X[fkey, i] = vector[i]
# 		else:
# 			X[fkey, i] = angles[i] ** 2

X = np.zeros((n, 3))
sigma = 1.0
for fkey, centroid in centroids.items():
	X[fkey,:] = centroids[fkey]
	X[fkey,2] = angles[fkey]

# A_dist = euclidean_distances(X, X)

# AN = np.zeros((n, 1))
# for fkey, angle in angles.items():
# 	AN[fkey] = angle

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.scatter(X[:,0], X[:,1], X[:,2])
# plt.show()

# ==========================================================================
# Heat up
# ==========================================================================

# # switch matrix from sparse to all
# A = A_dist

# # do gaussian heat
# print('heating up...')
# sigma = 1.0
# A = np.exp(-1.0 * (np.power(A, 2.0) / (2 * np.power(sigma, 2.0))))
# print('heated up')


# ==========================================================================
# Scipy KMeans - works okayish
# ==========================================================================

# n_clusters = 7

# clustering = KMeans(n_clusters=n_clusters, random_state=0)
# clustering.fit(X)

# centers = clustering.cluster_centers_
# labels = clustering.labels_

# print('coloring...')
# colors = {idx: i_to_rgb(label / (n_clusters - 1)) for idx, label in enumerate(labels)}
# print('colored')


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(X[:,0], X[:,1], X[:,2], c=labels)
# plt.show()

# ==========================================================================
# Scipy Spectral Clustering - doesn't work
# ==========================================================================

n_clusters = 7
print('spectral clustering...')  # has worked best so fa r(roundish clusters)
clustering = SpectralClustering(n_clusters=n_clusters, affinity="nearest_neighbors", assign_labels='kmeans')
print('fitting...')
clustering.fit(X)
print('fit')

# print('spectral clustering...')
# A = kneighbors_graph(X, n_neighbors=5).toarray()
# clustering = SpectralClustering(n_clusters=n_clusters, affinity="precomputed")
# print('fitting...')
# clustering.fit(A)
# print('fit')

# print('AgglomerativeClustering...')
# A = kneighbors_graph(X, n_neighbors=n_clusters, include_self=False)
# clustering = AgglomerativeClustering(n_clusters=n_clusters, connectivity=A)
# print('fitting...')
# clustering.fit(X)
# print('fit')

labels = clustering.labels_
print(labels)

print('coloring...')
colors = {idx: i_to_rgb(label / (n_clusters - 1)) for idx, label in enumerate(labels)}
print('colored')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X[:,0], X[:,1], X[:,2], c=labels)
plt.show()

# ==========================================================================
# Spectral Clustering Manual Mode - works okayish but with a non-sparse matrix
# ==========================================================================

# # find laplacian
# D = np.diag(A.sum(axis=1))
# L = D - A

# # find the eigenvalues and eigenvectors
# print('eigen valueing...')
# vals, vecs = np.linalg.eigh(L)
# # vals, vecs = eigsh(L, k=n_clusters, which='SM')
# print('eigen valued')

# # sort
# print('sorting....')
# vecs = vecs[:,np.argsort(vals)]
# vals = vals[np.argsort(vals)]
# print('sorted')

# # plot sorted eigenvalues
# plt.plot(vals)
# plt.show()

# # do kmeans
# print('kmeaning...')
# n_clusters = 7

# clustering = KMeans(n_clusters=n_clusters)
# # clustering.fit(vecs[:, 1:n_clusters+1])
# clustering.fit(vecs[:, :n_clusters])
# labels = clustering.labels_
# print('kmeant')

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(X[:,0], X[:,1], X[:,2], c=labels)
# plt.show()

# # color up
# print('coloring...')
# colors = {idx: i_to_rgb(label / (n_clusters - 1)) for idx, label in enumerate(labels)}
# print('colored')

# ==========================================================================
# Create Face Adjacency Network
# ==========================================================================

# network = Network()
# network.add_vertex(0)
# network.add_vertex(1)
# network.add_edge(0, 1)
# print(network)

# from compas.datastructures import mesh_connectivity_matrix
# from compas.datastructures import mesh_adjacency_matrix
# from compas.topology import face_adjacency

# mesh = Mesh.from_obj(compas.get('faces.obj'))
# print(mesh.number_of_faces())
# print(mesh.number_of_vertices())
# print(mesh.number_of_edges())

# mcm = mesh_connectivity_matrix(mesh)
# print(mcm.shape)

# mac = mesh_adjacency_matrix(mesh)  # this one?
# print(mac.shape)


# xyz, faces = mesh.to_vertices_and_faces()
# fa = face_adjacency(xyz, faces)  # ok!

# ==========================================================================
# Set up Plotter
# ==========================================================================

plotter = MeshPlotter(mesh, figsize=(12, 9))
plotter.draw_faces(facecolor=colors)
plotter.show()
