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

from streamlines.utilities import Utilities


ut = Utilities()

# ==========================================================================
# Globals
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

vector_tag_1 = 'ps_1_top'  # ps_1_top
vector_tag_2 = 'ps_2_top'  # ps_1_top
vector_tag = 'ps_12_top'  # ps_1_top
smooth_iters = 0

# ==========================================================================
# Import mesh
# ==========================================================================

mesh = Mesh()
mesh.load(HERE)
mesh_unify_cycles(mesh)

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
	# angle = angle_vectors([1.0, 0.0, 0.0], vector, deg=True)
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

anglemax = max(angles.values())
colors = {}
for idx, angle in angles.items():
	color = i_to_rgb(angle / anglemax)
	colors[idx] = color

# ==========================================================================
# Create Face Adjacency Network - keys from 0 to N!
# ==========================================================================

n = mesh.number_of_faces()

A = np.zeros((n, n))
for fkey in mesh.faces():
	nbrs = mesh.face_neighbors(fkey)
	for nbr in nbrs:
		if fkey == nbr:
			continue

		angle_diff = math.fabs(angles[fkey] - angles[nbr])
		A[fkey, nbr] = angle_diff

# ==========================================================================
# Peek in
# ==========================================================================

# peek into row values
# row = 10
# B = np.nonzero(A)
# rows, cols = B
# sel_rows = rows[rows == row]
# sel_cols = cols[rows == row]
# sel = A[sel_rows, sel_cols]
# print(sel_rows, sel_cols)
# print(sel)


# find max index
# a = np.argmax(A, axis=0)
# b = np.argmax(A, axis=1)
# c = np.argmax(A)

# print(A)
# print()
# print(A[:, b])
# print()
# print(A[a, :])
# print()
# print(c)
# print(A.flatten()[c])

# print('max', np.amax(A))

# ==========================================================================
# Normalization and reversal
# ==========================================================================

# normalization 0-1
A /= np.amax(A)

# reverse matrix so that large differences mean smaller weights
A = 1.0 - A

# ==========================================================================
# Angles as heights
# ==========================================================================

X = np.zeros((n, 3))
for fkey, centroid in centroids.items():
	X[fkey, :] = centroids[fkey]
	X[fkey, 2] = angles[fkey]

# ==========================================================================
# Heat up
# ==========================================================================

# print('heating up...')
# sigma = 1.0
# A = np.exp(-1.0 * (np.power(A, 2.0) / (2 * np.power(sigma, 2.0))))
# print('heated up')

# ==========================================================================
# Scipy Spectral Clustering - doesn't work
# ==========================================================================

n_clusters = 7

print('spectral clustering...')
clustering = SpectralClustering(n_clusters=n_clusters, affinity="nearest_neighbors")
print('fitting...')
clustering.fit(X)
print('fit')

labels = clustering.labels_

print('coloring...')
colors = {idx: i_to_rgb(label / (n_clusters - 1)) for idx, label in enumerate(labels)}
print('colored')

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# # ax.scatter(X[:,0], X[:,1], X[:,2], c=labels)
# ax.scatter(X[:,0], X[:,1], X[:,2])
# plt.show()

# ==========================================================================
# Spectral Clustering Manual Mode - works okayish but with a non-sparse matrix
# ==========================================================================

# n_clusters = 9

# # find laplacian
# D = np.diag(A.sum(axis=1))
# L = D - A
# L = csr_matrix(L)

# # normalize_laplacian
# # print('normalizing laplacian...')
# # print('normalized')
# # L = laplacian(A, normed=True)

# # find the eigenvalues and eigenvectors
# print('eigen valueing...')
# vals, vecs = eigsh(L, k=n_clusters + 1, which='LM')
# # vals, vecs = np.linalg.eigh(L)
# print('eigen valued')

# # sort
# print('sorting....')
# vecs = vecs[:,np.argsort(vals)]
# vals = vals[np.argsort(vals)]
# print('sorted')

# print('evecs', vecs.shape)
# print('evals', vals)

# # plot sorted eigenvalues
# plt.plot(vals)
# plt.show()

# # do kmeans
# print('kmeaning...')

# # 4-E) Create a low-dimensional representation of the data. 
# # Take the train_images data and multiply it by the top eigenvectors
# proj = np.matmul(L, vecs[:, :2])
# print('projection matrix shape', proj.shape)
# plt.scatter(x=proj[:, 0], y=proj[:, 1])
# plt.show()

# # cluster
# clustering = KMeans(n_clusters=n_clusters)
# clustering.fit(vecs[:, 1:n_clusters+1])
# # clustering.fit(vecs[:, :n_clusters])
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
# Set up Plotter
# ==========================================================================

plotter = MeshPlotter(mesh, figsize=(12, 9))
plotter.draw_faces(facecolor=colors)
plotter.show()
