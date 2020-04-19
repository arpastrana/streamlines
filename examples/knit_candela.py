from time import time
from functools import partial

from compas.geometry import scale_vector
from compas.geometry import add_vectors
from compas.geometry import normalize_vector

from compas.datastructures import Mesh
from compas.datastructures import mesh_transformed
from compas.geometry import Scale
# from compas_fofin.datastructures import Cablenet
from compas_plotters import MeshPlotter

from compas_viewers.meshviewer import MeshViewer
from compas_viewers.multimeshviewer import MeshObject
from compas_viewers.multimeshviewer import MultiMeshViewer


HERE = '/Users/arpj/code/brg/ITA19/modules/module1/Week2/knitcandela/data/knitcandela.fofin'
THERE = '/Users/arpj/Desktop/knitcandela.obj'


class CustomMesh(Mesh):
	def __init__(self):
		super(CustomMesh, self).__init__()

	def vertices_attributes(self, query):
		return self.get_vertices_attributes(query)


mesh = Mesh.from_json(HERE)
print(mesh.vertex_coordinates(mesh.get_any_vertex()))

# mesh = Cablenet.from_json(HERE)
# mesh.to_obj(THERE)
plotter = MeshPlotter(mesh, figsize=(12, 9))
plotter.draw_edges()
plotter.draw_faces()
plotter.show()

# meshes =[MeshObject(mesh_transformed(m, Scale([0.001] * 3))) for m in [mesh]]

# visualization
# viewer = MultiMeshViewer()
# viewer.meshes = meshes
# viewer = MeshViewer()
# viewer.mesh = mesh
# viewer.show()
