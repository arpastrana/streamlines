'''
'''

import compas

# 513 seconds / no check proximity during growth / all seeds
#  seconds / no check proximity during growth / all seeds / no kd trees / no check proximity at start

def trace(particle, *args, **kwargs):
	'''
	'''
	particle.grow('forth')
	particle.grow('back')


class Particle(object):
	'''
	'''
	def __init__(self, xyz):
		'''
		'''
		self.xyz = xyz


	def grow_back(self):
		'''
		'''
		pass

	def grow_forth(self):
		'''
		'''
		pass


def grow(particle):
	'''
	'''
	# check (proximity, length, outside)
	# on which face the particle is?
	# interpolate grow vector
	# translate (runge-kutta 4)

	pass


# current MEBARKI implementation
# input: seeds

# methods:
# get_new_streamline() -> make_streamline() ->
# -> grow() -> intersect() ->
# -> arrange() -> create_polyline()

# get_new_streamline(): creates seed Node() and at start checks proximity to other nodes
# make_streamline(): encompasses grow(), arrange() and create_polyline(). Makes Streamline().

# grow():

# 1. check if target length
# 2. query growth direction
# 3. query max_dist
# 4. check if too close to others dist > max_dist
# 5. interpolate vector to follow
# 6. correct vector to be on plane
# 7. new Node() displace by adding vector to previous Node() position
# 8. with intersect(), check intersections with host Mesh() edges
# 9. check if outside of host Mesh()
# 10. check for iterations > max_iterations
# 11. add new Node() to Streamline()


# intersect()
# 1. intersects line [prev Node(), new Node()] with all current face edges
# 2. if there is an intersection:
# 	2a. check if it hit a boundary edge
# 	2b. return intersection point as new Node()
# 3. otherwise:
#   3a. check if it is still withing current face, else it is outside host Mesh()
# 	3b. return new Node()


# arrange(): orients list of Node() according to 'forth' / 'back' flags
# create_polyline(): Makes Polyline() if # of nodes == 2.