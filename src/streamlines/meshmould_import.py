'''
An text interface to communicate with the Mesh Mould technology.
'''


__name__ = "Mesh Mould Interface"
__author__ = "Rafael Pastrana"
__version__ = "0.0.1"
__creation__ = "2018.09.16"
__date__ = "2018.09.16"



path = '/Users/arpj/Dropbox/Documentos/ETH Zuerich/ETH MAS Thesis/meshmould/'
name = 'new_stool_15deg_0'

filepath = path + name + '.txt'

with open(filepath, 'r') as f:
    data = f.readlines()

data = [x[1:-2] for x in data]
data = [x.split(', ') for x in data]
data = [x for x in data if len(x) > 2]

points = []

for line in data:
    polygon = []
    index = 0
    lim = 3

    for num in range(4):
        polygon.append(list(map(lambda x: float(x), line[index:index + lim])))
        index += lim

    points.append(polygon)
