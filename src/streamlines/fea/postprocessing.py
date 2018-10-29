import json
import math
import numpy as np
import matplotlib.pyplot as plt

filepath = '/Volumes/[C] My Boot Camp/Abaqus/slab_270818_01/'
name = 'slab_270818_01-results.json'

filepath_2 = '/Volumes/[C] My Boot Camp/Abaqus/slab_270818_01_ortho/'
name_2 = 'slab_270818_01_ortho-results.json'


with open(filepath+name) as f:
    data = json.load(f)

with open(filepath_2+name_2) as f:
    data_2 = json.load(f)


print('number of elements is {}'.format(len(data['step_load']['element']['rbfor'].keys())))

rb_forces = []
rb_forces_min = []
diameter = 8

area = math.pi * math.pow(diameter / 1000, 2) * 0.25
print('area is {}'.format(area))
for ekey, values in data['step_load']['element']['rbfor'].items():
    temp = []
    for ip, value in values.items():
        temp.append((math.fabs(float(value)) / area) * 10e-6)
    rb_forces.append(max(temp))
    rb_forces_min.append(min(temp))

print('number of rebar forces are {}'.format(len(rb_forces)))
print('max rebar force is {}'.format(max(rb_forces)))
print('average rebar force is {}'.format((sum(rb_forces)/len(rb_forces))))

array_rb_forces = np.array(sorted(rb_forces))
array_rb_forces_min = np.array(sorted(rb_forces_min))

print('---------')
print('number of elements is {}'.format(len(data['step_load']['element']['rbfor'].keys())))


rb_forces_2= []
rb_forces_min_2 = []
diameter_2 = 10
area_2 = math.pi * math.pow(diameter_2 / 1000, 2) * 0.25
print('area 2 is {}'.format(area_2))
for ekey, values in data_2['step_load']['element']['rbfor'].items():
    temp = []
    for ip, value in values.items():
        temp.append((math.fabs(float(value)) / area_2) * 10e-6)
    rb_forces_2.append(max(temp))
    rb_forces_min_2.append(min(temp))

print('number of rebar forces 2 are {}'.format(len(rb_forces_2)))
print('max rebar force 2 is {}'.format(max(rb_forces_2)))
print('average rebar 2 force is {}'.format((sum(rb_forces_2)/len(rb_forces_2))))

array_rb_forces_2 = np.array(sorted(rb_forces_2))
array_rb_forces_min_2 = np.array(sorted(rb_forces_min_2))

#plt.plot(sorted(rb_forces))

# red dashes, blue squares and green triangles
x_val = range(len(rb_forces))
plt.plot(x_val, array_rb_forces, 'r--',
         x_val, array_rb_forces_2, 'b--',
         )

# plt.plot(x_val, sorted(rb_forces), 'r--', x_val, sorted(rb_forces_min), 'b--')
plt.ylabel('rebar stress - MPa')
plt.show()
