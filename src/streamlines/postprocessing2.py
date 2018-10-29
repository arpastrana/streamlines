import json
import math
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

filepath = '/Volumes/BOOTCAMP/Abaqus/parapluie_180914_01_ortho/'
name = 'parapluie_180914_01_ortho-results.json'


filepath_2 = '/Volumes/BOOTCAMP/Abaqus/parapluie_180914_01_psmid/'
name_2 = 'parapluie_180914_01_psmid-results.json'


with open(filepath+name) as f:
    database_1 = json.load(f)

with open(filepath_2+name_2) as f:
    database_2 = json.load(f)


def get_data(database, field, mode='max', ab=True, diameter=10, mpa=True):
    data = []
    area = math.pi * math.pow(diameter / 1000, 2) * 0.25

    for ekey, values in database['step_load']['element'][field].items():
        temp = []
        values = [value for ip, value in values.items()]

        if mode not in ['first', 'last']:
            for value in values:
                value = float(value)
                if ab is True:
                    value = math.fabs(value)

                if field == 'rbfor':
                    value = (value / area)
                temp.append(value)

                if mode == 'max':
                    out = max(temp)
                elif mode == 'min':
                    out = min(temp)
        else:
            if mode == 'first':
                out = values[0]
            elif mode == 'last':
                out = values[-1]

            if ab is True:
                out = math.fabs(out)

        if mpa is True:
            out = out * 10e-6
        data.append(out)

    print('number of elements is {}'.format(len(database['step_load']['element'][field].keys())))
    print('number of values processed is {}'.format(len(data)))
    print('max field value {} is {}'.format(field, max(data)))
    print('average field value {} is {}'.format(field, (sum(data)/len(data))))

    return sorted(data)


field = 'rbfor'

values_1 = get_data(database_1, field, 'max', True, 8, True)
print('---------')
values_2 = get_data(database_2, field, 'max', True, 8, True)


# make dataframe
# d = {'xy': values_1, 'psmid': values_2}
# d = pd.DataFrame(d)
#
#
# sns.set(style="whitegrid")
# ax = sns.boxplot(data=d)

# red dashes, blue squares and green triangles
x_val = range(len(values_1))

plt.plot(x_val, values_1, 'r--',
         x_val, values_2, 'b--',
         )

plt.ylabel('rebar stress - mPa')
plt.show()
