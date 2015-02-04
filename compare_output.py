#!/bin/env python

import numpy as np

py_data = []

with open('single_out.txt') as py_out:
    for line in py_out:
        data = line.split(',')
        data = [float(d) for d in data]
        py_data.append(data)

cpp_data = []

with open('cpp_output.txt') as cpp_out:
    for line in cpp_out:
        data = line.split(',')
        data = [float(d) for d in data]
        cpp_data.append(data)

py_data = np.array(py_data)
cpp_data = np.array(cpp_data)
diff = py_data - cpp_data
percent = (diff / py_data) * 100

avg_pc = np.mean(percent, axis=0)
cols = ['RT', 'MZ', 'Amp', 'Min Score', 'AB', 'A0', 'B0', 'r1']

print "    AVERAGE DIFFERENCES    "
print "---------------------------"
for avg, col in zip(avg_pc, cols):
    print "{0:<10}: {1:>12.5f} %".format(col, avg)

np.savetxt('compare_out.txt', percent, delimiter=',', fmt='%5.5f')
