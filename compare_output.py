#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt

py_data = []

with open('data/python_single_out.txt') as py_out:
    for line in py_out:
        data = line.split(',')
        data = [float(d) for d in data]
        py_data.append(data)

cpp_data = []

with open('data/cpp_single_output.txt') as cpp_out:
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

py_mz = py_data[:, 1]
py_ms = py_data[:, 3]
cpp_mz = cpp_data[:, 1]
cpp_ms = cpp_data[:, 3]

dark2_3 = [(27, 158, 119), (217, 95, 2), (117, 112, 179)]

for i in range(len(dark2_3)):
    r, g, b = dark2_3[i]
    dark2_3[i] = (r / 255., g / 255., b / 255.)

plt.figure(figsize=(12, 9))
ax = plt.subplot(111)

ax.spines['top'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.spines['left'].set_visible(False)

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

plt.ylim(-0.2, 6)
plt.xlim(150, 154)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.xlabel('m/z value', fontsize=16)
plt.ylabel('Score', fontsize=16)
plt.title('Correlation scores for a single spectrum', fontsize=22)

plt.tick_params(axis='both', which='both', bottom='on', top='off',
                labelbottom='on', left='on', right='off', labelleft='on')

plt.plot(py_mz, py_ms, lw=2.5, color=dark2_3[0])
plt.plot(cpp_mz, cpp_ms, lw=2.5, color=dark2_3[2])
plt.savefig('comparison.png', bbox_inches='tight')

np.savetxt('compare_out.txt', percent, delimiter=',', fmt='%5.5f')
