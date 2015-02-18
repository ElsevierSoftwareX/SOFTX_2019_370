#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

py_data = []

with open('python_output_full.txt') as py_out:
    for line in py_out:
        data = line.split(',')
        data = [float(d) for d in data]
        py_data.append(data)

cpp_data = []

with open('output2.txt') as cpp_out:
    for line in cpp_out:
        data = line.split(',')
        data = [float(d) for d in data]
        cpp_data.append(data)

for rowi in xrange(len(py_data)):
    py_row = py_data[rowi]
    cpp_row = cpp_data[rowi]

    if py_row[3] == 0.0 and not cpp_row[3] == 0.0:
        print rowi, "PY ZERO <<<<<<<<<<<<<<<<<<<"
        print py_row
        print cpp_row
        print

    if not py_row[3] == 0.0 and cpp_row[3] == 0.0:
        print rowi, "CPP ZERO >>>>>>>>>>>>>>>>>>"
        print py_row
        print cpp_row
        print 

py_data = np.array(py_data)
cpp_data = np.array(cpp_data)

diff = py_data - cpp_data

avg_diff = np.mean(diff, axis=0)
cols = ['RT', 'MZ', 'Amp', 'Min Score', 'AB', 'A0', 'B0', 'r1']

print "     AVERAGE DIFFERENCES    "
print "----------------------------"
for avg, col in zip(avg_diff, cols):
    print "{0:<10}: {1:>12.5f}".format(col, avg)

print
print "            RMSE            "
print "----------------------------"
for coli in xrange(py_data.shape[1]):
    rmse = np.sqrt(np.mean((cpp_data[:, coli] - 
                                        py_data[:, coli]) ** 2))
    print "{0:<10}: {1:>16.10f}".format(cols[coli], rmse)

py_mz = py_data[:, 1]
py_ms = py_data[:, 3]
cpp_mz = cpp_data[:, 1]
cpp_ms = cpp_data[:, 3]

dark2_3 = [(27, 158, 119), (217, 95, 2), (117, 112, 179)]

for i in range(len(dark2_3)):
    r, g, b = dark2_3[i]
    dark2_3[i] = (r / 255., g / 255., b / 255.)
"""
plt.figure(figsize=(12, 9))
ax = plt.subplot(111)

ax.spines['top'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.spines['left'].set_visible(False)

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

plt.ylim(-0.2, 10)
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
plt.close()

#np.savetxt('compare_out.txt', diff, delimiter=',', fmt='%5.5f')

""" # HERE

"""

fig_amp = plt.figure()
ax = fig_amp.add_subplot(111, projection='3d')

x = py_data[:, 1]
y = cpp_data[:, 0]
z = py_data[:, 2]
#z[z<10000] = np.nan

ax.plot_trisurf(x, y, z, cmap=cm.YlOrRd, linewidths=0.00, 
                vmin=-300000, vmax=600000)

#ax.scatter(x, y, z, 
#           c=z,
#           cmap=cm.YlGnBu,
#           vmin=0,#-100000.0,
#           vmax=600000,
#           marker='o',
#           s=20,
#           linewidths=0.0,
#           alpha=1.0, 
#           label='Peaks'
#           )
 

ax.set_zlim3d(-10000, 600000)
ax.set_title("Example Twin Ion Peaks", size=22)
ax.set_ylabel("Retention Time", size=16)
ax.set_xlabel("m/z value", size=16)
ax.set_zlabel("Intensity (x1000)", size=16)
ax.view_init(azim=-56, elev=25)
zticks = ax.get_zticks().tolist()
zticks = [int(z) / 1000 for z in zticks]
ax.set_zticklabels(zticks[1:])

plt.savefig('twin_ion_peaks.png', bbox_inches='tight', 
            dpi=300)
plt.close()

"""
""" # HERE
fig_ms = plt.figure()
ax = fig_ms.add_subplot(111, projection='3d')

x = cpp_data[:, 0]
y = py_data[:, 1]
z = py_data[:, 3]

ax.scatter(x, y, z, c=dark2_3[1], marker='o')
"""
""" #HERE
fig_cpp = plt.figure()
ax = fig_cpp.add_subplot(111, projection='3d')

x = cpp_data[:, 1]
y = cpp_data[:, 0]
z = py_data[:, 3]
z[z==0] = np.nan
z2 = cpp_data[:, 3]
z2[z2==0] = np.nan

ax.scatter(x, y, z, 
           c=dark2_3[1],
           #cmap=cm.RdPu,
           #vmin=0.0,
           #vmax=100,
           marker='o',
           s=60,
           linewidths=0.0,
           alpha=0.20, 
           label='Python'
           )
           
ax.scatter(x, y, z2, 
           c=dark2_3[2],
           #cmap=cm.BuGn,
           marker='.',
           s=50,
           linewidths=0.0,
           alpha=0.9,
           label='C++')

ax.set_xlim3d(149.5, 154.5)
ax.set_zlim3d(-0.5, 10)

ax.set_title("Comparison of Correlation Scores", size=22)
ax.set_ylabel("Retention Time", size=16)
ax.set_xlabel("m/z value", size=16)
ax.set_zlabel("Correlation Score", size=16)
ax.view_init(azim=-56, elev=25)

plt.savefig('score_comparison.png', bbox_inches='tight', 
            dpi=300)

plt.close()


""" #HERE

fig_diff = plt.figure()
ax = fig_diff.add_subplot(111)#, projection='3d')


ax.spines['top'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.spines['left'].set_visible(False)

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tick_params(axis='both', which='both', bottom='on', top='off',
                labelbottom='on', left='on', right='off', labelleft='on')

xx = cpp_data[:, 1]
yy = cpp_data[:, 0]
zz = py_data[:, 3] - cpp_data[:, 3]
#z[z==0] = np.nan

#ax.plot_trisurf(x, y, z, cmap=cm.Spectral, linewidths=0.00, 
#                vmin=-2, vmax=2)

ax.axhline(0.0, linestyle='--', color='k', alpha=0.8)

ax.scatter(xx, zz,#y, z, 
           c=zz,
           cmap=cm.Spectral,
           #vmin=0.0,
           #vmax=100,
           marker='o',
           s=60,
           linewidths=0.0,
           alpha=0.8, 
           label='Difference'
           )
           
ax.set_xlim(149.5, 154.5)
ax.set_ylim(-2, 2)

ax.set_title("Difference in Correlation Scores", size=22)
#ax.set_ylabel("Retention Time", size=16)
ax.set_xlabel("m/z value", size=16)
ax.set_ylabel("Python - C++", size=16)
#ax.view_init(azim=-56, elev=25)

#plt.show()
plt.savefig('score_difference.png', bbox_inches='tight', 
            dpi=300)

plt.close()
