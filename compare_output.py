#!/bin/env python

## @package compare_output
# Program for comparing HiTIME output.
#
# Primarily designed for comparing output from the Python and C++ verions of
# HiTIME. The program identifies points that have a zero score in one output
# and not the other, calculates errors between the two outputs and produces
# a series of useful plots.

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

def main():

    # Read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('output1', help='path to first output file')
    parser.add_argument('output2', help='path to second output file')
    args = parser.parse_args()

    print "Loading data..."
    data1 = read_data(args.output1)
    data2 = read_data(args.output2)

    # Check output files have the same number of points, if not exit
    if not len(data1) == len(data2):
        print "Output files have different lengths and cannot be compared"
        print "Try producing full output"
        print "Ending..."
        sys.exit()

    print "Comparing zero score points..."
    compare_zeros(data1, data2)

    # Convert to numpy arrays for remaining processing
    data1 = np.array(data1)
    data2 = np.array(data2)

    print "Computing statistics..."
    calc_stats(data1, data2)

    # Brewer colours used for plotting
    dark2_3 = [(27, 158, 119), (217, 95, 2), (117, 112, 179)]

    # Convert from 0-255 to 0-1, required by matplotlib
    for i in range(len(dark2_3)):
        r, g, b = dark2_3[i]
        dark2_3[i] = (r / 255., g / 255., b / 255.)

    print "Creating plots..."
    plot_line_comp(data1, data2, dark2_3)
    plot_peaks(data1, data2)
    plot_score_comp(data2, data2, dark2_3)
    plot_score_diff(data1, data2)

    print "Done!"

## Read data from a HiTIME output file.
#
# @param datafile Path to a HiTIME output file.
#
# @return List of lists where each contains a line of data
def read_data(datafile):
    
    data_list = []

    with open(datafile) as df:
        for line in df:
            data = line.split(',')
            data = [float(d) for d in data]
            data_list.append(data)

    return data_list

## Indentify points that are zero in one data set and not in the other.
#
# @param data1 First list of lists containing HiTIME data
# @param data2 Second list of lists containing HiTIME data 
def compare_zeros(data1, data2):
   
    zero_count = 0

    for rowi in xrange(len(data1)):
        row1 = data1[rowi]
        row2 = data2[rowi]

        if row1[3] == 0.0 and not row2[3] == 0.0:
            print rowi, "FILE 1 ZERO <<<<<<<<<<<<<<<<<<<"
            print row1
            print row2
            print
            zero_count += 1

        if not row1[3] == 0.0 and row2[3] == 0.0:
            print rowi, "FILE 2 ZERO >>>>>>>>>>>>>>>>>>"
            print row1
            print row2
            print 
            zero_count += 1

    print "{0} Zero Mismatches Found".format(zero_count)
    print "These are points that have a zero score in one file and not the" \
          "other"

## Calculate difference and RMSE stats comparing two dataset. Results are 
# printed.
#
# @param data1 First numpy array containing HiTIME data
# @param data2 Second numpy array containing HiTIME data 
def calc_stats(data1, data2):

    # Find difference between values
    diff = data1 - data2

    # Find average differences
    avg_diff = np.mean(diff, axis=0)
    cols = ['RT', 'MZ', 'Amp', 'Min Score', 'AB', 'A0', 'B0', 'r1']

    print "     AVERAGE DIFFERENCES    "
    print "----------------------------"
    for avg, col in zip(avg_diff, cols):
        print "{0:<10}: {1:>16.10f}".format(col, avg)

    print
    print "            RMSE            "
    print "----------------------------"
    for coli in xrange(data1.shape[1]):
        rmse = np.sqrt(np.mean((data2[:, coli] - data1[:, coli]) ** 2))
        print "{0:<10}: {1:>16.10f}".format(cols[coli], rmse)


## Plot line chart showing difference between datasets. Data is compressed
# in the RT dimension. Result saved as 'line_comparision.png'.
#
# @param data1 First numpy array containing HiTIME data
# @param data2 Second numpy array containing HiTIME data
# @param colours List of (R, G, B) colours used for plotting 
def plot_line_comp(data1, data2, colours):

    # Extract data
    mz1 = data1[:, 1]
    ms1 = data1[:, 3]
    mz2 = data2[:, 1]
    ms2 = data2[:, 3]

    plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    #plt.ylim(-0.2, 10)
    #plt.xlim(150, 154)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.xlabel('m/z value', fontsize=16)
    plt.ylabel('Score', fontsize=16)
    plt.title('Correlation scores in m/z direction', fontsize=22)

    plt.tick_params(axis='both', which='both', bottom='on', top='off',
                    labelbottom='on', left='on', right='off', labelleft='on')

    plt.plot(mz1, ms1, lw=2.5, color=colours[0])
    plt.plot(mz2, ms2, lw=2.5, color=colours[2])
    
    plt.savefig('line_comparison.png', bbox_inches='tight')
    
    plt.close()


## Plot intensity peaks to demonstrate the original data. Results saved as
# 'intensity_peaks.png'.
#
# @param data1 First numpy array containing HiTIME data
# @param data2 Second numpy array containing HiTIME data 
def plot_peaks(data1, data2):

    fig_peaks = plt.figure()
    ax = fig_peaks.add_subplot(111, projection='3d')

    x = data1[:, 1]
    y = data2[:, 0]
    z = data1[:, 2]

    ax.plot_trisurf(x, y, z, cmap=cm.YlOrRd, linewidths=0.00, vmin=-300000, 
                    vmax=600000)

    #ax.set_zlim3d(-10000, 600000)
    ax.set_title("Intensity Peaks", size=22)
    ax.set_ylabel("Retention Time", size=16)
    ax.set_xlabel("m/z value", size=16)
    ax.set_zlabel("Intensity (x1000)", size=16)
    ax.view_init(azim=-56, elev=25)
    zticks = ax.get_zticks().tolist()
    zticks = [int(z) / 1000 for z in zticks]
    ax.set_zticklabels(zticks)

    plt.savefig('intensity_peaks.png', bbox_inches='tight', dpi=300)
    
    plt.close()


## Plot 3D scatter plot to compare datasets. Result saved as 
# 'score_comparision.png'.
#
# @param data1 First numpy array containing HiTIME data
# @param data2 Second numpy array containing HiTIME data
# @param colours List of (R, G, B) colours used for plotting 
def plot_score_comp(data1, data2, colours):

    fig_comp = plt.figure()
    ax = fig_comp.add_subplot(111, projection='3d')

    x = data2[:, 1]
    y = data2[:, 0]
    z = np.copy(data1[:, 3])
    z2 = np.copy(data2[:, 3])
    z[z==0] = np.nan
    z2[z2==0] = np.nan

    ax.scatter(x, y, z, 
               c=colours[1],
               marker='o',
               s=60,
               linewidths=0.0,
               alpha=0.20, 
              )
           
    ax.scatter(x, y, z2, 
               c=colours[2],
               marker='.',
               s=50,
               linewidths=0.0,
               alpha=0.9,
              )

    #ax.set_xlim3d(149.5, 154.5)
    #ax.set_zlim3d(-0.5, 10)

    ax.set_title("Comparison of Correlation Scores", size=22)
    ax.set_ylabel("Retention Time", size=16)
    ax.set_xlabel("m/z value", size=16)
    ax.set_zlabel("Correlation Score", size=16)
    ax.view_init(azim=-56, elev=25)

    plt.savefig('score_comparison.png', bbox_inches='tight', dpi=300)

    plt.close()


## Plot scatter plot showing difference between datasets. Data is compressed
# in the RT dimension and each dot show the difference between the datasets
# at that point. Result saved as 'score_difference.png'.
#
# @param data1 First numpy array containing HiTIME data
# @param data2 Second numpy array containing HiTIME data 
def plot_score_diff(data1, data2):

    fig_diff = plt.figure()
    ax = fig_diff.add_subplot(111)

    xx = data1[:, 1]
    yy = data2[:, 0]
    z1 = data1[:, 3] 
    z2 = data2[:, 3]
    zz = z1 - z2

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.tick_params(axis='both', which='both', bottom='on', top='off',
                    labelbottom='on', left='on', right='off', labelleft='on')


    ax.axhline(0.0, linestyle='--', color='k', alpha=0.8)

    ax.scatter(xx, zz, 
               c=zz,
               cmap=cm.Spectral,
               marker='o',
               s=60,
               linewidths=0.0,
               alpha=0.8, 
               label='Difference'
              )
           
    #ax.set_xlim(149.5, 154.5)
    #ax.set_ylim(-2, 2)

    ax.set_title("Difference in Correlation Scores", size=22)
    ax.set_xlabel("m/z value", size=16)
    ax.set_ylabel("Difference", size=16)

    plt.savefig('score_difference.png', bbox_inches='tight', dpi=300)

    plt.close()

if __name__ == "__main__":
    main()
