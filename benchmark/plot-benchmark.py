import csv
import numpy as np
from matplotlib import cycler
import matplotlib.pyplot as plt

import sys 

def plot_setup():
    plt.rcParams.update({'font.size': 10})
    colors = cycler('color',
                ['#EE6666', '#3388BB', '#9988DD',
                '#EECC55', '#88BB44', '#FFBBBB'])
    plt.rc('axes', facecolor='#E6E6E6', edgecolor='none',
        axisbelow=True, prop_cycle=colors)  
    plt.rc('grid', color='w', linestyle='solid')
    plt.rc('xtick', direction='out', color='black')
    plt.rc('ytick', direction='out', color='black')
    plt.rc('patch', edgecolor='#E6E6E6')
    plt.rc('lines', linewidth=2)

def find_versions(num_lines):
    versions = []
    for i in range(num_lines):
        line = fileAsList[i].split()
        if line[0] in versions :
            continue
        else:
            versions.append(line[0])
    return versions

def process_line(line):
    line = fileAsList[i].split()
    n = int(line[1])
    performance = float(line[4])
    return (n,performance)

def simplecount(filename):
    lines = 0
    for line in open(filename):
        lines += 1
    return lines

try:
    filename = sys.argv[-1]
    f = open(filename, 'r')
    fileAsList = f.readlines()

    num_lines = len(open(filename).readlines())
    versions = find_versions(num_lines)

    myDict = {} 

    for version in versions:
        myDict[version] = []

    for i in range(num_lines):
        line = fileAsList[i].split()
        version = line[0]
        n = int(line[1])
        performance = float(line[4])
        myDict[version].append((n,performance))

    plot_setup()

    for version in versions:
        points = myDict[version]
        x0 = [n for (n,p) in points]
        y0 = [p for (n,p) in points]
        plt.plot(x0,y0,'.-', linewidth=2,label=version)

    plt.gcf().text(0.125, 0.9, "Performance [Flops/Cycle] vs. input size", fontsize=10)
    plt.title('Performance-Plot',loc='left',y=1.06,fontsize=12)

    ax = plt.gca()
    ax.grid(which='major', axis='y')
    plt.ylim((0, 20))
    plt.yticks(fontsize=10)
    plt.axes().xaxis.set_label_coords(0.5, -0.15)
    plt.axes().tick_params(left= False)
    plt.xscale('log',basex=2)
    plt.xticks(fontsize=10, rotation=0)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True) # labels along the bottom edge are off
    plt.legend()
    plt.savefig(filename + ".pdf")
    plt.show()
except:
    print(filename + " can not be opened")
