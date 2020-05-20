import csv
import numpy as np
from matplotlib import cycler
import matplotlib.pyplot as plt

import sys 

import platform, subprocess


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

filename = sys.argv[-1]
try:
    f = open(filename, 'r')
except:
    print(filename + " can not be opened")

fileAsList = f.readlines()
num_lines = len(open(filename).readlines())

myDict = {} 
versions = []

for i in range(num_lines):
    line = fileAsList[i].split()
    version = str(line[0])
    if version not in versions:
      versions.append(version)
      myDict[version] = []
    n = float(line[1])
    diff = float(line[4])
    myDict[version].append((n,diff))

plot_setup()

for version in versions:
    points = myDict[version]
    x0 = [n for (n,p) in points]
    y0 = [p for (n,p) in points]
    plt.plot(x0,y0,'.-', linewidth=2,label=version)

plt.gcf().text(0.125, 0.9, "Exp fast error based on number of multiplications", fontsize=10)

ax = plt.gca()
ax.grid(which='major', axis='y')
plt.yticks(fontsize=10)
plt.axes().xaxis.set_label_coords(0.5, -0.15)
plt.axes().tick_params(left= False)
plt.xticks(fontsize=10, rotation=0)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True) # labels along the bottom edge are off
#plt.legend()
plt.savefig(filename + ".pdf")
plt.show()
