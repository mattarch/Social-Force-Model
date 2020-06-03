import csv
import numpy as np
from matplotlib import cycler
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
import sys 

import platform, subprocess

def get_processor_info():
    if platform.system() == "Windows":
        return ' AMD Ryzen 7 3700U with Radeon Vega Mobile Gfx'
    elif platform.system() == "Darwin":
        return "Intel(R) Core(TM) i7-9750H CPU @ 2.60GHz"
        # return " " + str(subprocess.check_output(['/usr/sbin/sysctl', "-n", "machdep.cpu.brand_string"]).strip().decode("utf-8"))
    elif platform.system() == "Linux":
        with open('/proc/cpuinfo') as f:
            for line in f:
                # Ignore the blank line separating the information between
                # details about two processing units
                if line.strip():
                    if line.rstrip('\n').startswith('model name'):
                        model_name = line.rstrip('\n').split(':')[1]
                        return model_name
    return "" 

def plot_setup():
    plt.rcParams.update({'font.size': 10})
    colors = cycler('color',
                ['#000000', '#2F4F4F', '#696969',
                '#A9A9A9', '#88BB44', '#FFBBBB',
                '#AA22FF', '#0492BA'])
    plt.rc('axes', facecolor='#E6E6E6', edgecolor='none',
        axisbelow=True, prop_cycle=colors)  
    plt.rc('grid', color='w', linestyle='solid',  linewidth=3)
    plt.rc('xtick', direction='out', color='black')
    plt.rc('ytick', direction='out', color='black')
    plt.rc('patch', edgecolor='#E6E6E6')
    plt.rc('lines', linewidth=10)

def find_versions(num_lines):
    versions = []
    for i in range(num_lines):
        line = fileAsList[i].split()
        if line[0] in versions :
            continue
        else:
            versions.append(line[0])
    return versions

def find_x_values_int(num_lines):
    values = []
    for i in range(num_lines):
        line = fileAsList[i].split()
        if int(line[1]) in values :
            continue
        else:
            values.append(int(line[1]))
    return values

def find_x_values(num_lines):
    values = []
    for i in range(num_lines):
        line = fileAsList[i].split()
        if line[1] in values :
            continue
        else:
            values.append(line[1])
    return values

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

filename = sys.argv[-1]
try:
    f = open(filename, 'r')
except:
    print(filename + " can not be opened")

markers = ['-o', '-^',  '-s', '-p', '-D', '-H', '-.', '-,']
marker_counter = -1
fileAsList = f.readlines()
num_lines = len(open(filename).readlines())
versions = find_versions(num_lines)
x_values_int = find_x_values_int(num_lines)
x_values_str = find_x_values(num_lines)

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
    marker_counter = marker_counter + 1
    points = myDict[version]
    x0 = [n for (n,p) in points]
    y0 = [p for (n,p) in points]
    plt.plot(x0,y0,markers[marker_counter % len(markers)], linewidth=2,label=version, markersize=4)

plt.gcf().text(0.125, 0.9, "Performance [Flops/Cycle] vs. input size", fontsize=14)
title = "Performance Plot on " + get_processor_info()
plt.title(title,loc='left',y=1.06,fontsize=16, weight='bold')

plt.text(512, 11.25, 'vectorize_5', weight='bold', fontsize=14)
plt.text(512, 1.25, 'stdc_2.5.1', weight='bold', fontsize=14, color='#2F4F4F')
plt.text(512, 0.60, 'basic version', weight='bold', fontsize=14, color='#696969')


ax = plt.gca()
ax.grid(which='major', axis='y')
plt.ylim((0, 12))
plt.yticks(fontsize=14)
plt.axes().xaxis.set_label_coords(0.5, -0.15)
plt.axes().tick_params(left= False)

# plt.yscale('log',basey=2)
plt.xscale('log',basex=2)
# plt.xticks(x_values_int, x_values_str, fontsize=10, rotation=0)
plt.xticks([32, 64, 128, 256, 512, 1024, 2048, 4096, 8192], ['32', '64', '128','256', '512', '1k', '2k', '4k', '8k'], fontsize=14, rotation=0)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True) # labels along the bottom edge are off
# plt.legend(frameon = False, loc = 'center right', fontsize = 14)
plt.savefig(filename + ".pdf")
plt.show()
