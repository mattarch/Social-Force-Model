import csv
import numpy as np
from matplotlib import cycler
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl

import sys

import platform
import subprocess


################################################################
# PARAMETERS
################################################################
FREQUENCY = 2.6
DRAM_BANDWIDTH = 17.3 / FREQUENCY
L3_BANDWIDTH = 37.7 / FREQUENCY
L2_BANDWIDTH = 69.9 / FREQUENCY
L1_BANDWIDTH = 218.5 / FREQUENCY
MAX_FLOPS_PER_CYCLE = 4
END_PLOT = 100
START_X = 0.001
END_X = END_PLOT * 10
STEP_SIZE = 100000

dram_bandwidth_bound = [(START_X, START_X * DRAM_BANDWIDTH),(END_X, END_X * DRAM_BANDWIDTH)]
l3_bandwidth_bound = [(START_X, START_X * L3_BANDWIDTH),(END_X, END_X * L3_BANDWIDTH)]
l2_bandwidth_bound = [(START_X, START_X * L2_BANDWIDTH),(END_X, END_X * L2_BANDWIDTH)]
l1_bandwidth_bound = [(START_X, START_X * L1_BANDWIDTH),(END_X, END_X * L1_BANDWIDTH)]

flops_bound = [(START_X, MAX_FLOPS_PER_CYCLE),(END_X, MAX_FLOPS_PER_CYCLE)]
simd_double_bound = [(START_X, 4 * MAX_FLOPS_PER_CYCLE),(END_X, 4 * MAX_FLOPS_PER_CYCLE)]
simd_float_bound = [(START_X, 8 * MAX_FLOPS_PER_CYCLE),(END_X , 8 * MAX_FLOPS_PER_CYCLE)]

memory_bound_str = "DRAM: y = "+str(round(DRAM_BANDWIDTH,2))+" x"
l3_bound_str = "L3: y = "+str(round(L3_BANDWIDTH,2))+" x"
l2_bound_str = "L2: y = "+str(round(L2_BANDWIDTH,2))+" x"
l1_bound_str = "L1: y = "+str(round(L1_BANDWIDTH,2))+" x"

scalar_flops_bound_str = "scalar FMA peak: "+str(MAX_FLOPS_PER_CYCLE)
simd_double_flops_bound_str = "single precision SIMD ADD peak: "+str(MAX_FLOPS_PER_CYCLE * 4)
simd_float_flops_bound_str = "single precision SIMD FMA peak: "+str(MAX_FLOPS_PER_CYCLE * 8)

################################################################
################################################################

def get_processor_info():
    return "Intel Core i7-9750H CPU"

def plot_setup_style():
    mpl.rcParams.update({'font.size': 10})
    colors = cycler('color',
                ['#000000', '#2F4F4F', '#696969',
                '#A9A9A9', '#88BB44', '#FFBBBB',
                '#AA22FF', '#0492BA'])
    mpl.rc('axes', facecolor='#E6E6E6', edgecolor='none',
        axisbelow=True, prop_cycle=colors)  
    mpl.rc('grid', color='w', linestyle='solid',  linewidth=3)
    mpl.rc('xtick', direction='out', color='black')
    mpl.rc('ytick', direction='out', color='black')
    mpl.rc('patch', edgecolor='#E6E6E6')
    mpl.rc('lines', linewidth=10)


def find_versions(num_lines):
    versions = []
    for i in range(num_lines):
        line = fileAsList[i].split()
        if line[0] in versions:
            continue
        else:
            versions.append(line[0])
    return versions


def process_line(line):
    line = fileAsList[i].split()
    n = int(line[1])
    performance = float(line[4])
    return (n, performance)


def simplecount(filename):
    lines = 0
    for line in open(filename):
        lines += 1
    return lines

def plot_rooflines():

    x_value = 0.035

    angle_l1_bandwidth = 37
    angle_l2_bandwidth = 37
    angle_l3_bandwidth = 39
    angle_dram_bandwith = 39

    x_text_vector_fma_peak_line = 0.64
    x_text_vector_add_peak_line = 0.62
    x_text_scalar_add_peak_line = 2.2

    # plot DRAM BANDWITH LINE
    x0 = [n for (n, p) in dram_bandwidth_bound]
    y0 = [p for (n, p) in dram_bandwidth_bound]
    plt.plot(x0, y0, '--', linewidth=1.1, color='grey')
    plt.text(x_value, DRAM_BANDWIDTH * x_value + 0.03 , memory_bound_str, {'ha': 'left', 'va': 'bottom'}, rotation=angle_dram_bandwith, fontsize=10,bbox=dict(facecolor='#E6E6E6', alpha=1,pad=0.0))

    # plot L3 BANDWITH LINE
    x0 = [n for (n, p) in l3_bandwidth_bound]
    y0 = [p for (n, p) in l3_bandwidth_bound]
    plt.plot(x0, y0, '--', linewidth=1.1, color='grey')
    plt.text(x_value, L3_BANDWIDTH * x_value + 0.06, l3_bound_str, {'ha': 'left', 'va': 'bottom'}, rotation=angle_l3_bandwidth,fontsize=10,bbox=dict(facecolor='#E6E6E6', alpha=1,pad=0.0))

    # plot L2 BANDWITH LINE
    x0 = [n for (n, p) in l2_bandwidth_bound]
    y0 = [p for (n, p) in l2_bandwidth_bound]
    plt.plot(x0, y0, '--', linewidth=1.1, color='grey')
    plt.text(x_value, L2_BANDWIDTH * x_value + 0.15, l2_bound_str, {'ha': 'left', 'va': 'bottom'}, rotation=angle_l2_bandwidth,fontsize=10,bbox=dict(facecolor='#E6E6E6', alpha=1,pad=0.0))

    # plot L1 BANDWITH LINE
    x0 = [n for (n, p) in l1_bandwidth_bound]
    y0 = [p for (n, p) in l1_bandwidth_bound]
    plt.plot(x0, y0, '-', linewidth=1.1, color='black')
    plt.text(x_value, L1_BANDWIDTH * x_value + 0.4, l1_bound_str, {'ha': 'left', 'va': 'bottom'}, rotation=angle_l1_bandwidth,fontsize=10,bbox=dict(facecolor='#E6E6E6', alpha=1,pad=0.0))


    # plot scalar add peak line
    x0 = [n for (n, p) in flops_bound]
    y0 = [p for (n, p) in flops_bound]
    plt.plot(x0, y0, '--', linewidth=1.1, color='grey')
    plt.text(x_text_scalar_add_peak_line, MAX_FLOPS_PER_CYCLE + 0.375, scalar_flops_bound_str,bbox=dict(facecolor='#E6E6E6', alpha=1,pad=0.0))

    # plot vector add peak line
    x0 = [n for (n, p) in simd_double_bound]
    y0 = [p for (n, p) in simd_double_bound]
    plt.plot(x0, y0, '--', linewidth=1.1, color='grey')
    plt.text(x_text_vector_add_peak_line, (MAX_FLOPS_PER_CYCLE*4) + 1.25, simd_double_flops_bound_str,bbox=dict(facecolor='#E6E6E6', alpha=1,pad=0.0) )

    # plot vector fma peak line
    x0 = [n for (n, p) in simd_float_bound]
    y0 = [p for (n, p) in simd_float_bound]
    plt.plot(x0, y0, '-', linewidth=1.1,color='black')
    plt.text(x_text_vector_fma_peak_line, (MAX_FLOPS_PER_CYCLE*8) + 2.5, simd_float_flops_bound_str,bbox=dict(facecolor='#E6E6E6', alpha=1,pad=0.0))

def myticks(x,pos):

    if x == 0: return "$0$"

    if(x >= 1):
        return r"${:4.0f}$".format(x)
    else:
        return r"{0}".format(x)

########################################################################################################

filename = sys.argv[-1]
try:
    f = open(filename, 'r')
except:
    print(filename + " can not be opened")

fileAsList = f.readlines()
num_lines = len(open(filename).readlines())
versions = find_versions(num_lines)

performance_all = {}
runtime_all = {}
names = {}
colors = {}
all_markers = {}

names["simplified_float"] = "straightforward float"
names["simplified_double"] = "straightforward double"
names["stdc_optv_2_4_double"] = "stdc_2_4"
names["vectorize_2_5_1"] = "vectorize_2_5_1"
names["stdc_optv_2_5_1_double"] = "stdc_2_5_1"
names["vectorize_1"] = "vectorize_1"
names["vectorize_2"] = "vectorize_2"
names["vectorize_3"] = "vectorize_3"
names["vectorize_5"] = "vectorize_4"

colors["simplified_double"] = '#2F4F4F'
colors["simplified_float"] = '#2F4F4F'
colors["stdc_optv_2_4_double"] = '#000000'
colors["stdc_optv_2_5_1_double"] = '#696969'
colors["vectorize_2_5_1"] = '#696969'
colors["vectorize_1"] = '#708090'
colors["vectorize_2"] = '#0492BA'
colors["vectorize_3"] = '#778899'
colors["vectorize_5"] = '#000000'

all_markers["simplified_double"] = '-^'
all_markers["simplified_float"] = '-^'
all_markers["stdc_optv_2_4_double"] = '-H'
all_markers["stdc_optv_2_5_1_double"] = '-s'
all_markers["vectorize_2_5_1"] = '-s'
all_markers["vectorize_1"] = '-p'
all_markers["vectorize_2"] = '-D'
all_markers["vectorize_3"] = '-H'
all_markers["vectorize_5"] = '-o'

for version in versions:
    performance_all[version] = []

for i in range(num_lines):
    line = fileAsList[i].split()
    version = line[0]
    intensity = float(line[2])  # I(n)
    performance = float(line[1]) / FREQUENCY
    performance_all[version].append((intensity, performance))

# plot setup
plt.figure(figsize=(6, 4.5))
plot_setup_style()
plt.gcf().text(0.085, 0.895,
               "Performance [Flop/Cycle] vs. Operational Intensity [Flop/Bytes]", fontsize=12)
title = "Roofline Plot on " + get_processor_info()
plt.title(title, loc='left', y=1.06, fontsize=12,weight='bold')

plot_rooflines()

# plot points
for version in versions:
    points = performance_all[version]
    x0 = [n for (n, p) in points]
    y0 = [p for (n, p) in points]
    plt.plot(x0, y0, '.-',color=colors[version], markersize=14, label=version)

# add labels to points
plt.text(1.45, 27.28/FREQUENCY, names["vectorize_5"], weight='bold', fontsize=10,color=colors["vectorize_5"],horizontalalignment='left',verticalalignment='center')
plt.text(0.49, 0.945/FREQUENCY, names["simplified_float"], weight='bold', fontsize=10,color=colors["simplified_float"],horizontalalignment='left',verticalalignment='center')
plt.text(0.32, 1.409/FREQUENCY, names["simplified_double"], weight='bold', fontsize=10,color=colors["simplified_double"],horizontalalignment='left',verticalalignment='center')
plt.text(0.16, 2.618/FREQUENCY, names["stdc_optv_2_5_1_double"], weight='bold', fontsize=10,color=colors["stdc_optv_2_5_1_double"],horizontalalignment='left',verticalalignment='center',bbox=dict(facecolor='#E6E6E6', alpha=1,pad=0.0))
plt.text(1.45, 19.615/FREQUENCY, names["vectorize_2_5_1"], weight='bold', fontsize=10,color=colors["vectorize_2_5_1"],horizontalalignment='left',verticalalignment='center',bbox=dict(facecolor='#E6E6E6', alpha=1,pad=0.0))

# plot axes
ax = plt.gca()
ax.grid(which='major', axis='y')
plt.ylim((1/(2**2), 50))
plt.xlim((1/(2**5), 10))
plt.yticks(fontsize=12)
plt.axes().tick_params(left=False)
plt.yscale('log',basey=2)
plt.xscale('log',basex=10)
plt.yticks([0.5,1,2,4,8,16,32],fontsize=12)
plt.xticks([0.1,1,10],fontsize=12)
plt.axes().xaxis.set_major_formatter(ticker.FuncFormatter(myticks))
plt.axes().yaxis.set_major_formatter(ticker.FuncFormatter(myticks))

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True)  # labels along the bottom edge are off

plt.tight_layout()
plt.savefig(filename + ".eps")
plt.show()
