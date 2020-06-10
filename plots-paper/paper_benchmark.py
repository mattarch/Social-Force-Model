import csv
import numpy as np
from matplotlib import cycler
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys 

import platform, subprocess

def get_processor_info():
    return "Intel Core i7-9750H CPU"
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

#start of code

filename = sys.argv[1]
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

# parse versions

performance_all = {}
runtime_all = {}
names = {}
colors = {}
all_markers = {}

names["simplified_float"] = "straightforward"
names["simplified_double"] = "straightforward"
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
    runtime_all[version] = []

for i in range(num_lines):
    line = fileAsList[i].split()
    version = line[0]
    n = int(line[1])
    performance = float(line[4])
    runtime = float(line[3])
    performance_all[version].append((n,performance))
    runtime_all[version].append((n,runtime))

def plot_title_performance_12():
    plt.gcf().text(0.076, 0.895, "Performance [Flop/Cycle] vs. input size", fontsize=12)
    title = "Performance Plot on " + get_processor_info()
    plt.title(title,loc='left',y=1.05,fontsize=12, weight='bold')

def plot_title_performance_4():
    plt.gcf().text(0.06, 0.895, "Performance [Flop/Cycle] vs. input size", fontsize=12)
    title = "Performance Plot on " + get_processor_info()
    plt.title(title,loc='left',y=1.05,fontsize=12, weight='bold')

def plot_title_speedup():
    plt.gcf().text(0.076, 0.895, "Speedup w.r.t. straightforward version vs. input size", fontsize=12)
    title = "Speedup Plot on " + get_processor_info()
    plt.title(title,loc='left',y=1.05,fontsize=12, weight='bold')

def plot_lines(data):
    for version in versions:
        points = data[version]
        x0 = [n for (n,p) in points]
        y0 = [p for (n,p) in points]
        plt.plot(x0,y0,all_markers[version], linewidth=2,label=version, markersize=4,color=colors[version])


if(filename == "paper_max_performance.txt"):
    plt.figure(figsize=(6, 4.5))
    plot_setup()
    plot_title_performance_12()
    plot_lines(performance_all)

    plt.text(3000, 10.4, names["vectorize_5"], weight='bold', fontsize=10,color=colors["vectorize_5"])
    plt.text(3000, 1.3, names["stdc_optv_2_5_1_double"], weight='bold', fontsize=10, color=colors["stdc_optv_2_5_1_double"])
    plt.text(3000, 0.60, names["simplified_float"], weight='bold', fontsize=10, color=colors["simplified_float"])

    plt.gca().grid(which='major', axis='y')
    plt.ylim((0, 12))
    plt.yticks(fontsize=12)
    plt.axes().xaxis.set_label_coords(0.5, -0.15)
    plt.axes().tick_params(left= False)
    #plt.yscale('log',basey=2)
    #plt.xscale('log',basex=2)
    plt.xticks([0,1000,2000,3000,4000,5000,6000], fontsize=12, rotation=0)

if(filename == "paper_std.txt"):
    plt.figure(figsize=(6, 4.5))
    plot_setup()
    plot_title_performance_4()
    plot_lines(performance_all)

    plt.text(3000, 1.2, names["stdc_optv_2_5_1_double"], weight='bold', fontsize=10,color=colors["stdc_optv_2_5_1_double"])
    plt.text(3000, 0.75, names["stdc_optv_2_4_double"], weight='bold', fontsize=10, color=colors["stdc_optv_2_4_double"])
    plt.text(3000, 0.38, names["simplified_double"], weight='bold', fontsize=10, color=colors["simplified_double"])

    plt.gca().grid(which='major', axis='y')
    plt.ylim((0, 4))
    plt.yticks(fontsize=12)
    plt.axes().xaxis.set_label_coords(0.5, -0.15)
    plt.axes().tick_params(left= False)
    #plt.yscale('log',basey=2)
    #plt.xscale('log',basex=2)
    plt.yticks([0,1,2,3,4], fontsize=12, rotation=0)

    plt.xticks([0,1000,2000,3000,4000,5000,6000], fontsize=12, rotation=0)

if(filename == "paper_simd.txt"):
    plt.figure(figsize=(6, 4.5))
    plot_setup()
    plot_title_performance_12()
    plot_lines(performance_all)

    plt.text(9000, 10.2, names["vectorize_5"], weight='bold', fontsize=10,color=colors["vectorize_5"])
    plt.text(9000, 8.9, names["vectorize_3"], weight='bold', fontsize=10,color=colors["vectorize_3"])
    plt.text(9000, 6.9, names["vectorize_2"], weight='bold', fontsize=10,color=colors["vectorize_2"])
    plt.text(9000, 6.20, names["vectorize_1"], weight='bold', fontsize=10,color=colors["vectorize_1"])

    plt.text(9000, 7.8, names["vectorize_2_5_1"], weight='bold', fontsize=10, color=colors["vectorize_2_5_1"])
    plt.text(9000, 0.35, names["simplified_float"], weight='bold', fontsize=10, color=colors["simplified_float"])
 

    plt.gca().grid(which='major', axis='y')
    plt.ylim((0, 12))
    plt.yticks(fontsize=12)
    plt.axes().xaxis.set_label_coords(0.5, -0.1)

    plt.axes().tick_params(left= False)
    #plt.yscale('log',basey=2)
    plt.xticks(fontsize=12)
    plt.xscale('log',basex=2)
    # plt.xlabel("n",fontsize=12)

    # def format_func(value, tick_number):
    #     if (value%2) == 0:
    #         return r"${0}$".format(2*tick_number+1)
    #     else:
    #         return r"$2^{tick_number}$"

    # plt.axes().xaxis.set_major_formatter(plt.FuncFormatter(format_func))

if(filename == "paper_speedup.txt"):
    plt.figure(figsize=(6, 4.5))
    plot_setup()
    plot_title_speedup()

    # normalize lines
    input_sizes_double,stf_double = zip(*runtime_all["simplified_double"])
    input_sizes_double = list(input_sizes_double)
    stf_double = list(stf_double)

    input_sizes_float,stf_float = zip(*runtime_all["simplified_float"])
    input_sizes_float = list(input_sizes_float)
    stf_float = list(stf_float)
    versions = [v for v in versions if v != "simplified_double"]
    for version in versions:
        num,version_values = zip(*runtime_all[version])
        version_values = list(version_values)
        num = list(num)
        for i in range(len(runtime_all["simplified_float"])):
            if(version != "stdc_optv_2_5_1_double"):
                version_values[i] = stf_float[i]/version_values[i]
            else:
                version_values[i] = stf_double[i]/version_values[i]
        runtime_all[version] = list(zip(num,version_values))

    plot_lines(runtime_all)


    plt.text(3000, 18.3, names["vectorize_5"], weight='bold', fontsize=10,color=colors["vectorize_5"])
    plt.text(3000, 2.5, names["stdc_optv_2_5_1_double"], weight='bold', fontsize=10, color=colors["stdc_optv_2_5_1_double"])
    plt.text(3000, 1.15, names["simplified_float"], weight='bold', fontsize=10, color=colors["simplified_float"])

    plt.gca().grid(which='major', axis='y')
    plt.ylim((0, 25))
    plt.yticks(fontsize=12)
    plt.axes().xaxis.set_label_coords(0.5, -0.1)

    plt.axes().tick_params(left= False)
    #plt.yscale('log',basey=10)
    plt.xticks(fontsize=12)
    #plt.yscale('log',basey=10)

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True) # labels along the bottom edge are off

#plt.legend(frameon = False, loc = 'center right', fontsize = 14)
plt.subplots_adjust(top=0.85)
plt.tight_layout()

plt.savefig(filename + ".eps")
plt.show()
