import matplotlib.pyplot as plt
import csv
import numpy as np

from matplotlib import cycler
import sys 

def plot_setup():
    plt.rcParams.update({'font.size': 10})
    colors = cycler('color',
                ['#EE6666', '#3388BB', '#9988DD',
                '#EECC55', '#88BB44', '#FFBBBB'])
    plt.rc('axes', facecolor='#E6E6E6', edgecolor='none',
        axisbelow=True, prop_cycle=colors)  
    plt.rc('grid', color='w', linestyle='solid')
    plt.rc('xtick', direction='out', color='gray')
    plt.rc('ytick', direction='out', color='gray')
    plt.rc('patch', edgecolor='#E6E6E6')
    plt.rc('lines', linewidth=2)

    plt.xticks(np.arange(0,4201 , step=200))
    plt.xticks(fontsize=10, rotation=0)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True) # labels along the bottom edge are off

    plt.yticks(np.arange(0,3.1 , step=0.5))
    plt.yticks(fontsize=10)

    plt.axes().xaxis.set_label_coords(0.5, -0.15)
    plt.axes().yaxis.set_label_coords(0.04,1.02)
    plt.axes().tick_params(left= False)
    plt.xscale('log',basex=2)


    plt.gcf().subplots_adjust(bottom=0.15)
    ax = plt.gca()
    ax.grid(which='major', axis='y')
    ax.yaxis.set_label_coords(0.18,1.02)
    plt.ylabel('[flops/cycle] vs input size',rotation=0, fontsize=10)
    plt.title('Performance-Plot',loc='left',y=1.08)
    plt.show()
    plt.savefig(filename + ".pdf")


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



filename = sys.argv[-1] # file to open
with open(filename, 'r') as f:
    fileAsList = f.readlines()


num_lines = len(fileAsList)

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
    

for version in versions:
    points = myDict[version]
    x0 = [n for (n,p) in points]
    y0 = [p for (n,p) in points]
    plt.plot(x0,y0,'.-',color='r', linewidth=2,)
    #plt.text(200,0.35,'-O0',fontsize=10,color='r')

plot_setup()
