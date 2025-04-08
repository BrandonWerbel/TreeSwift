#!/usr/bin/env python3
MEMORY_DEN = 1024*1024

# load data
from pickle import load
from sys import argv
if len(argv) != 2:
    print("USAGE: %s <data.pkl>"%argv[0]); exit(1)
if argv[1].lower().endswith('.gz'):
    from gzip import open as gopen
    data = load(gopen(argv[1]))
else:
    data = load(open(argv[1]))

print(data)

# set up seaborn
from matplotlib import rcParams
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set_context("paper", rc={"font.size":12,"axes.titlesize":16,"axes.labelsize":14,"legend.fontsize":10,"xtick.labelsize":10,"ytick.labelsize":10})
sns.set_style("ticks")
rcParams['font.family'] = 'serif'
T = ['DendroPy','Bio.Phylo','ETE Toolkit','TreeSwift', 'NetworkX']
name_to_key = {'DendroPy':'dendropy', 'Bio.Phylo':'biophylo', 'TreeSwift':'treeswift', 'ETE Toolkit':'ete3', 'NetworkX' : 'networkx'}
pal = {'TreeSwift':'black', 'Bio.Phylo':'red', 'DendroPy':'blue', 'ETE Toolkit':'green', 'NetworkX' : 'purple'}
linestyle = {'TreeSwift':'-', 'DendroPy':':', 'Bio.Phylo':'--', 'ETE Toolkit':'-.', 'NetworkX':'solid'}
meta = {
    'name': {
        'height': "Height",
        'average_length': "Average Length",
        'distance': "Distance",
        'treeness': "Treeness",
    }
}
N = None
for m in data:
    if N is None or len(data[m]) > len(N):
        N = sorted(data[m].keys())

# create figures
for m in sorted(data.keys()):
    fig = plt.figure()
    used = []
    for t in T:
        x = []; y = []
        for n in N:
            if n not in data[m]:
                continue
            if name_to_key[t] in data[m][n]:
                x += [n]*10
                y += data[m][n][name_to_key[t]]
            else:
                continue
        if len(x) != 0:
            x = np.array(x); y = np.array(y)
            if m == 'memory':
                y = y/float(MEMORY_DEN)
            sns.pointplot(x=x, y=y, color=pal[t], linestyles=linestyle[t])
            used.append(t)
    plt.yscale('log')
    plt.title(meta['name'][m])
    plt.xlabel("Number of Leaves")
    plt.ylabel("Execution Time (seconds)")
    handles = [Line2D([0],[0],color=pal[t],label=t,linestyle=linestyle[t],linewidth=1.5) for t in used]
    legend = plt.legend(handles=handles,bbox_to_anchor=(0.005, 0.995), loc=2, borderaxespad=0., frameon=True)
    fig.savefig('%s.pdf'%m.lower(), format='pdf', bbox_inches='tight')
