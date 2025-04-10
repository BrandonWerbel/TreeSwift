# this patch makes the version of DendroPy we're using work in Python 3.10
#!sed -i '2i import collections; from collections.abc import MutableMapping; collections.MutableMapping = MutableMapping' TreeSwift-Paper/scripts/time.py

# we want 10 replicates per measurement
NUM_REPS = 10

# set the "number of leaves" values we want to run for each task based on Figure 1 of the paper
# we're only going to go up to 10,000 leaves for the sake of time, but feel free to go all the way up to 1,000,000 if you'd like by modifying the lists below
NUM_LEAVES = {
    'height':          [100, 1000, 10000],
    'average_length':  [100, 1000, 10000],
    'distance':        [100, 1000, 10000],
    'treeness':        [100, 1000, 10000],
}

# some tools are particularly inefficient or lack certain functionality, so let's use the same max "number of leaves" for each tool from Figure 1 of the paper
MAX_NUM_LEAVES = {
    'height': {
        'ete3' : 1000
    },
    'distance': {
        'dendropy' : 1000
    }
}

# run the benchmark experiment
import subprocess
data = dict()
for task in NUM_LEAVES:
    data[task] = dict()
    print("=== Running task: %s ===" % task)
    for X in NUM_LEAVES[task]:
        data[task][X] = dict()
        print("- X = %d leaves." % X, end='')
        for tool in ['dendropy', 'treeswift', 'ete3', 'networkx', 'biophylo']:
            if task in MAX_NUM_LEAVES and tool in MAX_NUM_LEAVES[task] and X > MAX_NUM_LEAVES[task][tool]:
                continue # skip this tool if X > max num leaves we want to run
            data[task][X][tool] = list()
            print(" %s" % tool, end='')
            command = ['python3', 'TreeSwift/scripts/time.py', 'TreeSwift/data/tree_n%d.tre.gz' % X, tool, task]
            for replicate in range(1, NUM_REPS + 1):
                print('.', end='')
                print(command)
                measurement = subprocess.check_output(command).decode().strip()
                if measurement != 'NA':
                    data[task][X][tool].append(float(measurement))
        print()


# save `data` to a file called `data.pkl.gz`
import gzip
import pickle
with gzip.open('data.pkl.gz', 'wb') as f:
    pickle.dump(data, f)