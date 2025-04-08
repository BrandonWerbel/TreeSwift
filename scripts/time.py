#!/usr/bin/env python3
from Bio import Phylo
from itertools import combinations
from treeswift import read_tree_newick
import dendropy
import ete3
import numpy
import networkx
from statistics import mean
NA = "NA" # what to print when function is not implemented

# main code
from io import StringIO
from sys import argv
from time import time
if len(argv) != 4 or argv[2] not in {'treeswift','dendropy','biophylo','ete3', 'networkx'}:
    print("USAGE: %s <tree_file> <treeswift_or_dendropy_or_biophylo_or_ete3> <task>"%argv[0]); exit(1)
if argv[1].lower().endswith('.gz'):
    from gzip import open as gopen
    treestr = gopen(argv[1]).read().decode().strip()
else:
    treestr = open(argv[1]).read().strip()
treeio = StringIO(treestr) # for Bio.Phylo

def avg_tree_length(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        tree.max_distance_from_root()
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        max(tree.depths().values())
        t_end = time()
    elif m == 'networkx':
        tree = Phylo.read(treeio, 'newick')
        G = Phylo.to_networkx(tree)
        t_start = time()
        max(networkx.shortest_path_length(G, source=tree.root).values())
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        tree.height()
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        max(tree.get_distance(leaf) for leaf in tree.iter_leaves())
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

def tree_height(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        tree.mean((edge.length if edge.length is not None else 0) for edge in tree.postorder_edge_iter())
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        max(tree.depths().values())
        t_end = time()
    elif m == 'networkx':
        tree = Phylo.read(treeio, 'newick')
        G = Phylo.to_networkx(tree)
        t_start = time()
        max(networkx.shortest_path_length(G, source=tree.root).values())
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        tree.height()
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        max(tree.get_distance(leaf) for leaf in tree.iter_leaves())
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

def treeness(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        
        t_end = time()
    elif m == 'networkx':
        tree = Phylo.read(treeio, 'newick')
        G = Phylo.to_networkx(tree)
        t_start = time()
        
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

TASKS = {
    'height':tree_height,
    'average_length':avg_tree_length,
    'distance':distance,
    'treeness':treeness
}

# run
if argv[3] not in TASKS:
    print("Invalid task. Valid options: %s"%', '.join(sorted(TASKS.keys()))); exit(1)
print(TASKS[argv[3]](argv[2]))
