#!/usr/bin/env python3

import collections; from collections.abc import MutableMapping; collections.MutableMapping = MutableMapping
from Bio import Phylo
from itertools import combinations
from treeswift import read_tree_newick
import dendropy
import ete3
import numpy
import networkx
from statistics import mean
import random
NA = "NA" # what to print when function is not implemented
#test
# main code
from io import StringIO
from sys import argv
from time import time
import re
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
        numpy.average([0 if edge.length is None else edge.length for edge in tree.postorder_edge_iter()])
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        numpy.average([clade.branch_length for clade in tree.find_clades() if clade.branch_length])
        t_end = time()
    elif m == 'networkx':
        tree = Phylo.read(treeio, 'newick')
        G = Phylo.to_networkx(tree)
        t_start = time()
        numpy.average([data["weight"] for _, _, data in G.edges(data=True)])
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        tree.avg_branch_length(terminal=True, internal=True)
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        numpy.average([node.dist for node in tree.traverse()])
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

def tree_height(m):
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

def treeness(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        total = sum(edge.length for edge in tree.edges() if edge.length is not None)
        internal = sum(edge.length for edge in tree.internal_edges() if edge.length is not None)
        treeness_val = internal / total if total else 0
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        total = sum(clade.branch_length for clade in tree.find_clades() if clade.branch_length)
        internal = sum(clade.branch_length for clade in tree.get_nonterminals() if clade.branch_length)
        treeness_val = internal / total if total else 0
        t_end = time()
    elif m == 'networkx':
        tree = Phylo.read(treeio, 'newick')
        G = Phylo.to_networkx(tree)
        t_start = time()
        total = sum(G[u][v].get('weight', 1.0) for u, v in G.edges())
        leaves = {n for n in G if G.degree[n] == 1 and n != tree.root}
        internal = sum(G[u][v].get('weight', 1.0) for u, v in G.edges() if u not in leaves and v not in leaves)
        treeness_val = internal / total if total else 0
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        total = sum(n.edge_length for n in tree.traverse_preorder() if n.edge_length is not None)
        internal = sum(n.edge_length for n in tree.traverse_internal() if n.edge_length is not None)
        treeness_val = internal / total if total else 0
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        total = sum(n.dist for n in tree.traverse())
        internal = sum(n.dist for n in tree.traverse() if not n.is_leaf())
        treeness_val = internal / total if total else 0
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

def distance(m):  # Alexis - this is complete I beilieve let me know if there are any issues
    matches = re.findall(r"L\d+", treestr)
    taxon1, taxon2 = random.choice(matches), random.choice(matches)
    
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        taxon1, taxon2 = random.choice(tree.taxon_namespace), random.choice(tree.taxon_namespace)    
        t_start = time()
        tree.phylogenetic_distance_matrix().distance(taxon1, taxon2)
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        node1, node2 = list(tree.find_clades(name=taxon1))[0], list(tree.find_clades(name=taxon2))[0]
        t_start = time()
        tree.distance(taxon1, taxon2)
        t_end = time()
    elif m == 'networkx':
        tree = Phylo.read(treeio, 'newick')
        G = Phylo.to_networkx(tree)
        node1 = next(node for node in G if taxon1 == str(node))
        node2 = next(node for node in G if taxon2 == str(node))
        t_start = time()
        networkx.shortest_path_length(G, source=node1, target=node2)
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        node1 = tree.find_node(taxon1, True, True)
        node2 = tree.find_node(taxon2, True, True)
        t_start = time()
        tree.distance_between(node1, node2)
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr, format=1)
        node1 = tree.search_nodes(name=taxon1)[0]
        node2 = tree.search_nodes(name=taxon2)[0]
        t_start = time()
        node1.get_distance(node2)
        t_end = time()
    else:
        assert False, "Invalid tool: %s" % m
    return t_end - t_start



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
