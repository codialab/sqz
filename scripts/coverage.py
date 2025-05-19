#!/usr/bin/env python

from sys import stdout, stderr, exit
from itertools import repeat, chain
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
import re

#
# third party libraries
#
import numpy as np
import networkx as nx

def parse_gfa(data):

    nodes = list()
    wlines = list()
    qlines = list()
    zlines = list()
    for line in data:
        line_type = line[0]
        if line_type == 'S':
            nodes.append(line[2:line.find('\t', 2)])
        elif line_type == 'P':
            repl = {'+': '>', '-': '<'}
            splits1 = line.split('\t')
            splits2 = splits1[2].split(',')
            wlines.append((splits1[1], list(map(lambda x: (repl[x[-1]], x[:-1]), splits2))))
        elif line_type == 'W':
            splits1 = line.split('\t')
            splits2 = re.split(r'(<|>)', splits1[6][:-1])
            wlines.append(('#'.join(splits1[1:4]), list(zip(splits2[1::2],splits2[2::2]))))
        elif line_type == 'Q':
            splits1 = line.split('\t')
            splits2 = re.split(r'(<|>)', splits1[2][:-1])
            qlines.append((splits1[1], list(zip(splits2[1::2],splits2[2::2]))))
        elif line_type == 'Z':
            splits1 = line.split('\t')
            splits2 = re.split(r'(<|>)', splits1[6][:-1])
            zlines.append(('#'.join(splits1[1:4]), list(zip(splits2[1::2],splits2[2::2]))))
    return (nodes, wlines, qlines, zlines)

def build_dag(qlines):

    G = nx.DiGraph()
    for (qid, children) in qlines:
        _, cids = zip(*children)
        G.add_edges_from(list(zip(repeat(qid), cids)))
    return G

def propagate_counts(C, G, v2id):
    for u in nx.topological_sort(G):
        uid = v2id[u]
        for v in G.neighbors(u):
            C[v2id[v]] |= C[uid]

def main(nodes, wlines, qlines, zlines, out):

    coverage = None
    v2id = dict(zip(chain(nodes, map(lambda x: x[0], qlines)), range(len(nodes)+len(qlines))))
    if len(zlines):
        C = np.zeros(shape=(len(v2id), len(zlines)), dtype=bool)
        for i, (_, seq) in enumerate(zlines):
            for _, v in seq:
                C[v2id[v], i] |= True
        G = build_dag(qlines)
        propagate_counts(C, G, v2id)
        coverage = C.sum(axis=1)
    elif len(wlines):
        coverage = np.zeros(shape=len(v2id), dtype=np.int64)
        for _, seq in wlines:
            C = np.zeros(shape=len(v2id), dtype=bool)
            for _, v in seq:
                C[v2id[v]] |= True
            coverage += C.astype(np.int64)

    for i in range(len(nodes)):
        print(f'{nodes[i]}\t{coverage[i]}', file=out)

if __name__ == '__main__':

    description = '''count number of nodes in P/W/Z lines'''
    parser = ArgumentParser(formatter_class=ADHF, description=description)
    parser.add_argument('gfa_file', type=open, help='GFA file')
    args = parser.parse_args()

    out = stdout
    (nodes, wlines, qlines, zlines) = parse_gfa(args.gfa_file)

    main(nodes, wlines, qlines, zlines, out)

