#!/usr/bin/python

import argparse
import numpy as np
import scipy.stats as st
import subprocess
import os
from enum import Enum

benchmark = ["deriche", "fdtd-2d", "heat-3d", "jacobi-1d", "jacobi-2d", "adi", "2mm", "3mm", "mvt", "atax", "bicg", "gemm", "gemver", "gesummv", "stencil"]
# benchmark = ["jacobi-1d","2mm", "atax", "gemver", "fdtd-2d"]

class ReuseType(Enum):
    INTER_CHUNK_REUSE = 0
    SRC_NEIGHBOR_REUSE = 1  # neighbor reuse was found in source's neighbors
    SINK_NEIGHBOR_REUSE = 2 # neighbor reuse was found in sink's neighbors
    SCALE_REUSE = 3
    SRC_FOLD_REUSE = 4
    SINK_FOLD_REUSE = 5

class PLUMReference(object):
    """docstring for PLUMReference"""
    def __init__(self, name, idx):
        super(PLUMReference, self).__init__()
        self.name = name
        self.idx = idx

    def idx_to_tuple(self):
        return tuple(map(int, self.idx[1:-1].split(', ')))

    def is_equal(self, pr):
        print(self.name, pr.name, self.idx, pr.idx)
        return (self.name == pr.name) and (self.idx == pr.idx)
        

class ReuseStat(object):
    """docstring for ReuseStat"""
    def __init__(self, src, sink, ris, rip, rtype):
        super(ReuseStat, self).__init__()
        self.src = src
        self.sink = sink
        self.ris = ris
        self.rip = rip
        self.type = rtype

    def is_equal(self, rs):
        return self.src.is_equal(rs.src)


# load ri stat from file
# ri stat has the following format:
#   src_ref, src_idx, sink_ref, sink_idx, ris, rip, type 
def load_stat(fname):
    stats = list()
    with open(fname) as f:
        for l in f:
            l = l.rstrip()
            l = list(filter(None, l.split(";")))
            rs = ReuseStat(PLUMReference(l[0],l[1]), 
                            PLUMReference(l[2],l[3]),
                            int(l[4]),
                            int(l[5]),
                            int(l[6])
                            )
            stats.append(rs)
    return stats

# compare two list of stat running in different threads
# two reuse pair are equivalent when they have the same src.
# (src1.ref == src2.ref && src1.idx == src2.idx)
def is_equal_reuse_pair(stat1, stat2):
    return stat1.src.is_equal(stat2.src)

# check if the reuse type changes
def is_type_switch(stat1, stat2):
    return stat1.rtype == stat2.rtype

# compute the rip difference
def get_delta_reuse(stat1, stat2):
    return stat2.rip - stat1.rip

def run_command(command):
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return str(stdout.decode())

if __name__ == '__main__':
    tune_result = dict()
    for prog in benchmark:
        tuning_thread = run_command("./bin/{}_tuning".format(prog)) 
        tune_result[prog] = tuning_thread
    for p in tune_result:
        print(f"{p}\t{tune_result[p]}")
