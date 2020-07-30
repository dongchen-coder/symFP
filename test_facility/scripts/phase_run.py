#!usr/bin/python

import os
import numpy as np
import subprocess
import sys
import argparse

def dump_stats(fname):
    with open(fname) as f:
        for ln, line in enumerate(f):
            if line.startswith("Total"):
                print(line)
    f.close()

if __name__ == '__main__':
    for root, dirs, files in os.walk("./bin"):
        for exe in files:
            if "phase" in exe:
                print(os.path.join(root, exe))
                os.system(os.path.join(root, exe) + " > " + exe + ".txt")
                # os.system("./bin/"+exe+" > "+exe+".txt")
