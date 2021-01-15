#!usr/bin/python

import os
import numpy as np
import subprocess
import sys
import argparse
'''
BENCHMARK = [
	"2mm",
	"3mm",
	"mvt",
	"atax",
	"bicg",
	"gemm",
	"gemver",
        "deriche",
	"gesummv",
	"stencil",
        "fdtd-2d",
        "heat-3d",
        "adi",
        "jacobi-1d",
        "jacobi-2d"
]
'''
BENCHMARK = [#"adi", "atax", "deriche", "fdtd-2d", "gesummv", "gemver",
        #"jacobi-2d", "stencil-7p"
        "2mm"]

TESTBENCH = [
	"test1",
	"test2"
]

def dump_stats(fname):
    with open(fname) as f:
        for ln, line in enumerate(f):
            if line.startswith("Total"):
                print(line)
    f.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--all", action='store_true', help='run on all benchmarks')
	parser.add_argument('--prog', default='2mm', help='set the benchmark name to be tuned')
	parser.add_argument('--chk', type=int, default=4, help='set the chunk size')
	parser.add_argument('--exp', default='ps', help='set the interleaving pattern plum will run')
	parser.add_argument('--test', action='store_true', help='run on test benchmarks')
	args = parser.parse_args()
	if args.all:
		try:
			# output = subprocess.check_output(['make', 'check_gen', 'CHK_SIZE='+str(args.chk)], stderr=subprocess.PIPE)
			for p in BENCHMARK:
				print("Start %s" % p)
				# os.system("./bin/"+p+"_"+args.exp+" > "+p+".txt")
				os.system("./bin/"+p+"_"+args.exp+"> ./ss_data/"+p+"-chk-"+str(args.chk)+"-plum-"+args.exp+".data")
				#dump_stats("./ss_data/"+p+"-chk-"+str(args.chk)+"-plum-"+args.exp+".data")
		except subprocess.CalledProcessError as e:
			print("Compile Error %s", e)
	elif args.test:
		try:
			# output = subprocess.check_output(['make', 'test_parallel_gen', 'CHK_SIZE='+str(args.chk)], stderr=subprocess.PIPE)
			for p in TESTBENCH:
				print("Start %s" % p)
				os.system("./bin/"+p+"_"+args.exp+"> ./ss_data/"+p+"-chk-"+str(args.chk)+"-plum-"+args.exp+".data")
				dump_stats("./ss_data/"+p+"-chk-"+str(args.chk)+"-plum-"+args.exp+".data")
		except subprocess.CalledProcessError as e:
			print("Compile Error %s", e)
	else:
		try:
			# output = subprocess.check_output(['make', 'check_gen', 'CHK_SIZE='+str(args.chk)], stderr=subprocess.PIPE)
			p = args.prog
			print("Start %s" % p)
			os.system("./bin/"+p+"_"+args.exp+"> ./ss_data/"+p+"-chk-"+str(args.chk)+"-plum-"+args.exp+".data")
			dump_stats("./ss_data/"+p+"-chk-"+str(args.chk)+"-plum-"+args.exp+".data")
		except subprocess.CalledProcessError as e:
			print("Compile Error %s", e)
