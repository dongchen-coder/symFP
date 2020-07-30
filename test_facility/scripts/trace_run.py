#!usr/bin/python

import os
import numpy as np
import subprocess
import sys
import argparse

BENCHMARK = [
	"2mm",
	"3mm",
	"mvt",
	"atax",
	"bicg",
	"gemm",
	"gemver",
	"gesummv",
	"convolution2d",
	"convolution3d",
	"syrk",
	"syr2k",
	"stencil",
	"doitgen"
]

TESTBENCH = [
	"test1",
	"test2"
]

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--all", action='store_true', help='run on all benchmarks')
	parser.add_argument('--prog', default='2mm', help='set the benchmark name to be tuned')
	parser.add_argument('--test', action='store_true', help='run on test benchmarks')
	args = parser.parse_args()
	if args.all:
		for p in BENCHMARK:
			print("Start %s" % p)
			os.system("./bin/"+p+"_trace > ./ss_data/"+p+"-trace.data")
	if args.test:
		for p in TESTBENCH:
			print("Start %s" % p)
			os.system("./bin/"+p+"_trace > ./ss_data/"+p+"-trace.data")
	else:
		p = args.prog
		print("Start %s" % p)
		os.system("./bin/"+p+"_trace > ./ss_data/"+p+"-trace.data")







