#!usr/bin/python3

import os
import numpy as np
import testlib
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
	"deriche",
	"gesummv",
	"stencil",
	"fdtd-2d",
	"heat-3d",
	"adi",
	"jacobi-1d",
	"jacobi-2d"
]

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	pparser.add_argument('-p', '--prog', action='store', dest='benchmarks',
					type=str, nargs='*', default=BENCHMARK,
					help='set the program to be run, if not set, run all benchmarks')
	parser.add_argument('-t', '--tgroup', action='store', dest='thread_cnts',
					type=int, nargs='*', default=[2,4,6,8,10],
					help='set the set of threds to be run, if not set, run the default [2,4,6,8,10] setting')
	parser.add_argument('-e', '--epoch', type=int, default=5,
					help='set the times each program will run')
	parser.add_argument('--exp', default='ps', help='set the interleaving pattern plum will run')
	args = parser.parse_args()
	for p in sorted(args.benchmarks):
		print(f"Start {p}")