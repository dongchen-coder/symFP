#!/usr/bin/python
import os
import subprocess
import numpy as np

BIN_DIR = "./bin"
BENCHMARK_LIST = "2mm 3mm adi atax bicg cholesky correlation covariance deriche doitgen durbin fdtd_2d floyd_warshall gemm gemver gesummv gramschmidt heat_3d jacobi_1d jacobi_2d lu ludcmp mvt nussinov seidel_2d symm syr2d syrk trisolv trmm".split(" ")
# BENCHMARK_LIST = "atax bicg".split(" ")

EPOCH = 3

def load_oh(content):
	trace_oh_ms = 0.0
	carl_oh_ms = 0.0
	for ln, line in enumerate(content):
		line = line.rstrip()
		if line.startswith("[CARL] Overhead in"):
			carl_oh_ms = float(filter(None, line.split(" "))[-1])
		elif line.startswith("[SPS] Overhead in"):
			trace_oh_ms = float(filter(None, line.split(" "))[-1])
		# if line.startswith("CARL Lease Assignment"):
		# 	collect_carl_oh = True
		# if line.startswith("Overhead in"):
		# 	if collect_carl_oh:
		# 		carl_oh_ms = float(filter(None, line.split(" "))[-1])
		# 	else:
		# 		trace_oh_ms = float(filter(None, line.split(" "))[-1])

	return (trace_oh_ms, carl_oh_ms)


def compute_mean_overhead(prog, trace, carl, total):
	trace_avg_oh = np.mean(trace) / 1000000
	carl_avg_oh = np.mean(carl) / 1000000
	total_avg_oh = np.mean(total) / 1000000
	print("[%s] Average Overhead for SPS is %f (%f), CARL is %f (%f), TOTAL is %f (%f)" %  (prog, trace_avg_oh, trace_avg_oh/60, carl_avg_oh, carl_avg_oh/60, total_avg_oh, total_avg_oh/60))

if __name__ == '__main__':
	for p in BENCHMARK_LIST:
		trace_time_cost = list()
		carl_time_cost = list()
		total_time_cost = list()
		for i in xrange(0, EPOCH):
			proc = subprocess.Popen([BIN_DIR+"/"+p+"_ref_osl_trace"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdout, stderr = proc.communicate()
			(trace_oh, carl_oh) = load_oh(filter(None, stdout.split("\n")))
            # {oh}
			trace_time_cost.append(trace_oh)
			carl_time_cost.append(carl_oh)
			total_time_cost.append(trace_oh + carl_oh)
		compute_mean_overhead(p, trace_time_cost, carl_time_cost, total_time_cost)

