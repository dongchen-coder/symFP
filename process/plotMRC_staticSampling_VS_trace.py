
import os
import matplotlib.pyplot as plt

path = "../test/polyBench/result/"
files = os.listdir(path)

ss_c = {}
ss_mr = {}

for f in files:
	
	FileName = f.replace("_staticSampling_result.txt", "")

	if "result.txt" in f:
	
		File = open(path + f, "r")

		ss_c[FileName] = []
		ss_mr[FileName] = []

		mrFlag = 0
		for line in File:
			lineTmp = line.replace("\n","").split(" ")
			if (lineTmp[0] == "miss"):
				mrFlag = 1
				continue
			if (mrFlag == 0):
				continue
			if (lineTmp[0].isdigit()):
				ss_c[FileName].append(float(lineTmp[0]) * 8 / 1024)
				ss_mr[FileName].append(float(lineTmp[1]))
		
				if (float(lineTmp[1]) == 0):
					break

path = "../trace/polyBench/result/"
files = os.listdir(path)

trace_c = {}
trace_mr = {}

for f in files:

	FileName = f.replace("_trace_result.txt", "")
	
	if "result.txt" in f:
		
		File = open(path + f, "r")

		trace_c[FileName] = []
		trace_mr[FileName] = []

		mrFlag = 0
		for line in File:
			lineTmp = line.replace("\n","").split(" ")
			if (lineTmp[0] == "miss"):
				mrFlag = 1
				continue
			if (mrFlag == 0):
				continue
			if (lineTmp[0].isdigit()):
				trace_c[FileName].append(float(lineTmp[0]) * 8 / 1024)
				trace_mr[FileName].append(float(lineTmp[1]))

				if (float(lineTmp[1]) == 0):
					break


for key in trace_c.keys():
	if (key in ss_c.keys() ):

		plt.plot(ss_c[key], ss_mr[key], label = 'Static sampling')
		plt.plot(trace_c[key], trace_mr[key], label = 'Trace')
		plt.xscale('log')
		plt.tick_params(axis='both', which='major', labelsize=18)
		plt.xlabel('Cache Size (KB)', fontsize = 18)
		plt.ylabel('Miss Ratio', fontsize = 18)
		plt.legend(fontsize = 16)
		plt.savefig(key + "_Fig.pdf")		
		plt.clf()	
		
