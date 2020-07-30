#!usr/local/bin/python
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *
import numpy as np
import argparse
import subprocess

# BENCHMARK = ["stencil", "bicg", "2mm", "syrk"]
BENCHMARK = [
	"stencil",
	"2mm",
	"3mm",
	"bicg",
	"atax",
	"gemm",
	"convolution2d",
	"convolution3d",
	"gemver",
	"mvt",
	"gesummv",
	"syrk",
	"syr2k"
]

# fake_output = '''
# Start to dump reuse time histogram for A0
# 1 10
# 2 5
# Start to dump reuse time histogram for TMP0
# 2 3
# 5 10
# 20 512
# Start to dump reuse time histogram fro B3
# 1 50
# 46 7
# '''

RTCOLOR = ["#4472CA", "#5E7CE2", "#92B4F4", "#CFDEE7"]

BIN_DIR = "./instr_code/bin"
HOME_DIR = "./instr_code/data"
PLOT_DIR = "./instr_code/data"

def load_rt(l, delimiter=" "):
	rthist = dict()
	start_load = False
	refName = ""
	for ln, line in enumerate(l):
		line = line.rstrip()
		if line.startswith("Start to dump reuse time histogram"):
			if ln == 0 or refName != line.split(" ")[-1]:
				start_load = True
				refName = line.split(" ")[-1]
				rthist[refName] = dict()
				continue
		if start_load:
			rt = float(list(filter(None, line.split(delimiter)))[0])
			rtcnt = float(list(filter(None, line.split(delimiter)))[1].split(" ")[0])
			rthist[refName][rt] = [rtcnt]
	sumrtcnt = 0
	for ref in rthist:
		sumrtcnt += sum([rthist[ref][r][0] for r in rthist[ref]])
	for ref in rthist:
		for rt in rthist[ref]:
			if len(rthist[ref][rt]) == 2:
				rthist[ref][rt] = [float(rthist[ref][rt][0]) / sumrtcnt, float(rthist[ref][rt][1]) / sumrtcnt]
			else:
				rthist[ref][rt] = [float(rthist[ref][rt][0]) / sumrtcnt]
	return rthist


def plot_rtHisto(prog, rthist):
	refcnt = len(rthist)
	assert refcnt != 0, 'rthist is empty'
	if refcnt > 1:
		col = 3
		row = refcnt / col if (refcnt % col == 0) else (refcnt / col) + 1
	else:
		row = 1
		col = 1
	# fig = figure(figsize=(12,8))
	fig_dir = PLOT_DIR + "/" + prog + "-histograms"
	os.makedirs(fig_dir)
	for i, ref in enumerate(rthist):
		fig, ax = subplots()
		# ax = fig.add_subplot(row,col,i+1)
		#avg
		x = sorted(rthist[ref].keys())
		y = [rthist[ref][e][0] for e in x]
		index = np.arange(len(y))
		opacity = 1.0
		bar_width = 0.45
		# plot the reuse interval histo in one subfigure
		rect = ax.bar(index, y, bar_width, alpha=opacity, color=RTCOLOR[0])
		ax.set_title(ref, fontsize = 15, y = 0.90)
		ax.set_xlabel('Reuse Interval', fontsize=15)
		ax.set_ylabel('Reference Count', fontsize=15)
		ax.set_xticks(np.arange(len(rthist[ref])) + bar_width / 2)
		ax.set_xticklabels(tuple([str(i) for i in x]))
		ax.xaxis.set_tick_params(rotation=45)
		# if i == 0:
		# 	fig.legend(ncol=len(PLUM_EXPLISTS), mode="expand", loc = "upper center", fontsize=20)
		title(prog+"-"+ref)
		grid(True)
		tight_layout(pad=0.4, w_pad=0.4, h_pad=-1.2)
		subplots_adjust(left=0.13, right = 0.98, top = 0.90, bottom = 0.21, hspace=0.5)
		savefig(fig_dir + "/" + ref + "-rid.pdf")

def main(is_all=False, prog=None):
	# output = subprocess.check_output(['make', 'instr_code'], stderr=subprocess.PIPE)
	if is_all:
		for prog in BENCHMARK:
			try:
				print("collect from " + prog)
				proc = subprocess.Popen([BIN_DIR+"/"+prog], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				stdout, stderr = proc.communicate()
				data = load_rt(filter(None, stdout.split("\n")))
				plot_rtHisto(prog, data)
			except subprocess.CalledProcessError, e:
				print("Compile Error for benchmark %s" % prog)
	else:
		try:
			print("collect from " + prog)
			proc = subprocess.Popen([BIN_DIR+"/"+prog], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdout, stderr = proc.communicate()
			data = load_rt(filter(None, stdout.split("\n")))
			# data = load_rt(filter(None, fake_output.split("\n")))
			plot_rtHisto(prog, data)
		except subprocess.CalledProcessError, e:
			print("Compile Error for benchmark %s" % prog)
		


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--all", action='store_true', help='run on all benchmarks')
	parser.add_argument('--prog', default='2mm', help='set the benchmark name to be tuned')
	args = parser.parse_args()
	main(is_all=args.all, prog=args.prog)
