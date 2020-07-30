
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

P = {0:1}
DATA_DIR = "ss_data"
BENCHMARK_LIST = [
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

EXP_LIST = [
	("Baseline", "-chk-4-avg"),
	("Acc UI", "-chk-4-plum-ps"),
	("Acc UI Uniform Smoothing", "-chk-4-plum-ps-uismoothing"),
	("Acc UI Gaussian Smoothing", "-chk-4-plum-ps-ndsmoothing")
]

def rth_to_rid(fname, prog):
	rid = dict()
	rih = load_rt(DATA_DIR+"/"+fname+".data", prog=prog, delimiter=",")
	n = sum(rih.values())
	accumulate_num_RT = 0
	for k in sorted(rih.keys(), reverse=True):
		P[k] = accumulate_num_RT / n
		rid[k] = rih[k] / n
		accumulate_num_RT += rih[k]
	return rid

def mrc_to_rdh(fname, prog):
	rdh = dict()
	mrc = load_mrc(DATA_DIR+"/"+fname, prog=prog, delimiter=",")
	prev_c = 0
	for i, c in enumerate(sorted(mrc.keys())):
		if i > 0:
			# print("rd[{}] is {} - {} = {}".format(c, mrc[prev_c], mrc[c], mrc[prev_c] - mrc[c]))
			rdh[c] = mrc[prev_c] - mrc[c]
			prev_c = c
	rdh = dict(filter(lambda elem: elem[1] != 0.0, rdh.items()))
	return rdh

def load_rt(filename, prog=None, delimiter=" "):
	rthist = dict()
	f = open(filename)
	start_load = False
	with open(filename) as f:
		for ln, line in enumerate(f):
			line = line.rstrip()
			if line == ("mean-std rid for " + prog) or line == "Start to dump reuse time histogram":
				start_load = True
				continue
			if line == ("mean-std mrc for " + prog) or line == "miss ratio":
				start_load = False
			if start_load:
				rt = float(list(filter(None, line.split(delimiter)))[0])
				rtcnt = float(list(filter(None, line.split(delimiter)))[1])
				rthist[rt] = rtcnt
	f.close()
	return rthist


def load_mrc(filename, prog=None, delimiter=" "):
	mrc = dict()
	f = open(filename)
	start_load = False
	with open(filename) as f:
		for ln, line in enumerate(f):
			line = line.rstrip()
			if line == ("mean-std mrc for " + prog) or line == "miss ratio":
				start_load = True
				continue
			if start_load:
				csize = int(list(filter(None, line.split(delimiter)))[0])
				mr = float(list(filter(None, line.split(delimiter)))[1])
				mrc[csize] = mr
	f.close()
	return mrc


def aet_cal(c):
	c_sim = 0
	aet = 0;
	prev_aet = 0;
	max_RT = max(P.keys())
	while (c_sim < c and aet <= max_RT):
		if aet in P:
			c_sim += P[aet]
			# print("[rt = %d] Add %f in c_sim" % (aet, P[aet]))
			prev_aet = aet
		else:
			c_sim += P[prev_aet]
			# print("[rt = %d] Add %f in c_sim" % (prev_aet, P[prev_aet]))
		aet += 1
	return prev_aet

def rid_acc_cal(rid1, rid2):
	for k in rid1:
		if k not in rid2:
			rid2[k] = 0.0
	for k in rid2:
		if k not in rid1:
			rid1[k] = 0.0
	diff = 0.0
	for k in rid1:
		diff += abs(rid1[k] - rid2[k])
	return 1 - (diff / 2)
	

def plot_dict(dict_g1, dict_g2, labels, xlabel, ylabel, title, savefn):
	keys_g1 = dict_g1.keys()
	keys_g2 = dict_g2.keys()
	x_groups = sorted(list(set(keys_g1).union(set(keys_g2))))
	values_g1 = list()
	values_g2 = list()
	for k in x_groups:
		if k in dict_g1:
			values_g1 += [dict_g1[k]]
		else:
			values_g1 += [0.0]
	for k in x_groups:
		if k in dict_g2:
			values_g2 += [dict_g2[k]]
		else:
			values_g2 += [0.0]
	x = np.arange(len(x_groups))  # the label locations
	width = 0.5  # the width of the bars

	fig, ax = plt.subplots(figsize=(10,6))
	rects1 = ax.bar(x - width/2, values_g1, width, label=labels[0], color="#DF485E")
	rects2 = ax.bar(x + width/2, values_g2, width, label=labels[1], color="#6B4E55")

	# Add some text for labels, title and custom x-axis tick labels, etc.
	ax.set_ylabel(ylabel, fontsize=14)
	ax.set_xlabel(xlabel, fontsize=14)
	ax.set_title(title, fontsize=18)
	ax.set_xticks(x)
	ax.set_xticklabels(x_groups, rotation=90)
	ax.set_ylim([0,0.05])
	ax.legend()
	def autolabel(rects):
	    """Attach a text label above each bar in *rects*, displaying its height."""
	    for rect in rects:
	        height = rect.get_height()
	        if height == 0.0:
	        	continue
	        if height >= 0.05:
	        	ax.text(rect.get_x(), 0.05*0.95, str(round((rect.get_height()), 4)), rotation=90, fontweight='bold')
	        else:
	        	ax.text(rect.get_x(), rect.get_height()*0.95, str(round((rect.get_height()), 4)), rotation=90, fontweight='bold')
	autolabel(rects1)
	autolabel(rects2)
	fig.tight_layout()
	plt.savefig(savefn+".pdf", format='pdf')

def main(prog):
	print("[%s]" % prog)
	rid_g1 = rth_to_rid(prog+EXP_LIST[0][1], prog)
	aet_baseline = aet_cal(args.csize)
	print("[BASELINE] AET for cache size {} is {}".format(args.csize, aet_baseline))
	print("miss ratio for {} is {}".format(aet_baseline, P[aet_baseline]))
	miss_ratio_baseline = P[aet_baseline]
	P.clear()
	P[0] = 1
	for i, exp in enumerate(EXP_LIST):
		if i == 0: continue
		rid_g2 = rth_to_rid(prog+EXP_LIST[i][1], prog)
		aet_plum = aet_cal(args.csize)
		print(f"[{EXP_LIST[i][1]}] AET for cache size {args.csize} is {aet_plum}")
		print("miss ratio for {} is {}".format(aet_plum, P[aet_plum]))
		miss_ratio_plum = P[aet_plum]
		print(f"[{EXP_LIST[i][0]}] Reuse Interval Distribution Accuracty is %.4f" % (rid_acc_cal(rid_g1, rid_g2)))
		P.clear()
		P[0] = 1
		# plot_dict(rid_g1, rid_g2, labels=["BASELINE", EXP_LIST[i][0]], xlabel="Reuse Interval", ylabel="Probability", title="[%s] Reuse Interval Distribution" % prog, savefn="%s_baseline_plum-uismoothing_rid" % prog)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Process some integers.')
	parser.add_argument('--csize', type=int,
                    help='cache size')
	parser.add_argument('--all', action='store_true', help='run all benchmarks')
	parser.add_argument('--prog', type=str,
                    help='benchmark name')
	args = parser.parse_args()
	if args.all:
		for p in BENCHMARK_LIST:
			main(p)
			P = {0:1}
	else:
		main(args.prog)
