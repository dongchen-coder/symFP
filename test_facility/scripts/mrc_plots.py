from matplotlib.pyplot import *
matplotlib.use('Agg')
import numpy as np
import argparse

# BENCHMARK = ["stencil", "bicg", "2mm", "syrk"]
NEIGHBOR_BENCHMARK = [
	"stencil", 
	"convolution2d", 
	"convolution3d"
]

SCALE_BENCHMARK = [
	"2mm", 
	"3mm", 
	"gemm",
	"syr2k"
]

MIX_BENCHMARK = [
	"atax", 
	"bicg", 
  	"mvt", 
	"gesummv",
	"gemver", 
	"syrk"
]

EXTRA_BENCH = [
	"stencil", 
	"3mm",
	"convolution2d", 
	"convolution3d",
	"gemver",
	"syrk",
	"syr2k"
]

BENCHMARK = NEIGHBOR_BENCHMARK + SCALE_BENCHMARK + MIX_BENCHMARK

YLIM = {
	"2mm": 0.3, 
	"3mm": 0.3, 
	"gemm": 0.3,
	"stencil": 0.08, 
	"convolution2d": 0.08, 
	"convolution3d": 0.1, 
	"atax": 0.08, 
	"bicg": 0.08, 
  	"mvt": 0.2, 
	"gesummv": 0.04,
	"gemver": 0.15, 
	"syrk": 0.1,
	"syr2k": 0.1
}
# RTCOLOR = ["#0c120c", "#c20114", "#6d7275", "#c7d6d5", "#ecebf3"]
# RTCOLOR = ["#513543", "#D983A6", "#568259", ""]
RTCOLOR = ["#DE858A", "#0096FF", "#332C4A", "#FF9300"]
STYLE = ['-', '-.', '--', ':', "-"]
HOME_DIR = "./ss_data"
PLOT_DIR = "."
MARK = ["o", "x", "s", "^", "+"]

PLUM_EXPLISTS = [
	("chk-4-avg-sample", "Baseline"),
	("chk-4-plum-model", "PLUM"),
	# ("chk-4-plum-model", "Model"),
	("chk-4-plum-smoothing", "Smoothing")
]

def load_capacity_data(filename):
	cs = 0.0
	with open(filename) as f:
		for ln, line in enumerate(f):
			line = line.rstrip()
			if line.startswith("Number of Cache Line"):
				cs = float(list(filter(None, line.split(":")))[1])
				break

	f.close()
	return cs

def read_reuse_interval_histogram(fname, delimiter=" "):
	print(f"Read rih/mrc from {fname}")
	rih = dict()
	mrc = dict()
	start_load_rid = False
	start_load_mrc = False
	with open(fname) as f:
		for line in f:
		    line = line.rstrip()
		    if line.startswith("mean-std rid for") or line.startswith("Start to dump reuse time"):
		    	start_load_rid = True
		    	continue
		    if line.startswith("mean-std mrc for") or line.startswith("miss ratio"):
		        start_load_rid = False
		        start_load_mrc = True
		        continue
		    if start_load_rid:
		        rt = int(list(filter(None, line.split(delimiter)))[0])
		        rtcnt = float(list(filter(None, line.split(delimiter)))[1])
		        rih[rt] = rtcnt
		    if start_load_mrc:
		    	csize = int(list(filter(None, line.split(delimiter)))[0])
		    	mr = float(list(filter(None, line.split(delimiter)))[1])
		    	mrc[csize] = mr
	s = sum(rih.values())
	for ri in rih:
		rih[ri] = rih[ri] / s
	return rih, mrc


def plot_mr(fname):
	imgcnt = len(BENCHMARK)
	if imgcnt > 1:
		col = 3
		row = imgcnt / col if (imgcnt % col == 0) else (imgcnt / col) + 1
	else:
		row = 1
		col = 1
	fig = figure(figsize=(12,8))
	for i, prog in enumerate(BENCHMARK):
		# load the MRC data
		ax = fig.add_subplot(row,col,i+1)
		axes = list()
		data_occupancy = load_capacity_data(HOME_DIR+"/"+prog+"-chk-4-avg.data")
		color = RTCOLOR
		for j in range(0, len(PLUM_EXPLISTS)):
			exp = PLUM_EXPLISTS[j][0]
			if exp == "chk-4-plum-ss": continue
			# data: [[list_of_cachesize], [list_of_miss_ratio]]
			rih, mrc = read_reuse_interval_histogram(HOME_DIR+"/"+prog+"-"+exp+".data", delimiter=", ")
			cache_size = sorted(mrc.keys())
			miss_ratio = [mrc[c] for c in cache_size]
			# print(cache_size)
			# print(miss_ratio)
			assert len(cache_size) == len(miss_ratio), "Length not match"
			# plot the miss ratio in one subfigure
			axes.append(plot(cache_size, miss_ratio, label=PLUM_EXPLISTS[j][1], linestyle=STYLE[j], linewidth=2, color=color[j], marker=MARK[j], markersize=4.0))
			# if len(data) == 3:
			# 	axes.append(errorbar(cache_size, miss_ratio, mr_std, label=PLUM_EXPLISTS[j][1], fmt=STYLE[j], linewidth=2, color=color[j], marker=MARK[j], markersize=4.0))
			# else:
			# 	axes.append(plot(cache_size, miss_ratio, label=PLUM_EXPLISTS[j][1], linestyle=STYLE[j], linewidth=2, color=color[j], marker=MARK[j], markersize=4.0))

		# print(len(axes))
		if i == 0:
			fig.legend(ncol=len(PLUM_EXPLISTS), mode="expand", loc = "upper center", fontsize=20)

		# print(data_occupancy)
		ax.axvline(data_occupancy, color='r', lw=2.0, linestyle="--")
		ax.set_title(prog, fontsize = 18, y = 0.7, fontweight='bold')
		ax.set_xlabel('Cache Size')
		ax.set_ylabel('Miss Ratio')
		# ax.set_xscale('linear')
		ax.set_xscale('log')
		# ax.legend(fontsize=20)
		# grid(True)
		ylim((0, 1.0))
		# ylim((0,YLIM[prog]))
		# xlim((10, 1000))
		subplots_adjust(left=None, bottom=None, right=None, top=0.90, wspace=0.26, hspace=0.46)
	tight_layout(pad=0.4, w_pad=0.4, h_pad=-1.2)
	# subplots_adjust(left=0.08, right = 0.98, top = 0.90, bottom = 0.1, hspace=0.5)
	subplots_adjust(left=0.08, right = 0.98, top = 0.91, bottom = 0.1, hspace=0.5)
	# savefig(PLOT_DIR + "/2mm-sample.pdf")
	savefig(PLOT_DIR + "/" + fname)
	# show()


def plot_rt(fname, ctrl, ref):
	imgcnt = len(BENCHMARK)
	if imgcnt > 1:
		row = imgcnt / 3 if (imgcnt % 3 == 0) else (imgcnt / 3) + 1
		col = 3
	else:
		row = 1
		col = 1
	fig = figure(figsize=(12,8))
	for i, prog in enumerate(BENCHMARK):
		ax = fig.add_subplot(row,col,i+1)
		rih_ctrl, mrc_ctrl = read_reuse_interval_histogram(HOME_DIR+"/"+prog+"-"+ctrl[0]+".data", delimiter=", ")
		rih_ref, mrc_ref = read_reuse_interval_histogram(HOME_DIR+"/"+prog+"-"+ref[0]+".data", delimiter=", ")
		rtset = set(rih_ctrl.keys()) | set(rih_ref.keys())
		for rt in rtset:
			if rt not in rih_ref:
				rih_ref[rt] = 0.0
			elif rt not in rih_ctrl:
				rih_ctrl[rt] = 0.0
		assert len(rih_ctrl) == len(rih_ref), f"Lenght Error for RT. {len(rih_ctrl)}, {len(rih_ref)}"
		x = sorted(rtset)
		y1 = [rih_ctrl[e] for e in x]
		y2 = [rih_ref[e] for e in x]
		index = np.arange(len(y1))
		opacity = 1.0
		bar_width = 0.45
		# plot the reuse interval histo in one subfigure
		rect1 = ax.bar(index-bar_width/2, y1, bar_width, label=ctrl, alpha=opacity, color=RTCOLOR[0])
		
		rect2 = ax.bar(index+bar_width/2, y2, bar_width, label=ref, alpha=opacity, color=RTCOLOR[1])
		# if i == 0:
		# 	fig.legend(ncol=len(PLUM_EXPLISTS), mode="expand", loc = "upper center")
		ax.set_title(prog, fontsize = 15, y = 0.90)
		ax.set_xlabel('Reuse Interval', fontsize=15)
		ax.set_ylabel('Reference Count', fontsize=15)
		ax.set_xticks(np.arange(len(rih_ctrl)))
		ax.set_xticklabels(tuple([str(i) for i in x]))
		ax.xaxis.set_tick_params(rotation=45)
		if i == 0:
			fig.legend(ncol=len(PLUM_EXPLISTS), mode="expand", loc = "upper center", fontsize=20)
		grid(True)
	tight_layout(pad=0.4, w_pad=0.4, h_pad=-1.2)
	subplots_adjust(left=0.13, right = 0.98, top = 0.90, bottom = 0.21, hspace=0.5)
	savefig(PLOT_DIR + "/"+fname)

def cal_accuracy_rt(rth_ctrl, rth_ref):
	rtset = set(rth_ctrl.keys()) | set(rth_ref.keys())
	for rt in rtset:
		if rt not in rth_ref:
			rth_ref[rt] = 0.0
		elif rt not in rth_ctrl:
			rth_ctrl[rt] = 0.0
	cnt = [0.0] * 4
	acc = [0.0] * 4
	for rt in rtset:
		err = abs(rth_ctrl[rt] - rth_ref[rt])
		cnt[3] += err
		if rt < 10:
			cnt[0] += err
		elif rt >= 10 and rt < 1000:
			cnt[1] += err
		else:
			cnt[2] += err
	for i in range(len(cnt)):
		acc[i] = (1 - (cnt[i] / 2)) * 100
	return acc


def monotone_check(mrc):
	prev_c = -1
	for c in sorted(mrc.keys()):
		if prev_c == -1:
			prev_c = c
			continue
		if mrc[c] > mrc[prev_c]:
			return False
		prev_c = c
	return True

def make_up_mr(mrc, c):
	mr_lb = 1.0
	mr_ub = 0.0
	for cs in sorted(mrc.keys()):
		if cs < c:
			mr_lb = mrc[cs]
		if cs > c:
			break
	return mr_lb


def cal_accuracy_mrc(bench, control, reference):
	print(f"{control} vs. {reference}")
	mean_err = [[], [], [], []]
	# fig = figure(figsize=(8,6))
	for p, prog in enumerate(bench):
		# load the MRC data
		mrdata = dict()
		cset = set()
		mrdata[control[0]] = dict()
		mrdata[reference[0]] = dict()

		rih_ctrl, mrc_ctrl = read_reuse_interval_histogram(HOME_DIR+"/"+prog+"-"+control[0]+".data", delimiter=", ")
		rih_ref, mrc_ref = read_reuse_interval_histogram(HOME_DIR+"/"+prog+"-"+reference[0]+".data", delimiter=", ")
		cset = set(mrc_ctrl.keys()) | set(mrc_ref.keys())
		for c in cset:
			if c not in mrc_ctrl:
				# cache size miss in data1
				mrc_ctrl[c] = make_up_mr(mrc_ctrl, c)
			elif c not in mrc_ref:
				# cache size miss in data2
				mrc_ref[c] = make_up_mr(mrc_ref, c)
		assert len(mrc_ctrl) == len(cset), 'Length Error (%d vs %d)' % ( len(mrc_ctrl), len(cset)) 
		assert len(mrc_ref) == len(cset), 'Length Error (%d vs %d)' % ( len(mrc_ref), len(cset)) 
		# print("MRC for " + PLUM_EXPLISTS[1][0])
		# for i in sorted(cset):
		# 	print("({},{})\t({},{})".format(i, mrdata[PLUM_EXPLISTS[0][0]][i],i, mrdata[PLUM_EXPLISTS[1][0]][i]))
		num = [0,0,0,0]
		cnt = [0.0,0.0,0.0,0.0]
		accuracy = [0.0, 0.0, 0.0, 0.0]
		max_err = 0
		for c in sorted(cset):
			err = abs(mrc_ctrl[c] - mrc_ref[c])
			num[3] += 1
			cnt[3] += err
			max_err = max(max_err, err)
			if c < 10:
				num[0] += 1
				cnt[0] += err
			elif c >= 10 and c < 1000:
				num[1] += 1
				# print c, abs(mrdata['sls'][c] - mrdata['avg'][c])
				cnt[1] += err
			else:
				num[2] += 1
				cnt[2] += err
		# print "%s\t%.4f\t%.4f" % (prog, 1 - cnt / len(mrdata['sls']), max_err)
		assert sum(num[0:3]) == num[3], 'Sum Num Error (%d vs %d)' % ( sum(num[0:3]), num[3]) 
		# assert sum(cnt[0:3]) == cnt[3], 'Sum Err Error (%.4f vs %.4f)' % ( sum(cnt[0:3]), cnt[3]) 
		for i in range(0, len(num)):
			if num[i] == 0:
				continue
			accuracy[i] = (1 - cnt[i] / num[i]) * 100.0
			mean_err[i] += [accuracy[i]]
		rtacc = cal_accuracy_rt(rih_ctrl, rih_ref)
		print ("[RT] %s\t%.2f\t%.2f\t%.2f\t%.2f" % (prog, rtacc[0], rtacc[1], rtacc[2], rtacc[3]))
		# print("Reuse Portion:")
		# print ("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (prog, num[0]/num[3], num[1]/num[3], num[2]/num[3], num[3], max_err))
		# print("Miss Ratio Acc:")
		print ("[MR] %s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (prog, accuracy[0], accuracy[1], accuracy[2], accuracy[3], max_err))
		# print("%s\t%.4f\t%.4f" % (prog, 1 - csnt[3] / num[3],max_err))
		
	print ("Mean: %.2f\t%.2f\t%.2f\t%.2f" % (np.mean(mean_err[0]), np.mean(mean_err[1]), np.mean(mean_err[2]), np.mean(mean_err[3])))



if __name__ == '__main__':
	'''
	{
		expname: [[list_of_cachesize], [list_of_miss_ratio]]
	}
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument('--img', default='all', help='set the figure name to be stored')
	args = parser.parse_args()
	plot_mr(f"{args.img}.pdf")
	for i, exp in enumerate(PLUM_EXPLISTS):
		if i == 0: continue
		# plot_rt(f"{args.img}-rih-{exp[1]}.pdf", PLUM_EXPLISTS[0], exp)
		# print("Neighbor Bench:")
		cal_accuracy_mrc(BENCHMARK, PLUM_EXPLISTS[0], exp)
	'''
        for i, exp in enumerate(PLUM_EXPLISTS):
		if i == 0: continue
		print("Scale Bench:")
		cal_accuracy_mrc(SCALE_BENCHMARK, PLUM_EXPLISTS[0], exp)
	for i, exp in enumerate(PLUM_EXPLISTS):
		if i == 0: continue
		print("Mix Bench:")
		cal_accuracy_mrc(MIX_BENCHMARK, PLUM_EXPLISTS[0], exp)
	# plot_rt(PLUM_EXPLISTS[1], PLUM_EXPLISTS[2])
        '''
