import matplotlib.pyplot as plt
import os

lru_result_path = './poly_trace_lru_2/'
carl_result_path = './ref_trace_osl_result/'
#carl_result_low_path = './carl_onlyRueses/'	
#carl_result_upper_path = './carl_64B_result/'
#carl_result_random_path = './carl_randomLeaseAssignment/'
opt_result_path = './poly_trace_opt_mrc/'
rll_maxNumOfCLSize_path = "./poly_trace_rll_maxCLSize/"


names = ["2mm", "3mm", "adi", "atax", "bicg", "cholesky", "correlation", "covariance", "deriche", "doitgen", "durbin", "fdtd_2d", "floyd_warshall", "gemm", "gemver", "gesummv", "gramschmidt", "heat_3d", "jacobi_1d", "jacobi_2d", "lu", "ludcmp", "mvt", "nussinov", "seidel_2d", "symm", "syr2d", "syrk", "trisolv", "trmm"]

yMax = {"mvt" : 0.3, "3mm" : 0.4, "bicg" : 0.1, "gemm" : 0.5, "heat_3d" : 0.2, "trmm" : 0.8, "covariance" : 1, "durbin" : 0.2, "doitgen" : 1, "deriche" : 0.4, "symm" : 0.7, "syr2d" : 0.6, "jacobi_2d" : 0.2, "adi" : 0.3, "gramschmidt" : 1, "gesummv" : 0.1, "correlation" : 0.8, "floyd_warshall" : 0.1, "atax" : 0.1, "ludcmp" : 0.8, "seidel_2d" : 0.1, "trisolv" : 0.1, "gemver" : 0.2, "fdtd_2d" : 0.2, "2mm" : 0.5, "jacobi_1d" : 0.1, "cholesky" : 0.2, "lu" : 0.5, "syrk" : 0.1, "nussinov" : 0.4}

fontSizeSet = 20

def processLRU(path, name, cls):

	lru_c = []
	lru_rd = []
	lru_mr = []
	
	lru_total_rd_cnt = 0
	coldMiss = 0

	#content = open(path + 'lru_64B_' + name + '_mrc.txt', 'r')
	content = open(path + name + "_trace_result.txt", 'r')
	for line in content:
		if ('Bin' in line):
			lineList = line.split()
			if ('to' not in line):
				lru_rd.append([int(lineList[4]), int(lineList[7])])
				lru_total_rd_cnt += int(lineList[7])
			else:
				lru_rd.append([int(lineList[6]), int(lineList[9])])
				lru_total_rd_cnt += int(lineList[9])
		elif ("Total number of data blocks is" in line):
			lineList = line.split()
			coldMiss = int(lineList[7])
	content.close()

	lru_total_hits = 0
	for item in lru_rd:
		lru_c.append(float(item[0]+1) * cls / 1024)
		lru_total_hits += item[1]
		lru_mr.append(1 - (float(lru_total_hits) / lru_total_rd_cnt) + float(coldMiss)/lru_total_rd_cnt)	
	
	return lru_c, lru_mr

def processCARL(path, name, cls):
	
	carl_c = []
	carl_mr = []
	total_reuse = 0

	content = open(path + name + '_ref_osl_trace_result.txt', 'r')
	for line in content:
		if ('Assign lease' in line and 'avg cache size' in line and 'miss ratio' in line):
			linetmp = line.split()
			carl_c.append(float(linetmp[9]) * cls / 1024)
			carl_mr.append(float(linetmp[12]))
		if ('Ref' in line and 'RI' in line and 'CNT' in line):
			linetmp = line.split()
			total_reuse += int(linetmp[5])
	content.close()

	return carl_c, carl_mr

def processOPT(path, name, cls):
	opt_c = []
	opt_mr = []
	
	content = open(path + "opt_" + name + ".cpp.o_mrc", 'r')
	for line in content:
		lineList = line.split()
		if (len(lineList) == 2):
			if (lineList[0].isdigit()):
				opt_c.append(float(lineList[0]))
				opt_mr.append(float(lineList[1]))

	while(len(opt_mr) > 2):
		if (opt_mr[-1] == opt_mr[-2]):
			del opt_mr[-1]
			del opt_c[-1]
		else:
			break
	return opt_c, opt_mr

def processMaxCLSizes(path, name, cls, cacheSizes):
	maxCLSizes_rll = []
	for i in range(len(cacheSizes)):
		if (os.path.exists(path + name + "_lease_" + str(i) + "_rll_maxCLSize.txt")):
			content = open(rll_maxNumOfCLSize_path + name + "_lease_" + str(i) + "_rll_maxCLSize.txt")
			for line in content:
				maxCLSizes_rll.append(float(line) * 64 / 1024)
			content.close()
		else:
			print i,

	return maxCLSizes_rll

'''
def processMaxCappedCLSizes(path, name, cls):
	maxCLSizes_capped_rll = []
	missRatios_capped_rll = []
	content = open(rll_cappedMaxNumOfCLSize_path + "all_capped_0_99.txt", 'r')
	capped_flag = False
	for line in content:
		if (capped_flag):
			maxCLSizes_capped_rll = eval(line)
			break
		if (name in line):
			capped_flag = True
	content.close()
	return 
'''

def processProg(name, cls, fig, fig_cnt):
	'''
	carl_64B_c_low = []
	carl_64B_mr_low = []
	carl_64B_c_upper = []
	carl_64B_mr_upper = []
	carl_64B_c_random = []
	carl_64B_mr_random = []
	if (cls == 64):
		carl_64B_c_low, carl_64B_mr_low = processCARL(carl_result_low_path, name, cls)
		carl_64B_c_upper, carl_64B_mr_upper = processCARL(carl_result_upper_path, name, cls)
		carl_64B_c_random, carl_64B_mr_random = processCARL(carl_result_random_path, name, cls)
	'''
	carl_64B_c = []
	carl_64B_mr = []
	if (cls == 64):
		carl_64B_c, carl_64B_mr = processCARL(carl_result_path, name, cls)

	lru_64B_c = []
	lru_64B_mr = []
	if (cls == 64):
		lru_64B_c, lru_64B_mr = processLRU(lru_result_path, name, cls)

	opt_64B_c = []
	opt_64B_mr = []
	if (cls == 64):
		opt_64B_c, opt_64B_mr = processOPT(opt_result_path, name, cls)

	maxCLSizes_rll = []
	if (cls == 64):
		maxCLSizes_rll = processMaxCLSizes(rll_maxNumOfCLSize_path, name, cls, carl_64B_c);

	ax = fig.add_subplot(10, 3, fig_cnt)
	ax.set_title(name, fontsize = fontSizeSet, y = 0.8, x = 0.6)
	ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0), fontsize = fontSizeSet)

	l1, l2, l3, l4, l5, l6 = None, None, None, None, None, None

	l1 = ax.plot(lru_64B_c, lru_64B_mr, color = '#5e3c99', linewidth=2.5)[0]
	
	#l3 = ax.plot(carl_64B_c_low, carl_64B_mr_low, color = '#b2abd2', linewidth=1.5)[0]
	#l4 = ax.plot(carl_64B_c_upper, carl_64B_mr_upper, color = 'g', linewidth=1.5)[0]
	#l5 = ax.plot(carl_64B_c_random, carl_64B_mr_random, color = 'b', linewidth=1.5)[0]
	l3 = ax.plot(carl_64B_c, carl_64B_mr, color = '#b2abd2', linewidth=1.5)[0]	

	l2 = ax.plot(opt_64B_c, opt_64B_mr, color = '#e66101', linewidth=1)[0]
	
	if (len(maxCLSizes_rll) == len(carl_64B_mr)):
		l4 = ax.plot(maxCLSizes_rll, carl_64B_mr, "3" ,color = '#fdb863', linewidth=2, markersize=5, label = "RL-MAX")[0]
	else:
		print "does not match", len(maxCLSizes_rll), len(carl_64B_mr)
	#for index in range(len(missRatios_rll)):
	#	ax.plot([maxCLSizes_capped_rll[index], maxCLSizes_rll[index]], [missRatios_rll[index], missRatios_rll[index]], color='sandybrown', linewidth=0.5, linestyle='-.')

	if (fig_cnt == 29):
		ax.set_xlabel('Cache Size for LRU / Averaged Cache Size for CARL (KB)', fontsize = fontSizeSet)
	if (fig_cnt == 13):
		ax.set_ylabel('Miss Ratio', fontsize = fontSizeSet)
	ax.set_xscale('symlog')
	
	ax.set_ylim(0,yMax[name])
	
	ax.yaxis.set_tick_params(labelsize=fontSizeSet)
	ax.xaxis.set_tick_params(labelsize=fontSizeSet)

	if (fig_cnt == 1):
		fig.legend([l1, l2, l3, l4], ["LRU", "OPT", "CARL", "RL-MAX"], ncol=5, mode="expand", loc = "upper center", fontsize = fontSizeSet)
		#fig.legend([l5, l6], ["LRU-512B", "CARL-512B"], ncol=2, mode="expand", loc = "upper center", fontsize = fontSizeSet)	
		#fig.legend([l2, l5, l4, l6], ["LRU-64B", "LRU-512B", "CARL-64B", "CARL-512B"], ncol=4, mode="expand", loc = "upper center", fontsize = fontSizeSet)

	fig_cnt += 1
	return fig, fig_cnt

fig = plt.figure(figsize=(15, 30))
fig_cnt = 1

cls = 64
for name in names:
	fig, fig_cnt = processProg(name, cls, fig, fig_cnt)
plt.subplots_adjust(left=0.08, right = 0.98, top = 0.96, bottom = 0.04)
#plt.show()
plt.savefig('all_yCapped.pdf')

