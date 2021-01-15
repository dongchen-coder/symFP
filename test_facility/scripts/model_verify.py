import numpy as np
import argparse
import subprocess
import os
BENCHMARK = [
        "deriche",
        "fdtd-2d",
        "adi",
        "jacobi-1d",
        "jacobi-2d",
        "heat-3d",
        "stencil", 
        "2mm", 
        "3mm", 
        "bicg", 
        "atax", 
        "gemm",
        "gemver", 
        "mvt", 
        "gesummv"
]
'''
BENCHMARK = [
        "2mm", 
        "3mm", 
        "bicg", 
        "atax", 
        "gemm",
        "gemver", 
        "mvt", 
        "gesummv",
        "stencil"
]
BENCHMARK = [
        "deriche",
        "fdtd-2d",
        "adi",
        "jacobi-1d",
        "jacobi-2d",
        "heat-3d",
]
'''
BIN_DIR = "./bin"

def load_rt(l, delimiter=" "):
    rthist = dict()
    start_load = False
    time = 0.0
    for ln, line in enumerate(l):
        line = line.rstrip()
        if line.startswith("Overhead in"):
            time = float(list(filter(None, line.split(" ")))[-1]) 
        if line.startswith("miss ratio"):
            start_load = False
        if line.startswith("Start to dump reuse time histogram"):
            start_load = True
            continue
        if start_load:
            rt = float(list(filter(None, line.split(delimiter)))[0])
            rtcnt = float(list(filter(None, line.split(delimiter)))[1])
            rthist[rt] = rtcnt
    sumrtcnt = sum([rthist[r] for r in rthist])
    for rt in rthist:
        rthist[rt] = rthist[rt] / sumrtcnt
    print(time)
    return rthist, sumrtcnt, time


def cal_accuracy(rthist_1, rthist_2):
    for rt in rthist_1:
        if rt not in rthist_2:
            rthist_2[rt] = 0.0
    for rt in rthist_2:
        if rt not in rthist_1:
            rthist_1[rt] = 0.0
    cnt = 0.0
    for rt in rthist_1:
        cnt += abs(rthist_1[rt] - rthist_2[rt])
    return float((1 - (cnt / 2))) * 100



def main(is_all=False, prog=None, epochs=5):
    if is_all:
        for prog in BENCHMARK:
            plum_time_list = list()
            model_time_list = list()
            try:
                print("collect from " + prog)
                for i in range(epochs):
                    proc = subprocess.Popen([BIN_DIR+"/"+prog+"_plum"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    stdout, stderr = proc.communicate()
                    plum_data, plum_sum, plum_time = load_rt(list(filter(None, \
                        stdout.decode().split("\n"))), delimiter=",")
                    plum_time_list.append(plum_time)
                    proc = subprocess.Popen([BIN_DIR+"/"+prog+"_model"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    stdout, stderr = proc.communicate()
                    model_data, model_sum, model_time = \
                    load_rt(list(filter(None, stdout.decode().split("\n"))),
                        delimiter=",")
                    model_time_list.append(model_time)
                    print("%s\t%.4f%%\t[%.2f vs. %.2f]" % (prog, cal_accuracy(plum_data, model_data), plum_sum, model_sum))
                print("%s SPEEDUP \t%.2fx [%.2f vs. %.2f]" % (prog,
                    np.mean(plum_time_list) / np.mean(model_time_list),
                    np.mean(plum_time_list), np.mean(model_time_list)))
            except subprocess.CalledProcessError:
                print("Compile Error for benchmark %s" % prog)
    else:
        plum_time_list = list()
        model_time_list = list()
        try:
            print("collect from " + prog)
            for i in range(epochs):
                proc = subprocess.Popen([BIN_DIR+"/"+prog+"_plum"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                plum_data, plum_sum, plum_time = load_rt(list(filter(None, \
                    stdout.decode().split("\n"))), delimiter=",")
                plum_time_list.append(plum_time)
                proc = subprocess.Popen([BIN_DIR+"/"+prog+"_model"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                model_data, model_sum, model_time = load_rt(list(filter(None, \
                    stdout.decode().split("\n"))),
                        delimiter=",")
                model_time_list.append(model_time)
                print("%s\t%.4f%%\t[%.2f vs. %.2f]" % (prog, cal_accuracy(plum_data, model_data), plum_sum, model_sum))
            print("%s SPEEDUP \t%.2fx [%.2f vs. %.2f]" % (prog,
                    np.mean(plum_time_list) / np.mean(model_time_list),
                    np.mean(plum_time_list), np.mean(model_time_list)))
        except subprocess.CalledProcessError:
            print("Compile Error for benchmark %s" % prog)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action='store_true', help='run all benchmarks')
    parser.add_argument("--prog", default='2mm', help='set the benchmark to run')
    parser.add_argument("--epoch", type=int, default=3, help='set the times each benchmark to run')
    args = parser.parse_args()
    main(is_all=args.all, prog=args.prog, epochs=args.epoch)
