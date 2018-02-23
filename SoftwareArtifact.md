
# File Descriptin
This repository contains following files or directories:
```
|---CMakeList.txt
|---Dockerfile
|---sps
|---test_facility
    |---bc
    |---ir
    |---polyBench
    |---polyBench_trace
    |---ss_bin
    |---ss_code
    |---ss_result
    |---trace_bin
    |---trace_result
|---test_run
    |---overhead
    |---parallel
    |---precision
```
#### Dockerfile
This file was used to build the docker image, this image contains ubuntu-16.04 and llvm-4.0.0.
#### CMakeList.txt
This file was used to configure the running environment of all benchmarks including generating the Makefile for SPS compiler passes.
#### sps/
This directory contains all the source code for SPS compiler passes.
#### test_facility/
This directory contains Poly Benchmark suite(for SPS and Trace analysis), the binary code file and IR file of these benchmarks, the source code and executable file for newly-generated static sampling code and trace analysis code,  and the test result after static sampling methods and trace analysis methods respectivley. 
#### test_run/
This directory contains all the script for compiling and running SPS/trace analysis or drawing figures, together wtih all test results.

# Quick Start

### Install Docker
Docker can be downloaded from [Docker Installation Guide](https://docs.docker.com/docker-for-mac/install/). Currently, it supports Mac, Windows, Linux(Ubuntu, Debian, CentOS, Fedora, Binaries). After installation, you can use `docker --version` command to check whether the docker was installed correctly. A correct installation will show the version of docker on the command line, i.e. `Docker version 17.12.0-ce, build c97c6d6`
 
### Build Docker Image
This image will occupy 2.33 GB on your machine, and the building process will last for about 30-40 minutes.
```bash
# run this command under that path contains Dockerfile
$ docker build -t sps-image ./sps_pldi18_aec
```

### Run Docker Container
Before running the Docker container, please set the memory contraint for docker container from 2GB to 3GB. In Mac or Windows, this can be easily done in Preferences/Advanced Tab after clicking the icon shown in status bar.

To launch the docker container, typing the following command on your command-line tools. 
```bash
# under the directory that contains sps_pldi18_aec/
$ docker run -it --memory-swap -1 -v $PWD/sps_pldi18_aec:/sps_pldi18_aec --name sps sps-image /bin/bash
```
After launching the docker container, you can see something like `#root@2bd62a73d523` on your command line. Notes that the string after `#root@` is a sequence of random number/letter, it may not be the same as the example we give here. 

### Compile LLVM compiler pass
```bash
# go to the sps_pldi18_aec/ directory
$ cd  sps_pldi18_aec

# create a build directory
$ mkdir build

# configure the compiler passes in build/ directory using cmake
$ cmake ..

# compile the compiler passes
$ make
```

### Precision Test (Fig. 5)
Here is the guide for running the precision test. This test result corresponds to Figure 5 in the paper. In this test, we compare the miss ratio curve after doing the SPS and trace analysis.
```bash
# go to the sps_pldi18_aec/test_run/ directory 
# run the script for static parallel sampling locality analysis and trace analysis. Notes that the running time for this two script will last around XXX minutes.
$ sh run_ss.sh
$ sh run_trace.sh

# Then run the python code showing the analysis result diagram, which comparing the miss ratio curve between our SPS method and Trace Analysis

# Notes that this code should be run out of Docker container.
$ python plotMRC_staticSampling_VS_trace_cl.py
```

### Overhead Test (Fig. 6)
Here is the guide for running the overhead test. This test result corresponds to Figure 6 in the paper. In this test, we compare the running time when doing the SPS and trace analysis towards the benchmark programs.
```bash
# go to the sps_pldi18_aec/test_run/ directory 
```

### Parallel Test (Fig. 7)
Here is the guide for running the overhead test. This test result corresponds to Figure 7 in the paper. In this test, we compare the running time when doing the SPS and trace analysis towards the benchmark programs.
```bash
# go to the sps_pldi18_aec/test_run/ directory 
```

# Explanation

This provided version is adjusted from the newest version of the tool we are developping. The implementation provided is extracted according to the experiments in the paper and compromised due to docker vm environment (regix).

The reason we gives the binary code(.bc) file for all benchmarks is that we doesn't install Clang in our docker container. Installing clang and llvm will consumes lots of time and space (fails a lot when building the docker file installing Clang).

When running our SPS analysis, the miss ratio curve may have a bit difference than what we showed in the paper because of the following reasons:
- We limit cache size to at most 100,000 when calculating the miss ratio curve due to the memory limitation of the Docker container. This modification will XXX the overhead. Here is the function we modified later(the outer for-loop).
```C++
for (uint64_t c = 0; c <= max_RT && c <= 100000; c++) {
        while (sum_P < c && t <= max_RT) {
            if (P.find(t) != P.end()) {
                sum_P += P[t];
                prev_t = t;
            } else {
                sum_P += P[prev_t];
            }
            t++;
        }
        MR[c] = P[prev_t];
}
```
- We use random sampling as one of the sampling methods. This randomness leads to 


