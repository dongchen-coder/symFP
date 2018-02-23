
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
This directory contains all the compiler pass for this paper
#### test_facility/
This directory contains Poly Benchmark suite, the binary code file and IR file of these benchmarks after compiling, the source code and executable file for newly-generated static sampling code and trace analysis code,  and the test result after static sampling methods and trace analysis methods respectivley. 
#### test_run/
This directory contains all the 

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
To launch the docker container, typing the following command on your command-line tools. 
```bash
# under the directory that contains sps_pldi18_aec/
$ docker run -it -v $PWD/sps_pldi18_aec:/sps_pldi18_aec --name sps sps-image /bin/bash
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

### Test
```bash
# go to the sps_pldi18_aec/test_run/ directory 
# run the script for static parallel sampling locality analysis and trace analysis. Notes that the running time for this two script will last around XXX minutes.
$ sh run_ss.sh
$ sh run_trace.sh

# Then run the python code showing the analysis result diagram, which comparing the miss ratio curve between our SPS method and Trace Analysis
$ python plotMRC_staticSampling_VS_trace_cl.py
```

# Explanation
When running our SPS analysis, the miss ratio curve may be a bit different than what we showed in the paper for the following reasons:
- We limit the miss ratio calculation due to the memory limitation of the Docker container.
- We use random sampling as one of the sampling methods.


