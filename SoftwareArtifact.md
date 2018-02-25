# What to expect

Here we provide the artifact that contains codes, testing shell scripts and plotting python scripts for Static Parallel Sampling (SPS). We extracted the minimal code needed from the newest version of the Static Parallel Sampling tools and provided it with a Dockerfile to regenerate the testing environment needed. 

Here we provide the artifact that can support the results in Fig 5-7 in the paper. (Miss ratio curve, overhead and parallel execution). Fig 5 and 6 can be re-plotted and the difference with the original submission is explained at the end. For Fig 7, we show our generated parallel sampling code as the support.

# Setup
Our artifact contains following components:
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
- **Dockerfile**: This file is used to build the docker image, this image contains ubuntu-16.04 and llvm-4.0.0. Clang-4.0.0 is not included as building Dockerfile with clang tool often fails due to limited resources.
- **CMakeList.txt**: Cmake setup for SPS which are implemented as LLVM analysis passes.
- **sps/**: The source code for SPS. It first provides one array index analysis pass (idxAnalysis.cpp, idxAnalysis.hpp) and one loop analysis pass (loopAnalysis.cpp, loopAnalysis.hpp) to generate the tree representation of the program. Then it provides one generation pass (ssCodeGen_ref.cpp, ssCodeGen_ref.hpp) to generate static sampling code from the extracted tree representation. At last, it provides one wrapper pass static sampling pass (sps.cpp) and one header file to define sampling rate and enable parallel sampling by macros.
- **test_facility/**: It contains Poly Benchmarks, fft and makefiles to assist auto scripts in test_run to generate our results. Here is the detail description for each directory:
    - **bc**: The binary code(.bc) file of all benchmarks, which will be used for SPS analysis.
    - **ir**: The IR for all benchmarks.
    - **polyBench & polyBench_trace**: The benchmark program suite for SPS and trace analysis respectively.
    - **ss_bin**: The executable program for SPS.
    - **ss_code**: The auto-generated program by SPS that do the sampling and calculate the miss ratio for each memory reference in parallel.
    - **ss_result**: The SPS analysis result.
    - **trace_bin**: The executable program for trace analysis.
    - **trace_result**: The trace analysis result.
- **test_run/**: It contains scripts for compiling and running SPS/trace analysis and python scripts that output the diagrams shown in paper.

To reproduce the results in the paper, you need to install Docker. Our virtual machine will occupy 2.33GB and the building time for our docker image will last for about an hour. 


# Setup

### Install Docker
Docker can be downloaded from [Docker Installation Guide](https://docs.docker.com/docker-for-mac/install/). Currently, it supports Mac, Windows, Linux(Ubuntu, Debian, CentOS, Fedora, Binaries). After installation, you can use `docker --version` command to check whether the docker was installed correctly. A correct installation will show the version of docker on the command line, i.e. `Docker version 17.12.0-ce, build c97c6d6`. We've tested our artifact on MacOS Serria and Ubuntu 16.04. Please set the memory constraint for docker container from 2GB to 3GB. In Mac or Windows, this can be easily done in Preferences/Advanced Tab after clicking the icon shown in status bar, in Linux, you have to do nothing.
 
### Build Docker Image
This image will occupy 2.33 GB on your machine, and the building process will last for about an hour. But this build is totally automated, you can just leave it alone and do whatever stuffs you like. 
```bash
# run this command under that path contains Dockerfile
$ docker build -t sps-image ./sps_pldi18_aec
```

### Run Docker Container
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
### Useful Docker command
Here we also provide some useful command that you may use when 'play' with our artifact. For more command, please check [Docker Command Reference](https://docs.docker.com/engine/reference/run/)
```bash
# list all docker containers
$ docker ps -a

# list all docker images
$ docker image ls

# exit from docker container
$ exit

# start an existed docker container
$ docker start -ai $(CONTAINER NAME)

# delete the docker container
$ docker container rm $(CONTAINER NAME)

# remove the docker image
$ docker rmi $(IMAGE NAME)
```
# Reproducing the result
Running the script for each test will take only a few minutes because in each script we will change the sampling ratio for each benchmark and recompile them. Overall, the full evaluation process will take no more than 20 min. The running order is not limited. Notes that all the python script should be run out of docker container.

### Precision Test (Fig. 5)
Here is the guide for running the precision test. This test result corresponds to Figure 5 in the paper. In this test, we compare the miss ratio curve after doing the SPS and trace analysis.
```bash
# go to the sps_pldi18_aec/test_run/precision directory 
$ cd /sps_pldi_aec/test_run/precision
$ sh run_ss.sh
$ sh run_trace.sh

# Then run the python code showing the analysis result diagram, which comparing the miss ratio curve between our SPS method and Trace Analysis
# To run this script, you have to exit the docker container and then go to the same directory
$ python plotMRC_staticSampling_VS_trace_cl.py
```

### Overhead Test (Fig. 6)
Here is the guide for running the overhead test. This test result corresponds to Figure 6 in the paper. In this test, we compare the running time when doing the SPS and trace analysis towards the benchmark programs.
```bash
# go to the sps_pldi18_aec/test_run/overhead directory 
$ cd /sps_pldi_aec/test_run/overhead
$ sh time_ss.sh
$ sh time_trace.sh

# Then run the python code showing the performance diagram, which comparing the runnint time between the static sampling and trace analysis.
# To run this script, you have to exit the docker container and then go to the same directory
$ python plotRT_staticSampling_VS_trace_cl.py
```

### Parallel Test (Fig. 7)
Here is the guide for running the overhead test. This test result corresponds to Figure 7 in the paper. In this test, we run the static parallel sampling analysis for all memory reference in parallel using C++ thread. The code can be checked in /test/facility/ss_code.
```bash
# go to the sps_pldi18_aec/test_run/parallel directory 
$ cd /sps_pldi_aec/test_run/parallel
$ sh time_pp.sh
```

# Explanation

This provided version is adjusted from the newest version of the tool we are developing. The implementation provided is extracted according to the experiments in the paper. 
Some of the codes are compromised due to docker vm environment:
- regix used in loopAnalysis.cpp are replace by matchCond(), matchInc() as using regix will generate error when using gcc.
- We limit cache size to at most 100,000 when calculating the miss ratio curve due to the memory limitation of the Docker container (running of the tests will be killed). 
    ```C++
    // In function RTtoMR_AET():
    for (uint64_t c = 0; c <= max_RT && c <= 100000; c++) { // original code does not have "&& c <= 100000"
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

The reason we gives the binary code(.bc) file for all benchmarks is that we doesn't install Clang in our docker container. Installing clang and llvm will consumes lots of time and space (fails a lot when building the docker file installing Clang).

The results in submitted version of the paper are performed on MacOS with clang compiler and 1.4 GHz Intel Core i5 with 4 GB 1600 MHz DDR3. The docker provides linux with g++ compiler may affect the results. Different hardware platform may affect the results.

### Fig. 5 (Precision)
When running our SPS analysis, the miss ratio curve may have a little difference than what we showed in the paper because of the following reasons:

- We limit cache size to at most 100,000 when calculating the miss ratio curve due to the memory limitation of the Docker container. 

- We use random sampling as the sampling methods. The final output may change in different environment as random() function may generate different sampling iteration points.

The newly generated result should show the same conclusion as the submitted version does.

### Fig. 6 (Overhead)

For code generation time, the difference is introduced by compiling llvm-4.0.0 using g++ in docker vm instead of using clang in MacOS. It will change the running time of opt binary and generated sps pass lib. But the order of magnitude of the time is the same.

For tracing time, the difference will be affected by different compilers used and different running environment. Limiting cache size by adding "c <= 100000" will change the performance. But the order of magnitude of the time is the same.

As we further improved the code after submission and limited the cache size by adding "c <= 100000". The overhead measured for static sampling are significantly smaller (can be barely seen from the plots). We will update the newest better results to the paper.

### Fig. 7 (Parallel)

As we improved the performance. We didn't show the parallel effect here but instead the gernerated static parallel sampling with C++ thread constructs can be checked in /test/facility/ss_code.



