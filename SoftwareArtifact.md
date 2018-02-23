# Quick Start

### Install Docker
Docker can be downloaded from [Docker Installation Guide](https://docs.docker.com/docker-for-mac/install/). Currently, it supports Mac, Windows, Linux(Ubuntu, Debian, CentOS, Fedora, Binaries). After installation, you can use `docker --version` command to check whether the docker was installed correctly. A correct installation will show the version of docker on the command line, i.e. `Docker version 17.12.0-ce, build c97c6d6`
 
### Build Docker Image
This image will occupy 2.33 GB on your machine, and the building process will last for a couple minutes.
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

### Test
```bash
# go to the sps_pldi18_aec/test directory 
$ make
# or if 
```


