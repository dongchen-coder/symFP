
## Install Docker
Docker can be downloaded from [Docker Installation Guid](https://docs.docker.com/docker-for-mac/install/). Currently, it supports Mac, Windows, Linux(Ubuntu, Debian, CentOS, Fedora, Binaries).
 
## Build Docker Image
This image requires about 2.33 GB in your machine.
```bash
$ docker build -t sps-image ./sps_pldi18_aec
```

## Run Docker Container
```bash
$ docker run -it -v $PWD/sps_pldi18_aec:/sps_pldi18_aec --name sps sps-image /bin/bash
```

## Compile and Test
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

