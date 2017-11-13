# ATAC-seq QC Pipeline Usage

To start using the ATAC-seq QC pipeline, there are three steps described below:
1. [Install Docker](#1-install-docker)
2. [Configure Docker image](#2-configure-docker-image)
3. [Ready to run](#3-ready-to-run)

## 1. Install Docker
Docker is a great way to run our QC pipeline, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

First, install docker on your system : [Docker Installation Instructions](https://docs.docker.com/engine/installation/).

If you are working on a server of Ubuntu, first make sure the essential packages are installed:
```bash
$ sudo apt-get update

$ sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    software-properties-common
```

Then, you can download the configuration file from Docker:
```bash
$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

$ sudo add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) \
    stable"
```

Makes sure that you have the permission to run Docker:
```bash
$ sudo usermod -a -G docker your_username
```

Common commands in Docker:
```bash
# check what images you have now
$ docker images

# check the current running containers
$ docker ps

# kill current running containers and remove it
$ docker stop $(docker ps -a -q)
$ docker rm $(docker ps -a -q)

# remove unused images, make sure stop and remove all dependent containers
$ docker rmi <repository>:<tag>
```

## 2. Configure Docker image
Our images are available in [Docker Hub](https://hub.docker.com/r/zhanglab/atac-seq/). There are 3 ways to configure the required image:
#### 1) Pull from Docker Hub
```bash
$ docker pull zhanglab/atac-seq:full
```
`full` contain all available reference genomes including mm9, mm10, hg19, hg38 and danRer10, and it is so large that you make take some time to do it.

If you are only interested in mm10, you can download `mm10` instead:
```bash
$ docker pull zhanglab/atac-seq:mm10
```

#### 2) Pull `base` from Docker Hub, and build `full` by yourself
```bash
$ docker pull zhanglab/atac-seq:base

# download Dockerfile
$ mkdir docker
$ cd docker/
$ curl http://brc.wustl.edu/SPACE/chengl/Docker/Dockerfile_full -o Dockerfile

# build atac-seq:full
$ docker build -t zhanglab/atac-seq:full .

$ cd ..
$ rm -r docker/
```

#### 3）Downlaod `base` from our sever, and build `full` by youself
In case that you cannot connect to Docker Hub, or the speed is low, you can download `base` from our server:
```bash
# build atac-seq:base
$ curl -fsSL http://brc.wustl.edu/SPACE/chengl/Docker/atac-seq-base.tar.gz -o ./atac-seq-base.tar.gz
$ docker load -i atac-seq-base.tar.gz
$ rm ./atac-seq-base.tar.gz

# download Dockerfile
$ mkdir docker
$ cd docker/
$ curl http://brc.wustl.edu/SPACE/chengl/Docker/Dockerfile_full -o Dockerfile

# build atac-seq:full
$ docker build -t zhanglab/atac-seq:full .
```








