# ATAC-seq QC Pipeline Usage

To start using the ATAC-seq QC pipeline, there are three steps described below:
1. [Install Docker](#1-install-docker)
2. [Configure Docker image](#2-configure-docker-image)
3. [Ready to run](#3-ready-to-run)

## 1. Install Docker
Docker is a great way to run our QC pipeline, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

First, install docker on your system : [Docker Installation Instructions](https://docs.docker.com/engine/installation/).

If you are working on a server of Ubuntu, use the script provided by Docker to simplify the installation:
```bash
curl -fsSL get.docker.com -o get-docker.sh
sudo sh get-docker.sh --mirror Aliyun
rm get-docker.sh
```

If you are working on Ubuntu, first make sure the essential packages are installed:
```bash
sudo apt-get update

sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    software-properties-common
```

Then, you can download the configuration file from Docker:
```bash
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

sudo add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) \
    stable"
sudo apt install docker.io
```

Makes sure that you have the permission to run Docker:
```bash
sudo usermod -a -G docker your_username
```

Common commands in Docker:
```bash
# check what images you have now
docker images

# check the current running containers
docker ps

# kill current running containers and remove it
docker stop $(docker ps -a -q)
docker rm $(docker ps -a -q)

# remove unused images, make sure stop and remove all dependent containers
docker rmi <repository>:<tag>
```

## 2. Configure Docker image
Our images are available in [Docker Hub](https://hub.docker.com/r/zhanglab/atac-seq/). There are 3 ways to configure the required image:
#### 1) Pull from Docker Hub
```bash
docker pull zhanglab/atac-seq:full
```
`full` contain all available reference genomes including mm9, mm10, hg19, hg38 and danRer10, and it is so large that you make take some time to do it.

If you are only interested in mm10, you can download `mm10` instead:
```bash
docker pull zhanglab/atac-seq:mm10
```

#### 2) Pull `base` from Docker Hub, and build `full` by yourself
```bash
docker pull zhanglab/atac-seq:base

# download Dockerfile
mkdir docker
cd docker/
curl http://brc.wustl.edu/SPACE/chengl/Docker/Dockerfile_full -o Dockerfile

# build atac-seq:full
docker build -t zhanglab/atac-seq:full .

cd ..
rm -r docker/
```

#### 3）Downlaod `base` from our sever, and build `full` by youself
In case that you cannot connect to Docker Hub, or the speed is low, you can download `base` from our server:
```bash
# build atac-seq:base
curl -fsSL http://brc.wustl.edu/SPACE/chengl/Docker/atac-seq-base.tar.gz -o ./atac-seq-base.tar.gz
docker load -i atac-seq-base.tar.gz
rm ./atac-seq-base.tar.gz

# download Dockerfile
mkdir docker
cd docker/
curl http://brc.wustl.edu/SPACE/chengl/Docker/Dockerfile_full -o Dockerfile

# build atac-seq:full
docker build -t zhanglab/atac-seq:full .

cd ..
rm -r docker/
```

If you want to have `mm10` instead, simply replace `Dockerfile_full` by `Dockerfile_mm10`.

## 3. Ready to run

```bash
# create a container called "container_name", as what ever you want, running in the background
docker run --name container_name -d -it zhanglab/atac-seq:full

# copy input files into container atac
docker cp input_file container_name:/data/input_file
# if it is paired-end data
docker cp input_file2 container_name:/data/input_file2

# run our QC pipeline
docker exec -i container_name bash /atac_seq/pipe_code/atac_pipe_v1.sh \
    -o input_file -p input_file2 \
    -g input_genome -r input_type
```
Options of our pipe:\
`-o` and `-p`: specify the input files, accepting format of `.sra`, `.fastq` and `.fastq.gz`, use `-p` for the paired-end data.\
`-g`: specify the reference genome, including `mm9`, `mm10`, `hg19`, `hg38` and `danRer10`.\
`-r`: specify the type of input files, including `PE` and `SE`.\
`-t`: specify the number of threads to be used, defaulted as `24`.\
`-m`: specify the marker, and can be named by user.\
`-h`: for help.

```bash
# copy results from container to current working directory
docker cp container_name:/data/ data

# remember to stop and remove background container
docker stop container_name
docker rm container_name
```

What if you want to run a bunch of files? Here is an example of how you can do it:
```bash
date
# create a container
docker run --name container_name -d -it zhanglab/atac-seq:full

# run loop in the current directory which contains all the files to be run
for file in `ls *_1.fastq.gz`
do
    name=`echo ${file%_1.fastq.gz}`
    echo "processing $name"
    docker cp $name'_1.fastq.gz' container_name:/data/$name'_1.fastq.gz'
    docker cp $name'_2.fastq.gz' container_name:/data/$name'_2.fastq.gz'

    docker exec -i container_name bash /atac_seq/pipe_code/atac_pipe_v1.sh \
       -o input_file -p input_file2 \
       -g input_genome -r input_type
done

docker cp container_name:/data/  data
docker stop container_name
docker rm container_name
```
Now you will have a folder called `data/` which contains several sub-folders of processed files.

That's it! Now you are ready to use our pipeline. Have fun!
