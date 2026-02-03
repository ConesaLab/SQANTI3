### Table of contents:

* <a href="#dependencies">SQANTI3 dependencies</a>

* <a href="#install">Getting ready to use SQANTI3</a>

    * <a href="#install0">Installing and updating Anaconda</a>
    * <a href="#install1">Downloading SQANTI3</a>
    * <a href="#install2">Creating the conda environment</a>
    * <a href="#install3">Installing gtfToGenePred</a>

* <a href="#docker">Using SQANTI3 docker container</a>

  
____________________________________________________

<a id="dependencies"></a>

## Dependencies

Non-comprehensive list of the main dependencies for SQANTI3:

### General

* Perl
* Minimap2 
* Python (3.7)
* R (>= 3.4.0)
* kallisto
* samtools
* STAR
* uLTRA
* deSALT
* pip

### Python-related libraries

* bx-python
* BioPython
* BCBioGFF
* cython
* NumPy
* pysam
* pybedtools
* psutil
* pandas
* scipy

### R-related libraries

* R packages for `sqanti3_qc.py` and `sqanti3_filter.py` (installed when creating the conda environment): 

### External scripts

* We have downloaded gtfToGenePred from [UCSC utilities](https://hgdownload.soe.ucsc.edu/admin/exe/) and gave execution permissions, so it is no longer necessary to do it after downloading SQANTI3.


<a id="install"></a>

## Getting ready to use SQANTI3

We recommend using Anaconda to substantially facilitate installation of all Python dependencies. 
Probably you already have Anaconda installed because you use 
[BioConda IsoSeq(3)](https://github.com/PacificBiosciences/IsoSeq). 
Please, follow the steps here to ensure an error-free installation. 
All the dependencies will be installed automatically in a conda environment. 
The installation will be done just once and it usually takes less than 10 minutes to be installed. This approach has been tested successfully macOS and Linux-based operating systems. When the environment has been entirely built, 
you just need to activate the conda environment of SQANTI3 and run it!


<a id="install0"></a>

### 0. Installing and updating Anaconda

Make sure you have installed Anaconda. If so, you may add it to your PATH variable and update it -if necessary- as follows:

```
export PATH=$HOME/anacondaPy37/bin:$PATH
conda -V
conda update conda
```

If you have not installed Anaconda, the generic installer for Linux OS can be found [here](http://docs.continuum.io/anaconda/install/#linux-install). Note that only Linux and Mac machines are currently supported.


<a id="install1"></a>

### 1. Downloading SQANTI3

Next, download the [latest stable SQANTI3 release](https://github.com/ConesaLab/SQANTI3/releases/tag/v5.4). 
As of April 23rd 2025, the **current version is v5.4**.

For general users, we recommend downloading the SQANTI3 repository as follows:

```
wget https://github.com/ConesaLab/SQANTI3/releases/download/5.4/SQANTI3_v5.4.zip
mkdir sqanti3
unzip SQANTI3_v5.4.zip -d sqanti3
```

If you have intentions of contributing to the development of SQANTI3, please clone the developer version. 
This option will set up a git repository within your SQANTI folder and is **NOT** recommended for general users.
Contributors outside the main development team are welcome submit a **pull-request** after performing the changes in their forked repositories.

```
 git clone https://github.com/ConesaLab/SQANTI3.git
```

If you are interested in a previous version of SQANTI3, see the [version history](Version-history.md) in our wiki site, which contains a detailed account of all changes introduced for versions >=5.0.

<a id="install2"></a>

### 2. Creating the conda environment

To use SQANTI3, you will need to move into the SQANTI3 folder that you just downloaded and create a virtual 
environment including all required packages. All you need to do is run the `SQANTI3.conda_env.yml` script that 
you can find in the main SQANTI folder: 

```
conda env create -f SQANTI3.conda_env.yml
conda activate sqanti3
```

This script contains all the information required to install the SQANTI3 dependencies. 
As the environment creation progresses, you will need to type `y` when prompted to proceed with the installation. 
Note that you may change the name of the environment using the `-n` argument to the `conda env create` command. 
By default, the name of the environment will be **SQANTI3.env**.

Once you have activated the virtual environment, you should see your prompt changing to something like this:

```
(sqanti3)$
```

**Environment management after SQANTI3 updates**

If you had already installed SQANTI3 and have updated to a new version, creating a fresh environment may be required in **major and minor releases** (i.e. version bumps from 5.**1**.x to 5.**2**.x), but not in the case of **patches** (i.e. version bumps from 5.**1**.x to 5.**1**.y).
However, if you experience errors running SQANTI3 after downloading a new version, we advise creating the conda environment again to discard dependency version issues.

<a id="install3"></a>

### 3. Installing gtfToGenePred
  
For **Linux** users, the right version of gtfToGenePred is already in the `utilities` folder. However, if you are running 
SQANTI3 on **MacOS**, you must remove the `src/utilities/gtfToGenePred` file and download the Mac version [here](https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/).  

If you are not working on a Linux system, SQANTI requires manual installation of [gtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html).
This tool seems to have some issues with Python 3.7 (or openssl) when installed via conda.
At this point, the easiest solution is to download it from the [UCSC website](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) 
and add it to the `SQANTI3/src/utilities` folder you will need to make the file executable) or to your `PATH` variable. 

```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -P <path_to>/SQANTI3/utilities/
chmod +x <path_to>/SQANTI3/src/utilities/gtfToGenePred 
```


<a id="docker"></a>

## Using SQANTI3 docker container

A docker container is available from version 5.2.2 onwards in dockerhub

> docker pull anaconesalab/sqanti3

On the other hand, the dockerfile is available in the master branch of the git repository which you can use to create a Docker image of SQANTI3, not yet publicly available. The steps needed to use it are:

* **Build the image**

In the same folder as your dockerfile (preferably, but not necessarily, with only the dockerfile), run
```
docker build -t SQANTI3 -f Dockerfile .
```

* Execute the container: all three modules of SQANTI3 are hosted inside the container and can be executed:
```
docker run -it SQANTI3 sqanti3_qc.py -h 
docker run -it SQANTI3 sqanti3_filter.py -h 
docker run -it SQANTI3 sqanti3_rescue.py -h 
```

In order to make docker be aware of your files, input files and output folder must be bound to a path in your filesystem, as otherwise, docker won't be able to read them. The path /data2 inside the container is free to be mapped to your project folder, and it's the default working directory inside the container.

```
docker run -it -v <my-project-path>:/data2 SQANTI3 sqanti3_qc.py <sqanti_qc_args>
```

### Docker miniFAQ

* **I'm having trouble to build the docker image on my Mac**

Building the image on Mac systems may require adding the "--platform linux/amd64" parameter when building the image.

* **Can I use a singularity image instead of docker?** 

Once the docker image is built, you can turn it into a singularity image either from dockerhub or from the installed local version:
```
singularity build sqanti3.sif docker://anaconesalab/sqanti3
singularity build sqanti3.sif docker-daemon://sqanti3
```
However you build the singularity image, you can use "singularity run sqanti3.sif <command>" exactly as you would with the docker container. Note that using "singularity exec" will not load the conda environment by default

* **How can I make Docker generated output to be owned by my user and not root?**

This is the default behaviour of docker. The most straightforward way to solve it is by telling docker to execute the container with a specific user with -u parameter

```
docker run -it -u $(id -u):$(id -g) sqanti3 <sqanti3_qc.py/sqanti3_filter.py/sqanti3_rescue.py> 
```
This command will execute the container with your current user, and thus, the resulting files will be owned by your user and group