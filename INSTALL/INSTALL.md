# Full automated installation 


### Mamba/Conda

If you want to avoid potential conflicting versions or do not have root access on your device, you can use **conda** or **mamba** to install dependencies.

I recommend **mamba** for linux:

```
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge-pypy3-Linux-x86_64.sh
bash Miniforge-pypy3-Linux-x86_64.sh
#note: see here for other verions of mamba: https://github.com/conda-forge/miniforge
```

### Cloning the git

To get the workflow on your device:

```bash
git clone https://github.com/QuentinRougemont/EASYstrata
```

then create the following environnement:
```sh
cd EASYstrata 
mamba env create -f INSTALL/annotation_env.yml   

#to install to a specific directory:
#use mamba env create -p your_path/ -f annotation_env.yml 

#for busco:
mamba env create -f INSTALL/busco_env.yml

#for repeatmodeller:
mamba env create -f INSTALL/repeatmodeler.yml

#Note: I often observed bugs when installing all the software above in a single environnement

```
if you prefer, you can install all dependencies one by one, this takes more time but may help solving bugs. See [here](https://github.com/QuentinRougemont/EASYstrata/blob/main/.infos/install.md) to have a list of all dependencies.


### Non-conda dependencies

a number of software are not on conda (or not up-to-date) and need to be installed manually

here is an exemple procedure : 

## attempt to automatically install all softwares in one go: 

```bash
#all at once: 
conda activate superannot #wil activate some uefull dependancies

#install GENESPACE:
Rscript -e  'devtools::install_github("jtlovell/GENESPACE")'

#installing all other dependencies:
./INSTALL/dependencies.sh 2>&1 |tee log_installation.txt 

WARNING:
Augustus and Bamtools are the two tools that are most likely to fail 
please check that it works and that all path were correctly appended to your repository
```


all PATH to the new softwares should have been automatically append to your .bashrc



# ALTERNATIVELY: install each non-conda software one by one: 

below are attempts to install braker and a number of compulsory dependencies in a more or less automated fashion:

first activate the braker_environnment above to get basic dependencies 

```
conda activate superannot #wil activate some uefull dependancies

mkdir softs ; cd softs
```

then insall all tools one by one: 

### BRAKER 

```
command='braker.pl'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    git clone https://github.com/Gaius-Augustus/BRAKER
    cd BRAKER/scripts 
    chmod a+x *.pl *.py
    path=$(pwd)
    echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../../
fi
```

### ProtHint 

```
command='prothint.py'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    wget https://github.com/gatech-genemark/ProtHint/releases/download/v2.6.0/ProtHint-2.6.0.tar.gz 
    tar zxvf ProtHint-2.6.0.tar.gz
    cd ProtHint-2.6.0/bin 
    #then add to ~/.bashrc
    protpath=$(pwd)
    echo -e "\n#Path to $command\nexport PROTHINT_PATH=:$protpath" >> ~/.bashrc 
    echo -e "\nexport PATH=\$PATH:$protpath" >> ~/.bashrc
    source ~/.bashrc  
    cd ../../
fi
```

###  Diamond

```
command='diamond'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    mkdir diamond ; cd diamond
    wget https://github.com/bbuchfink/diamond/releases/download/v2.1.1/diamond-linux64.tar.gz
    tar zxvf diamond-linux64.tar.gz
    #then add to ~/.bashrc
    path=$(pwd)
    echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../ 

fi
```

###  Cdbfasta

```
command='cdbfasta'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    #direct install: 
    git clone https://github.com/gpertea/cdbfasta.git
    cd cdbfasta
    make all 
    if [ $? -eq 0 ]; then
        echo $command installation worked successfully
        cdbpath=$(pwd)
        echo -e "\n#Path to $command\nexport CDBTOOLS_PATH=$cdbpath" >> ~/.bashrc 
        echo -e "\nexport PATH=\$PATH:$cdbpath" >> ~/.bashrc
        source ~/.bashrc  
	    cd ../
    fi
fi
```

###  bamtools

```
command='bamtools'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    #direct install: 
    #bamtools : https://github.com/pezmaster31/bamtools
    git clone https://github.com/pezmaster31/bamtools
    cd bamtools
    mkdir build
    cd build
    path=$(pwd)
    cmake -DCMAKE_INSTALL_PREFIX=$(pwd) ..
    make -j8 
    make install
    if [ $? -eq 0 ]; then
        echo $command installation worked successfully
	cd $(find . -name "include" )
	bt_inc=$(pwd)
	cd ../../
	cd $(find . -name "lib*" |head -n 1 ) #can be lib/lib64
	bt_lib=$(pwd)

	cd ../src
	#$(find . -name "bamtools" )
        path=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
        source ~/.bashrc  
        cd ../../../
    else
       echo installation failed\nmake sur to have git, make and cmake
       exit 1
    fi
fi

```

note1: bamtools path need to be stored and declared in the common.mk file for augustus.  

This is what I attempt to do in the above two variable **$bt_inc** and **$bc_lib**   

Depending on your cluster this may fail...  

note2: many other dependencies may be required for you to have bamtools to run properly.  
in case of bug see details [here](https://github.com/pezmaster31/bamtools/wiki) 

### htslib 

```
command='htsfile'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    #download the latest htslib to recover the path for Augustus::
    wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2
    bzip2 -d htslib-1.18.tar.bz2
    tar xf htslib-1.18.tar
    cd htslib-1.18/

    ./configure --prefix=$(pwd)
    make
    make install
    htpath=$(pwd)

    cd ../
else
    htcmd=$(command -v "$command")
    htpath=$(echo $htcmd |sed 's/bin\/htsfile/include\/htslib/')
fi
```

### Augustus

```
command='augustus'
if command -v $command &> /dev/null && $commnand --version 2>&1 |grep -q "3.5.0" ; 
then 
   echo  "augustus v3.5.0" installed
else
    echo "$command could not be found"
    echo "will attempt a manual install" 
    echo "you may encounter several error and have to comment/uncomment or change path in "common.mk" especially without root privilege"
    git clone https://github.com/Gaius-Augustus/Augustus.git
    cd Augustus
    sed -i.bkp 's/COMPGENEPRED = true/COMPGENEPRED = false/g' common.mk
    sed -i 's/ZIPINPUT = true/ZIPINPUT = false\nBOOST = false/g' common.mk
    #bt=$(command -v bamtools)
    sed -i 's/#INCLUDE_PATH_BAMTOOLS/INCLUDE_PATH_BAMTOOLS/g' common.mk
    sed -i 's/#LIBRARY_PATH_BAMTOOLS/LIBRARY_PATH_BAMTOOLS/g' common.mk
    sed -i 's/#INCLUDE_PATH_HTSLIB/INCLUDE_PATH_HTSLIB/g' common.mk
    sed -i 's/#LIBRARY_PATH_HTSLIB/LIBRARY_PATH_HTSLIB/g' common.mk
    sed -i 's|usr/include/bamtools|baminclude|g' common.mk
    sed -i "s#baminclude#$bt_inc/bamtools#" common.mk
    sed -i "33s#usr/lib/x86_64-linux-gnu#$bt_lib#g" common.mk

    #same with HTSLIB:
    #sed -i "s#usr/include/htslib#$htpath/include/htslib#g" common.mk
    sed -i "s#usr/include/htslib#$htpath#g" common.mk

    n=$(grep -n "LIBRARY_PATH_HTSLIB" common.mk |awk -F":" '{print $1}' )
    sed -i  "${n}s#usr/lib/x86_64-linux-gnu#$htpath/lib#g" common.mk
    #attempt to make augustus:
    make augustus
    #if all was succesffull: ""
    if [ $? -eq 0 ]; then
        echo $command installation worked successfully
        augustuspath=$(pwd)
	echo -e "#path to AUGUSTUS:" >> ~/.bashrc 
        echo -e "export AUGUSTUS_CONFIG_PATH=$augustuspath/config " >> ~/.bashrc
        echo -e "export AUGUSTUS_BIN_PATH=$augustuspath/bin/ " >> ~/.bashrc
	echo -e "export PATH=\$PATH:$augustuspath/bin" >> ~/.bashrc
        echo -e "export AUGUSTUS_SCRIPTS_PATH=/$augustuspath/augustus_scripts " >> ~/.bashrc
	source ~/.bashrc
	cd ./auxprogs/joingenes
	make -j8
        cd ../bam2hints
	make
	cd ../bam2wig
	make
	cd ../filterBam
	make
	cd ../../../
    else
	    echo "installation of augustus failed, look at the logs ! "
	    echo "see details here: https://github.com/Gaius-Augustus/Augustus/blob/master/docs/INSTALL.md"
    fi
fi
```

### genemark 

```
command="gms2hints.pl" 
if ! command -v $command &> /dev/null
    then
    echo "$command could not be found"
    echo "will try a manual installation through git"
    git clone https://github.com/gatech-genemark/GeneMark-ETP/
    cd GeneMark-ETP/
    cd bin/gmes
    gmarkpath=$(pwd)
    echo -e "export GENEMARK_PATH=$gmarkpath/ " >> ~/.bashrc
    echo -e "\nexport PATH=\$PATH:$gmarkpath" >> ~/.bashrc
    source ~/.bashrc  
    cd ../../../ 
fi
```


### TSEBRA available [here](https://github.com/Gaius-Augustus/TSEBRA)

```
## direct install: 
command='tsebra.py'
if ! command -v $command &> /dev/null
    then
    echo "$command could not be found"
    echo "will try a manual installation through git"
    git clone https://github.com/Gaius-Augustus/TSEBRA 
    cd TSEBRA/bin
    tsebrapath=$(pwd)
    echo -e "\n#Path to $command\nexport TSEBRA_PATH=$tsebrapath" >> ~/.bashrc 
    echo -e "\n#Path to $command\nexport PATH=\$PATH:$tsebrapath" >> ~/.bashrc 

    source ~/.bashrc  
    cd ../../
fi 
```

### GSNAP for read mappig available [here](http://research-pub.gene.com/gmap/)

```
#direct install :
command='gsnap'
if ! command -v $command &> /dev/null
    then
    echo "$command could not be found"
    echo "will try a manual installation through git"
    wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2023-10-10.v2.tar.gz
    tar zxf gmap-gsnap-2023-10-10.v2.tar.gz
    cd gmap-2023-10-10/
    path=$(pwd)
    ./configure --prefix=$path  
    make -j8
    make install
    if [ $? -eq 0 ]; then
        echo $command installation worked successfully
	cd bin/
        path=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
        source ~/.bashrc  
        cd ../../
fi
```

## gffread  to extract CDS from fasta [see](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)

```
command='gffread'
if ! command -v $command &> /dev/null
    then
    echo "$command could not be found"
    echo "will try a manual installation through git"
    git clone https://github.com/gpertea/gffread
    cd gffread
    make release
    #if command was successfull then add to path:
    if [ $? -eq 0 ]; then
        echo $command installation worked successfully
        path=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
        source ~/.bashrc  
        cd ../
fi
```

## OrthoFinder

```
command='orthofinder'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through wget"
    wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.5/OrthoFinder.tar.gz
    tar zxf OrthoFinder.tar.gz
    cd OrthoFinder
    #if command was successfull then add to path:
    path=$(pwd)
    echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    
    #also add diamond, fastme and mcl which are present within orthofinder:
    cd bin/
    path=$(pwd)
    echo -e "\n#Path to diamond\n export PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../../
fi
```

## MCScanX : 

```
command='MCScanX'
if ! command -v $command &> /dev/null
then
   #direct install: 
   git clone https://github.com/wyp1125/MCScanX
   cd MCScanX ; make -j 8
   #if command was successfull then add to path:
   if [ $? -eq 0 ]; then
        echo $command installation worked successfully
	MCScanpath=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$MCScanpath" >> ~/.bashrc 
        source ~/.bashrc  
        cd ../

    else
       echo installation failed\nmake sur to have make and git
    fi
fi
```

## muscle
 
```
command='muscle'
if ! command -v $command &> /dev/null
then
   #direct install: 
   mkdir muscle ; cd muscle
   #wget https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.linux_intel64
   #ln -s muscle5.1.linux_intel64 muscle
   #if using muscle5 then modify the -in option into -align in the perl code of translocator
   wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
   tar zxf muscle3.8.31_i86linux64.tar.gz
   ln -s muscle3.8.31_i86linux64 muscle
   chmod +x muscle
   path=$(pwd)
   echo -e "\n#Path to $command\n export PATH=\$PATH:$path" >> ~/.bashrc 
   source ~/.bashrc  
   cd ../
fi
```

##  yn00/paml 

```
command='yn00'
if ! command -v $command &> /dev/null
then
   #direct install: 
   wget https://github.com/abacus-gene/paml/releases/download/4.10.7/paml-4.10.7-linux-X86_64.tgz
   tar zxvf paml-4.10.7-linux-X86_64.tgz
   cd paml-4.10.7/bin
   path=$(pwd)
   echo -e "\n#Path to $command\n export PATH=\$PATH:$path" >> ~/.bashrc 
   source ~/.bashrc  
   cd ../../
fi
```


## translatorx_vLocal.pl 

```
command='translatorx_vLocal.pl'
if ! command -v $command &> /dev/null
then
   #direct install: 
    mkdir translatorx ; cd translatorx
    wget http://161.111.160.230/cgi-bin/translatorx_vLocal.pl
    chmod +x translatorx_vLocal.pl
    path=$(pwd)
    echo -e "\n#Path to translatorx\n export PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../
fi
```


## subread - featureCounts: 
```
command='featureCounts' 
if ! command -v $command &> /dev/null
then
    #direct install: 
    wget https://sourceforge.net/projects/subread/files/subread-2.0.8/subread-2.0.8-Linux-x86_64.tar.gz/download -O subread-2.0.8-Linux-x86_64.tar.gz
    tar zxvf subread-2.0.8-Linux-x86_64.tar.gz 
    cd subread-2.0.8-Linux-x86_64/bin
    path=$(pwd)
    echo -e "\n#Path to translatorx\n export PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../../
fi
```

## bedtools :

```
command='bedtools'
if ! command -v $command &> /dev/null
then
    #direct install: 
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
    mv bedtools.static bedtools
    chmod a+x bedtools
    path=$(pwd)
    echo -e "\n#Path to bedtools  export PATH=\$PATH:$path" >> ~/.bashrc
    source ~/.bashrc  
fi
```
