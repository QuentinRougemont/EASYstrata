#!/bin/bash
#PURPOSE: script to run braker
#AUTHOR: QR
#Date updated: 10-03-2023

#------------- EXTERNAL VARIABLE FROM CONFIG FILE -------------- #
source ../config/config
source ../config/cpu_mem
#------------- CONDA ACTIVATION  -------------- #
eval "$(conda shell.bash hook)"
conda activate superannot

#--- start of setting path ---- " 
CDB_PATH
TSEBR_PATH
PROTH_PATH
GMARK_PATH
AUGCO_PATH
AUGBI_PATH
AUGSC_PATH
#--- end of setting path ---- " 

#set -e

current_command=$BASH_COMMAND
last_command=""

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

#  ---- external data required arguments --------------- #
if (( $# < 4 )) ; then
    echo "USAGE: $0 <reference_genome> <species> <RNAseq>(YES/NO) <fungus>(YES/NO) <bamlist> (optional) <NCPUS>(optional)" 
    echo -e "Expecting the following parameters:\n
        1 - the reference genome\n
        2 - a species name\n
        3 - a string YES/NO for wether RNAseq should be considered or not\n\n
        4 - a string YES/NO stating wether data are from fungus or not \n\n" 
        exit 1
else
    genome=$1 #the reference genome here! "genome.wholemask.fa"
    species=$2
    RNAseq=$3 #YES/NO
    fungus=$4
    bamlist=$5
    NCPUS=$6
fi

NCPUS="$NCPUS_BRAKER"

TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)

#----------------- BAM data for RNAseq ------------------------ #
#previously we used a concatenated and filtered bam.
#however, providing a list of all bam gives similar results with a significant gain in time:
if [ -z "$bamlist" ] ; then
    alnBAM=$(echo 04_mapped/*sorted.bam |sed 's/ /,/g' )
else

    #assuming list of bam already exist !
    #these will be soft linked in 04_mapped
    if [ ! -d "04_mapped" ] ; then
           mkdir 04_mapped
    fi

    cd 04_mapped || exit 

    #for i in $(cat "$bamlist" ) ; do
    while read -r line 
    do 
        ln -s "$line" . 
    done < "$bamlist" 
    cd ../ ; 

    alnBAM=$(echo 04_mapped/*sorted.bam |sed 's/ /,/g' )
fi

#----------------OrthoDB and Other Protein data -------------- #
target="$orthoDBspecies"

if [ -z "${orthoDBspecies}" ]; 
then
    echo "no orthoDB species provided"
        if [ -z "${RelatedProt}" ] ; 
        then
            echo -e "\nno related protein"
            echo -e "WARNING - you have not provided any external data\n"
            echo -e "braker will try running with RNAseq only\n\n" 
        else
            echo -e "related protein is $relatProt\n\n"
            relatProt="$RelatedProt"
        fi
else
    clades=("Metazoa" 
        "Vertebrata" 
        "Viridiplantae" 
        "Arthropoda" 
        "Eukaryota" 
        "Fungi" 
        "Alveolata" 
        "Stramenopiles")
    if [[ ${clades[*]} =~ $target ]]
    then
        if [ -d odb11 ] ; then
            rm -rf odb11 
            mkdir odb11 
        else 
           mkdir odb11
        fi

        cd odb11 || exit
        for file in "$target".fa.*
        do 
            if [ -f "$file" ] 
            then
                echo "warning file $target.fa already present "
                echo "please verify if this is the file that you need"
        #       exit 1 
            else
                wget -q https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/"${target}".fa.gz
                gunzip "${target}".fa.gz
                cd ../ 
                if [ -z ${RelatedProt+x} ] ; then
                    echo "no related protein"
                    cp odb11/"${target}".fa  relatProt.fa
                else
                    cat "$RelatedProt"  odb11/"${target}".fa > relatProt.fa
                    relatProt="relatProt.fa"
                fi
            fi
        done 
    else
        echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        echo -e "error the clade name you provided is not in the orthoDB list !!\n"
        echo -e "please check the clade name"
    	echo -e "this should be one among (casee sensitive):
	    [Metazoa Vertebrata Viridiplantae Arthropoda Eukaryota Fungi Alveolata]"
        echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        exit 1
    fi
fi

## --------- step 1 : BRAKER WITH RNA SEQ  ---------  ##
## Note on braker: we found that running on RNAseq & proteins separately and combining afterwards 
## manually with TSEBRA give higher BUSCO scores. 
## in addition manually setting parameter for TSEBRA config file may prevent the removal of 
## biologically important genes, a case we have observed on all out dataset.
## Therefore we choose to only implement this approach
## for user having trouble installing braker and all it's dependencies it should still be possible
## to modify the code lightly and run through singularity.

#Run braker with RNAseq 

output="braker.gtf"

if [[ $RNAseq = "YES" ]]
then
    wd=06_braker/rnaseq
    if [ -f "$wd"/"$output" ]
    then
        echo "file $output RNAseq already exist will skip the run"
    else
        mkdir -p $wd

        if [[ $fungus = "YES" ]]
        then
            echo -e "------ \n running braker on rnaseq data \n -------"
            echo -e "------ \n data are from fungus \n -------"
                braker.pl --species="$species"_"$TIME"_rnaseq --species="$species"_"$TIME"_rnaseq --fungus \
                    --genome="$genome" --threads="$NCPUS"  --softmasking --bam="$alnBAM" --workingdir=$wd 
        else
            echo -e "------ \n running braker on rnaseq data \n -------"
                braker.pl --species="$species"_"$TIME"_rnaseq --species="$species"_"$TIME"_rnaseq\
                    --genome="$genome" --threads="$NCPUS"  --softmasking --bam="$alnBAM" --workingdir=$wd 
        fi
    fi 
fi

##  --------- step 2 : BRAKER WITH REFERENCE DATABASE USING THREE ROUNDS --------- ## 
if [[ -z "$relatProt" ]]
then
    echo "no related protein - only RNAseq were used"
    exit
fi

#prepare architecture:
FOLDER1=06_braker/round1_braker_on_refprot #_$TIME
FOLDER2=06_braker/round2_braker_on_refprot #_$TIME
FOLDER3=06_braker/round3_braker_on_refprot #_$TIME 
FOLDER4=06_braker/round4_braker_on_refprot #_$TIME
FOLDER5=06_braker/round5_braker_on_refprot #_$TIME

mkdir -p $FOLDER1 $FOLDER2 $FOLDER3 $FOLDER4 $FOLDER5 2>/dev/null

echo -e "\n\n----------- round 1 ------------\n\n" 
echo AUGUSTUS_SCRIPTS_PATH is "$AUGUSTUS_SCRIPTS_PATH" 
echo AUGUSTUS_BINS_PATH is "$AUGUSTUS_BIN_PATH"
echo AUGUSTUS_CONFIG_PATH is "$AUGUSTUS_CONFIG_PATH" 


wd=${FOLDER1}

if [ -f "$wd"/"$output" ]
then
    echo "file $output round 1 already exist will skip the run"
else
   rm -rf "${wd:?}"/*
   if [[ $fungus = "YES" ]]
   then
       braker.pl --species="$species"_"$TIME"_round1  --genome="$genome" --threads="$NCPUS" \
           --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus  
   else
       braker.pl --species="$species"_"$TIME"_round1  --genome="$genome" --threads="$NCPUS" \
           --softmasking --prot_seq=$relatProt --workingdir=$wd 
   fi
fi

echo "----------- round 2 ------------" 
wd=${FOLDER2}

if [ -f "$wd"/"$output" ]
then
    echo "file $output round 1 already exist will skip the run"
else
    rm -rf "${wd:?}"/*
    if [[ $fungus = "YES" ]]
    then
        braker.pl --species="$species"_"$TIME"_round2 --genome="$genome" --threads="$NCPUS" \
            --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus #--hints=${FOLDER1}/hintsfile.gff 
    else
        braker.pl --species="$species"_"$TIME"_round2 --genome="$genome" --threads="$NCPUS" \
            --softmasking --prot_seq=$relatProt --workingdir=$wd #--hints=${FOLDER1}/hintsfile.gff 
    fi
fi

echo "----------- round 3 ------------" 
wd=${FOLDER3}

if [ -f "$wd"/"$output" ]
then
    echo "file $output round 1 already exist will skip the run"
else
    rm -rf "${wd:?}"/*
    if [[ $fungus = "YES" ]]
    then
        braker.pl --species="$species"_"$TIME"_round3 --genome="$genome" --threads="$NCPUS" \
            --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus #--hints=${FOLDER2}/hintsfile.gff 
    else
        braker.pl --species="$species"_"$TIME"_round3 --genome="$genome" --threads="$NCPUS" \
            --softmasking --prot_seq=$relatProt --workingdir=$wd #--hints=${FOLDER2}/hintsfile.gff 
    fi
fi

echo "----------- round 4 ------------" 
wd=${FOLDER4}
if [ -f "$wd"/"$output" ]
then
    echo "file $output round 1 already exist will skip the run"
else
    rm -rf "${wd:?}"/*
    if [[ $fungus = "YES" ]]
    then
        braker.pl --species="$species"_"$TIME"_round4 --genome="$genome" --threads="$NCPUS" \
            --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus #--hints=${FOLDER3}/hintsfile.gff 
    else
        braker.pl --species="$species"_"$TIME"_round4 --genome="$genome" --threads="$NCPUS" \
            --softmasking --prot_seq=$relatProt --workingdir=$wd #--hints=${FOLDER3}/hintsfile.gff 
    fi
fi

echo "----------- round 5 ------------" 
wd=${FOLDER5}

if [ -f "$wd"/"$output" ]
then
    echo "file $output round 1 already exist will skip the run"
else
    rm -rf "${wd:?}"/*
    if [[ $fungus = "YES" ]]
    then
        braker.pl --species="$species"_"$TIME"_round5 --genome="$genome" --threads="$NCPUS" \
            --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus #--hints=${FOLDER4}/hintsfile.gff 
    else
        braker.pl --species="$species"_"$TIME"_round5 --genome="$genome" --threads="$NCPUS" \
            --softmasking --prot_seq=$relatProt --workingdir=$wd #--hints=${FOLDER4}/hintsfile.gff 
    fi
fi
echo -e "\n${BLU}-----------------------------\n
    \t ALL BRAKER REPLICATES finished\n
    ------------------------------------${NC}\n\n" 

