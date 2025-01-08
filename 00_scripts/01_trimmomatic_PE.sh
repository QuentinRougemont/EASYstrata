#!/bin/bash
#date: 14-11-22
#Author:QR
#script to run trimmomatic
source config/cpu_mem

#required arguments: name of read1.fastq.gz and read2.fastq.gz
if [ $# -ne 2 ]; then
    echo "USAGE: $0 R1.fastq.gz R2.fastq.gz file "
    echo "Expecting fastq files from read1 and read2 as input"
    echo "Extension should be '*fastq.gz'"
    exit 1
else
    file1=$1
    file2=$2 
    echo "fastq file1 is : ${file1}"
    echo "fastq file2 is : ${file2}"
    echo " "
fi

#base=$(basename $file1)
base=$(basename ${file1%[_.]*[._]f**gz} )

base=$(basename ${file1%R1.*f**gz} )
echo base is $base

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
LOG_FOLDER="log_files"

#create folder if not existent:
mkdir $LOG_FOLDER 2>/dev/null
mkdir 02_trimmed  2>/dev/null

ADAPTERFILE="Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"
NCPU="$NCPUS_TRIMMO"

for file in 02_trimmed/*gz 
do
    if [ ! -s "${file}" ] 
    then
    
    java -jar -Xmx10G Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -threads "$NCPU" \
        -phred33 \
        "$file1" \
        "$file2" \
        02_trimmed/"$base"_R1.paired.fastq.gz \
        02_trimmed/"$base"_R1.single.fastq.gz \
        02_trimmed/"$base"_R2.paired.fastq.gz \
        02_trimmed/"$base"_R2.single.fastq.gz \
        ILLUMINACLIP:"$ADAPTERFILE":2:30:10:2:True \
        HEADCROP:9 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:30:30 \
        MINLEN:36 2> $LOG_FOLDER/log.trimmomatic.pe."$base"."$TIMESTAMP"
        
    else
        echo "file $base already present"
    fi
done

