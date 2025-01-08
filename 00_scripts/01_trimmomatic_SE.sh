#!/bin/bash
source config/cpu_mem

#required arguments: name of read1.fastq.gz
if [ $# -ne 1  ]; then
    echo "USAGE: $0 R1.fastq.gz file "
    echo "Expecting a fastq file from read1 as input"
    echo "Extension should be '*fastq.gz'"
    exit 1
else
    file=$1
    echo "fastq file is : ${file}"
    echo " "
    NCPU=$2 #number of CPU (optional)
fi

#file : read1.fastq.gz located in 01_raw folder 
base=$(basename "$file")

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
LOG_FOLDER="log_files"

#create folder if not existent:
mkdir $LOG_FOLDER 2>/dev/null
mkdir 02_trimmed  2>/dev/null

ADAPTERFILE="Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"

NCPU="$NCPUS_TRIMMO"

java -jar -Xmx10G Trimmomatic-0.39/trimmomatic-0.39.jar SE \
        -threads "$NCPU" \
        -phred33 \
        "$file" \
        02_trimmed/"$base"_R1.fastq.gz \
        ILLUMINACLIP:"$ADAPTERFILE":2:30:10 \
        HEADCROP:9 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 2> $LOG_FOLDER/log.trimmomatic.pe."$base"."$TIMESTAMP"

