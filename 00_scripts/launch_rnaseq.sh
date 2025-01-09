#!/bin/bash 

source ../config/config
source ../config/colors

echo rnaseqlist is $RNAseqlist 
#microscript to run all rnaseq steps from read trimming to read mapping and count
#RNAseq=YES

TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)
LOG_FOLDER="LOGS"

#create folder if not existent:
mkdir $LOG_FOLDER 2>/dev/null

haplotype=$1
genome=03_genome/"$haplotype".fa

#launch gmap :
../00_scripts/02_gmap.sh "$genome" 

ncol=$(awk '{print NF}'  "$RNAseqlist" |uniq)

if [[ $ncol = 2 ]] ; then
    echo -e "\nrunning gsnap assuming reads are Paired-End\n" 

    #launch gsnap - samtools and read count:
    for read1 in ../02_trimmed/*R1.paired.fastq.gz  ; do
        [ -e "$read1" ] || continue 
        ../00_scripts/03_gsnap_PE.sh "$genome" "$read1" 2>&1 |tee "$LOG_FOLDER"/gsnap_"$(basename "$read1")"_"$TIME".log
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR GSNAP FAILED"
            exit 1
        fi

    done 

else
    #assuming SE:
    #launch gsnap - samtools and read count:
    echo -e "\nrunning gsnap assuming reads are Single-End\n"

    for read1 in ../02_trimmed/*R1.fastq.gz ; do
        [ -e "$read1" ] || continue 
        ../00_scripts/03_gsnap_SE.sh "$genome" "$read1" 2>&1 |tee "$LOG_FOLDER"/gsnap_"$(basename "$read1")"_"$TIME".log
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR GSNAP FAILED"
            exit 1
        fi

    done 

fi

rm -rf Trimmomatic-0.39* 
