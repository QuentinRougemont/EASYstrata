#!/bin/bash
#author: QR
#script to run gsnap
#input : fastq and genome
#output: bamm file 
if [ $# -ne 2  ]; then
    echo "USAGE: $0 reference_genome trimmed_fastq_file"
    echo "Expecting the name of the reference genome and the name of the fastq_file (read1 only)"
    echo "Extension should be '*fastq.gz'"
    exit 1
else
    genome=$1
    fq=$2
    echo -e "reference genome is $genome \n"
    echo "fastq file is : ${fq}"
    echo " "
fi
source ../config/cpu_mem
# Global variables
DATAOUTPUT="04_mapped/"
DATAINPUT="02_trimmed"
mkdir -p $DATAOUTPUT 2>/dev/null

NCPUS=$NCPUS_GSNAP

# For genome
#check:
genome=$(basename "$genome" )
GENOMEFOLDER="03_genome"
GENOME=gmap_"${genome%.fa**}"
platform="Illumina"

input=$(basename "$fq")
base=${input%_R1.fastq.gz}

# Align reads
echo "Aligning $base"
#gsnap  -t "$NCPUS" -A sam \
gsnap --gunzip -t "$NCPUS" -A sam \
         -M 2 -n 10 -N 1 \
         -w 200000 --pairmax-rna=200000 \
         -E 1 -B 2 \
         --clip-overlap \
        --dir="$GENOMEFOLDER" -d "$GENOME" \
        --split-output="$DATAOUTPUT"/"$base" \
        --read-group-id="$base" \
        --read-group-platform="$platform" \
        "$DATAINPUT"/"$base"_R1.fastq.gz

# concatenate sam
samtools view -b "$DATAOUTPUT"/"$base".uniq >"$DATAOUTPUT"/"$base".concordant_uniq.bam

#sorting bam
echo "Creating sorted bam for $base"
samtools sort "$DATAOUTPUT"/"$base".concordant_uniq.bam -o "$DATAOUTPUT"/"$base".sorted.bam
samtools index "$DATAOUTPUT"/"$base".sorted.bam
# Clean up
echo "Removing ""$TMP""/""$base"".bam"

rm $DATAOUTPUT/"$base".concordant*
rm $DATAOUTPUT/"$base".halfmapping*
rm $DATAOUTPUT/"$base".nomapping*
#rm $DATAOUTPUT/"$base".paired*
rm $DATAOUTPUT/"$base".unpaired*

