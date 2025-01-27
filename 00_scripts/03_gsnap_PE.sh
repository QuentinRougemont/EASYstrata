#!/bin/bash
source ../config/cpu_mem
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

# Global variables

DATAOUTPUT="04_mapped/"
DATAINPUT="../02_trimmed"
mkdir -p "$DATAOUTPUT" 2>/dev/null

NCPUS=$NCPUS_GSNAP


# For genome
#check:
genome=$(basename "$genome" )
GENOMEFOLDER="03_genome"
GENOME=gmap_"${genome%.fa**}"
platform="Illumina"

input=$(basename "$fq")
base=${input%_R1.paired.fastq.gz}

minsize=150999000 #any bam smaller than this (150 mb) will be ignored and recreated
bamfile="$DATAOUTPUT"/"$base".sorted.bam
if [ -s $bamfile ]
then
    filesize=$(wc -c <$bamfile )
else
    filesize=0
fi
echo filesize is "$filesize"

if [ "$filesize" -lt "$minsize" ]
then
    echo "running gsnap"

     # Align reads
     echo "Aligning $base"
     gsnap --gunzip -t "$NCPUS" -A sam \
         -M 2 -n 10 -N 1 \
         -w 200000 --pairmax-rna=200000 \
         -E 1 -B 2 \
         --clip-overlap \
         --dir="$GENOMEFOLDER" -d "$GENOME" \
         --split-output="$DATAOUTPUT"/"$base" \
         --read-group-id="$base" \
         --read-group-platform="$platform" \
         "$DATAINPUT"/"$base"_R1.paired.fastq.gz "$DATAINPUT"/"$base"_R2.paired.fastq.gz 
     
     
     # concatenate sam
         samtools view -b "$DATAOUTPUT"/"$base".concordant_uniq >"$DATAOUTPUT"/"$base".concordant_uniq.bam
     # name sorting bam
         echo "Creating sorted bam for $base"
         samtools sort "$DATAOUTPUT"/"$base".concordant_uniq.bam -o "$DATAOUTPUT"/"$base".sorted.bam
         samtools index "$DATAOUTPUT"/"$base".sorted.bam
     # Clean up
         echo "Removing ""$TMP""/""$base"".bam"
     
         rm $DATAOUTPUT/"$base".concordant*
         rm $DATAOUTPUT/"$base".halfmapping*
         rm $DATAOUTPUT/"$base".nomapping*
     #   rm $DATAOUTPUT/"$base".paired*
         rm $DATAOUTPUT/"$base".unpaired*
     
     #counting the number of mapped reads :
     cd "$DATAOUTPUT" || exit
     
     samtools view -c "$base".sorted.bam |awk -v var="$base" 'END {print var"\t"$1}' > comptage_brute."${base}".txt
     samtools view -F 0x4 "$base".sorted.bam | cut -f 1 | sort | uniq | wc -l |\
                 awk -v var="$base" 'END {print var"\t"$1}' > comptage_F04."${base}".txt ;
     
     samtools view "$base".sorted.bam |cut -f 3-5|uniq |awk -v var="$base" '{print var"\t"$0}' |gzip > mapq."${base}".txt.gz
     samtools depth "$base".sorted.bam |gzip > "$base".dp.gz  
     
     #plot depth along the genome:
     Rscript ../../00_scripts/Rscripts/plot_dp.R "$base".dp.gz
     
     #plot mapq along the genome: 
     Rscript  ../../00_scripts/Rscripts/plot_mapq.R mapq."$base".txt.gz
 else
     echo "$bamfile already ok"
 fi
