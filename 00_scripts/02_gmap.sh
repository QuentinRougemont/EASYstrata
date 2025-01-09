#!/bin/bash

#script to build genome for gsnap
if [ $# -ne 1  ]; then
    echo "USAGE: $0 fasta file "
    echo "Expecting a fasta file from the reference genome -- located in '03_genome' folder "
    echo "Extension should be '.fa/.fas/.fasta/.fa'"
    echo "File should be uncompressed"
    exit 1
else
    FASTA=$1
    echo "fasta file is : ${FASTA}"
    echo " "
fi

#source config/config

GENOMEFOLDER="03_genome/"
base=$(basename "$FASTA" )
GENOME=gmap_"${base%.fa**}"  # "gmap_genome_MspDse1_assembly.fa"

if [ ! -s $GENOMEFOLDER/"$GENOME"/"$GENOME".sarray16 ] 
then
    #run gmap
    gmap_build --dir="$GENOMEFOLDER" "$FASTA" -d "$GENOME"
else
    echo "gmap index ok"
fi
