#!/bin/bash
#Author: QR
#Date: 11-2022
#script to detect repeated sequences

#--- EXTERNAL VARIABLE ---- #
#input : reference genome and a database name (e.g. the species name)
#output: several files
if [ $# -ne 2  ]; then
    echo "USAGE: $0 reference_genome database name"
    echo "Expecting the name of the reference genome and a basename for database building (e.g. 'myfavoritespecies')"
    exit 1
else
    genome=$1
    database=$2
    echo -e "reference genome is $genome \n"
    echo "database is : ${database}"
    echo " "
fi


base=$(basename $genome)

#--- USUAL GENERIQ STUFF:  ---#
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="99_log_files"
#create log folder
mkdir $LOG_FOLDER 2>/dev/null
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

## ----- step 1 -- run repeatmodeler ------ ##
sed 's/ [0-9A-Za-z=-]*//g' $genome > ${genome%.fa}.simpl.fa

genome=${genome%.fa}.simpl.fa
base=$(basename $genome )

BuildDatabase -name $database -engine ncbi $genome 2>&1 | tee $LOG_FOLDER/buildDatabase.$base.$TIMESTAMP.log
RepeatModeler -pa 24 -engine ncbi -database $database 2>&1 | tee $LOG_FOLDER/repeatmodeler_$base.$TIMESTAMP.log



## ----- step 2 -- run repeatmasker ----- ##
## ----- step2.1: based on fugrep ----- ##
FOLDER1="${base}"_mask_fugrep.$TIMESTAMP
mkdir $FOLDER1
lib1=/home/quentin/database/fugrep.ref #change according to your repeat lib
RepeatMasker -pa 24 -e ncbi -lib $lib1 -xsmall -dir "$FOLDER1" $genome 2>&1 | tee $LOG_FOLDER/repeatmasker_fugrep_$base.$TIMESTAMP.log

## ----- step2.2: based on fngrep ----- ##
# database 2:
lib2=/home/quentin/database/fngrep.aa.fasta #edit this if necessary!
fnbase=$(basename $lib2)

FOLDER2="${base}"_mask_fngrep.$TIMESTAMP
mkdir $FOLDER2
RepeatMasker -pa 24 -e ncbi -lib $lib2 -xsmall -dir "$FOLDER2" "$FOLDER1"/$base.masked 2>&1 |\
	tee $LOG_FOLDER/repeatmasker_fngrep_$base.$TIMESTAMP.log


## ----- step2.3: based on de-novo repeat + repbase ----- ##
lib3=/home/quentin/database/repbase20.05_aaSeq_cleaned_TE.fa
lib3base=$(basename $lib3)
FOLDER3="${base}"_mask_"$base"_"$lib3base"."$TIMESTAMP"
mkdir "$FOLDER3"

cat $database-families.fa $lib3 > $base.repbase.fa
RepeatMasker -pa 24 -e ncbi -lib $base.repbase.fa -xsmall -dir "$FOLDER3" "$FOLDER2"/"$base".masked.masked 2>&1 |\
	tee $LOG_FOLDER/repeatmasker_$base.repbase20.05.$base.$TIMESTAMP.log


## ----- step 2.4: based on online data ----- ## 
#online database
FOLDER4="${base}"_mask_fungi.$TIMESTAMP
mkdir "$FOLDER4"
RepeatMasker -pa 24 -e ncbi -species fungi -xsmall -dir "$FOLDER4"   "$FOLDER3"/"$base".masked.masked.masked 2>&1 | \
	tee $LOG_FOLDER/repeatmasker_fungi.$base.$TIMESTAMP.log
