#!/bin/bash

#Author: QR
#Date: 11-2022
#script to detect repeated sequences

#------------ EXTERNAL VARIABLE FROM CONFIG FILE -------------- #
source ../config/config
source ../config/cpu_mem

#------------- CONDA ACTIVATION  -------------- #
eval "$(conda shell.bash hook)"
conda activate repeatmodeler_env


echo -e "TEdatabase is ""$TEdatabase"" "
echo -e "NCBI species is ""$ncbi_species"" "
echo -e "genome is ""$genome"" "

#------------- CHECK PARAMETERS -------------- #
if [ $# -ne 3  ]; then
    echo "USAGE: $0 reference_genome database name rm_unknwon(yes/no)"
    echo -e "Expecting the following parameters:\n
          1 - the reference genome\n
          2 - a basename for database building (e.g. 'myfavoritespecies')\n
          3 - a string YES/NO about wether unknwon repeat should be removed\n\n"
    exit 1
else
    genome=$1
    database=$2
    rm_unknown=$3
    echo -e "reference genome is $genome \n"
    echo -e "database is : ${database}\n"
    echo -e "rm_unknown option is set to $rm_unknown"
    echo " "
fi


base=$(basename "$genome")
#--------------DECLARE THE USUAL GENERIQ STUFF: -----------------#
TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)
LOG_FOLDER="LOGS"
#create log folder
mkdir $LOG_FOLDER 2>/dev/null
NCPUS="$NCPUS_REPEATEMODELER"

# ----- check compression of fasta  ------ ##
#check compression
if file --mime-type "$genome" | grep -q gzip$; then
   echo "$genome is gzipped"
   gunzip "$genome"
   genome=${genome%.gz}
else
   echo "$genome is not gzipped"
   genome=$genome
   echo "$genome" 
fi


if file --mime-type "$TEdatabase" | grep -q gzip$; then
   echo "$TEdatabase is gzipped"
   gunzip "$TEdatabase"
   TEdatabase=${TEdatabase%.gz}
else
   echo "$TEdatabase is not gzipped"
   TEdatabase=${TEdatabase%.gz} #just in case
fi

base=$(basename "$genome" )

mkdir 05_TE 2>/dev/null
cd 05_TE || exit

#--------------STEP1 : RUN REPEATMODELER  -----------------------#

##build db:
BuildDatabase -name "$database" -engine ncbi ../"$genome" 2>&1 |\
    tee ../$LOG_FOLDER/buildDatabase."$base"."$TIME".log

#de novo TE annotations:
RepeatModeler -threads "$NCPUS" -engine ncbi -database "$database" 2>&1 |\
    tee ../$LOG_FOLDER/repeatmodeler_"$base"."$TIME".log 
if [[  "${PIPESTATUS[0]}" -ne 0 ]]
then
    echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo -e "\terror repeatmodeler failed"
    echo -e "please check file LOGS/repeatmodeler_*log"
    echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    exit 1
fi


#--------------STEP2 : RUN REPEATMASKER -------------------------#

# ROUND1 BASED ON DATABASE : 
FOLDER1=FOLDER1_"${base}"_mask."$TIME"
mkdir "$FOLDER1"
lib1="$TEdatabase" 
RepeatMasker -pa "$NCPUS" -e ncbi -lib "$lib1" -noint -xsmall -dir "$FOLDER1" ../"$genome" 2>&1 |\
    tee ../"$LOG_FOLDER"/F1_repeatmasker_"$base"."$TIME".log 
if [[  "${PIPESTATUS[0]}" -ne 0 ]]
then
    echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo -e "\terror repeatmasker failed"
    echo -e "please check file LOGS/F1_repeatmasker_*log"
    echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    exit 1
fi
awk '{print $5"\t"$6"\t"$7"\t"$11}' "$FOLDER1"/"$database".fa.out |sed '1,3d' > Round1.bed

# Based on de-novo repeat + database:
FOLDER2=FOLDER2_"${base}"_mask."$TIME"

# test if we keep Unknwon repeat or not
## without Unknwon repeat ##
if [[ "$rm_unknown" = "YES" ]]
then
    echo "removing Unknown TE Repeats ..."
    awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "$database"-families.fa |\
    sed -e '/Unknown/,+1d' |\
    cat "$lib1" - > "$base".no_unknown.repbase.fa
    libcat="$base".no_unknown.repbase.fa
else
    #with known repeat
    echo "keep all candidate TEs... "
    cat "$database"-families.fa "$lib1" > "$base".repbase.fa
    libcat="$base".repbase.fa
fi

#ROUND2:
#run repeatmasker:
RepeatMasker -pa "$NCPUS" -e ncbi -lib "$libcat" -xsmall -dir "$FOLDER2" "$FOLDER1"/"$base".masked 2>&1 |\
    tee ../$LOG_FOLDER/F2_repeatmasker_"$base"."$TIME".log  
if [[  "${PIPESTATUS[0]}" -ne 0 ]]
then
    echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo -e "\terror repeatmasker failed"
    echo -e "please check file LOGS/F2_repeatmasker_*log"
    echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    exit 1
fi
awk '{print $5"\t"$6"\t"$7"\t"$11}' "$FOLDER2"/"$database".fa.masked.out |sed '1,3d' > Round2.bed

#
FOLDER3=FOLDER3_"${base}"_mask.$TIME
mkdir "$FOLDER3"

#run repeatmasker:
RepeatMasker -pa "$NCPUS" -e ncbi -species  "${ncbi_species}" -xsmall -dir "$FOLDER3"   "$FOLDER2"/"$base".masked.masked 2>&1 | \
    tee ../$LOG_FOLDER/F3_repeatmasker_"$base"."$TIME".log  ||\
if [[  "${PIPESTATUS[0]}" -ne 0 ]]
then
    echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo -e "\terror repeatmasker failed"
    echo -e "please check file LOGS/F3_repeatmasker_*log"
    echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    exit 1
fi
awk '{print $5"\t"$6"\t"$7"\t"$11}' "$FOLDER3"/"$database".fa.masked.masked.out |sed '1,3d' > Round3.bed

#FOLDER4=FOLDER4_"${base}"_mask.$TIME
#mkdir "$FOLDER4"

#run repeatmasker:
#awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "$database"-families.fa |\
#    grep -A1 Unknown |sed '/^--$/d' > unknown.fa

cat Round*bed >> ../03_genome/raw."$database".TE.bed
#make overly simplistic filtered bed: 
grep -v "Simple\|Low\|Unsp" ../raw."$database".TE.bed > ../03_genome/filtered."$database".TE.bed

if [[ $rm_unknown = "YES" ]]
then
   bedtools maskfasta -soft \
    -fi ../"$genome"  \
    -bed ../03_genome/filtered."$database".TE.bed \
    -fo ../03_genome/genome.wholemask.fa
else
   bedtools maskfasta -soft \
    -fi ../"$genome"  \
    -bed ../03_genome/raw."$database".TE.bed \
    -fo ../03_genome/genome.wholemask.fa
fi

#RepeatMasker -pa 18 -e ncbi -lib unknown.fa  -dir "$FOLDER4" "$FOLDER3"/"$base".masked.masked.masked 2>&1 | \
#    tee ../$LOG_FOLDER/F4_repeatmasker_"$base"."$TIME".log  ||\
#if [[  "${PIPESTATUS[0]}" -ne 0 ]]
#then
#    echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
#    echo -e "\terror repeatmasker failed"
#    echo -e "please check file LOGS/F3_repeatmasker_*log"
#    echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
#    exit 1
#fi
#awk '{print $5"\t"$6"\t"$7"\t"$11}' "$FOLDER4"/"$database".fa.masked.masked.masked.out |sed '1,3d' > Round4.bed
#

#cd ../03_genome || exit


