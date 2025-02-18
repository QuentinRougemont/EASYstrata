#!/bin/bash 
source config/config
source config/colors

echo rnaseqlist is $RNAseqlist 
#launch RNAseq trimming 

TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)
LOG_FOLDER="LOGS"

#create folder if not existent:
mkdir $LOG_FOLDER 2>/dev/null

if [ -f Trimmomatic-0.39/trimmomatic-0.39.jar ] ; then 
    echo "found trimmomatic jar " ; 
else 
    echo -e "trimmomatic jar not found\nwill attempt to download" 
    wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip 
    if ! unzip Trimmomatic-0.39.zip
    then
        echo -e "${RED} ERROR! failed downloading trimmomatic - check your internet connexion  \n${NC}"
        exit 1
    else
        echo -e "trimmomatic download successfull\n"
fi

fi
echo -e "\n\ntrimming read for RNAseq\n" 

#in case there would be space instead of tab:
awk -v OFS="\t" '$1=$1' "$RNAseqlist" > tmp
mv tmp "$RNAseqlist" 
rm tmp
ncol=$(awk '{print NF}'  "$RNAseqlist" |uniq)

if [[ $ncol = 2 ]] ; then
    #assuming PE:
    echo -e "\nrunning trimmomatic assuming reads are Paired-End\n" 
    while IFS=$'\t' read -r -a read ; 
    #while IFS=  read -r -a read ; 
    do 
        ./00_scripts/01_trimmomatic_PE.sh "${read[0]}" "${read[1]}"  
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR TRIMMOMATIC FAILED"
            exit 1
        fi

    done < "$RNAseqlist" #  file1file2.tmp 

    if [ $? -eq 0 ]; then
        echo "trimmomatic complete "
        if [ ! -e read_count.txt ] ; then
            echo -e "\ncounting the number of retained reads\n" 
           ./00_scripts/utility_scripts/count_read_fastq.sh 02_trimmed/*gz > read_count.txt
        fi

    else
        echo -e "\n${RED}#ERROR : Runnning trimmomatic failed. please check your input files${NC}"
        exit 1
    fi

else
    #assuming SE:
    echo "running trimmomatic" 
    while IFS=$'\t' read -r -a read ; 
    do 
        ./00_scripts/01_trimmomatic_SE.sh "${read[0]}" #2>&1 |tee "$LOG_FOLDER"/trimmo_"$(basename ${read[0]})"_log  
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR TRIMMOMATIC FAILED"
            exit 1
        fi

    done < "$RNAseqlist" #  file1file2.tmp 
    
    if [ $? -eq 0 ]; then
        echo trimmomatic complete
        if [ ! -e read_count.txt ] ; then
            echo -e "\ncounting the number of retained reads\n" 
           ./00_scripts/utility_scripts/count_read_fastq.sh 02_trimmed/*gz > read_count.txt
        fi
    else
        echo -e "\n${RED} ERROR : Runnning trimmomatic failed. please check your input files ${NC}"
        exit 1
    fi

fi

rm -rf Trimmomatic-0.39* 
