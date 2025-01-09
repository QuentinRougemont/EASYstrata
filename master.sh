#!/bin/bash

#master file to run the pipeline :
#Date: 11-2023
#Author: QR
source config/colors 
source ./config/config 

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo " "
   echo "Usage: $0 [-h]"
   echo "-h|--help: Print this Help."
   echo -e "Main script to launch all subsequent analyses\n"
   echo -e "${RED}----ALL RECQUIRED PARAMETERS MUST BE PROVIDED IN THE CONFIG 
   FILE (see ./config/config) ----${NC}\n"
   echo -e "use cases:

    ${bolditalic}${BLU}**** options **** : ${NC}${normal} \n
    -o 1 : perform whole inference where: \n
            step 1 : repeatmodeler + genome annotation + quality checks ; \n
            step 2 : GeneSpace  Gene Synteny + Whole Genome Synteny; \n 
            step 3 : Ds + related analyses  ;
    -o 2 : perform step 1 & 2 where :\n
            step 1 : repeatmodeler + genome annotation + quality checks ; \n
            step 2 : GeneSpace/Synteny 
    
    -o 3 :  perform step 2 & 3 where :\n
            step 2 : GeneSpace  Gene Synteny + Whole Genome Synteny; \n 
            step 3 : Ds + related analyses ;

    -o 4:  perform step 3 only where :\n
            step 3 : Ds + related analyses ;

    -o 5: perform step2 only where : 
            step 2: GeneSapce/Gene Synteny + Whole Genome Synteny
            
    -o 6 : perform step1 only where :
           step  1 : repeatmodeler + genome annotation + quality checks ; \n
    
    -o 7 : perfrom only the changepoint analyis of evolutionary strata

    All the details of the dat MUST be provided in the config file) \n
    see details in .infos/input_data_info.md
     \n
        ${bolditalic} - Analyses that will be performed: ${normal} 
        genome annotation, single copy orthologs inference, genespace, Dsplot, 
        Rideogram, dotplot, circos, changepoint analysis
    
        "
}

############################################################
# Process the input options.                               #
############################################################


while [ $# -gt 0 ] ; do
  case $1 in
    -h  | --help ) Help ; exit 2 ;;
    -o  | --option ) option="$2" ; 
        echo -e "options activated : ***${option}*** \n" >&2;;
   esac
   shift
done 

if [ -z "$option" ] ; then
    echo "Error no option provided ! I don't know what to do"
    Help
    exit 2
fi

###########################################################
# a few minor chek here:
###########################################################
# verify that at least a genome is provided - 
# otherwise no work can be done at all:
if  [ -z "${genome1}" ]  && [ -z "${genome2}" ] ; 
then 
    echo "Error! provide the genome of at least one species"
    Help
    exit 2
else 
    if [ -n "${genome1}" ] ; 
    then
        b1=$(basename "${genome1%.fa*}" )
        if [[ "$b1" =~ [^a-zA-Z0-9.] ]] ; 
        then 
            echo "error only alphanumeric character allowed in genomeIDs" 
            echo "see readme.md on github"
        else 
            echo "genome 1 is $genome1" 
        fi
    fi
    if [ -n "${genome2}" ] ;
    then
        b1=$(basename "${genome2%.fa*}" )
        if [[ "$b1" =~ [^a-zA-Z0-9.] ]] ; 
        then 
            echo "error only alphanumeric character allowed in genomeIDs" 
            echo "see readme.md on github"
        else 
            echo "genome 2 is $genome2" 
        fi
    fi
fi

# if no TE database then exit
# if no lineage for busco then exit and ask for it 
if [ -z "$TEdatabase" ] && [[ "$annotate" = YES ]] ;
then
    echo "Error NO Database for TE is provided "
    echo "I will not be able to identify and mask TEs prior to the genome annotation step "
    echo "this is very bad" 
    exit
    Help
fi


# if no lineage for busco then exit and ask for it 
if [ -z "${busco_lineage}" ]  && [[ "$annotate" = YES ]] ;
then
    echo "Error NO lineage provided for busco analyses" 
    echo "I will not be able to assess the quality of the annotation runs, 
    which is compulsory for these anaylses"
    exit
    Help
fi

########################################
## Global variables
########################################

echo "------------------------------------------------------------"
echo "-----check all variables from the configuration file  ------"
echo "  ancestral genome ${ancestral_genome}  "
echo "  ancestral gff is ${ancestral_gff} "
echo "*** fasta for genome1 is ${genome1} **** "
echo "*** fasta for genome2 is ${genome2} **** "
echo "*** haplotype1 is ${haplotype1} ***"
echo "*** haplotype2 is ${haplotype2}"
echo "*** shoulbe genome be annotated ? $annotate" 
echo "*** is RNAseq provided? $rnaseq ***"
echo "*** RNAseq list is: ${RNAseqlist} ****" 
echo "*** Alternatively, BAM list might be : ${bamlist1} & ${bamlist2}****"
echo "*** should TE be annotated ? $annotateTE ****"
echo "*** TE database is : $TEdatabase ***"
echo "*** lineage for busco analyses will be:  $busco_lineage ***"
echo "*** is species a fungus? $fungus ****"
echo -e "------------------------------------------------------------\n\n"



############################################################
# Generate architecture:
############################################################

echo "Generate Architecture"
mkdir -p haplo1/03_genome 

if [[ -z "${haplotype1}" ]] ; then
    haplotype1="$(basename "${genome1%.fa**}")"
fi

# ----- check compression of fasta  ------ ##

#eval "$(conda shell.bash hook)"
#conda activate superannot

#check compression
if file --mime-type "$genome1" | grep -q gzip$; then
   echo "$genome1 is gzipped"
   gunzip "$genome1"
   genome1=${genome1%.gz}
   cd haplo1/03_genome || exit 1 
       cp "$genome1" "$haplotype1".fa
   cd ../../
else
   echo "$genome1 is not gzipped"
   genome1=$genome1
   cd haplo1/03_genome || exit 1 
       cp "$genome1" "$haplotype1".fa
   cd ../../
fi

#--- handling genome 2 -----
if [[ -n "$genome2" ]] ; then
    
    if [[ -n "${haplotype2}" ]] && [[ -n "${genome2}" ]]; then 
        mkdir -p haplo2/03_genome ; 
    fi
    if [[ -z "${haplotype2}" ]] && [[ -n "${genome2}" ]]; then 
        haplotype2="$(basename "${genome2%.fa**}")"
        mkdir -p haplo2/03_genome ; 
    fi
    
    #check genome compression:
    # ----- check compression of fasta  ------ ##
    #check compression
    if file --mime-type "$genome2" | grep -q gzip$; then
    echo "$genome2 is gzipped"
    gunzip "$genome2"
    genome2="${genome2%.gz}"
    cd haplo2/03_genome  || exit 1
        cp "$genome2" "$haplotype2".fa
    cd ../../
    else
    echo "$genome2 is not gzipped"
    cd haplo2/03_genome || exit 1 
        cp "$genome2" "$haplotype2".fa
    
    genome2=$genome2
    cd ../../
    fi
fi

#conda deactivate

if [[ -n "${ancestral_genome}" ]] ; then 
    mkdir ancestral_sp  2>/dev/null
    #test if ancestral gtf is also provided in the config file:
    if [ -z "${ancestral_gff}" ] ; then 
        echo "error ! you provided an ancestral genome but no corresponding annotation (gff file)"
        echo "please correct that"
    fi

fi

#test if RNAseq info are provided or not:
if [ -z "${rnaseq}" ] 
then 
    echo -e "no info on RNAseq data\nwe assume no data are available"
    rnaseq="NO"
fi

if [[ $interpro = "YES" ]]
then
    echo "attempting to download and install inter-pro - this will take some time" 
    ./00_scripts/get_interpro.sh
fi

#artchitecture should be OK to proceed

############################################################
# Process the input options.                               #
############################################################

mkdir LOGS 2>/dev/null

echo -e "\n\ntesting use cases and checking inputs\n\n"


################################################################################
# ------------------ Section for option 1 (whole workflow) --------------------#
################################################################################

if [ "$option" == 1 ]; then 

    echo "----------------------------------------------------------------"
    echo "      the whole pipeline will be launched                       " 
    echo "          checking config files settings                        "
    echo "----------------------------------------------------------------"

    opt="synteny_and_Ds"

    #download uniprot for later checks
    ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 


    #if both species are provided without RNAseq:
    if [ -n "${genome1}" ] && [ -n "${genome2}" ]  && [[ $rnaseq = "NO" ]]  && [ -z "$ancestral_genome" ]  ; then
        echo "we will perform TE detection - genome annotation - Ds computation and plots"
        echo "genome are $genome1 and $genome2 "
        opt="synteny_and_Ds"
       
    
        cd haplo1 || exit 1
        #partie en erreur à débuguer:
        if [ -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n-------------------------------------------"
            echo "      running TE detection and gene prediction    "
            echo -e "---------------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa \
                -s "$haplotype1" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
       
        cd haplo2 || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n-------------------------------------------"
            echo "      running TE detection and gene prediction    "
            echo -e "---------------------------------------------\n"
            
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa \
                -s "$haplotype2" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../
        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
        
        #if all is OK then run GeneSpace - paml etc :
        #TO DO: modifiy the script RunGeneSpace etc to handle case with/without ancestral species
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

    #step1: 
    elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "YES" ]]  && [ -n "${RNAseqlist}" ]  && [ -z "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ] ; then
    
        echo "we will perform TE detection - genome annotation with RNAseq - Ds computation and plots"
        echo "genomes are ${genome1} and ${genome2}"
        
        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        ./00_scripts/launch_trimm.sh 2>&1 |tee LOGS/log_trimmomatic
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror trimmomatic failed"
                echo -e "please check file LOGS/log_trimmomatic"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
        fi

        cd haplo1/ || exit 1 
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            ../00_scripts/launch_rnaseq.sh "${haplotype1}"   2>&1 |\
                tee ../LOGS/log_rna_haplo1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi
            
            echo -e "\n----------------------------------------\n"
            echo -e "launching next steps....\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
       
        cd haplo2/ || exit 1 
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            ../00_scripts/launch_rnaseq.sh "${haplotype2}"   2>&1 |\
                tee ../LOGS/log_rna_haplo2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/"$haplotype2".fa  \
                -s "$haplotype2" \
                -r YES \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi

        #then run GeneSpace etc :
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

    #step1: 
    elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "YES" ]]  && [ -n "${bamlist1}" ] && [ -n "${bamlist2}" ] && [ -z "$ancestral_genome" ] ; then
    
        echo -e "\nwe will perform TE detection - genome annotation with RNAseq - 
        dS computation and plots - bam files already provided"

        echo -e "genomes are ${genome1} and ${genome2}\n\n"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        cd haplo1/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES \
                -b YES \
                -b "$bamlist1" 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../
       
        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
        
        cd haplo2/ || exit 1 
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa  \
                -s "$haplotype2" \
                -r YES \
                -m YES \
                -f YES \
                -b "$bamlist2" 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

           cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
        
        #then run GeneSpace etc :
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log__GeneSpace_and_Co"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

    #step1: 
    elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "NO" ]] && [ -n "$ancestral_genome" ]  ; then
        echo "we will perform all analyses with annotations performed without rnaseq "
        echo "genomes are ${genome1} and ${genome2}"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        cd haplo1 || exit 1 

        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        cd haplo2  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa \
                -s "$haplotype2" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../
       
        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        #then run GeneSpace etc :
        #modifiy the script RunGeneSpace etc to handle case with/without ancestral species
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -a "$ancestral_genome" \
            -g "$ancestral_gff" \
            -c "$scaffold" \
            -o "$opt"  2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi


    #step1: 
    elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "YES" ]] && [ -n "${bamlist1}" ] && [ -n "${bamlist2}" ] && [ -n "$ancestral_genome" ]  ; then
        echo -e "\n\nwe will perform all analyses with annotations performed with rnaseq - 
        list of bam provided "
        echo "genomes are ${genome1} and ${genome2}\n\n"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        cd haplo1  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES \
                -b "$bamlist1"  2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        cd haplo2  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa \
                -s "$haplotype2" \
                -r NO \
                -m YES \
                -f YES \
                -b "$bamlist2" 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
                echo "please copy your genome here"    
            help
            exit 1
        fi
        
        #then run GeneSpace etc :
        #modifiy the script RunGeneSpace etc to handle case with/without ancestral species
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -a "$ancestral_genome" \
            -g "$ancestral_gff" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

 
    #step1: 
    elif [ -n "${genome1}" ] && [ -n "${genome2}" ]  && [[ $rnaseq = "YES" ]] && [ -n "${RNAseqlist}" ]  && [ -n "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ]
    then
        echo "we will perform all analyses including annotation with rnaseq"
        echo "genomes are ${genome1} and ${genome2}"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        ./00_scripts/launch_trimm.sh 2>&1 |tee LOGS/log_trimmomatic
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror trimmomatic failed"
                echo -e "please check file LOGS/log_trimmomatic"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
        fi

        cd haplo1/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            ../00_scripts/launch_rnaseq.sh \
                "${haplotype1}"   2>&1 |\
                tee ../LOGS/log_rna_haplo1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES 2>&1 \
                | tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        cd haplo2/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            ../00_scripts/launch_rnaseq.sh \
                "${haplotype2}"   2>&1 |\
                tee ../LOGS/log_rna_haplo1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n\n"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/"$haplotype2".fa \
                -s "$haplotype2" \
                -r YES \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        #then run GeneSpace etc :
        #modifiy the script RunGeneSpace etc to handle case with/without ancestral species
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -a "$ancestral_genome" \
            -g "$ancestral_gff" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

    ################################################################################################
    #handle case where only genome 1 is provided with full option requested: 
    ################################################################################################
    
    #if one single genome are provided without RNAseq:
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ]  && [[ $rnaseq = "NO" ]]  && [ -z "$ancestral_genome" ]  ; then
        echo -e "!!! WARNING !!! only one genome provided\n\tnot ideal for single copy ortholog inference"
        echo -e "make sure to provide information on the target haplotype in scaffold.txt" 
        echo "we will perform TE detection - genome annotation - Ds computation and plots"
        echo "genome is $genome1"
        opt="synteny_and_Ds"

        cd haplo1 || exit 1

        if [ -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n-------------------------------------------"
            echo "      running TE detection and gene prediction    "
            echo -e "---------------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa \
                -s "$haplotype1" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi

        #handle haplotype2 now - we assumme that haplotype2 is present in the genome of "$haplotype1"
        #we will extract proteins for $haplotype2 from the whole $haplotype1
        #we will remove them from $haplotype1 to create a separate dataset for genespace etc:
        mkdir -p haplo2/03_genome

        #linearise genome of haplo1  and extract scaffold to study:
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
          haplo1/03_genome/"$haplotype1".fa | grep -A1 "$haplotype2" > haplo2/03_genome/"$haplotype2".fa
        
        #
        #extract the gene from the gff :
        mkdir -p haplo2/08_best_run
        grep "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf > haplo2/08_best_run/"$haplotype2".final.gtf
        
        #extract the corresponding protein: 
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
        haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa |\
        grep -A1 "$haplotype2"  \
        > haplo2/08_best_run/"$haplotype2"_prot.final.clean.fa
        
        #then we must create a sub file from haplotype1 by excluding haplotype2 from within it:
        cp haplo1/08_best_run/"$haplotype1".final.gtf haplo1/08_best_run/"$haplotype1".final.gtf.bkp
        cp haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        
        grep -v "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf.bkp > haplo1/08_best_run/"$haplotype1".final.gtf
        #remove seq:
        input=haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
            else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$input" \
            |\
            grep -v "$haplotype2" |tr "\t" "\n"  > haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa

        #if all is OK then run GeneSpace - paml etc :
        #TO DO: modifiy the script RunGeneSpace etc to handle case with/without ancestral species
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

    #step1: 
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ] && [[ $rnaseq = "YES" ]]  && [ -n "${RNAseqlist}" ]  && [ -z "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ] ; then
        echo -e "!!! WARNING !!! only one genome provided\n\tnot ideal for single copy ortholog inference"
        echo -e "make sure to provide information on the target haplotype in scaffold.txt" 

        echo "we will perform TE detection - genome annotation with RNAseq - Ds computation and plots"
        echo "genomes are ${genome1} and ${genome2}"
        
        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        ./00_scripts/launch_trimm.sh 2>&1 |tee LOGS/log_trimmomatic
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror trimmomatic failed"
                echo -e "please check file LOGS/log_trimmomatic"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
        fi

        cd haplo1/ || exit 1 
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            ../00_scripts/launch_rnaseq.sh "${haplotype1}"   2>&1 |\
                tee ../LOGS/log_rna_haplo1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi
            
            echo -e "\n----------------------------------------\n"
            echo -e "launching next steps....\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
       

        #handle haplotype2 now - we assumme that haplotype2 is present in the genome of "$haplotype1"
        #we will extract proteins for $haplotype2 from the whole $haplotype1
        #we will remove them from $haplotype1 to create a separate dataset for genespace etc:
        mkdir -p haplo2/03_genome

        #linearise genome of haplo1  and extract scaffold to study:
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
          haplo1/03_genome/"$haplotype1".fa | grep -A1 "$haplotype2" > haplo2/03_genome/"$haplotype2".fa
        
        #
        #extract the gene from the gff :
        mkdir -p haplo2/08_best_run
        grep "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf > haplo2/08_best_run/"$haplotype2".final.gtf
        
        #extract the corresponding protein: 
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
        haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa |\
        grep -A1 "$haplotype2"  \
        > haplo2/08_best_run/"$haplotype2"_prot.final.clean.fa

        #then we must create a sub file from haplotype1 by excluding haplotype2 from within it:
        cp haplo1/08_best_run/"$haplotype1".final.gtf haplo1/08_best_run/"$haplotype1".final.gtf.bkp
        cp haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        
        grep -v "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf.bkp > haplo1/08_best_run/"$haplotype1".final.gtf
        #remove seq:
        input=haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
            else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$input" \
            |\
            grep -v "$haplotype2" |tr "\t" "\n" > haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa


        #then run GeneSpace etc :
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

    #step1: 
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ] && [[ $rnaseq = "YES" ]]  && [ -n "${bamlist1}" ] && [ -z "${bamlist2}" ] && [ -z "$ancestral_genome" ] ; then
        echo -e "!!! WARNING !!! only one genome provided\n\tnot ideal for single copy ortholog inference"
        echo -e "make sure to provide information on the target haplotype in scaffold.txt" 

        echo -e "\nwe will perform TE detection - genome annotation with RNAseq - 
        dS computation and plots - bam files already provided"

        echo -e "genomes are ${genome1} and ${genome2}\n\n"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        cd haplo1/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES \
                -b YES \
                -b "$bamlist1" 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../
       
        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
        
        #handle haplotype2 now - we assumme that haplotype2 is present in the genome of "$haplotype1"
        #we will extract proteins for $haplotype2 from the whole $haplotype1
        #we will remove them from $haplotype1 to create a separate dataset for genespace etc:
        mkdir -p haplo2/03_genome

        #linearise genome of haplo1  and extract scaffold to study:
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
          haplo1/03_genome/"$haplotype1".fa | grep -A1 "$haplotype2" > haplo2/03_genome/"$haplotype2".fa
        
        #
        #extract the gene from the gff :
        mkdir -p haplo2/08_best_run
        grep "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf > haplo2/08_best_run/"$haplotype2".final.gtf
        
        #extract the corresponding protein: 
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
        haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa |\
        grep -A1 "$haplotype2"  \
        > haplo2/08_best_run/"$haplotype2"_prot.final.clean.fa

        #then we must create a sub file from haplotype1 by excluding haplotype2 from within it:
        cp haplo1/08_best_run/"$haplotype1".final.gtf haplo1/08_best_run/"$haplotype1".final.gtf.bkp
        cp haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        
        grep -v "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf.bkp > haplo1/08_best_run/"$haplotype1".final.gtf
        #remove seq:
        input=haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
            else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$input" \
            |\
            grep -v "$haplotype2" |tr "\t" "\n" > haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa


        #then run GeneSpace etc :
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log__GeneSpace_and_Co"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

    #step1: 
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ] && [[ $rnaseq = "NO" ]] && [ -n "$ancestral_genome" ]  ; then
        echo -e "!!! WARNING !!! only one genome provided\n\tnot ideal for single copy ortholog inference"
        echo -e "make sure to provide information on the target haplotype in scaffold.txt" 

        echo "we will perform all analyses with annotations performed without rnaseq "
        echo "genomes are ${genome1} and ${genome2}"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        cd haplo1 || exit 1 

        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        #handle haplotype2 now - we assumme that haplotype2 is present in the genome of "$haplotype1"
        #we will extract proteins for $haplotype2 from the whole $haplotype1
        #we will remove them from $haplotype1 to create a separate dataset for genespace etc:
        mkdir -p haplo2/03_genome

        #linearise genome of haplo1  and extract scaffold to study:
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
          haplo1/03_genome/"$haplotype1".fa | grep -A1 "$haplotype2" > haplo2/03_genome/"$haplotype2".fa
        
        #
        #extract the gene from the gff :
        mkdir -p haplo2/08_best_run
        grep "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf > haplo2/08_best_run/"$haplotype2".final.gtf
        
        #extract the corresponding protein: 
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
        haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa |\
        grep -A1 "$haplotype2"  \
        > haplo2/08_best_run/"$haplotype2"_prot.final.clean.fa

        #then we must create a sub file from haplotype1 by excluding haplotype2 from within it:
        cp haplo1/08_best_run/"$haplotype1".final.gtf haplo1/08_best_run/"$haplotype1".final.gtf.bkp
        cp haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        
        grep -v "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf.bkp > haplo1/08_best_run/"$haplotype1".final.gtf
        #remove seq:
        input=haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
            else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$input" \
            |\
            grep -v "$haplotype2" |tr "\t" "\n" >  haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa

        
        #then run GeneSpace etc :
        #modifiy the script RunGeneSpace etc to handle case with/without ancestral species
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -a "$ancestral_genome" \
            -g "$ancestral_gff" \
            -c "$scaffold" \
            -o "$opt"  2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi


    #step1: 
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ] && [[ $rnaseq = "YES" ]] && [ -n "${bamlist1}" ] && [ -z "${bamlist2}" ] && [ -n "$ancestral_genome" ]  ; then
        echo -e "!!! WARNING !!! only one genome provided\n\tnot ideal for single copy ortholog inference"
        echo -e "make sure to provide information on the target haplotype in scaffold.txt" 

        echo -e "\n\nwe will perform all analyses with annotations performed with rnaseq - 
        list of bam provided "
        echo "genomes are ${genome1} and ${genome2}\n\n"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        cd haplo1  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES \
                -b "$bamlist1"  2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        #handle haplotype2 now - we assumme that haplotype2 is present in the genome of "$haplotype1"
        #we will extract proteins for $haplotype2 from the whole $haplotype1
        #we will remove them from $haplotype1 to create a separate dataset for genespace etc:
        mkdir -p haplo2/03_genome

        #linearise genome of haplo1  and extract scaffold to study:
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
          haplo1/03_genome/"$haplotype1".fa | grep -A1 "$haplotype2" > haplo2/03_genome/"$haplotype2".fa
        
        #
        #extract the gene from the gff :
        mkdir -p haplo2/08_best_run
        grep "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf > haplo2/08_best_run/"$haplotype2".final.gtf
        
        #extract the corresponding protein: 
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
        haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa |\
        grep -A1 "$haplotype2"  \
        > haplo2/08_best_run/"$haplotype2"_prot.final.clean.fa

        #then we must create a sub file from haplotype1 by excluding haplotype2 from within it:
        cp haplo1/08_best_run/"$haplotype1".final.gtf haplo1/08_best_run/"$haplotype1".final.gtf.bkp
        cp haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        
        grep -v "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf.bkp > haplo1/08_best_run/"$haplotype1".final.gtf
        #remove seq:
        input=haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
            else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$input" \
            |\
            grep -v "$haplotype2" |tr "\t" "\n" >  haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa


        #then run GeneSpace etc :
        #modifiy the script RunGeneSpace etc to handle case with/without ancestral species
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -a "$ancestral_genome" \
            -g "$ancestral_gff" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

 
    #step1: 
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ]  && [[ $rnaseq = "YES" ]] && [ -n "${RNAseqlist}" ]  && [ -n "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ]
    then
        echo -e "!!! WARNING !!! only one genome provided\n\tnot ideal for single copy ortholog inference"
        echo -e "make sure to provide information on the target haplotype in scaffold.txt" 

        echo "we will perform all analyses including annotation with rnaseq"
        echo "genomes are ${genome1} and ${genome2}"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        ./00_scripts/launch_trimm.sh 2>&1 |tee LOGS/log_trimmomatic
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror trimmomatic failed"
                echo -e "please check file LOGS/log_trimmomatic"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
        fi

        cd haplo1/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            ../00_scripts/launch_rnaseq.sh \
                "${haplotype1}"   2>&1 |\
                tee ../LOGS/log_rna_haplo1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES 2>&1 \
                | tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        

        #handle haplotype2 now - we assumme that haplotype2 is present in the genome of "$haplotype1"
        #we will extract proteins for $haplotype2 from the whole $haplotype1
        #we will remove them from $haplotype1 to create a separate dataset for genespace etc:
        mkdir -p haplo2/03_genome

        #linearise genome of haplo1 and extract scaffold to study:
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
          haplo1/03_genome/"$haplotype1".fa | grep -A1 "$haplotype2" > haplo2/03_genome/"$haplotype2".fa
        
        #
        #extract the gene from the gff :
        mkdir -p haplo2/08_best_run
        grep "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf > haplo2/08_best_run/"$haplotype2".final.gtf
        
        #extract the corresponding protein: 
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
        haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa |\
        grep -A1 "$haplotype2"  \
        > haplo2/08_best_run/"$haplotype2"_prot.final.clean.fa

        #then we must create a sub file from haplotype1 by excluding haplotype2 from within it:
        cp haplo1/08_best_run/"$haplotype1".final.gtf haplo1/08_best_run/"$haplotype1".final.gtf.bkp
        cp haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        
        grep -v "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf.bkp > haplo1/08_best_run/"$haplotype1".final.gtf
        #remove seq:
        input=haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
            else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$input" \
            |\
            grep -v "$haplotype2" |tr "\t" "\n" > haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa


        #then run GeneSpace etc :
        #modifiy the script RunGeneSpace etc to handle case with/without ancestral species
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -a "$ancestral_genome" \
            -g "$ancestral_gff" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

    fi
fi


################################################################################
# ------section for option 2 : TE + gene prediction + GeneSpace/Synteny--------#
################################################################################

if [ "$option" == 2 ] ; 
then 
    echo "----------------------------------------------------------------"
    echo "  TE + gene prediction + GeneSpace/Synteny will be launched     " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"

    #the test case are the same as option 1 for the gene prediction part
    
    #then we only need to do the geneSpace - so we must insert an option in script 11
    # to avoid running script 12 from within it. 
    # we still want minimap output for synteny so minimap should really be 


    #if both species are provided without RNAseq:
    if [ -n "${genome1}" ] && [ -n "${genome2}" ]  && [[ $rnaseq = "NO" ]]  && [ -z "$ancestral_genome" ]  ; then
        echo "genomes are $genome1 and $genome2 "
        
        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 
    
        cd haplo1  || exit 1
        #partie en erreur à débuguer:
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa \
                -s "$haplotype1" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
                cd ../
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
        
        cd haplo2  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa  \
                -s "$haplotype2" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
        
        #if all is OK then run GeneSpace / Synteny and optionally Ds (paml) etc :
        opt="synteny_only"
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi


    #step1: 
    elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "YES" ]]  && [ -n "${RNAseqlist}" ]  && [ -z "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ] ; then
        echo "genomes are ${genome1} and ${genome2}"
        
        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        ./00_scripts/launch_trimm.sh 2>&1 |tee LOGS/log_trimmomatic
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror trimmomatic failed"
                echo -e "please check file LOGS/log_trimmomatic"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
        fi

        cd haplo1/ || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            ../00_scripts/launch_rnaseq.sh \
                "${haplotype1}"   2>&1 |\
                tee ../LOGS/log_rna_haplo1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

                cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
        
        cd haplo2/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            ../00_scripts/launch_rnaseq.sh \
                "${haplotype2}"   2>&1 |\
                tee ../LOGS/log_rna_haplo2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa  \
                -s "$haplotype2" \
                -r YES \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

                cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi

        #if all is OK then run GeneSpace / Synteny and optionally Ds (paml) etc :
        opt="synteny_only"
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -c "$scaffold"  \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi


    #step1: 
    elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "YES" ]]  && [ -n "${bamlist1}" ] && [ -n "${bamlist2}" ] && [ -z "$ancestral_genome" ] ; then
    
        echo "we will perform TE detection - genome annotation with RNAseq - 
        Ds computation and plots - bam files already provided"
        echo "genomes are ${genome1} and ${genome2}"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        cd haplo1/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES \
                -b YES \
                -b "$bamlist1" 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

                cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
        
        cd haplo2/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa \
                -s "$haplotype2" \
                -r YES \
                -m YES \
                -f YES \
                -b "$bamlist2" 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

                cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
        
        #if all is OK then run GeneSpace / Synteny and optionally Ds (paml) etc :
        opt="synteny_only"
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -c "$scaffold" \
            -o "$opt" 2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi


    #step1: 
    elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "NO" ]] && [ -n "$ancestral_genome" ]  ; then
        echo "we will perform analyses with annotations performed without rnaseq "
        echo "genomes are ${genome1} and ${genome2}"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        cd haplo1  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        cd haplo2  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa  \
                -s "$haplotype2" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        #if all is OK then run GeneSpace / Synteny and optionally Ds (paml) etc :
        opt="synteny_only"
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -a "$ancestral_genome" \
            -g "$ancestral_gff" \
            -c "$scaffold" \
            -o "$opt"  2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi


    #step1: 
    elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "YES" ]] && [ -n "${bamlist1}" ] && [ -n "${bamlist2}" ] && [ -n "$ancestral_genome" ]  ; then
        echo "we will perform analyses with annotations performed with rnaseq - 
        list of bam provided "
        echo "genomes are ${genome1} and ${genome2}"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        cd haplo1  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES \
                -b "$bamlist1"  2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        cd haplo2  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa  \
                -s "$haplotype2" \
                -r NO \
                -m YES \
                -f YES \
                -b "$bamlist2" 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            
            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        #if all is OK then run GeneSpace / Synteny and optionally Ds (paml) etc :
        opt="synteny_only"
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -a "$ancestral_genome" \
            -g "$ancestral_gff" \
            -c "$scaffold" \
            -o "$opt"  2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi


    #step1: 
    elif [ -n "${genome1}" ] && [ -n "${genome2}" ]  && [[ $rnaseq = "YES" ]] && [ -n "${RNAseqlist}" ]  && [ -n "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ]   ; then
        echo "we will perform all analyses including annotation with rnaseq"
        echo "genomes are ${genome1} and ${genome2}"

        #download uniprot for later checks
        ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

        ./00_scripts/launch_trimm.sh 2>&1 |tee LOGS/log_trimmomatic
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror trimmomatic failed"
                echo -e "please check file LOGS/log_trimmomatic"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
        fi

        cd haplo1/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            ../00_scripts/launch_rnaseq.sh \
                "${haplotype1}"   2>&1 |\
                tee ../LOGS/log_rna_haplo1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            echo -e "----------------------------------------\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        cd haplo2/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            ../00_scripts/launch_rnaseq.sh \
                "${haplotype2}"   2>&1 |\
                tee ../LOGS/log_rna_haplo1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa \
                -s "$haplotype2" \
                -r YES \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            help
            exit 1
        fi
        
        #if all is OK then run GeneSpace / Synteny and optionally Ds (paml) etc :
        opt="synteny_only"
        echo -e "\n\t\t~~~~~~~~~~~~\n\tlaunching GeneSpace & other analysis\n\t\t~~~~~~~~~~~\n"
        ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
            -s1 "$haplotype1" \
            -s2 "$haplotype2" \
            -a "$ancestral_genome" \
            -g "$ancestral_gff" \
            -c "$scaffold" \
            -o "$opt"  2>&1 |\
            tee LOGS/log_GeneSpace_and_Co
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo -e "error some steps failed during GeneSpace/Ds analyses"
            echo -e "please check file LOGS/log_GeneSpace"
            echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
            exit 1
        fi

    fi    # launched earlier


fi


################################################################################
# ------section for option 3 : GeneSpace/Synteny + Ds analyses ----------------#
################################################################################


################################################################################
# ----------------- section for option 4 : Ds analyses  -----------------------#
################################################################################

################################################################################
# -------------- section for option 5 : GeneSpace/Synteny analyses ------------#
################################################################################

if [ "$option" == 3 ] || [ "$option" == 4 ] || [ "$option" == 5 ]; then 

        
    #test if the architecture and data are already present from a previous incomplete run  :
    if [ -f haplo1/08_best_run/"$haplotype1"_prot.fa ] && [ -f haplo1/03_genome/"$haplotype1".fa ]  ; then
        echo "genome1 already present ---- cleaned protein file for genome1 already present "

        if [ -f haplo2/08_best_run/"$haplotype2"_prot.fa ]  && [ -f haplo2/03_genome/"$haplotype2".fa ] ; then
            echo "genome1 already present ---- cleaned protein file for genome1 already present "
            echo "it seems all data for GeneSpace/Synteny and Ds are present, 
            will try running from here"
            
        fi
        
    elif [ -n "${gtf1}" ] && [ -n "${gtf2}" ] && [ -n "${genome1}" ] && [ -n "${genome2}" ] ; then
        
        #else we expect them to be provided in the config file 
        echo "gtf and genomes file were provided" 
        echo "we will run geneSpace, compute Ds and other kind of analyses"
        echo -e "\----------------------------------------\n"
        mkdir -p haplo1/08_best_run haplo1/03_genome 
        mkdir -p haplo2/08_best_run haplo2/03_genome

        
        #check compression
        if file --mime-type "$gtf1" | grep -q gzip$; then
           echo "$gtf1 is gzipped"
           gunzip "$gtf1"
           gtf1=${gtf1%.gz}
        else
           echo "$gtf1 is not gzipped"
        fi
        
        if file --mime-type "$gtf2" | grep -q gzip$; then
           echo "$gtf2 is gzipped"
           gunzip "$gtf2"
           gtf2=${gtf1%.gz}
        else
           echo "$gtf2 is not gzipped"
        fi

        cp "$gtf1" haplo1/08_best_run/"${haplotype1}".final.gtf
        cp "$gtf2" haplo2/08_best_run/"${haplotype2}".final.gtf
        cp "$genome1" haplo1/03_genome/
        cp "$genome2" haplo2/03_genome/


        if ! gffread -g "$genome1" -w haplo1/08_best_run/"$haplotype1".spliced_cds.fa "$gtf1" 
        then 
            echo "error failed to extract cds from genome $genome1 and gtf $gtf1"
            echo "please check input file synchronisation"
            exit
        else
            transeq -sequence haplo1/08_best_run/"$haplotype1".spliced_cds.fa \
                    -outseq haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa

        fi     
        if ! gffread -g "$genome2" -w haplo2/08_best_run/"$haplotype2".spliced_cds.fa  "$gtf2"
        then
            echo "error failed to extract cds from genome $genome2 and gtf $gtf2"
            echo "please check input file synchronisation"
            exit
        else
            transeq -sequence haplo2/08_best_run/"$haplotype2".spliced_cds.fa \
                    -outseq haplo2/08_best_run/"$haplotype2"_prot.final.clean.fa
        fi 
        
        #if [ -z "$ancestral_genome" ] ; then
        #    echo "no ancestral species"
        #    #leave the variable empty
        #    #do nothing 
        #else [ -n "$ancestral_genome" ] ; then
        #    echo "ancestral species existant"
        #fi
    fi
fi

#        echo ancestral_genome is "$ancestral_genome"
        
if [ "$option" == 3 ] ; then 

    echo "----------------------------------------------------------------"
    echo "       GeneSpace/Synteny + Ds analyses will be launched         " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------" 

    opt="synteny_and_Ds"
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
        -s1 "$haplotype1" \
        -s2 "$haplotype2" \
        -a "$ancestral_genome" \
        -g "$ancestral_gff" \
        -c "$scaffold" \
        -o "$opt" 2>&1 |\
        tee LOGS/log_GeneSpace_and_Co
    if [[  "${PIPESTATUS[0]}" -ne 0 ]]
    then
        echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        echo -e "error some steps failed during GeneSpace/Ds analyses"
        echo -e "please check file LOGS/log_GeneSpace"
        echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
        exit 1
    fi

elif [ "$option" == 4 ] ; then 
            
    echo "----------------------------------------------------------------"
    echo "         only Ds + associated analyses will be launched         " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"

    opt="Ds_only"
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
        -s1 "$haplotype1" \
        -s2 "$haplotype2" \
        -a "$ancestral_genome" \
        -g "$ancestral_gff" \
        -c "$scaffold" \
        -o "$opt" 2>&1 |\
        tee LOGS/log_GeneSpace_and_Co
    if [[  "${PIPESTATUS[0]}" -ne 0 ]]
    then
        echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        echo -e "error some steps failed during GeneSpace/Ds analyses"
        echo -e "please check file LOGS/log_GeneSpace"
        echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
        exit 1
    fi


elif [ "$option" == 5 ] ; then 
        
    echo "----------------------------------------------------------------"
    echo "          only GeneSpace/Synteny  analyses will be launched     " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"

    opt="synteny_only"
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
        -s1 "$haplotype1" \
        -s2 "$haplotype2" \
        -a "$ancestral_genome" \
        -g "$ancestral_gff" \
        -c "$scaffold"  \
        -o "$opt" 2>&1 |\
        tee LOGS/log_GeneSpace_and_Co
    if [[  "${PIPESTATUS[0]}" -ne 0 ]]
    then
        echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        echo -e "error some steps failed during GeneSpace/Ds analyses"
        echo -e "please check file LOGS/log_GeneSpace"
        echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
        exit 1
    fi


elif [ "$option" == 7 ] ; then

    echo "----------------------------------------------------------------"
    echo "         only changepoint analyses will be launched         "
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"

    opt="changepoint"
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
        -s1 "$haplotype1" \
        -s2 "$haplotype2" \
        -a "$ancestral_genome" \
        -g "$ancestral_gff" \
        -c "$scaffold" \
        -o "$opt" 2>&1 |\
        tee LOGS/log_GeneSpace_and_Co
    if [[  "${PIPESTATUS[0]}" -ne 0 ]]
    then
        echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        echo -e "error some steps failed during GeneSpace/Ds analyses"
        echo -e "please check file LOGS/log_GeneSpace"
        echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
        exit 1
    fi

fi

################################################################################
# -------------- section for option 6: TE + gene Prediction only --------------#
################################################################################

if [ "$option" = 6 ]; then 

    echo "----------------------------------------------------------------"
    echo "          Only TE + gene prediction  will be performed          " 
    echo "              checking config files settings                    "
    echo "----------------------------------------------------------------"

    #download uniprot for later checks
    ./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

    # option = 6 - only genome 1 - No RNAseq : 
    if [ -n "${genome1}" ] && [ -z "${genome2}" ]  && [[ $rnaseq = "NO" ]]   && [ -z "$ancestral_genome" ] ; then
        cd haplo1  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "only the genome of one species was provided" 
            echo "genome is ${genome1} "
            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa \
                -s "$haplotype1" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
    
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ]  && [[ $rnaseq = "NO" ]]   && [ -n "$ancestral_genome" ] ; then
        cd haplo1  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "only the genome of one species was provided" 
            echo "genome is ${genome1} "
            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa  \
                -s "$haplotype1" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
        
    # option = 6 - only genome 1 - No RNAseq : 
    elif [ -z "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "NO" ]] && [ -z "$ancestral_genome" ]   ; then
        cd haplo2  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "only the genome of one species was provided" 
            echo "genome is ${genome2} "
            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa \
                -s "$haplotype2" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
    
    elif [ -z "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "NO" ]] && [ -n "$ancestral_genome" ]   ; then
        cd haplo2  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo "only the genome of one species was provided" 
            echo "genome is ${genome2} "
            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n\n"
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa \
                -s "$haplotype2" \
                -r NO \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi

    # option = 6 - only genome 2 & RNAseq : 
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ] && [[ $rnaseq = "YES" ]] && [ -n "${RNAseqlist}" ]  && [ -z "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ] ; then
        cd haplo1/  || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "only the genome of one species was provided with RNAseq data" 
            echo "genome is ${genome1} "

        ./00_scripts/launch_trimm.sh 2>&1 |tee LOGS/log_trimmomatic
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror trimmomatic failed"
                echo -e "please check file LOGS/log_trimmomatic"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
        fi
            
            ../00_scripts/launch_rnaseq.sh \
                "${haplotype1}" 2>&1 |\
                tee ../LOGS/log_rna_haplo1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            echo "running TE detection and gene prediction"
            echo -e "\----------------------------------------\n\n"

            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap1
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap1"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi

    # option = 6 - only genome 1 & RNAseq from bamfiles: 
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ]  && [[ $rnaseq = "YES" ]]  && [ -z "$ancestral_genome" ] && [ -n "$bamlist1" ] && [ -z "$bamlist2" ] ; then
        cd haplo1/ || exit 1
        if [  -n "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "only the genome of one species was provided along with a list of bam" 
            echo "genome is ${genome1} "
            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------"
            
            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype1".fa \
                -s "$haplotype1" \
                -r YES \
                -m YES \
                -f YES \
                -b "$bamlist1" 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
    
    
    # option = 6 - only genome 2 & RNAseq :
    elif [ -z "${genome1}" ] && [ -n "${genome2}" ]  && [[ $rnaseq = "YES" ]] && [ -n "${RNAseqlist}" ]  && [ -z "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ]  ; then
        cd haplo2/ || exit 1
        if [  "$(ls -A 03_genome/ --ignore=Readme )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "only the genome of one species was provided with RNAseq data" 
            echo "genome is ${genome2} "

        ./00_scripts/launch_trimm.sh 2>&1 |tee LOGS/log_trimmomatic
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror trimmomatic failed"
                echo -e "please check file LOGS/log_trimmomatic"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
        fi


            ../00_scripts/launch_rnaseq.sh \
                "${haplotype2}" 2>&1 |\
                tee ../LOGS/log_rna_haplo2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror RNA seq failed"
                echo -e "please check file LOGS/log__rna_haplo2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n\n"

            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa  \
                -s "$haplotype2" \
                -r YES \
                -m YES \
                -f YES 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
    
    # option = 6 - only genome 2 & RNAseq from bam file:
    elif [ -z "${genome1}" ] && [ -n "${genome2}" ]  && [[ $rnaseq = "YES" ]] && [ -z "$ancestral_genome" ]  && [ -z "$bamlist1" ] && [ -n "${bamlist2}" ]  ; then
        cd haplo2/ || exit 1
        if [ "$(ls -A 03_genome/ --ignore=Readme  )"  ] ; then
            echo -e "\n\n----------------------------------------"
            echo "only the genome of one species was provided along with a list of bam" 
            echo "genome is ${genome2} "
            echo "running TE detection and gene prediction"
            echo -e "----------------------------------------\n\n"

            ../00_scripts/launch_step05_to_08.sh \
                -g 03_genome/"$haplotype2".fa \
                -s "$haplotype2" \
                -r YES \
                -m YES \
                -f YES \
                -b "$bamlist2" 2>&1 |\
                tee ../LOGS/log_step05_08_hap2
            if [[  "${PIPESTATUS[0]}" -ne 0 ]]
            then
                echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                echo -e "\terror some steps failed"
                echo -e "please check file LOGS/log__step05_08_hap2"
                echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
                exit 1
            fi

            cd ../

        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"    
            Help
            exit 1
        fi
    
    fi

fi
