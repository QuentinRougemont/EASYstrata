#!/bin/bash
#Purpose:
#master script to prepare bed files, laucnch GeneSpace, run paml and launch downstream Rscript 
#will check existence of all dependencies
#Date: 2023
#Author: QR

# -- some colors for warnings in the terminal  --:
source config/colors

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "master script to: \n 1 - create bed files, \n 2 - launch GeneSpace, \n 3 - run paml and launch downstream Rscripts (Rideogram, plot paml, etc)"
   echo " "
   echo "Usage: $0 [-s1|-s2|-f|-a|-g|-h|]"
   echo "options:"
   echo " -h|--help: Print this Help."
   echo " -s1|--haplo1: the name of the first  focal haplotype "
   echo " -s2|--haplo2: the name of the second focal haplotype "
   echo " -a|--ancestral_genome: the name of the ancestral haplo to infer orthology and plot gene order"
   echo " -g|--ancestral_gff: the name of the ancestral gff associated with the ancestral genome"
   echo " -f|--folderpath: the path to the global folder containing haplo1 and haplo 2"
   echo " -c|--chromosome: a tab separated txt file listing the name of the reference species (e.g sp1), the corresponding set of chromosomes (e.g.: chrX , supergene, etc) and the orientation of the chromosome (N: Normal, R: Reverse) if their is more than one"
   echo " -o|--options : the type of analysis to be performed: either 'synteny_and_Ds' (GeneSpace+Minimap2+Ds+changepoint), 'synteny_only' (GeneSpace+Minimap2), 'Ds_only' (paml and changepoint)"
   echo " "
   echo "dependancies: orthofinder, mcscanx, GeneSpace, paml (yn00), Rideogram, translatorX minimap2"
}


############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -s1 | --haplo1) haplo1="$2" ; echo -e "haplotype 1 Name is ***${haplo1}*** \n" >&2;;
    -s2 | --haplo2) haplo2="$2" ; echo -e "haplotype 2 Name is ***${haplo2}*** \n" >&2;;
    -a  | --ancestral_genome) ancestral_genome="$2" ; 
        echo -e "ancestral haplo  Name is ***${ancestral_genome}*** \n" >&2;;
    -g  | --ancestral_gff) ancestral_gff="$2" ; 
        echo -e "ancestral gff  Name is ***${ancestral_gff}*** \n" >&2;;
    -f  | --folderpath  ) folderpath="$2"   ; 
        echo -e "global folder is  ${folderpath} \n" >&2;;
    -c  | --chromosome )  chromosome="$2"   ; 
        echo -e "target chromosome are ${chromosome} \n" >&2 ;; 
    -o  | --options ) options="$2" ; 
        echo -e "options for computation are ***${options}*** \n" >&2 ;;
    -h  | --help ) Help ; exit 2 ;;
   esac
   shift
done 

if [ -z "${haplo1}" ] || [ -z "${haplo2}" ] || [ -z "${chromosome}" ] || [ -z "${options}" ]  ; then
    Help
    exit 2
fi


scaffold=$chromosome
mkdir 02_results 2>/dev/null #ignore if already existent

#make ancestral species optional
if [ -n "${ancestral_genome}" ] ; then
    echo "ancestral_species is $ancestral_genome "
    echo "will attempt to extract the CDS and PROT from it "
    mkdir ancestral_sp/03_genome 2>/dev/null #note: this folder exist already 
    # ----- check compression of fasta  ------ ##
    #check compression
    if file --mime-type "$ancestral_genome" | grep -q gzip$; then
       echo "$ancestral_genome is gzipped"
       gunzip "$ancestral_genome"
       ancestral_genome="${ancestral_genome%.gz}"
    else
       #trim any eventual .gz extension from file
       ancestral_genome="${ancestral_genome%.gz}"
       echo "$ancestral_genome is not gzipped"

    fi
   
    if file --mime-type "$ancestral_gff" | grep -q gzip$; then
       echo "$ancestral_gff is gzipped"
       gunzip "$ancestral_gff"
       ancestral_gff="${ancestral_gff%.gz}"
    else
       echo "$ancestral_gff is not gzipped"
       #trim any eventual .gz extension from file
       ancestral_gff="${ancestral_gff%.gz}"

    fi

    cd ancestral_sp || exit 
    if [ -f ancestral_sp.fa ] ; then
        rm ancestral_sp.fa
    fi

    ln -s "${ancestral_genome}" ancestral_sp.fa ; samtools faidx ancestral_sp.fa ; cd ../
    gffread -g "${ancestral_genome}" -w ancestral_sp/ancestral_sp.spliced_cds.fa  "${ancestral_gff}" 
    transeq -sequence ancestral_sp/ancestral_sp.spliced_cds.fa \
        -outseq ancestral_sp/ancestral_sp_prot.fa
    awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10}' "$ancestral_gff" |\
        sed 's/"//g' > ancestral_sp/ancestral_sp.bed
    sed -i 's/_1 CDS=.*$//g'  ancestral_sp/ancestral_sp_prot.fa

fi


#------------------------------ step 1 prepare bed file for each haplo -------------------------------------#
#

#test options :
if [[ $options = "Ds_only" ]] ;  
then
    mkdir -p 02_results/paml 02_results/plots
elif [[ $options = "plots" ]] ; 
then
    mkdir -p 02_results/plots
elif [[ $options = "synteny_and_Ds" ]] ; 
then
    rm -rf genespace peptide 02_results/paml 02_results/plots  #2>/dev/null
    mkdir -p genespace/bed genespace/peptide 
    mkdir -p  02_results/paml 02_results/plots 
elif [[ $options = "synteny_only" ]] ; 
then
    rm -rf genespace/bed genespace/peptide 02_results/paml 02_results/plots #2>/dev/null
    mkdir -p genespace/bed genespace/peptide 
    mkdir 02_results/plots 
    mkdir 02_results/paml
elif [[ $options == "changepoint" ]] ;
then
        echo "only changepoint will be performed"
fi


# create bed
awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10}' haplo1/08_best_run/"$haplo1".final.gtf |\
    sed 's/"//g' |sed 's/;//g' > genespace/bed/"$haplo1".bed
awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10}' haplo2/08_best_run/"$haplo2".final.gtf |\
    sed 's/"//g' |sed 's/;//g'  > genespace/bed/"$haplo2".bed

# simplify the protein file to match the bed (i.e. remove the _1 inserted by transeq and the CDS length info):
sed 's/_1 CDS=.*$//g' haplo1/08_best_run/"$haplo1"_prot.final.clean.fa \
    > genespace/peptide/"$haplo1".fa
sed 's/_1 CDS=.*$//g' haplo2/08_best_run/"$haplo2"_prot.final.clean.fa \
    > genespace/peptide/"$haplo2".fa

#verify that IDs in bed and fasta file are matching - else exit  
grep ">" genespace/peptide/"$haplo1".fa |sed 's/>//g' > tmp1
grep ">" genespace/peptide/"$haplo2".fa |sed 's/>//g' > tmp2

check1=$(grep -Ff tmp1 genespace/bed/"$haplo1".bed |wc -l )
check2=$(grep -Ff tmp2 genespace/bed/"$haplo2".bed |wc -l )

echo -e "check2 size is $check2"
echo -e "check1 size is $check1"

bedsize1=$(wc -l genespace/bed/"$haplo1".bed |awk '{print $1}' )
bedsize2=$(wc -l genespace/bed/"$haplo2".bed |awk '{print $1}' )

echo -e "bedisze2  size is $bedsize2"
echo -e "bedisze1  size is $bedsize1"

#check that all is matching:
if [ "$bedsize1" = "$check1" ]
then
    echo "input1 is ok" 
    rm tmp1
else
    echo "input1 is not ok"
    echo "check your data"
    exit 2
fi

if [ "$bedsize2" = "$check2" ]
then
    echo "input2 is ok" 
    rm tmp2
else
    echo "input2 is not ok"
    echo "check your data"
    exit 2
fi

# -- handling ancestral haplo ------
# -- this part assumes that a bed and peptide file are existant for the ancestral haplo
# -- here we used a genome annotated with the same pipeline relying on braker 
if [ -n "${ancestral_genome}" ] ; then
    cd genespace/bed/ || exit 1
    ln -s ../../ancestral_sp/ancestral_sp.bed . 
    cd ../peptide || exit 1
    ln -s ../../ancestral_sp/ancestral_sp_prot.fa ancestral_sp.fa
    
    cd ../../
fi

#------------------------------ step 2 run GeneSpace ---------------------------------------------------------#

if [[ $options = "synteny_and_Ds" ]]  || [[ $options = "synteny_only" ]] ; then
    cd genespace  || exit 1

    
    MCScanpath=$(command -v MCScanX |xargs dirname )
    #just in case this is not already done:
    Rscript -e  'devtools::install_github("jtlovell/GENESPACE")'

    sed -i "s#mcpath#$MCScanpath#" ../00_scripts/Rscripts/01.run_geneSpace.R
    
    Rscript ../00_scripts/Rscripts/01.run_geneSpace.R || exit 1
    
    #plot genespace subspace of target chromosomes: 
    #a refaire en fonction de si ancestral species or not:
    echo scaffold is "$scaffold"
    ln -s "$scaffold" scaffold.txt
    
    echo -e "---------- making subplots using scaffold data ----------------"
    if [ -n "${ancestral_genome}" ] ; then
        Rscript ../00_scripts/Rscripts/02.plot_geneSpace.R ancestral_sp
    else
        Rscript ../00_scripts/Rscripts/02.plot_geneSpace.R "$haplo1"
    fi
    
    cd ../

    cp genespace/*pdf 02_results/

    echo -e "\n--------------------\n \tperform whole genome synteny \n------------------------\n" 
    echo -e "\n-------------------- running minimap  ------------------------\n\n" 
    
    minimap2 -cx asm5 \
        haplo1/03_genome/"$haplo1".fa \
        haplo2/03_genome/"$haplo2".fa \
        > 02_results/aln."$haplo1"_"$haplo2".paf || \
        { echo -e "${RED} ERROR! minimap2 failed - check your data\n${NC} " ; exit 1 ; }

    
    #/!\ to do: replace by the simple NO.file!
    if [ -n "${ancestral_genome}" ]; then
        #ancestral genome exist
        ./00_scripts/extract_singlecopy.sh -h1 "$haplo1" -h2 "$haplo2" -s "$scaffold" -a ancestral_sp
    else
        #ancestral genome not provided	
        ./00_scripts/extract_singlecopy.sh -h1 "$haplo1" -h2 "$haplo2" -s "$scaffold" 
    fi


    if [ -n "${ancestral_genome}" ] ; then
        echo -e "\n------- an ancestral genome was provided ------ "
        echo -e "running minimap for genome broad synteny plots  -------\n" 
    
        minimap2 -cx asm5 \
            ancestral_sp/ancestral_sp.fa \
            haplo2/03_genome/"$haplo2".fa \
            > 02_results/aln.ancestral_sp_"$haplo2".paf || \
        { echo -e "${RED} ERROR! minimap2 faield - check your data\n${NC} " ; exit 1 ; }
 
        minimap2 -cx asm5 \
            ancestral_sp/ancestral_sp.fa \
            haplo1/03_genome/"$haplo1".fa \
            > 02_results/aln.ancestral_sp_"$haplo1".paf  || \
        { echo -e "${RED} ERROR! minimap2 faield - check your data\n${NC} " ; exit 1 ; }

    
    
        #preparing scaffold to highlight in dotplot:
        awk '{gsub("_","\t",$0) ; print $2"_"$3"_"$4"\t"$6"_"$7}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>10 '  > 02_results/scaff.anc.haplo1.txt
        awk '{gsub("_","\t",$0) ; print $2"_"$3"_"$4"\t"$9"_"$10}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c\
            |awk '$1>10 ' > 02_results/scaff.anc.haplo2.txt
        awk '{gsub("_","\t",$0) ; print $6"_"$7"\t"$9"_"$10}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>10 ' \
            |sed -e 's/^    //g' -e 's/ /\t/g' > 02_results/scaff.haplo1.haplo2.txt 
        
        #do a clean-up in case there is false positive single copy orthologs:  
        grep -Ff <(cut -f 2 02_results/scaff.haplo1.haplo2.txt ) 02_results/paml/single.copy.orthologs \
            |grep -Ff <(cut -f3 02_results/scaff.haplo1.haplo2.txt ) - > 02_results/paml/single.copy.orthologs_cleaned 


        Rscript 00_scripts/Rscripts/dotplot_paf.R  02_results/aln."$haplo1"_"$haplo2".paf 
        Rscript 00_scripts/Rscripts/dotplot_paf.R  02_results/aln.ancestral_sp_"$haplo1".paf 
        Rscript 00_scripts/Rscripts/dotplot_paf.R  02_results/aln.ancestral_sp_"$haplo2".paf 
    
        Rscript 00_scripts/Rscripts/synteny_plot.R 02_results/aln.ancestral_sp_"$haplo1".paf \
             02_results/scaff.anc.haplo1.txt 
        Rscript 00_scripts/Rscripts/synteny_plot.R 02_results/aln.ancestral_sp_"$haplo2".paf \
            02_results/scaff.anc.haplo2.txt 
        Rscript 00_scripts/Rscripts/synteny_plot.R 02_results/aln."$haplo1"_"$haplo2".paf \
            02_results/scaff.haplo1.haplo2.txt 
    
    else 
        awk '{gsub("_","\t",$0) ; print $2"_"$3"\t"$5"_"$6}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>10 ' \
           |sed -e 's/^    //g' -e 's/ /\t/g'  > 02_results/scaff.haplo1.haplo2.txt
        
        #do a clean-up in case there is false positive single copy orthologs:  
        grep -Ff <(cut -f 2 02_results/scaff.haplo1.haplo2.txt ) 02_results/paml/single.copy.orthologs \
            |grep -Ff <(cut -f3 02_results/scaff.haplo1.haplo2.txt ) - > 02_results/paml/single.copy.orthologs_cleaned 

        #then run pafr to generate a whole genome dotplot and eventually dotplot for some target scaffold:
        Rscript 00_scripts/Rscripts/dotplot_paf.R  02_results/aln."$haplo1"_"$haplo2".paf 
        Rscript 00_scripts/Rscripts/synteny_plot.R 02_results/aln."$haplo1"_"$haplo2".paf \
            02_results/scaff.haplo1.haplo2.txt 
    fi
fi


#optional - to be optimize: 
if [[ $options = "Ds_only" ]] ; then

    if [ -n "${ancestral_genome}" ]; then
        #ancestral genome exist
        ./00_scripts/extract_singlecopy.sh -h1 "$haplo1" -h2 "$haplo2" -s "$scaffold" -a ancestral_sp
        #preparing scaffold to highlight in dotplot:
        awk '{gsub("_","\t",$0) ; print $2"_"$3"_"$4"\t"$6"_"$7}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>4 '  > 02_results/scaff.anc.haplo1.txt
        awk '{gsub("_","\t",$0) ; print $2"_"$3"_"$4"\t"$9"_"$10}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c\
            |awk '$1>4 ' > 02_results/scaff.anc.haplo2.txt
        awk '{gsub("_","\t",$0) ; print $6"_"$7"\t"$9"_"$10}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>4 ' \
            |sed -e 's/^    //g' -e 's/ /\t/g' > 02_results/scaff.haplo1.haplo2.txt

        #do a clean-up in case there is false positive single copy orthologs:  
        grep -Ff <(cut -f 2 02_results/scaff.haplo1.haplo2.txt ) 02_results/paml/single.copy.orthologs \
            |grep -Ff <(cut -f3 02_results/scaff.haplo1.haplo2.txt ) - > 02_results/paml/single.copy.orthologs_cleaned


    else
        #ancestral genome not provided  
        ./00_scripts/extract_singlecopy.sh -h1 "$haplo1" -h2 "$haplo2" -s "$scaffold"
        awk '{gsub("_","\t",$0) ; print $2"_"$3"\t"$5"_"$6}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>4 ' \
           |sed -e 's/^    //g' -e 's/ /\t/g'  > 02_results/scaff.haplo1.haplo2.txt

        #do a clean-up in case there is false positive single copy orthologs:  
        grep -Ff <(cut -f 2 02_results/scaff.haplo1.haplo2.txt ) 02_results/paml/single.copy.orthologs \
            |grep -Ff <(cut -f3 02_results/scaff.haplo1.haplo2.txt ) - > 02_results/paml/single.copy.orthologs_cleaned


    fi
fi

if [ -n "${ancestral_genome}" ]; then
    cut  -f3 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo1".txt
    cut  -f4 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo2".txt
else
    cut  -f2 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo1".txt
    cut  -f3 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo2".txt
fi
#------------------------------ step 3 run paml  -------------------------------------------------------------#

if [[ $options = "synteny_and_Ds" ]] || [[ $options = "Ds_only" ]] ; then 
    #echo haplo1 is "$haplo1"
    #echo haplo2 is "$haplo2"
    echo -e  "\n${BLU}----------------------\npreparing data for paml\n-----------------------${NC}\n"
   
    
    if [ -n "${ancestral_genome}" ]; then
        #ancestral genome exist
        ./00_scripts/12_command_line.paml.sh \
        -h1 "$haplo1" \
        -h2 "$haplo2" \
        -s "$scaffold" \
        -a ancestral_sp || { echo -e "${RED} ERROR! paml failed - check your data\n${NC} " ; exit 1 ; }
    else
        #ancestral genome not provided	
        ./00_scripts/12_command_line.paml.sh \
        -h1 "$haplo1" \
        -h2 "$haplo2" \
        -s "$scaffold" || { echo -e "${RED} ERROR! paml failed - check your data\n${NC} " ; exit 1 ; }
    fi
    
fi

if [[ $options = "synteny_and_Ds" ]]  || [[ $options = "Ds_only" ]] || [[ $options = "plots" ]] ; then 

    pamlsize=$(wc -l 02_results/paml/results_YN.txt |awk '{print $1}' ) 
    scpo=$(wc -l 02_results/paml/single.copy.orthologs_cleaned |awk '{print $1}' )
    
    echo -e "there is $pamlsize results for PAML \n"
    echo -e "there is $scpo single copy orthologs \n" 
    
    #just in case:
    #sed -i 's/  */\t/g' paml/single.copy.orthologs
    
    #----------------------------------- step4 -- plot paml results  -----------------------------------------#
    mkdir  02_results/plots/ 2>/dev/null
    
    if [ -n "${ancestral_genome}" ]; then
        echo "using ancestral genome"
        if ! Rscript ./00_scripts/Rscripts/03.plot_paml.R "$haplo1" "$haplo2" "$scaffold" ancestral_sp 
        then
            echo -e "\nERROR plotting paml results failed\n"
            echo -e "\nplease check input files and logs!!\n\n"
            exit 1
        fi
    else
        if ! Rscript ./00_scripts/Rscripts/03.plot_paml.R "$haplo1" "$haplo2" "$scaffold"
        then
            echo -e "\nERROR plotting paml results failed\n"
            echo -e "\nplease check input files and logs!!\n\n"
            exit 1
        fi
    fi
       
    # ---------------------------------- step5 -- plot ideogram -----------------------------------------------#
    #this part has been done elsewhere and should be removed:
    pathN0="genespace/orthofinder/Results_*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
    haplo=$(head -n1 $pathN0 |awk '{print $7}')

    if [ -n "${ancestral_genome}" ]
    then
     ancestral=$(head -n1 ancestral_sp/ancestral_sp.fa.fai \
            |cut -f1 \
            |awk '{gsub("_","\t",$0) ; print $1}')
        
    echo "ancestral genome ID is $ancestral " 
      
       if [[ $haplo1 == $haplo ]] ;
       then
             echo " " ;
           awk -v var1="$haplo1" -v var2="$haplo2" -v var3="$ancestral" 'NF==6 && $4 ~ var1 && $5 ~ var2 && $6 ~ var3 ' $pathN0 \
           | grep -Ff <(awk '{print $2}' "$scaffold") - > orthologues
       sed -i -e "s/\r//g" orthologues
       else
           echo "not egal - reversing haplotype names to match columns"
           awk -v var1="$haplo2" -v var2="$haplo1" -v var3="$ancestral" 'NF==6 && $4 ~ var1 && $5 ~ var2 && $6 ~ var3 ' $pathN0 \
           | grep -Ff <(awk '{print $2}' "$scaffold") - |awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6}' > orthologues
           sed -i -e "s/\r//g" orthologues
      fi

    else
       #echo "assuming no ancestral species:" 
       if [[ $haplo1 == $haplo ]] ;
       then
           echo " " ;
           awk -v var1="$haplo1" -v var2="$haplo2" 'NF==5 && $4 ~ var1 && $5 ~ var2 ' $pathN0 \
           | grep -Ff <(awk '{print $2}' "$scaffold") - > orthologues
           sed -i -e "s/\r//g" orthologues
       else
           echo "not egal - reversing haplotype names to match columns"
           awk -v var1="$haplo2" -v var2="$haplo1" 'NF==5 && $4 ~ var1 && $5 ~ var2 ' $pathN0 \
           | grep -Ff <(awk '{print $2}' "$scaffold") - |awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' > orthologues
           sed -i -e "s/\r//g" orthologues
       fi
    fi

    #creating different synteny table 
    #note: we already have that with the joint bed from the ideogram 
    #this is redundant 

    if [ -n "${ancestral_genome}" ]
    then
        echo "inferring synteny with ancestral species: "
        join  -1 6 -2 4 <(sort -k6,6 orthologues)  \
                        <(sort -k4,4 genespace/bed/ancestral_sp.bed ) \
            | sed 's/ /\t/g' \
            | join -1 5 -2 4 <(sort -k5,5 -) \
                           <(sort -k4,4 genespace/bed/"$haplo1".bed )  \
            |awk 'NR==1 {print "HOG\tOG\tN0\tchrom1\tGene1\tstart1\tend1\tchrom2\tGene2\tstart2\tend2"}
                    {print $3"\t"$4"\t"$5"\t"$7"\t"$2"\t"$8"\t"$9"\t"$10"\t"$1"\t"$11"\t"$12}' \
            > 02_results/synteny_ancestral_sp_"$haplo1".txt
       
        if [ -s 02_results/synteny_ancestral_sp_"$haplo1".txt ]
        then
            size1=(wc -l   02_results/synteny_ancestral_sp_"$haplo1".txt)
        else
            echo "synteny file between ancestral species and $haplo1 is empty"
            echo "please check your data"
            exit 1
        fi
    
        join  -1 6 -2 4 <(sort -k6,6 orthologues)  \
                <(sort -k4,4 genespace/bed/ancestral_sp.bed ) \
                | sed 's/ /\t/g' \
                |join -1 6 -2 4 <(sort -k6,6 -) \
                                <(sort -k4,4 genespace/bed/"$haplo2".bed ) \
                |awk 'NR==1 {print "HOG\tOG\tN0\tchrom1\tGene1\tstart1\tend1\tchrom2\tGene2\tstart2\tend2"}
                    {print $3"\t"$4"\t"$5"\t"$7"\t"$2"\t"$8"\t"$9"\t"$10"\t"$1"\t"$11"\t"$12}' \
                >  02_results/synteny_ancestral_sp_"$haplo2".txt
    
         if [ -s 02_results/synteny_ancestral_sp_"$haplo2".txt ]
         then
             size2=(wc -l   02_results/synteny_ancestral_sp_"$haplo2".txt)
         else
             echo "synteny file between ancestral species and $haplo2 is empty"
             echo "please check your data"
             exit 1
         fi
       
           echo -e "number of lines in synteny file ancestral_sp vs $haplo1 is $size1"
           echo -e "number of lines in synteny file ancestral_sp vs $haplo2 is $size2"
     
        join  -1 4 -2 4 <(sort -k4,4 orthologues)  \
                        <(sort -k4,4 genespace/bed/"$haplo1".bed ) \
            | sed 's/ /\t/g' \
            | join -1 5 -2 4 <(sort -k5,5 -) \
                           <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
            |awk 'NR==1 {print "HOG\tOG\tN0\tchrom1\tGene1\tstart1\tend1\tchrom2\tGene2\tstart2\tend2"}
                    {print $3"\t"$4"\t"$5"\t"$7"\t"$2"\t"$8"\t"$9"\t"$10"\t"$1"\t"$11"\t"$12}' \
            > 02_results/synteny_"$haplo1"_"$haplo2".txt
    
         if [ -s 02_results/synteny_"$haplo1"_"$haplo2".txt ] ;
         then
             size3=(wc -l   02_results/synteny_"$haplo1"_"$haplo2".txt)
         else
             echo "synteny file between $haplo1 and $haplo2 is empty"
             echo "please check your data"
             exit 1
         fi
                   
         echo -e "number of lines in synteny file $haplo1 vs $haplo2 is $size3"

     else
        echo "no ancestral species assumed "
        echo "inferring synteny between $haplo1 and $haplo2"


       join  -1 4 -2 4 <(sort -k4,4 orthologues)  \
                       <(sort -k4,4 genespace/bed/"$haplo1".bed ) \
            | sed 's/ /\t/g' \
            | join -1 5 -2 4 <(sort -k5,5 -) \
                           <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
            |awk 'NR==1 {print "HOG\tOG\tN0\tchrom1\tGene1\tstart1\tend1\tchrom2\tGene2\tstart2\tend2"}
                    {print $3"\t"$4"\t"$5"\t"$6"\t"$2"\t"$7"\t"$8"\t"$9"\t"$1"\t"$10"\t"$11}'   \
            > 02_results/synteny_"$haplo1"_"$haplo2".txt

         if [ -s 02_results/synteny_"$haplo1"_"$haplo2".txt ] ;
         then
             size3=(wc -l 02_results/synteny_"$haplo1"_"$haplo2".txt)
         else
             echo "synteny file between $haplo1 and $haplo2 is empty"
             echo "please check your data"
             exit 1
         fi
                   
         echo -e "number of lines in synteny file $haplo1 vs $haplo2 is $size3"

    fi

    #/!\ chromosomes should be reconstructed on the fly from the N0.tsv file
    file="02_results/paml/single.copy.orthologs"
    if [ -n "${ancestral_genome}" ]
    then
        awk '{gsub("_g[0-9]*.t[1-9]","\t",$0); print $2"\t"$3"\t"$4}' $file \
           |sort \
           |uniq -c\
           |awk -v var1="$ancestral" -v var2="$haplo1" -v var3="$haplo2" '$1>5 {print var1"\t"$2"\n"var2"\t"$3"\n"var3"\t"$4}' \
           |sort \
           |uniq  \
            > 02_results/chromosomes.txt
    else
        awk '{gsub("_g[0-9]*.t[1-9]","\t",$0); print $2"\t"$3"\t"$4}' $file \
            |sort \
            |uniq -c\
            |awk -v var1="$haplo1" -v var2="$haplo2" '$1>5 {print var1"\t"$2"\n"var2"\t"$3}' \
            |sort \
            |uniq  \
            > 02_results/chromosomes.txt

    fi
    chromosomes="02_results/chromosomes.txt"
    awk '{print $0"\tN"}' $chromosomes > 02_results/chromosomes_orientation.txt
    scafforientation="02_results/chromosomes_orientation.txt"

    #test if previous step was successfull else plot or exit with high levels of pain
    #take advantage of samtools to get length of genome
    samtools faidx haplo1/03_genome/"$haplo1".fa 
    samtools faidx haplo2/03_genome/"$haplo2".fa
    
    eval "$(conda shell.bash hook)"
    conda activate superannot
    #conda deactivate 
    awk '{print $1"\t"$3"\t"$4}' 02_results/paml/single.copy.orthologs_cleaned > 02_results/sco

    source config/config #to get infos on scaffold orientation 
    if [  -n "${links}" ] ; then    
        #links were provided and will be colored
        if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
            -c 02_results/sco \
            -i genespace/bed/"$haplo1".bed \
            -j genespace/bed/"$haplo2".bed  \
            -f haplo1/03_genome/"$haplo1".fa.fai \
            -g haplo2/03_genome/"$haplo2".fa.fai \
            -l "$links" \
            -s "$scafforientation"
        then
                echo -e "\nERROR: ideograms failed /!\ \n
                please check logs and input data\n" 
                exit 1
        fi
    else
        #no links were provided
        if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
            -c 02_results/sco  \
            -i genespace/bed/"$haplo1".bed  \
            -j genespace/bed/"$haplo2".bed  \
            -f haplo1/03_genome/"$haplo1".fa.fai \
            -g haplo2/03_genome/"$haplo2".fa.fai \
            -s "$scafforientation"
        then
                echo -e "\nERROR: ideograms failed /!\ \n
                please check logs and input data\n" 
                exit 1
        fi
    fi

    if [  -n "${ancestral_genome}" ] ; then
        echo -e "ancestral genome was provided for inference" 
        #we will make an ideogram with it 
        awk '{print $1"\t"$2"\t"$3}' 02_results/paml/single.copy.orthologs_cleaned > 02_results/sco_anc	
        if [  -n "${links}" ] ; then    
            if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco_anc \
                -i genespace/bed/ancestral_sp.bed \
                -j genespace/bed/"$haplo1".bed  \
                -f "${ancestral_genome}".fai \
                -g haplo1/03_genome/"$haplo1".fa.fai \
                -l "$links" \
                -s "$scafforientation"
            then
                    echo -e "\nERROR: ideograms failed /!\ \n
                    please check logs and input data\n" 
                    exit 1
            fi

        else
            if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco_anc \
                -i genespace/bed/ancestral_sp.bed  \
                -j genespace/bed/"$haplo1".bed \
                -f "${ancestral_genome}".fai \
                -g haplo1/03_genome/"$haplo1".fa.fai \
                -s "$scafforientation"
            then
                echo -e "\nERROR: ideograms failed /!\ \n
                please check logs and input data\n" 
                exit 1
            fi
        fi

    fi
    
    # ---------------------------------- step6 -- create circos plot --------------------------------
    echo -e "\n~~~~~~~~~~~~~~~contstructing circos plots ~~~~~~~~~~~~~~~~~~~"

    #check if a TEfile exist for genome1 and genome2:
    if [ -s haplo1/03_genome/"$haplo1".TE.bed ] ;
    then 
        genome1TE="haplo1/03_genome/"$haplo1".TE.bed"
        annotateTE = "YES" 
    elif [ -n "$TEgenome1" ] ;
    then 
        genome1TE="$TEgenome1"
        annotateTE = "YES" 
    else
        annotateTE="NO"
    fi
    
    if [ -s haplo2/03_genome/"$haplo2".TE.bed ] ;
    then 
        genome2TE="haplo2/03_genome/"$haplo2".TE.bed"
        annotateTE = "YES" 
    elif [ -n "$TEgenome2" ] ;
    then 
        genome2TE="$TEgenome2"
        annotateTE = "YES" 
    else
        annotateTE="NO"
    fi


    if [ -n "${ancestral_genome}" ] ; then

        echo "ancestral genome was provided" 
        ancestralTE="$ancestralTE" #to be provided in config file #ignored for now
  
        if [ -n "$links" ] ; then
            echo "links file were provided"
            if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$ancestral" -p "$haplo1" \
                -c "$chromosomes" \
                -y 02_results/synteny_ancestral_sp_"$haplo1".txt \
                -f ancestral_sp/ancestral_sp.fa.fai \
                -g haplo1/03_genome/"$haplo1".fa.fai \
                -i genespace/bed/ancestral_sp.bed  \
                -j genespace/bed/"$haplo1".bed  \
                -l "$links"
            then
                echo -e "\nERROR: circos plots failed /!\ \n
                please check logs and input data\n" 
                exit 1
            fi
            if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$ancestral" -p "$haplo2" \
                -c "$chromosomes" \
                -y 02_results/synteny_ancestral_sp_"$haplo2".txt \
                -f ancestral_sp/ancestral_sp.fa.fai \
                -g haplo2/03_genome/"$haplo2".fa.fai \
                -i genespace/bed/ancestral_sp.bed  \
                -j genespace/bed/"$haplo2".bed  \
                -l "$links"
            then
                echo -e "\nERROR: circos plots failed /!\ \n
                please check logs and input data\n" 
                exit 1
            fi
            if [[ $annotateTE = "YES" ]] ; then
                echo "assuming TE bed file exist"
                echo "assuming links"
                #bed file of TE should exist:
                if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplo1" -p "$haplo2" \
                    -c "$chromosomes" \
                    -y 02_results/synteny_"$haplo1"_"$haplo2".txt  \
                    -f haplo1/03_genome/"$haplo1".fa.fai \
                    -g haplo2/03_genome/"$haplo2".fa.fai \
                    -i genespace/bed/"$haplo1".bed  \
                    -j genespace/bed/"$haplo2".bed  \
                    -t "$genome1TE" \
                    -u "$genome2TE" \
                    -l "$links"
                then
                    echo -e "\nERROR: circos plots failed /!\ \n
                    please check logs and input data\n" 
                    exit 1
                fi
            else #assume no TE: 
                echo "assuming no TE bed files"
                echo "assuming links"
                if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplo1" -p "$haplo2" \
                    -c "$chromosomes" \
                    -y 02_results/synteny_"$haplo1"_"$haplo2".txt  \
                    -f haplo1/03_genome/"$haplo1".fa.fai \
                    -g haplo2/03_genome/"$haplo2".fa.fai \
                    -i genespace/bed/"$haplo1".bed  \
                    -j genespace/bed/"$haplo2".bed  \
                    -l "$links"
                then
                    echo -e "\nERROR: circos plots failed /!\ \n
                    please check logs and input data\n" 
                    exit 1
                fi
            fi

        else 
       	    echo "assuming no links"
            echo "assuming no TE"
            echo "plotting between ancestral genome and haplotype1"
            if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$ancestral" -p "$haplo1" \
                -c "$chromosomes" \
                -y 02_results/synteny_ancestral_sp_"$haplo1".txt \
                -f ancestral_sp/ancestral_sp.fa.fai \
                -g haplo1/03_genome/"$haplo1".fa.fai \
                -i genespace/bed/ancestral_sp.bed  \
                -j genespace/bed/"$haplo1".bed  
            then
                echo -e "\nERROR: circos plots failed /!\ \n
                please check logs and input data\n" 
                exit 1
            fi
        fi
            echo "assuming no links"
            echo "assuming no TE"
            echo "plotting between ancestral genome and haplotype2"
            if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$ancestral" -p "$haplo2" \
                -c "$chromosomes" \
                -y 02_results/synteny_ancestral_sp_"$haplo2".txt \
                -f ancestral_sp/ancestral_sp.fa.fai \
                -g haplo2/03_genome/"$haplo2".fa.fai \
                -i genespace/bed/ancestral_sp.bed  \
                -j genespace/bed/"$haplo2".bed  
            then
                echo -e "\nERROR: circos plots failed /!\ \n
                please check logs and input data\n" 
                exit 1
            fi

            if [[ $annotateTE = "YES" ]] ; then
                #bed file of TE should exist:
                echo "assuming no TE"
                echo "plotting between haplotype1  and haplotype2"

               if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplo1" -p "$haplo2" \
                   -c "$chromosomes" \
                   -y 02_results/synteny_"$haplo1"_"$haplo2".txt  \
                   -f haplo1/03_genome/"$haplo1".fa.fai \
                   -g haplo2/03_genome/"$haplo2".fa.fai \
                   -i genespace/bed/"$haplo1".bed  \
                   -j genespace/bed/"$haplo2".bed  \
                   -t "$genome1TE" \
                   -u "$genome2TE"
               then
                   echo -e "\nERROR: circos plots failed /!\ \n
                   please check logs and input data\n" 
                   exit 1
               fi
           else
               if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplo1" -p "$haplo2" \
                   -c "$chromosomes" \
                   -y 02_results/synteny_"$haplo1"_"$haplo2".txt  \
                   -f haplo1/03_genome/"$haplo1".fa.fai \
                   -g haplo2/03_genome/"$haplo2".fa.fai \
                   -i genespace/bed/"$haplo1".bed  \
                   -j genespace/bed/"$haplo2".bed  
               then
                   echo -e "\nERROR: circos plots failed /!\ \n
                   please check logs and input data\n" 
                   exit 1
               fi
            fi

    else
        echo "no ancestral genome" 
        echo "plotting between haplotype1  and haplotype2"
        if [ -n "$links" ] ; then
            echo "links file were provided"
            if [[ $annotateTE = "YES" ]] ; then
                echo "TE bed file are assumed to exist"
                if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplo1" -p "$haplo2" \
                    -c "$chromosomes" \
                    -y 02_results/synteny_"$haplo1"_"$haplo2".txt  \
                    -f haplo1/03_genome/"$haplo1".fa.fai \
                    -g haplo2/03_genome/"$haplo2".fa.fai \
                    -i genespace/bed/"$haplo1".bed  \
                    -j genespace/bed/"$haplo2".bed  \
                    -t "$genome1TE" \
                    -u "$genome2TE" \
                    -l "$links"
                then
                    echo -e "\nERROR: circos plots failed /!\ \n
                    please check logs and input data\n" 
                    exit 1
                fi
            else
                echo "no TE bed file are assumed to exist"
                if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplo1" -p "$haplo2" \
                    -c "$chromosomes" \
                    -y 02_results/synteny_"$haplo1"_"$haplo2".txt  \
                    -f haplo1/03_genome/"$haplo1".fa.fai \
                    -g haplo2/03_genome/"$haplo2".fa.fai \
                    -i genespace/bed/"$haplo1".bed  \
                    -j genespace/bed/"$haplo2".bed  \
                    -l "$links"
                then
                    echo -e "\nERROR: circos plots failed /!\ \n
                    please check logs and input data\n" 
                    exit 1
                fi
            fi
       else
            echo "no links file were provided"
            if [[ $annotateTE = "YES" ]] ; then
                echo "TE bed file are assumed to exist"
                if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplo1" -p "$haplo2" \
                    -c "$chromosomes" \
                    -y 02_results/synteny_"$haplo1"_"$haplo2".txt  \
                    -f haplo1/03_genome/"$haplo1".fa.fai \
                    -g haplo2/03_genome/"$haplo2".fa.fai \
                    -i genespace/bed/"$haplo1".bed  \
                    -j genespace/bed/"$haplo2".bed  \
                    -t "$genome1TE" \
                    -u "$genome2TE" 
                then
                    echo -e "\nERROR: circos plots failed /!\ \n
                    please check logs and input data\n" 
                    exit 1
                fi
            else
                echo "no TE bed file are assumed to exist"
                if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplo1" -p "$haplo2" \
                    -c  "$chromosomes" \
                    -y  02_results/synteny_"$haplo1"_"$haplo2".txt  \
                    -f  haplo1/03_genome/"$haplo1".fa.fai \
                    -g haplo2/03_genome/"$haplo2".fa.fai \
                    -i genespace/bed/"$haplo1".bed  \
                    -j genespace/bed/"$haplo2".bed  
                then
                    echo -e "\nERROR: circos plots failed /!\ \n
                    please check logs and input data\n" 
                    exit 1
                fi
            fi

        fi    
    fi
    #
    if [ $? -eq 0 ]; then
        echo -e  "\n${BLU}------------------\ncircos plot worked successfully------------------${NC}\n"
    else
        echo -e "\n${RED}-------------------\nERROR: circos plots failed /!\ \n
        PLEASE CHECK PACKAGES AND INTPUT DATA------------------${NC}\n"
        exit 1
    fi
    #
    ##---------------------------------- step7 -- run minimap between the genomes -----------------------------#
    #run minimap on the genome 
    #assumption : each genome MUST BE located in folder 03-genome
    
    #------------------------ step 8 -- model comparison -------------------------------------------------#
    echo -e "\n\n~~~~~~~~~~~\n\trunning changepoint\n~~~~~~~~~~~~~~~\n"
    if [[ -d 02_results/modelcomp ]]
    then
        echo -e "WARNING directory modelcomp already exists! check its content first
        Do you wish to remove it?\n
        the data will be lost\n"
        select yn in "Yes" "No"; do
            case $yn in
                Yes ) rm -rf 02_results/modelcomp/ ; if [ -n "$ancestral_genome" ] ; then Rscript 00_scripts/Rscripts/06.MCP_model_comp.R YES else Rscript 00_scripts/Rscripts/06.MCP_model_comp.R NO ; fi  || \
            { echo -e "${RED} ERROR! changepoint failed - check your data\n${NC} " ; exit 1 ; } ;
                break;;
                No ) ; exit;;
            esac
        done
    else
        mkdir 02_results/modelcomp/
       if [ -n "${ancestral_genome}" ] ; then
          Rscript 00_scripts/Rscripts/06.MCP_model_comp.R YES || \
          { echo -e "${RED} ERROR! changepoint failed - check your data\n${NC} " ; exit 1 ; }
       else
          Rscript 00_scripts/Rscripts/06.MCP_model_comp.R NO || \
          { echo -e "${RED} ERROR! changepoint failed - check your data\n${NC} " ; exit 1 ; }
       fi
    fi 
elif [[ $options = "changepoint" ]]
    then

    #eval "$(conda shell.bash hook)"
    #conda activate superannot

    #eventually check its existence - ask user if he wants to remove it
    mkdir 02_results/modelcomp/ 2>/dev/null
    if [[ -d 02_results/modelcomp ]]
    then
        echo -e "WARNING directory modelcomp already exists! check its content first
        Do you wish to remove it?\n
        the data will be lost\n"
        select yn in "Yes" "No"; do
            case $yn in
                Yes ) rm -rf 02_results/modelcomp/ ; if [ -n "$ancestral_genome" ] ; then Rscript 00_scripts/Rscripts/06.MCP_model_comp.R YES else Rscript 00_scripts/Rscripts/06.MCP_model_comp.R NO ; fi  || \
            { echo -e "${RED} ERROR! changepoint failed - check your data\n${NC} " ; exit 1 ; } ;
                break;;
                No ) exit;;
            esac
        done
    else
       if [ -n "${ancestral_genome}" ] ; then
          Rscript 00_scripts/Rscripts/06.MCP_model_comp.R YES || \
          { echo -e "${RED} ERROR! changepoint failed - check your data\n${NC} " ; exit 1 ; }
       else
           Rscript 00_scripts/Rscripts/06.MCP_model_comp.R NO || \
          { echo -e "${RED} ERROR! changepoint failed - check your data\n${NC} " ; exit 1 ; }
       fi

    fi

fi 


for links in 02_results/modelcomp/noprior/classif.s*haplo1.haplo2 ; 
do 
    if [ -n "$scafforientation" ] ; then
        echo "particular orientation provided"
        Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco \
                -i genespace/bed/"$haplo1".bed \
                -j genespace/bed/"$haplo2".bed \
                -f haplo1/03_genome/"$haplo1".fa.fai \
                -g haplo2/03_genome/"$haplo2".fa.fai \
                -l "$links" \
                -s "$scafforientation" 
    else
        Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco \
                -i genespace/bed/"$haplo1".bed \
                -j genespace/bed/"$haplo2".bed \
                -f haplo1/03_genome/"$haplo1".fa.fai \
                -g haplo2/03_genome/"$haplo2".fa.fai \
                -l "$links" 
    fi
done

if [ -n "${ancestral_genome}" ] ; then
   for links in 02_results/modelcomp/noprior/classif.s*ancestral.haplo1 ; 
   do 
      if [ -n "$scafforientation" ] ; then
        echo "particular orientation provided"
        Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco_anc \
                -i genespace/bed/ancestral_sp.bed \
                -j genespace/bed/"$haplo1".bed  \
                -f "$ancestral_genome".fai \
                -g haplo1/03_genome/"$haplo1".fa.fai \
                -l "$links" \
                -s "$scafforientation" 
     else
         Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco_anc \
                -i genespace/bed/ancestral_sp.bed \
                -j genespace/bed/"$haplo1".bed  \
                -f "$ancestral_genome".fai \
                -g haplo1/03_genome/"$haplo1".fa.fai \
                -l "$links" 
     fi
   done
fi

## for fun we now make circos plot with links consisting of the strata and colored by their ds Values:
#preparer des bed file pour faire des circos plots:
mkdir 02_results/bed 2>/dev/null
if [ -n "${ancestral_genome}" ] ;
then
    cut -f 1-3,20 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.4strata.bed
    cut -f 1-3,21 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.5strata.bed
    cut -f 1-3,22 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.6strata.bed
    cut -f 1-3,23 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.7strata.bed
    cut -f 1-3,24 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.8strata.bed
    cut -f 1-3,25 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.9strata.bed

    #haplotype1 bed:
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,20 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo1".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo1".4strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,21 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo1".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo1".5strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,22 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo1".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo1".6strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,23 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo1".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo1".7strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,24 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo1".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo1".8strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,25 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo1".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo1".9strata.bed

    #haplotype2 bed:
      join -1 1 -2 4 <(sort -k1,1 <(cut -f12,20 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".4strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,21 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".5strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,22 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".6strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,23 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".7strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,24 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".8strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,25 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".9strata.bed

    cat 02_results/bed/ancestralspecies.4strata.bed 02_results/bed/"$haplo1".4strata.bed > 02_results/bed/ancestralspecies_"$haplo1".4strata.bed
    cat 02_results/bed/ancestralspecies.5strata.bed 02_results/bed/"$haplo1".5strata.bed > 02_results/bed/ancestralspecies_"$haplo1".5strata.bed
    cat 02_results/bed/ancestralspecies.6strata.bed 02_results/bed/"$haplo1".6strata.bed > 02_results/bed/ancestralspecies_"$haplo1".6strata.bed
    cat 02_results/bed/ancestralspecies.7strata.bed 02_results/bed/"$haplo1".7strata.bed > 02_results/bed/ancestralspecies_"$haplo1".7strata.bed
    cat 02_results/bed/ancestralspecies.8strata.bed 02_results/bed/"$haplo1".8strata.bed > 02_results/bed/ancestralspecies_"$haplo1".8strata.bed
    cat 02_results/bed/ancestralspecies.9strata.bed 02_results/bed/"$haplo1".9strata.bed > 02_results/bed/ancestralspecies_"$haplo1".9strata.bed

    cat 02_results/bed/ancestralspecies.4strata.bed 02_results/bed/"$haplo2".4strata.bed > 02_results/bed/ancestralspecies_"$haplo2".4strata.bed
    cat 02_results/bed/ancestralspecies.5strata.bed 02_results/bed/"$haplo2".5strata.bed > 02_results/bed/ancestralspecies_"$haplo2".5strata.bed
    cat 02_results/bed/ancestralspecies.6strata.bed 02_results/bed/"$haplo2".6strata.bed > 02_results/bed/ancestralspecies_"$haplo2".6strata.bed
    cat 02_results/bed/ancestralspecies.7strata.bed 02_results/bed/"$haplo2".7strata.bed > 02_results/bed/ancestralspecies_"$haplo2".7strata.bed
    cat 02_results/bed/ancestralspecies.8strata.bed 02_results/bed/"$haplo2".8strata.bed > 02_results/bed/ancestralspecies_"$haplo2".8strata.bed
    cat 02_results/bed/ancestralspecies.9strata.bed 02_results/bed/"$haplo2".9strata.bed > 02_results/bed/ancestralspecies_"$haplo2".9strata.bed


else
    #tester si pas d'ancestral
    #haplotype1 bed:
    cut -f1-3,18 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplo1".4strata.bed
    cut -f1-3,19 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplo1".5strata.bed
    cut -f1-3,20 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplo1".6strata.bed
    cut -f1-3,21 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplo1".7strata.bed
    cut -f1-3,22 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplo1".8strata.bed
    cut -f1-3,23 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplo1".9strata.bed
   #haplotype2 bed:
      join -1 1 -2 4 <(sort -k1,1 <(cut -f11,18 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".4strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f11,19 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".5strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f11,20 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".6strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f11,21 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".7strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f11,22 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".8strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f11,23 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 genespace/bed/"$haplo2".bed )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplo2".9strata.bed
fi

cat 02_results/bed/"$haplo1".4strata.bed 02_results/bed/"$haplo2".4strata.bed> 02_results/bed/"$haplo1"."$haplo2".4strata.bed
cat 02_results/bed/"$haplo1".5strata.bed 02_results/bed/"$haplo2".5strata.bed> 02_results/bed/"$haplo1"."$haplo2".5strata.bed
cat 02_results/bed/"$haplo1".6strata.bed 02_results/bed/"$haplo2".6strata.bed> 02_results/bed/"$haplo1"."$haplo2".6strata.bed
cat 02_results/bed/"$haplo1".7strata.bed 02_results/bed/"$haplo2".7strata.bed> 02_results/bed/"$haplo1"."$haplo2".7strata.bed
cat 02_results/bed/"$haplo1".8strata.bed 02_results/bed/"$haplo2".8strata.bed> 02_results/bed/"$haplo1"."$haplo2".8strata.bed
cat 02_results/bed/"$haplo1".9strata.bed 02_results/bed/"$haplo2".9strata.bed> 02_results/bed/"$haplo1"."$haplo2".9strata.bed

if [ -n "${ancestral_genome}" ] ; then
    echo ancestral genome provided 
    for links in 02_results/bed/ancestralspecies$haplo1.*.strata.bed ; do
        if [ -n "${annotateTE}" ] ; then
            echo TE bed provided
            Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$ancestral" -p "$haplo1" \
               -c "$chromosomes" \
               -y 02_results/synteny_ancestral_sp_"$haplo1".txt \
               -f ancestral_sp/ancestral_sp.fa.fai \
               -g haplo1/03_genome/"$haplo1".fa.fai \
               -i genespace/bed/ancestral_sp.bed  \
               -j genespace/bed/"$haplo1".bed  \
               -t "$genome1TE" \
               -u "$genome2TE" \
               -l "$links"
        else 
            echo assuming noTE
            Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$ancestral" -p "$haplo1" \
               -c "$chromosomes" \
               -y 02_results/synteny_ancestral_sp_"$haplo1".txt \
               -f ancestral_sp/ancestral_sp.fa.fai \
               -g haplo1/03_genome/"$haplo1".fa.fai \
               -i genespace/bed/ancestral_sp.bed  \
               -j genespace/bed/"$haplo1".bed  \
               -t "$genome1TE" \
               -u "$genome2TE" \
               -l "$links"
        fi
    done
else
    echo no ancestral genome provided 
fi
#performing haplo1 vs haplo2 comparisons :
for links in 02_results/bed/"$haplo1"."$haplo2".*strata.bed ; do
    if [ -n "${annotateTE}" ] ; then
        echo TE bed provided
        Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplo1" -p "$haplo2" \
           -c "$chromosomes" \
           -y 02_results/synteny_"$haplo1"_"$haplo2".txt \
           -f haplo1/03_genome/"$haplo1".fa.fai\
           -g haplo2/03_genome/"$haplo2".fa.fai \
           -i genespace/bed/"$haplo1".bed  \
           -j genespace/bed/"$haplo2".bed  \
           -t "$genome1TE" \
           -u "$genome2TE" \
           -l "$links"
    else 
        echo assuming noTE
        Rscript 00_scripts/Rscripts/05_plot_circos.R \
           -s "$haplo1"  \
           -p "$haplo2" \
           -c "$chromosomes" \
           -y 02_results/synteny_"$haplo1"_"$haplo2".txt \
           -f haplo1/03_genome/"$haplo1".fa.fai\
           -g haplo2/03_genome/"$haplo2".fa.fai \
           -i genespace/bed/"$haplo1".bed  \
           -j genespace/bed/"$haplo2".bed  \
           -l "$links"
    fi
done

#now we will do the same discretisation of dS for plotting in ideogram: 
if [  -n "${ancestral_genome}" ] ; then
    echo -e "ancestral genome was provided for inference" 
    if [ -n $scafforientation ] ; then
        echo -e "particular orientation will be used"
        #we will make an ideogram with it 
        if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco_anc \
                -i genespace/bed/ancestral_sp.bed \
                -j genespace/bed/"$haplo1".bed  \
                -d 02_results/dS.values.forchangepoint.txt \
                -f "$ancestral_genome".fai \
                -g haplo1/03_genome/"$haplo1".fa.fai \
                -s "$scafforientation" \
        then
             echo -e "\nERROR: ideograms failed /!\ \n
             please check logs and input data\n" 
             exit 1
        fi
    else
        if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco_anc \
                -i genespace/bed/ancestral_sp.bed \
                -j genespace/bed/"$haplo1".bed  \
                -d 02_results/dS.values.forchangepoint.txt \
                -f "$ancestral_genome".fai \
                -g haplo1/03_genome/"$haplo1".fa.fai \
                -s "$scafforientation" 
        then
             echo -e "\nERROR: ideograms failed /!\ \n
             please check logs and input data\n" 
             exit 1
        fi
    fi
fi

if [ -n $scafforientation ] ; then
    echo -e "particular orientation will be used"

    if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
        -c 02_results/sco \
        -i genespace/bed/"$haplo1".bed \
        -j genespace/bed/"$haplo2".bed  \
        -d 02_results/dS.values.forchangepoint.txt \
        -f haplo1/03_genome/"$haplo1".fa.fai \
        -g haplo2/03_genome/"$haplo2".fa.fai \
        -s "$scafforientation" \
    then
        echo -e "\nERROR: ideograms failed /!\ \n
        please check logs and input data\n" 
        exit 1
    fi
else
    if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
        -c 02_results/sco \
        -i genespace/bed/"$haplo1".bed \
        -j genespace/bed/"$haplo2".bed  \
        -d 02_results/dS.values.forchangepoint.txt \
        -f haplo1/03_genome/"$haplo1".fa.fai \
        -g haplo2/03_genome/"$haplo2".fa.fai 
fi
