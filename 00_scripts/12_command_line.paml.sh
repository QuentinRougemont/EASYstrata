#!/bin/bash

#Author: QR
#date: 2024
#purpose: compute dS, dN, and their SE based on paml using yn00 model
#it will perform some checks and if all is ok will run muscle from TranslatorX first

#run: ./12.command_line.paml.sh -h1 haplo1.cds.fa \
    #-h2 haplo2.cds.fa \
    #-s scaffold.txt \
    #-a ancestral_genome 2>&1 |tee log

# -- some colors for warnings --:
RED='\033[0;31m'
NC='\033[0m' # No Color


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "master script to: \n run translatorX and paml "
   echo " "
   echo "Usage: $0 [-h1|-h2|-s|-a|-h|]"
   echo "options:"
   echo " -h|--help: Print this Help."
   echo " -h1|--haplo1: the name of the first  focal haplotype "
   echo " -h2|--haplo2: the name of the second focal haplotype "
   echo " -a|--ancestral_genome: the name of the ancestral haplo to infer orthology and plot gene order"
   echo " -s|--scaffold: the name of scaffold of interest"
   echo " "
   echo "dependancies: paml (yn00), translatorX "
}

############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -h1 | --haplo1) haplo1="$2" ; echo -e "haplotype 1 Name is ***${haplo1}*** \n" >&2;;
    -h2 | --haplo2) haplo2="$2" ; echo -e "haplotype 2 Name is ***${haplo2}*** \n" >&2;;
    -a  | --ancestral_genome) ancestral_genome="$2" ; 
        echo -e "ancestral haplo  Name is ***${ancestral_genome}*** \n" >&2;;
    -s  | --scaffold) scaffold="$2"  ; echo -e "target scaffold is ${scaffold} \n" >&2;;
    -h  | --help ) Help ; exit 2 ;;
   esac
   shift
done 

if [ -z "${haplo1}" ] && [ -z "${haplo2}" ]  && [ -z "${scaffold}" ]    ; then
    Help
    exit 2
fi

#test if ancestral genome is provided or not:
if [ -n "$ancestral_genome" ] ; then
    echo -e the ancestral reference name "$ancestral_genome"  will be used
else
    echo "no ancestral reference genome provided"
    echo "the genome 1 will be used"
fi

#------------------------------ step 1 prepare input files  -------------------------------------#
#cds file:
cdsfile1=haplo1/08_best_run/$haplo1.spliced_cds.fa
cdsfile2=haplo2/08_best_run/$haplo2.spliced_cds.fa

#remove the CDS length info that is introduced by gffread:
#sed -i 's/ CDS=.*$//g' $cdsfile1
#sed -i 's/ CDS=.*$//g' $cdsfile2
#
##-- get single copy orthologs from orthofinder ---
##-- criteria: we want 1:1:1 orthologs between the ancestralhaplo:haplo1:haplo2
##single copy path:
#scopy=$(echo "genespace/orthofinder/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt" ) 
#
##remove the trailing^M from OrthoFinder:
#sed -i -e "s/\r//g" $scopy
#sed -i -e "s/\r//g" genespace/orthofinder/Results_*/Orthologues/*/*tsv
#
#
##first we test if an ancestral ref is provided an extract orthologs accordingly:
#if [ -n "$ancestral_genome" ] ; then
#    echo "using ancestral genome"
#    ancestral_vs_hap1=$(echo "genespace/orthofinder/Results_*/Orthologues/Orthologues_"$ancestral_genome"/"$ancestral_genome"__v__"$haplo1".tsv ")
#    ancestral_vs_hap2=$(echo "genespace/orthofinder/Results_*/Orthologues/Orthologues_"$ancestral_genome"/"$ancestral_genome"__v__"$haplo2".tsv ")
#    #sed -i -e "s/\r//g" $ancestral_vs_hap1
#    #sed -i -e "s/\r//g" $ancestral_vs_hap2
#
#    paste <(grep -Ff "$(echo $scopy )" "$(echo $ancestral_vs_hap1 )" )  <(grep -Ff "$(echo $scopy )" "$(echo $ancestral_vs_hap2 )" )  |\
#       grep -Ff <(awk '{print $2}' $scaffold) - |\
#        awk '{ if ($1 == $4) { print $1"\t"$2"\t"$3"\t"$6; } else { print $0"\tdifference exitst -- error"; } }' > paml/single.copy.orthologs 
#        
#        #sed -i -e "s/\r//g"  paml/single.copy.orthologs
#
#        cut  -f3 paml/single.copy.orthologs > paml/sco.$haplo1.txt
#        cut  -f4 paml/single.copy.orthologs > paml/sco.$haplo2.txt
#
##if not we extract orthologs like this:
#else
#    echo "no ancestral genome"
#    hap1_vs_hap2=$(echo "genespace/orthofinder/Results_*/Orthologues/Orthologues_"$haplo1"/"$haplo1"__v__"$haplo2".tsv ")
#    sed -i -e "s/\r//g"  paml/single.copy.orthologs
#
#    paste <(grep -Ff "$(echo $scopy )" "$(echo $hap1_vs_hap2)" ) |\
#        grep -f <(cut -f 2 $scaffold ) - > paml/single.copy.orthologs  
#       
#        #sed -i -e "s/\r//g"  paml/single.copy.orthologs
#
#        cut  -f2 paml/single.copy.orthologs > paml/sco.$haplo1.txt
#        cut  -f3 paml/single.copy.orthologs > paml/sco.$haplo2.txt
#
#fi

##---  linearise the cds file ---#
awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "${cdsfile2}"  > 02_results/paml/"$haplo2".linearised.cds
awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "${cdsfile1}"  > 02_results/paml/"$haplo1".linearised.cds


##---- recover the wanted sequences in the CDS file #
cd 02_results/paml || exit 1

#all is run from paml folder now

rm sorted* 2>/dev/null
while read -r pattern ; 
do 
    grep -w -A1 "$pattern" "$haplo2".linearised.cds  >> sorted."$haplo2".wanted_cds.fa ;
done < sco."$haplo2".txt 
#
while read -r pattern ; 
do 
    grep -w -A1 "$pattern" "$haplo1".linearised.cds >> sorted."$haplo1".wanted_cds.fa ; 
done < sco."$haplo1".txt 
#

### ---- then run paml and dnds.sh  ----- #
#paml does not like long fasta name so we shorten them for the sake of computation 
#to do here

#------ part specific to PAML --------------------- #
fasta1=sorted."$haplo1".wanted_cds.fa  
fasta2=sorted."$haplo2".wanted_cds.fa

f1=$(basename "$fasta1" )
f2=$(basename "$fasta2" )
newf1=${f1%.fa**}.nostopcodon.fasta
newf2=${f2%.fa**}.nostopcodon.fasta

##1 ----- remove stop codon -------
awk '{if($1 !~ /^>/) {if( substr($0, length($0)-2, length($0)) ~ /(TGA|TAG|TAA)/){ print substr($0, 1, length($0)-3)} else  {print $0 }}else{print $0}}' "${fasta1}" > "${newf1}" 
awk '{if($1 !~ /^>/) {if( substr($0, length($0)-2, length($0)) ~ /(TGA|TAG|TAA)/){ print substr($0, 1, length($0)-3)} else  {print $0 }}else{print $0}}' "${fasta2}" > "${newf2}"

##2 ------ split and cat pairwise sequence -------
#split 
rm -rf sequence_files 2>/dev/null
rm wanted_sequence 2>/dev/null
rm ID1 ID2 2>/dev/null

mkdir sequence_files

grep ">"  "$newf1" > ID1
grep ">"  "$newf2" > ID2

#check that all gene names are below 32 characters otherwise paml will fail!
awk '{print $0"\t"length }' ID1 |awk '$2>32 {print $1}' > long_geneID.hap1
awk '{print $0"\t"length }' ID2 |awk '$2>32 {print $1}' > long_geneID.hap2

#if some gene have length above we rename them using a structure of the type: gene$id
#id is a seq from 1 to n with n the max number of gene to rename

if [ -s long_geneID.hap1 ] ; then
    #the file is not empty; so we will rename the "long" genes:
    #0 - print some important warning as this may affect the user expectation: 
    nb_genes=$(wc -l long_geneID.hap1 |awk '{print $1}' )
    echo -e "${RED}!!! warning !!!\n some gene names are too long! ${NC} \n 
    a total of $nb_genes genes names will be renamed\n
    you'll find their name in the file:\n
    correspondance.table.hap1.txt\n\n" 


    #1 - create a correspondance table:
    #j=0 ; for i in $(cat long_geneID.hap1) ; do j=$(( "$j" + 1)) ; echo -e "$i\t>gene.$j"  >> correspondance.table.hap1.txt; done
    j=0 ; 
    while read -r line ; 
    do
        j=$(( "$j" + 1)) ; 
        echo -e "$line\t>hap1_gene.$j"  >> correspondance.table.hap1.txt; 
    done < long_geneID.hap1

    #2 - keep a copy:
    oldf1=$newf1.original.genename.fa
    cp "$newf1" "$oldf1"

    #3 - rename the fasta with awk
    awk 'NR==1 { next } FNR==NR { a[$1]=$2; next } $1 in a { $1=a[$1] }1' \
        correspondance.table.hap1.txt "${oldf1}" > "${newf1}" 

    #4 - re-extract the corrected ID for paml to succeed!
    grep ">"  "$newf1" > ID1

fi

#do the same for haplotype2: 
if [ -s long_geneID.hap2 ] ; then
    #the file is not empty; so we will rename the "long" genes:
    #0 - print some important warning as this may affect the user expectation: 
    nb_genes=$(wc -l long_geneID.hap2 |awk '{print $1}' )
    echo -e "${RED}!!! warning !!!\n some gene names are too long! ${NC} \n 
    a total of $nb_genes genes names will be renamed\n
    you'll find their name in the file:\n
    correspondance.table.hap2.txt\n\n" 

    #1 - create a correspondance table:
    j=0 
    while read -r line ; 
    do
        j=$(( "$j" + 1)) ; 
        echo -e "$line\t>hap2_gene.$j"  >> correspondance.table.hap2.txt; 
    done < long_geneID.hap2
 
    #2 - keep a copy:
    oldf2=$newf2.original.genename.fa
    cp "$newf2" "$oldf2"

    #3 - rename the fasta with awk
    awk 'NR==1 { next } FNR==NR { a[$1]=$2; next } $1 in a { $1=a[$1] }1' \
        correspondance.table.hap2.txt "${oldf2}" > "${newf2}"

    #4 - re-extract the corrected ID for paml to succeed!
    grep ">"  "$newf2" > ID2

fi

paste ID1 ID2 |sed 's/>//g' > wanted_sequence

#---- this is the only interesting part that is actually making stuff: 
#---  create architecture run muscle and paml : 

while IFS=$'\t' read -r -a line
do
    mkdir sequence_files/tmp."${line[0]}".vs."${line[1]}"
    
    grep -w -A1 "${line[0]}"  "$newf1" > sequence_files/tmp."${line[0]}".vs."${line[1]}"/sequence.fasta
    grep -w -A1 "${line[1]}"  "$newf2" >> sequence_files/tmp."${line[0]}".vs."${line[1]}"/sequence.fasta
    
    #run muscle from within translatorX, so we also have gblocks output and the html files:
    translatorx_vLocal.pl -i sequence_files/tmp."${line[0]}".vs."${line[1]}"/sequence.fasta \
        -o sequence_files/tmp."${line[0]}".vs."${line[1]}"/results 2>&1 |tee log.translator
    
    cp ../../config/yn00_template.ctl sequence_files/tmp."${line[0]}".vs."${line[1]}"/
    
    cd sequence_files/tmp."${line[0]}".vs."${line[1]}"/ || exit 1
    
    #configure paml:
    path=$(pwd)
    #echo "$path"
    sed -i "s|PATH|$path|g" yn00_template.ctl #"
    
    #run paml :
    yn00 yn00_template.ctl
    
        cd ../../
    
    #extract the dS, dN and SE from the output: 
    awk '/\+\-/ && !/(dS|SE)/ {split(FILENAME, a, "."); 
        print $(NF-2), $(NF), $(NF-5), $(NF-3),"'${line[0]}'","'${line[1]}'"}'  \
            sequence_files/tmp."${line[0]}".vs."${line[1]}"/out_yn00_orthogp \
            >  sequence_files/tmp."${line[0]}".vs."${line[1]}"/resultat_Yang_Nielsen_2000_method.orthogp.txt 


done < wanted_sequence 2>&1 |tee log.paml 

#we concatenate everyone to work with them in the next scripts:
rm results_YN.txt 2>/dev/null   
cat sequence_files/tmp.*/resultat_Yang_Nielsen_2000_method.orthogp.txt >> results_YN.txt
cp results_YN.txt results_YN.txt.bkp


if [ -e correspondance.table.hap1.txt ] && [ ! -e correspondance.table.hap2.txt ] ; then
   sed -i 's/>//g' correspondance.table.hap1.txt
   awk 'NR==FNR{a[$2]=$1;next}$5 in a{$5=a[$5]}1' correspondance.table.hap1.txt results_YN.txt > tmp
   mv tmp results_YN.txt
fi

if [ -e correspondance.table.hap2.txt ] && [ -e correspondance.table.hap1.txt ] ; then
   sed -i 's/>//g' correspondance.table.hap1.txt
   awk 'NR==FNR{a[$2]=$1;next}$5 in a{$5=a[$5]}1' correspondance.table.hap1.txt results_YN.txt > tmp

   sed -i 's/>//g' correspondance.table.hap2.txt
   awk 'NR==FNR{a[$2]=$1;next}$6 in a{$6=a[$6]}1' correspondance.table.hap2.txt tmp > results_YN.txt
fi

if [ -e correspondance.table.hap2.txt ] && [ ! -e correspondance.table.hap1.txt ]  ; then
   sed -i 's/>//g' correspondance.table.hap2.txt
   awk 'NR==FNR{a[$2]=$1;next}$6 in a{$6=a[$6]}1' correspondance.table.hap2.txt results_YN.txt > tmp
   mv tmp results_YN.txt
fi
