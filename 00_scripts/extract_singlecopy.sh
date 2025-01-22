#!/bin/bash
#Author: QR
#date: 2024

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "master script to: \n extract single copies"
   echo " "
   echo "Usage: $0 [-h1|-h2|-s|-a|-h|]"
   echo "options:"
   echo " -h|--help: Print this Help."
   echo " -h1|--haplo1: the name of the first  focal haplotype\t "
   echo " -h2|--haplo2: the name of the second focal haplotype\t "
   echo " -a |--ancestral_genome: the name of the ancestral haplo to infer orthology and plot gene order"
   echo " -s |--scaffold: the name of scaffold of interest"
   echo " "
}
############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -h1 | --haplo1) haplo1="$2" ; echo -e "haplotype 1 Name is ***${haplo1}*** \n" >&2;;
    -h2 | --haplo2) haplo2="$2" ; echo -e "haplotype 2 Name is ***${haplo2}*** \n" >&2;;
    -a  | --ancestral_genome) ancestral_genome="$2" ; echo -e "ancestral haplo  Name is ***${ancestral_genome}*** \n" >&2;;
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
    echo -e "the ancestral reference name $ancestral_genome will be used"
else
    echo "no ancestral reference genome provided"
    echo "the genome 1 will be used"
fi

#------------------------------ step 1 prepare input files  -------------------------------------#
#cds file:
cdsfile1=haplo1/08_best_run/$haplo1.spliced_cds.fa
cdsfile2=haplo2/08_best_run/$haplo2.spliced_cds.fa

#remove the CDS length info that is introduced by gffread:
sed -i 's/ CDS=.*$//g' "$cdsfile1"
sed -i 's/ CDS=.*$//g' "$cdsfile2"

#-- get single copy orthologs from orthofinder ---
#-- criteria: we want 1:1:1 orthologs between the ancestralhaplo:haplo1:haplo2
scopy="genespace/orthofinder/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt" 

#remove the trailing^M from OrthoFinder:
sed -i -e "s/\r//g" $scopy
sed -i -e "s/\r//g" genespace/orthofinder/Results_*/Orthologues/*/*tsv

mkdir 02_results/paml/ 2>/dev/null
#first we test if an ancestral ref is provided an extract orthologs accordingly:
if [ -n "$ancestral_genome" ] ; then
    echo "using ancestral genome"
    ancestral_vs_hap1="genespace/orthofinder/Results_*/Orthologues/Orthologues_""$ancestral_genome""/""$ancestral_genome""__v__""$haplo1"".tsv"
    ancestral_vs_hap2="genespace/orthofinder/Results_*/Orthologues/Orthologues_""$ancestral_genome""/""$ancestral_genome""__v__""$haplo2"".tsv"

    paste <(grep -Ff $scopy $ancestral_vs_hap1 )  <(grep -Ff $scopy $ancestral_vs_hap2 )  |\
	    grep -Ff <(awk '{print $2}' "$scaffold") - |\
        awk '{ if ($1 == $4) { print $1"\t"$2"\t"$3"\t"$6; } else { print $0"\tdifference exitst -- error"; } }' > 02_results/paml/single.copy.orthologs 
        

#        cut  -f3 paml/single.copy.orthologs > paml/sco."$haplo1".txt
#        cut  -f4 paml/single.copy.orthologs > paml/sco."$haplo2".txt

#if not we extract orthologs like this:
else
    echo "no ancestral genome"
    hap1_vs_hap2="genespace/orthofinder/Results_*/Orthologues/Orthologues_""$haplo1""/""$haplo1""__v__""$haplo2"".tsv"

    paste <(grep -Ff $scopy $hap1_vs_hap2 ) |\
       grep -f <(cut -f 2 "$scaffold" ) - > 02_results/paml/single.copy.orthologs
    sed -i -e "s/\r//g"  02_results/paml/single.copy.orthologs


fi
