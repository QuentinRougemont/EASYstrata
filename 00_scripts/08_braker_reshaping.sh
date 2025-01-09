#!/bin/bash

#Purpose:
#script to rename the scaffold in the gtf, 
#rename the gtf, create a synchronised genome
#extract longest protein
#Date: 2023
#Author: QR

#set -e

current_command=$BASH_COMMAND
last_command=""

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "script to post-process braker gtf:\n 
       1 - run tsebra \n 
       2 - rename genes id \n
       3 - find longest transcript \n 
       4 - create gtf containing non overlapping transcript \n
       5 - extract protein and CDS fasta"
   echo " "
   echo "Usage: $0 [-s|r|]"
   echo "options:"
   echo "-h|--help: Print this Help."
   echo "-s|--haplo: the name of the focal haplo\t will be used to rename the 
   genes"
   echo "-r|--RNAseq: a float YES/NO stating whether RNAseq was used to 
   annotate the genome "
   echo "-g|--genome: the genome name (full path)"
   echo " "
   echo " "
   echo "dependancies: TSEBRA, samtools, gffread, transeq busco "
}


############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -s | --haplo )  haplo="$2" ; 
    echo -e "haplotype Name is ***${haplo}*** \n" >&2;;
    -g | --genome ) genome="$2"  ;
    echo -e "genome Name is ***${genome}*** \n" >&2;;
    -r | --rnaseq ) RNAseq="$2"  ;
    echo -e "annotation was performed with RNAseq ? ${RNAseq} \n" >&2;;
    -h | --help ) Help ; exit 2 ;;
   esac
   shift
done 


if [ -z "${haplo}" ] || [ -z "${RNAseq}" ] || [ -z "${genome}" ] ; then
    Help
    exit 2 
fi

mkdir -p 08_best_run/01_haplo_cds 2>/dev/null
mkdir 08_best_run/02_haplo_prot 2>/dev/null

#------------------------------ step 1 find best run --------------------------#
echo -e  "\n-----------------------------------------------------------------"
echo -e "\nfinding best run \n" 
echo -e  "-----------------------------------------------------------------\n"
source ../config/config
cd 06_braker

#there is a bit of variance between braker run train on a protein database,
#therefore we perform a few replicate run and choose the 'best run' where 
#best is defined as the run providing 
#the highest score in terms of busco completness 
#we sort on k3 (highest completeness), k11 (lowest amount of missing), 
#k9 (lowest amount of fragmented genes)  
best_round=$(grep "C:" round*_braker_on_refprot/busco_augustus*/short_summary.specific.*.busco_augustus*.txt |\
    sed -e 's/%//g' -e 's/\[/,/g' -e 's/]//g' -e 's/:/\t/g' -e  's/,/\t/g' |\
    LC_ALL=C sort -k3nr -k11n -k9n -k7n -k5  |head -n 1 |cut -d "/" -f 1 ) 
#alternative sorting: 
#LC_ALL=C sort  -k11 -k3 -n -k5 -k7 -k9  |tail -n 1 |cut -d "/" -f 1 ) 

echo -e "best_round is $best_round\n------------------------------------------"

cd ../

#optionally make a report:
#------------- CONDA ACTIVATION  -------------- #
eval "$(conda shell.bash hook)"
conda activate superannot

python3 ../00_scripts/utility_scripts/generateReport.py \
    06_braker/"$best_round"/braker.gtf \
    06_braker/"$best_round"/hintsfile.gff  \
    08_best_run/report_"$haplo"_"$best_round".pdf

conda deactivate
cp  08_best_run/report_"$haplo"_"$best_round".pdf ../02_results/

#------------------------------ step 2 finding runs ---------------------------#

if [[ $RNAseq = "YES" ]]
then

    #make a report on rnaseq:
    #------------- CONDA ACTIVATION  -------------- #
    eval "$(conda shell.bash hook)"
    conda activate superannot

    python3 ../00_scripts/utility_scripts/generateReport.py \
        06_braker/rnaseq/braker.gtf \
        06_braker/rnaseq/hintsfile.gff  \
        08_best_run/report_"$haplo"_rnaseq.pdf
    
    conda deactivate

cp  08_best_run/report_"$haplo"_rnaseq.pdf ../02_results

    echo -e "\nrunning tsebra\n" 
    
    #2 -- run tsebra
    #test if default.cfg can be found:
    tsebraconf=default.cfg
    cp ../config/$tsebraconf .
    if [ -f $tsebraconf ] ; then 
        echo " $tsebraconf exist";  cp "$tsebraconf" . 
    else 
        echo "FATAL ERROR: no config file for tsebra"
        exit
    fi
    
    #then run tsebra:
    ../00_scripts/09_tsebra.sh "$haplo" "$best_round"
    
    cd 08_best_run
    
    ln -s ../07-tsebra_results/"$haplo".combined.gtf "$haplo".tmp.gtf
    cd ../
    file=08_best_run/$haplo.tmp.gtf
    nb_genes=$(awk '$3=="gene" ' "$file" |wc -l)
    echo -e "\nthere is $nb_genes genes after running tsebra\n" 

else 
    #2 -- copy best run based on database in a final folder
    echo "running on protein only - not running tsebra" 
    
    #mkdir 08_best_run
    cp 06_braker/"$best_round"/braker.gtf 08_best_run/"$haplo".tmp.gtf
    file=08_best_run/$haplo.tmp.gtf
    
    nb_genes=$(awk '$3=="gene" ' "$file" |wc -l)
        echo -e "\nthere is $nb_genes genes in the best protein run\n"
fi

#--------------------------- step 3 ----------------------------------------#
echo -e  "\n----------------------------------------------------------------"
echo -e  "\n-----------renaming and fixing braker output now------------\n" 
echo -e  "--------------------------------------------------------------\n"

#then this is a part common to both RNAseq + Proteins or Proteins only:
rename_gtf.py --gtf "${file}" --prefix "${haplo}" \
    --translation_tab translation_tab."$haplo" \
    --out 08_best_run/"${haplo}".renamed.gtf 

# Fix TSEBRA output (source: https://github.com/Gaius-Augustus/BRAKER/issues/457 )
# Fix lack of gene_id and transcript_id tags in gtf file column 9
../00_scripts/utility_scripts/Fix_Augustus_gtf.pl \
    08_best_run/"${haplo}".renamed.gtf \
    > 08_best_run/"${haplo}".renamed.fixed.gtf

#---------------------------- step 4 ----------------------------------------#
#rename the gene/cds/transcript/exons/ ID in the gtf so they contain 
#the chromosome and haplo name (easier to parse)
#1 - declare the gtf:
cd 08_best_run/
gtf=${haplo}.renamed.fixed.gtf #current gtf

#2 - renaming with a simple command: 
../../00_scripts/utility_scripts/01.recode_braker_output.py "${gtf}" "${haplo}"
cd ../

gtf=${haplo}.IDchecked.gtf

#------------------- step 5 extracting prot and cds from new files  ----------#
echo -e  "\n-----------------------------------------------------------------"
echo "extract protein and cds from the renamed gtf" 
echo -e  "-----------------------------------------------------------------\n"

gtffull=08_best_run/$gtf

#extract the cds and protein from it:
output="$haplo"_cds.fa
echo output cds is "$output"
gffread -w 08_best_run/01_haplo_cds/"$output" \
        -g "${genome}" "$gtffull"

#then convert also the file to its cds:
echo "translate CDS into amino acid "
transeq -sequence 08_best_run/01_haplo_cds/"$output" \
        -outseq 08_best_run/02_haplo_prot/"$haplo".prot

echo -e "\n-----------------------------------------------------------------"
echo "extract longest transcript" 
echo -e  "-----------------------------------------------------------------\n"

cd 08_best_run/02_haplo_prot
#assumption :transcript finishes by ".t1, .t2, .t3 the dot (.) is the delimiter

awk '/^>/ {if (seqlen){print seqlen}; 
    printf(">%s\t",substr($0,2)) ;seqlen=0;next; } 
    { seqlen += length($0)}END{print seqlen}' "${haplo}".prot |\
    awk -F ".t[0-9]_1 " '{print $1"\t"$0}'  |\
    awk '$4>max[$1]{max[$1]=$4; row[$1]=$2} END{for (i in row) print row[i]}' \
        > longest.transcript.tmp

#linearize file so that the next command will work:
awk '$0~/^>/{if(NR>1){
    print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
        "$haplo".prot > "$haplo".prot.lin.fasta


#create a list of wanted busco:
echo busco_lineage are $busco_lineage # from config file

table=../../06_braker/"$best_round"/busco_augustus/run_"$busco_lineage"/full_table.tsv
#recover some busco gene that are lost based on gene length:

awk '$2 =="Complete" || $2 =="Duplicated" {print $1"\t"$2"\t"$3}' "$table" \
	| grep -Ff <(cut -d "_" -f3 longest.transcript.tmp) - \
    | grep -vf - <(awk '$2 =="Complete" || $2 =="Duplicated" {print $1"\t"$2"\t"$3}' "$table") \
    | grep "Complete" \
    | awk '{print $3}' \
    | cat longest.transcript.tmp - > all.transcripts

grep -A1 -Ff all.transcripts "$haplo".prot.lin.fasta |\
    sed '/^--/d' > \
    "$haplo".longest_transcript.fa


source ../../../config/config

eval "$(conda shell.bash hook)"
conda activate busco571
busco -c8 -o busco_check -i "$haplo".longest_transcript.fa -l "$busco_lineage" -m protein -f  || \
          { echo -e "${RED} ERROR! busco failed - check your data\n${NC} " ; exit 1 ; }

#now we add the dup:
grep -Ff busco_check/run_"$busco_lineage"/missing_busco_list.tsv "$table" \
	|grep -v "Missing\|#" \
	|cut -f3 \
	|grep -Ff - "$haplo".prot.lin.fasta \
	|awk '{gsub(/^>/,""); print $1}' \
	| cat - all.transcripts > all.transcripts2

#"
grep -A1 -Ff all.transcripts2 "$haplo".prot.lin.fasta \
    sed '/^--/d' > \
    "$haplo".longest_transcript.fa

busco -c8 -o busco_check2 -i "$haplo".longest_transcript.fa -l "$busco_lineage" -m protein -f  || \
          { echo -e "${RED} ERROR! busco failed - check your data\n${NC} " ; exit 1 ; }

if [[ $RNAseq = "YES" ]]
then
    buscorna=../../06_braker/rnaseq/busco_augustus/run_"$busco_lineage"/full_table.tsv
    awk '$2 =="Complete" || $2 =="Duplicated" {print $1"\t"$2"\t"$3}'  "$buscorna" > busco.rna
    grep -Ff longest.transcript.tmp busco.rna > busco.longest.transcript
    grep -vf busco.longest.transcript buso.rna |grep "Complete" |awk '{print $3}' > list.of.missing_busco
    cat longest.transcript.tmp list.of.missing_busco > all.ids
    grep -A1 -Ff all.ids "$haplo".prot.lin.fasta |sed '/--/d' > "$haplo".longest_transcript.fa
    busco -c8 -o busco_check3 -i "$haplo".longest_transcript.fa -l "$busco_lineage" -m protein -f || \
          { echo -e "${RED} ERROR! busco failed - check your data\n${NC} " ; exit 1 ; }
fi 

cd ../../

#~~~~~~~~~ step 6 : cleaning the GTF based on non-overlapping transcript ~~~~~~~"
echo -e "\n-----------------------------------------------------------------"
echo "remove redundant Gene in the CDS and Protein fasta files"
echo -e "-----------------------------------------------------------------\n"

cd 08_best_run

#look for putative fragented gene (this should be zero):
cut -f 9 "$gtf" |\
    sed 's/^gene_id//g' |\
    sed 's/transcript_id//g' |\
    cut -d ";" -f 1 |\
    sed 's/.t[1-9]//g' |\
    sort |uniq -c |\
    awk '$1==1 {print}'  > fragmented_gene.txt 

loss=$(wc -l fragmented_gene.txt |awk '{print $1}' )
echo " there is $loss fragmented gene"


#  declare full path to input:
protpath=02_haplo_prot
prot="$protpath"/"$haplo".longest_transcript.fa
echo -e "there is $(wc -l "$gtf" |awk '{print $1}') lines in ""$gtf"" " 

# subset our gtf to keep only the cds in the cds files!
#we keep the cds:
#now the genes:
#grep -f wanted.gene.tmp  <(awk '$3=="gene" ' "$gtf" )  > p3
#ideally I want to sort on the gene, then transcript, then CDS, then exon as well
#concatenate gene and cds/etc:
#cat p1 p3 |LC_ALL=C sort -k1,1 -k4,4n -k5,5n > "$haplo".longest_transcript.gtf


cat <( grep -Ff <(grep ">" "$prot" \
    |sed 's/>//g' \
    |sed 's/_1 CDS=.*//g'  ) "$gtf" ;
 grep -f <(grep ">" "$prot" \
    |sed 's/>//g' \
    |sed 's/_1 CDS=.*//g' \
    |sed 's/.t[0-9]$//g' ) <(awk '$3=="gene" ' "$gtf" )) \
    |sort -k1,1 -k4,4n -k5,5n > "$haplo".longest_transcript.gtf

echo -e "there is $(wc -l "$haplo".longest_transcript.gtf |\
    awk '{print $1}' ) lines in ""$haplo"".longest_transcript.gtf" 

#declare new gtf for new work:
gtf="$haplo".longest_transcript.gtf
gtf2="gtf.tmp"                                      
gtf3="$haplo".no_overlap.gtf
gtf4="$haplo".final.gtf
# now removing fully overlapping CDS with different IDs (<<1% of the data)
# rule  adopted here: 
# we removed any CDS with identidical CDS chr-start or identical CDS chr-end
# identical CHR-START:
awk '$3=="transcript" {print $1"_"$4"_"$5"\t"$1"_"$4"\t"$1"_"$5"\t"$0}'  "$gtf" > tmp
awk 'NR == FNR {count[$2]++; next} count[$2]>0 {print $2"\t"$7"\t"$8"\t"$13"\t"$8-$7}' tmp tmp |\
    awk '$5>max[$1]{max[$1]=$5; row[$1]=$4} END{for (i in row) print row[i]}' > longest.to.keep.tmp                   
grep -Ff longest.to.keep.tmp "$gtf"  > "$gtf2"                                                               

# identical CHR-END:
awk '$3=="transcript" {print $1"_"$4"_"$5"\t"$1"_"$4"\t"$1"_"$5"\t"$0}'  "$gtf2" > tmp2
awk 'NR == FNR {count[$3]++; next} count[$3]>0 {print $3"\t"$7"\t"$8"\t"$13"\t"$8-$7}' tmp2 tmp2 |\
    awk '$5>max[$1]{max[$1]=$5; row[$1]=$4} END{for (i in row) print row[i]}' > longest.to.keep.tmp2

grep -Ff longest.to.keep.tmp2 "$gtf2"  > "$gtf3"

grep -Ff <( awk '$3=="transcript" {print $10} ' "$gtf3" |sed 's/.t[1-9]//') <(awk '$3=="gene" ' "gtf" ) \
    |cat - "$gtf3" |LC_ALL=C sort -k1,1 -k4,4n -k5,5n  > "$gtf4"

echo -e "there is $(wc -l "$gtf4"|awk '{print $1}' ) lines in ""$gtf4"" (final gtf)"

#~~~~~~~~~~~~~~~~~~~~~~ step 7 : re-extracting the proteins ~~~~~~~~~~~~~~~~~~~#
echo -e "\n-----------------------------------------------------------------"
echo "re-extracting protein and CDS from the final non-redundan gtf" 
echo -e "-----------------------------------------------------------------\n"

gffread -w "$haplo".spliced_cds.fa -g ../"${genome}" "$gtf4" 
echo "translate CDS into amino acid "
transeq -sequence "$haplo".spliced_cds.fa \
    -outseq "$haplo"_prot.final.fa
transeq -clean -sequence "$haplo".spliced_cds.fa \
    -outseq "$haplo"_prot.final.clean.fa #for interproscan and other pipelines

echo -e "there is $( grep -c ">" "$haplo"_prot.final.fa |\
    awk '{print $1}' ) total protein corresponding to a single longest transcript in the final files"

#rm p1 p3
rm ./*tmp*
rm ./*renamed.*gtf
rm ./*IDchecked.gtf
#now run busco to validate. the score should be close to the initial score 
#insert call to busco here

source ../../config/config

eval "$(conda shell.bash hook)"
conda activate busco571
busco -c8 -o busco_final -i "$haplo"_prot.final.fa -l "$busco_lineage" -m protein -f   || \
          { echo -e "${RED} ERROR! busco failed - check your data\n${NC} " ; exit 1 ; }

#then launch quality check on the final dataset: 
chmod +x ../../00_scripts/quality.check.sh

#note: maybe this could be an option
echo -e "running quality checks now "
../../00_scripts/quality.check.sh -s "$haplo"  || \
          { echo -e "${RED} ERROR! quality check failed! - check your data\n${NC} " ; exit 1 ; }

# copy things : 
cp "$gtf4" ../../02_results
cp "$haplo"_prot.final.clean.fa ../../02_results
cp "$haplo".spliced_cds.fa ../../02_results
cp busco_final/short_summary*.txt ../../02_results/busco."$haplo".txt

#to do: copy other stuff to 02_results general folder (quality, busco summary, etc)

# checking if busco score is ok:
score=$(grep "C:" busco_final/short_summary.*.txt  |cut -d ":" -f 2 |sed 's/%.*//' )
target=85

if [ "$(echo "$score < $target" | bc)" -eq 1 ] ;
then
    echo "stop !! busco score too low for further analyses !!" ;
    echo "please check your data "
    exit 1
fi
