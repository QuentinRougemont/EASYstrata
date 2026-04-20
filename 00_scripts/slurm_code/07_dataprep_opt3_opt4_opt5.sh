#!/bin/bash
##SBATCH --account=youraccount
#SBATCH --time=0:50:00
#SBATCH --job-name=dataprep
#SBATCH --output=log_prep-%J.out
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#activate conda
#source /local/env/envconda3.sh 
#activate the env
mamba activate superannot/

#give expected number of arguments from config file:
source config/config
source config/cpu_mem

#PURPOSE: 
#very simple code to prepare the data:#

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
           sed -i -E '/^gtf1/ s/.gz//g' config/config
        else
           echo "$gtf1 is not gzipped"
        fi
        
        #check if gtf or gff:
        if file --mime-type "$gtf1" | grep -q gff; then
            echo file is gff ;
            echo converting into gtf
            gffread "$gtf1" -T -o tmp
            gtf1=${gtf1%.gff*}.gtf
            mv tmp "$gtf1"
            sed -i -E '/gtf1/ s/.gff./.gtf/' config/config
        else
            echo " "
        fi

        if file --mime-type "$gtf2" | grep -q gzip$; then
           echo "$gtf2 is gzipped"
           gunzip "$gtf2"
           gtf2=${gtf2%.gz}
           sed -i -E '/^gtf2/ s/.gz//g' config/config
        else
           echo "$gtf2 is not gzipped"
        fi

        #check if gtf or gff:
        if file --mime-type "$gtf2" | grep -q gff; then
            echo file is gff ;
            echo converting into gtf
            gffread "$gtf2" -T -o tmp
            gtf2=${gtf2%.gff*}.gtf
            mv tmp "$gtf2"
            sed -i -E '/gtf2/ s/.gff./.gtf/' config/config
        else
            echo " "
        fi

        cp "$gtf1" haplo1/08_best_run/"${haplotype1}".final.gtf
        cp "$gtf2" haplo2/08_best_run/"${haplotype2}".final.gtf
        cp "$genome1" haplo1/03_genome/"${haplotype1}".fa
        cp "$genome2" haplo2/03_genome/"${haplotype2}".fa
        if ! gffread -g "$genome1" -x haplo1/08_best_run/"$haplotype1".spliced_cds.fa "$gtf1" 
        then 
            echo "error failed to extract cds from genome $genome1 and gtf $gtf1"
            echo "please check input file synchronisation"
            exit
        else
            gffread -g "$genome1" -y haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa "$gtf1"
            sed -i '/^>/!s/\./*/g' haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa 
        fi     
        if ! gffread -g "$genome2" -x haplo2/08_best_run/"$haplotype2".spliced_cds.fa  "$gtf2"
        then
            echo "error failed to extract cds from genome $genome2 and gtf $gtf2"
            echo "please check input file synchronisation"
            exit
        else
            gffread -g "$genome2" -y haplo2/08_best_run/"$haplotype2"_prot.final.clean.fa  "$gtf2"
            sed -i '/^>/!s/\./*/g' haplo2/08_best_run/"$haplotype2"_prot.final.clean.fa 

        fi 
    elif [ -n "${gtf1}" ] && [ -n "${genome1}" ] ; then
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
        #check if gtf or gff:
        if file --mime-type "$gtf1" | grep -q gff; then
            echo file is gff ;
            echo converting into gtf
            gffread "$gtf1" -T -o tmp
            gtf1=${gtf1%.gff*}.gtf
            mv tmp "$gtf1"
            sed -i -E '/gtf1/ s/.gff./.gtf/' config/config
        else
            echo " "
        fi
        cp "$gtf1" haplo1/08_best_run/"${haplotype1}".final.gtf
        cp "$genome1" haplo1/03_genome/"${haplotype1}".fa

        if ! gffread -g "$genome1" -x haplo1/08_best_run/"$haplotype1".spliced_cds.fa "$gtf1" 
        then 
            echo "error failed to extract cds from genome $genome1 and gtf $gtf1"
            echo "please check input file synchronisation"
            exit
        else
            gffread -g "$genome1" -y haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa "$gtf1"
            sed -i '/^>/!s/\./*/g' haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa 
            #transeq -sequence haplo1/08_best_run/"$haplotype1".spliced_cds.fa \
            #        -outseq haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa

        fi     
        #handle haplotype2 now - we assumme that haplotype2 is present in the genome of "$haplotype1"
        #we will extract proteins for $haplotype2 from the whole $haplotype1
        #we will remove them from $haplotype1 to create a separate dataset for genespace etc:
        mkdir -p haplo2/03_genome
        #linearise genome of haplo1  and extract scaffold to study:       
        if [ -s "haplo1/03_genome/genome.wholemask.fa" ] ; then
            awk '$0~/^>/{if(NR>1){
            print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
            haplo1/03_genome/genome.wholemask.fa | grep -A1 "$haplotype2" > haplo2/03_genome/"$haplotype2".fa
        else 
            awk '$0~/^>/{if(NR>1){
            print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
            haplo1/03_genome/"$haplotype1".fa | grep -A1 "$haplotype2" > haplo2/03_genome/"$haplotype2".fa
        fi
        #extract the gene from the gff :
        mkdir -p haplo2/08_best_run
        grep "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf > haplo2/08_best_run/"$haplotype2".final.gtf
        
        #extract the corresponding protein: 
        awk '$3=="transcript" {print $10}'  haplo2/08_best_run/"$haplotype2".final.gtf \
            |sed -e 's/"//g' -e 's/;//' > haplo2/"$haplotype2".ids

        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
            haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa |\
            grep -Ff  haplo2/"$haplotype2".ids -A1  \
             > haplo2/08_best_run/"$haplotype2"_prot.final.clean.fa

        #extract the corresponding CDS: 
        awk '$0~/^>/{if(NR>1){
        print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' \
        haplo1/08_best_run/"$haplotype1".spliced_cds.fa |\
        grep -Ff haplo2/"$haplotype2".ids -A1  \
        > haplo2/08_best_run/"$haplotype2".spliced_cds.fa
        
        #then we must create a sub file from haplotype1 by excluding haplotype2 from within it:
        cp haplo1/08_best_run/"$haplotype1".final.gtf haplo1/08_best_run/"$haplotype1".final.gtf.bkp
        cp haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        cp haplo1/08_best_run/"$haplotype1".spliced_cds.fa haplo1/08_best_run/"$haplotype1".spliced_cds.fa.bkp
        cp haplo1/03_genome/"$haplotype1".fa haplo1/03_genome/"$haplotype1".fa.bkp 
        
        grep -v "$haplotype2" haplo1/08_best_run/"$haplotype1".final.gtf.bkp > haplo1/08_best_run/"$haplotype1".final.gtf
        #remove seq:
        input=haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa.bkp
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
            else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$input" \
            |\
            grep -vFf haplo2/"$haplotype2".ids |tr "\t" "\n"  > haplo1/08_best_run/"$haplotype1"_prot.final.clean.fa
        
        #remove seq from CDS file:
        input=haplo1/08_best_run/"$haplotype1".spliced_cds.fa.bkp
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
            else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$input" \
            |\
            grep -vFf haplo2/"$haplotype2".ids |tr "\t" "\n"  > haplo1/08_best_run/"$haplotype1".spliced_cds.fa

        #remove from genome
        input=haplo1/03_genome/"$haplotype1".fa.bkp
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
            else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$input" \
            |\
            grep -v "$haplotype2" |tr "\t" "\n"  > haplo1/03_genome/"$haplotype1".fa
         # F#INALLY DEAL WITH ANY EXISTING TE BED FILE HERE:
         if [ -s "$TEgenome1" ] ;
         then
             #TE file provided for genome1 - will try to split into 2
             echo "processing TE file "
             grep "$haplotype2" "$TEgenome1" > haplo2/03_genome/filtered."$haplotype2".TE.bed
             grep -v "$haplotype1" "$TEgenome1" > haplo1/03_genome/filtered."$haplotype1".TE.bed 
         fi
    fi
