#!/bin/bash
#downlaod uniprot and index

#this will be used for later quality checks
if [ -f uniprot/uniprot_sprot.fasta ] ;
then
    echo uniprot already present 

    #check diamond
    if [ -f uniprot/uniprot_sprot.fasta.dmd ] ;
    then
        echo indexing uniprot
        cd uniprot
        diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot.fasta
        cd ../
    else #index:
        echo uniprot already indexed 
    fi

else
    mkdir uniprot 2>/dev/null
    cd uniprot || exit 

    echo -e "\n\n------downloading uniprot--------\n"
    if ! wget -q https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    then
       echo -e "download of uniprot failed\n
       check your internet connexion"   
       exit

    else
       echo -e "\n---- download sucessfull\n make diamond db for later checks ------\n"
       gunzip uniprot_sprot.fasta.gz

       #index:
       diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot.fasta

    fi
fi

