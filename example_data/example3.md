# Example 3: full workflow without RNAseq one single genome (containing X/Y or two MAT)and no ancestral genome (a.k.a worst case scenario)


## 1 download example data genome

data needed:
 - genome assembly for haplotype1 containing X/Y or two MAT (fasta)
 - database of known TE elements 
 - txt file with the name of the sex chromosome/MAT and their orientation
 - optionally: database of known genes in closely related species


assemblies and annotation are all deposited on zenodo and available here: 

to download: 

`wget /path/to/zenodo/....tar.gz` #(to be updated later)

extract the data in the example_data folder


## 2 - setting your config file

| option in config | description |
| --- | --- |
| *genome1* | ="/path/to/EASYstrata/example_data/fakegenome.fa.gz" |
| *haplotype1* | ="fakegenome" |
| \[*genome2*\] | ="" |
| \[*haplotype2*\] | ="fakegenome_A2" |
| annotate | ="YES" |
| *RelatedProt* | ="" |
| \[*RNAseqlist*\] | ="" |
| \[*bamlist1*\] | ="" |
| \[*bamlist2*\] | ="" |
| \[*orthoDBspecies*\] | ="Fungi" |
| *fungus* | ="YES" |
| *TEdatabase* | ="/path/to/EASYstrata/example_data/TE.fa.gz" |
| *ncbi_species* | ="fungi" |
| \[*gtf1*\] | =""  |
| \[*gtf2*\] | ="" |
| *busco_lineage* | ="basidyomicota_odb10" |
| *interpro* | ="NO" |
| \[*ancestral_genome*\] |  ="" |
| \[*ancestral_gff*\] | ="" |
| *TEancestral* | ="" |
| *scaffolds* | ="/path/to/EASYstrata/example_data/scaffold_example3.txt" |

click [here](/example2.config) to see the full config file example 

## 3 - creating some files: 

use readlink to set full path to text files

## 4 - launching the workflow : 


```./master.sh -o1 2>&1 |tee log```  

## 5 - launching the workflow :  

insert resulting plots and other stats here
