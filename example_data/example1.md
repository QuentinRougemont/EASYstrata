
# Example 1: full workflow with RNAseq and ancestral genome

## 1 download example data genome

data needed:
 - genome assembly for haplotype1 (fasta)
 - genome assembly for haplotype2 (fasta)
 - ancestral genome assembly (fasta - proxy for the gene order)
 - ancestral genome annotation (gft file)
 - RNAseq data (single end or paired-end)
 - database of known TE elements 
 - txt file with the name of the ancestral sex chromosome and their orientation
 - optionally: database of known genes in closely related species


assemblies and annotation are all deposited on zenodo and available here: 

to download: 

`wget /path/to/zenodo/....tar.gz` #(to be updated later)

## 2 - setting your config file


| option in config | description |
| --- | --- |
| *genome1* | ="/path/to/EASYstrata/example_data/genome1.hap1.fa.gz" |
| *haplotype1* | ="genome1.hap1" |
| \[*genome2*\] | ="/path/to/EASYstrata/example_data/genome1.hap2.fa.gz" |
| \[*haplotype2*\] | ="genome1.hap2" |
| annotate | ="YES" |
| *RelatedProt* | ="/path/to/EASYstrata/example_data/relatedprot.fa.gz" |
| \[*RNAseqlist*\] | ="/path/to/EASYstrata/example_data/rnaseq.list.txt" |
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
| \[*ancestral_genome*\] |  ="/path/to/EASYstrata/example_data/Mlag129.A1.fa.gz" |
| \[*ancestral_gff*\] | ="/path/to/EASYstrata/example_data/Mlag129.A1.gff" |
| *TEancestral* | ="/path/to/EASYstrata/example_data/Mlag129.A1.TE.bed" |
| *scaffolds* | ="/path/to/EASYstrata/example_data/scaffold.txt" |

click [here](/example1.config) to see the example

## 3 - creating some files list: 

if you have PE data we expect two column file with R1 in column1 and R2 in column2:

simply use a command like this: 

`readlink -f your_rna_folder/*gz |awk 'ORS=NR%2?FS:RS ' > rnaseq.list.PE.txt `  

other multi column link could be generated similarly with awk if needed

I like to use the `readlink -f` command to extract full path and add them to txt.files


## 3 - launching the workflow : 


```./master.sh -o1 2>&1 |tee log```


## 4 - results visualisation :

insert resulting plots and other stats here
