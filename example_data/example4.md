# Example 4: full workflow with already annotated genome(s) 

If you already have a gtf prediction and BED files of TE then you can start here!

Assuming an ancestral species as well  
 
Example taken from fungi as they are small 

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
| annotate | ="NO" |
| \[*orthoDBspecies*\] | ="Fungi" |
| *fungus* | ="YES" |
|| \[*gtf1*\] | ="/path/to/EASYstrata/example_data/genome1.hap1.gff.gz"  |
| \[*gtf2*\] | ="/path/to/EASYstrata/example_data/genome2.hap2.gff.gz" |
| *busco_lineage* | ="basidyomicota_odb10" |
| *interpro* | ="NO" |
| \[*ancestral_genome*\] |  ="/path/to/EASYstrata/example_data/Mlag129.A1.fa.gz" |
| \[*ancestral_gff*\] | ="/path/to/EASYstrata/example_data/Mlag129.A1.gff" |
| *TEancestral* | ="/path/to/EASYstrata/example_data/Mlag129.A1.TE.bed" |
| *scaffolds* | ="/path/to/EASYstrata/example_data/scaffold.txt" |

click [here](/example4.config) to see the example

## 4 - launching the workflow : 


```./master.sh -o3 2>&1 |tee log```


## 5 - results visualisation :

the expected results are the same as for example1 as this relies also on RNAseq and the ancestral genome for the analysis.
see examples plot and files from [example1](/example1.config)
