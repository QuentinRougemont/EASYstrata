# Example 4: full workflow with already annotated genome(s) 

If you already have a gtf prediction and BED files of TE then you can start here!

Assuming an ancestral species as well  
 
Example taken from fungi as they are small 

## 1 download example data genome

data needed:
 - genome assembly for haplotype1 (fasta: *genome1.hap1.fa.gz* in example data folder)
 - genome assembly for haplotype2 (fasta: *genome2.hap2.fa.gz* in example data folder)
 - genome annotation for haplotype1 (fasta: *genome1.hap1.gtf* in example data folder)
 - genome annotation for haplotype2 (fasta: *genome2.hap2.gtf* in example data folder)
 - ancestral genome assembly (fasta - proxy for the gene order: *Mlag129A1.fa.gz* in example data folder)
 - ancestral genome annotation (gff file: *Mlag129A1.gff.gz* )
 - txt file with the name of the ancestral sex chromosome and their orientation


### donwload the data 

```
cd example_data
#assemblies and annotation are all deposited on zenodo and available through wget:
wget https://zenodo.org/records/14716941/files/data.tar.gz?download=1
tar zxf data.tar.gz\?download\=1
mv data/* .
```

this should contain all the necessary data mentionned earlier 

genome1hap1 is the name of haplotype1

genome2hap2 is the name of thaplotype2

ancestral scaffold are provided directly in github example folder and consist in the following tab separated file (named **scaffold.txt**):  

Mlag129A1       Mlag129A1_contig_8      R
Mlag129A1       Mlag129A1_contig_11     N


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

click [here](example_data/example4.config) to see the example

## 4 - launching the workflow : 


```./master.sh -o3 2>&1 |tee log```


## 5 - results visualisation :

the expected results are the same as for example1 as this relies also on RNAseq and the ancestral genome for the analysis.
see examples plot and files from [example1](/example1.config)
