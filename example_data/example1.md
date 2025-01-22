
# Example 1: full workflow with RNAseq and ancestral genome

## 1 download example data genome

data needed:
 - genome assembly for haplotype1 (fasta: *genome1.hap1.fa.gz* in example data folder)
 - genome assembly for haplotype2 (fasta: *genome2.hap2.fa.gz* in example data folder)
 - ancestral genome assembly (fasta - proxy for the gene order: *Mlag129A1.fa.gz* in example data folder)
 - ancestral genome annotation (gff file: *Mlag129A1.gff.gz* )
 - RNAseq data (single end or paired-end: folder of fastq named *rnaseq* in example data folder)
 - database of known TE elements (fasta: *TE.fa.gz* in example data folder) 
 - txt file with the name of the ancestral sex chromosome and their orientation (see: *scaffold.txt*)
 - optional : database of additional external genes evidence (fasta: *relatedprot.fa.gz* in example data folder)
 
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

**WARNING** : replace /path/to by your true working path below.

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

click [here](example_data/example1.config) to see the example

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
