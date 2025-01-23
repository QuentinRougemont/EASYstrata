# Example 2: full workflow with mapped RNAseq and ancestral genome

## 1 download example data folder from zenodo (see below)

data needed:
 - genome assembly for haplotype1 (fasta: *genome1.hap1.fa.gz* in example data folder)
 - genome assembly for haplotype2 (fasta: *genome2.hap2.fa.gz* in example data folder)
 - ancestral genome assembly (fasta - proxy for the gene order: *Mlag129A1.fa.gz* in example data folder)
 - ancestral genome annotation (gff file: *Mlag129A1.gff.gz* )
 - bam file mapped on  A1 and A2 respectively (bam: *bam_on_genome1hap1* and *bam_on_genome2hap2* in example data folder) 
 - txt file with the name of the ancestral sex chromosome and their orientation (see: *scaffold.txt*)
 - optional : database of additional external genes evidence (fasta: *relatedprot.fa.gz* in example data folder)

### download the data 

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

bam files are deposited in the folder with the corresponding names (i.e. *bam_on_genome1hap1* *bam_on_genome2hap2* )  

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
| annotate | ="YES" |
| *RelatedProt* | ="/path/to/EASYstrata/example_data/relatedprot.fa.gz" |
| \[*RNAseqlist*\] | ="/path/to/EASYstrata/example_data/rnaseq.list.txt" |
| \[*bamlist1*\] | ="/path/to/EASYstrata/example_data/bam.aligned_on_genome1.txt" |
| \[*bamlist2*\] | ="/path/to/EASYstrata/example_data/bam.aligned_on_genome1.txt" |
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

click [here](example2.config) to see the full config file example 

## 3 - creating some files: 
list of bam: 

readlink -f \*bam > list_of_bam.txt

everything else should be rather straightforward.


## 4 - launching the workflow : 


```./master.sh -o1 2>&1 |tee log```


## 5 - results visualisation :

the expected results are the same as for example1 as this relies also on RNAseq and the ancestral genome for the analysis.
see examples plot and files from [example1](/example1.config)
