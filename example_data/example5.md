
# Example 5: full workflow with RNAseq and no ancestral genome

## WORKFLOW TESTED ON A 2.8 gb Genome of Silene latifolia 

## 1 download example data genome

data needed:
 - genome assembly for Silene latifolia
 - RNAseq data from Silene latifolia 
 - database of known TE elements 
 - txt file with the name of the ancestral sex chromosome and their orientation
 - optionally: database of known genes in closely related species

S. latifolia genome assembly is available on NCBI at the following link:
[insert_link_here](/path/to/link/)


RNAseq data are available at XXXXX and can be download from there 

 
a set of TE for plant can be obtained [here](https://trep-db.uzh.ch/):

for plants:
`wget https://trep-db.uzh.ch/downloads/trep-db_nr_Rel-19.fasta.gz`

then this should be added to the config file 

note that repeatexplorer could also work 

download with:
`wget http://repeatexplorer.org/repeatexplorer/wp-content/uploads/Viridiplantae_v3.0_ALL_protein-domains_repet_formated.fa`


## 2 - setting your config file


| option in config | description |
| --- | --- |
| *genome1* | ="/path/to/EASYstrata/example_data/Slatifolia.fa.gz" |
| *haplotype1* | ="Slatifolia" |
| \[*haplotype2*\] | ="Slatifolia_Y" |
| annotate | ="YES" |
| *RelatedProt* | ="/path/to/EASYstrata/example_data/relatedprot.fa.gz" |
| \[*RNAseqlist*\] | ="/path/to/EASYstrata/example_data/rnaseq.list.example5.txt" |
| \[*orthoDBspecies*\] | ="Viridiplantae" |
| *fungus* | ="NO" |
| *TEdatabase* | ="/path/to/EASYstrata/example_data/trep-db_nr_Rel-19.fasta.gz" |
| *ncbi_species* | ="eudicots" |
| *busco_lineage* | ="viridiplantae_odb10" |
| *interpro* | ="NO" |
| *scaffolds* | ="//path/to/EASYstrata/example_data/scaffold_example5.txt" |

click [here](/example5.config) to see the example


3 - launching the workflow : 


```./master.sh -o1 2>&1 |tee log```


4 - results visualisation :


see [example1](/example1.md) to see examples plots
