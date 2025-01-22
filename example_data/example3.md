# Example 3: full workflow without RNAseq one single genome (containing X/Y or two MAT)and no ancestral genome (a.k.a worst case scenario)


## 1 download example data genome

data needed:
 - genome assembly for haplotype1 containing X/Y or two MAT ( *fakegenome.fa.gz* in fasta format in the example folder)
 - database of known TE elements (*TE.fa.gz* in the example folder) 
 - txt file with the name of the sex chromosome/MAT and their orientation (*scaffold_example3.txt* in the example folder)
 - optionally: database of known genes in closely related species (*relatedprot.fa.gz* in the example folder)

### prepare the data 

```
cd example_data
#assemblies and annotation are all deposited on zenodo and available through wget:
wget https://zenodo.org/records/14716941/files/data.tar.gz?download=1
tar zxf data.tar.gz\?download\=1
mv data/* .
```

in this example we will work with the "fakegenome.fa.gz" it contains two MAT the A1 and A2 with following contig ID:
fakegenome_A1
fakegenome_A2 

the basename being fakegenome and matching the basename assembly, as expected.

scaffold are provided in the example_folder (named **scaffold_example3.txt**):  

fakegenome     fakegenome_A1  N



## 2 - setting your config file

based on the previous data it is very easy to set up the config file:

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

click [here](example_data/example3.config) to see the full config file example 

## 3 - creating some files: 

use readlink to set full path to text files

## 4 - launching the workflow : 


```./master.sh -o1 2>&1 |tee log```  

## 5 - launching the workflow :  

insert resulting plots and other stats here
