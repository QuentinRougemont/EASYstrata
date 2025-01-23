
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

click [here](example1.config) to see the example

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



![Fig4.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig4.png)

Figure 1: A) Synteny plot from GeneSpace showing gene synteny between ancestral species (ancestral_sp) and the two mating type of *Microbotryum lychnidis dioiciae 1064*  and B) Circos plot between the ancestral species and mating type A1 (left part) and cicros plot between mating type A1 and mating tpye A2. External links show the position of some major gene (red, green and light blue single link as well as the centromeres in violet). External density plot in lightblue display the gene density, the green density plot displays TE density. Red and darkblue interior links display single copy orthologs links.




![Fig5.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig5.png)

Figure 5: Ds plot and arrangements. A) dS values along the ancestral chromosomes. B) dS values along the ancestral gene order after returning the chromosomes and removing the large autosomal part on contig 8.
C) and D) arrangement as infered based on gene rank in mating type A1 and A2 respectively. 



![Fig6.A.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig6.A.png)

**Figure 6A:** Results from the changepoint analysis for 3  (panel A) to 8 changepoints (panel F) that will be automatically performed to infer evolutionary strata. Each changepoint panel displays the distribution of raw data (i.e., dS values as black dots) along with 25 draws from the joint posterior distribution (grey lines) and 95% highest density interval (red lines). Posterior distributions of the changepoints are shown in blue with one line for each chain. Note that in general the "strata" with zero dS value on the left most and rigth most side respectively will correspond to the PAR, not true evolutionary strata.  


it is important to check the convergence of the runs for each parameters : 
this will be perform automatically in our code resulting in these plots for each changepoint tested.

![Fig6.B.svg](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig6.B.svg)

**Figure 6B:** posterior fit of the model parameter and mixing of the chains. Here only the values inferred for the 7-changepoint model are shown - the one with highest support.
values for all other models are generated on the flye.


The MCP produced many usefull informations that will be extracted and automatically exported in *.txt* tables:

* `modelX.chpt.txt`:  X = number of changepoint tested (from 1 to 9). 

    These files is the output of the summary function from MCP. 
    
    It contains the following columns:
 
    1. name: name of the parameter (changepont and interval) 
    2. mean: mean value of dS and interval (gene order based) 
    3. lower/upper: lower and upper boundaries, 
    4. Rhat: is the Gelman-Rubin convergence diagnostic which is often taken to be acceptable if <1.1. 
    5. n.eff:  is the effective sample size computed using effectiveSize. Low effective sample sizes are also obvious as poor mixing in trace plots .

* `modelchoice.txt` : 
    This file contains info from the loo model choice operation 

    It contains the following columns:

    1. elpd_diff 
    2. se_diff 
    3. elpd_loo 
    4. se_elpd_loo 
    5. p_loo 
    6. se_p_loo looic 
    7. se_looic

* `weights.txt` : 

    This file contains the weights of each tested models
    higher weights indicates higher supports.


* `HypothesisXstrata.txt` :  X = number of changepoint tested (from 1 to 9). 

    These file contains results from hypothesis testing (BayesFactor and posterior probabilities aiming at testing difference among strata) 

    Here only difference in dS values among adjacent strata are tested when moving forward from the left to the right of the gene order. 

    The two directionalyty of differences are tested, i.e.: 

    "int_1 > int_2": the intercept is greater in strata 1 than 2. 
    "int_1 < int_2": the intercept is greater in strata 2 than 1. 
    
    This is repeated for all comparison of adjacent interval for 1 to 9 changepoints.


* `classif.sX.$haplo1.$haplo2` :   X = number of changepoint tested (from 1 to 9). 


    A three column file containing the assignment of single copy orthologs to a strata :
    1. column1: genes in $haplo1 

    2. column2: genes in $haplo2 

    3. column3: strata of appartenance 

    These file are use to automatically color links in Ideograms. 

* `df.txt` : a dataframe recaputilating all infos 


## other output : 

vio-boxplot with statistical tests. 

here's an example for the two best model in the studied species: 
 
![Fig7.svg](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig7.svg)

**Figure 7:** example violinboxplot for the two "most likely" models inferred by the loo analysis.
By default plots are constructed for all models (from 1 to 9 changepoints). Default statiscal test from the ggstats plot package are used
assuming parametric tests. 



dS colored by strata along the ancestral gene order:

![Fig8.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig8.png)

**Figure 8:** dS values plotted along the ancestral gene order for all possible models from three to eight changepoints  
each point is a gene dS value colored according to the strata of assignation. 

dS colored by strata along the ancestral genome:

automatically generated for each changepoint values: 
![Fig9.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig9.png)

**Figure 9:** dS values plotted along the ancestral genome for all possible models from three to eight changepoints  
each point is a gene dS value colored according to the strata of assignation

a posterior colored ideogram: 
automatically generated for each changepoint values: 
![Fig10.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig10.png)

**Figure10:**  example ideograms infered for the most likely models here. Links are colored according to their strata of appartenance. 



# to do: add circos plots colored by Ds quantile and strata - same with ideogram
