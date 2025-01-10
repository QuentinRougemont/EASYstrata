# genome annotation - Synteny - Ds computation - changepoint analysis - whole genome alignments 
====================================================================================

# TO DO: 
	fill this readme 

# Requirements

This software is suitable only linux-like systems  (Unfortunately not Windows or MAC) 

# Table of content 

   * [Purpose](#purpose)
   * [Installation](#installation)
   * [Before-launching-the-workflow](#before-launching-the-workflow)
   * [How to use](#how-to-use)
	   * [From genome annotation to strata inference](#from-genome-annotation-to-strata-inference)
	   * [Genome annotation only](#genome-annotation-only)
	   * [Synteny and strata inference from existing data](#synteny-and-strata-inference-from-existing-data)
	   * [Restart at any step](#restart-at-any-step) 
   * [Input data](#input-data)
   * [Example input data](#example-input-data)
   * [Details of the worfklow and results](#details-of-the-workflow-and-results)
   * [Output files](#output-files)


# Purpose:
##  sets of scripts to : 
[I - Perform TE and gene prediction](#I---Perform-TE-and-gene-prediction)

[II - Identify synteny blocks and rearragements](#II---Identify-synteny-blocks-and-rearragements)

[III - Plot dS along the genome](#III---Plot-dS-along-the-genome)

[IV - Perform changepoint analysis to identify evolutionary strata](#IV---Perform-changepoint-analysis-to-identify-evolutionary-strata)

<img src="https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig1.png" width = "490" heigth = "490">


Installation: 
-------

- [Installation instructions](INSTALL/INSTALL.md)



# Before Launching the workflow

After cloning the pipeline, please, work from within it to preserve the architecture  

We recommend that you clone the pipeline ***for each of your new project*** and work within it.   

Keep all project separated otherwise it will be difficult to recover your results.   


All options and full path to input files **must** be set in the **config** file provided in : `config/config` .

PLEASE, carefully read the user guide below before any attempte at running the workflow


# How to use: 

### In short : 

simply run:
```
./master.sh --help to see all options 
```

The script provides several options depending on what you want: 

-o 1: to perform all analyses (TE+gene prediction, GeneSpace, single copy orthogls inference between X/Y,  ds, evoluationary strata inference)
 
-o 2: to perform TE and gene prediction, as well as synteny analysis (no dS or Strata inference)

-o 3: to perform synteny analysis (from GeneSpace) and subsequent analysis 

-o 4: to perform only d*S* and subsequent analysis 

-o 5: to only do only perform GeneSpace/Gene Synteny + Whole Genome Synteny 

-o 6: to only predict TE and genes

-o 7: to perform only the changepoint analysis

-o 8: to only perform the plots after dS computations


## Detailed use cases

as described above the workflow is simply run with the command:
```
./master.sh -o X #with X an option from 1 to 7. 
```

each option and their requirement in the config file are described below :

an example config file is provided [here](https://github.com/QuentinRougemont/EASYstrata/blob/main/example_data/example.config)


### From genome annotation to strata inference

basically this means running the whole workflow.  

to do so run : 

```./master.sh -o1 2>&1 |tee log```

this will:

*	 Perform TE annotation using **repeatmodeller** (de-novo prediction) and **repeatMasker**
*	 Perform gene prediction on the softmasked genome using **BRAKER** (with or without RNAseq). 
*	 Evaluate the quality of the gene prediction (with **BUSCO** mostly)
*	 Run **GeneSpace** between your genomes (and eventual ancestral genome) to infer broad pattern of synteny including inference of single copy orthologs from **orthofinder**  
*	Run **paml** to estimate synonymous divergence between the sequences/region of interests
*	Perform various plots (circos plot, ds along the genome, ideogram, etc) 
*	Infer the most likely number of evolutionary strata using a changepoint analysis


### other usefull option: 

* **Genome annotation only**

to simply perform genome annotation (gene and TE) run :
```./master.sh -o 6 2>&1 |tee log```

NOTE: if you already have a softmasked genome simply set :  

**annotateTE="NO"** in the config file and this step will be skipped



### Synteny and strata inference from existing data

in case you already have a pair of genome (fasta format) along with their gene prediction (gff):    

```./master.sh -o 3 2>&1 |tee log```

This will perform all steps appart from the gene prediction.  

### dS strata inference from existing data

in case you already have a pair of genome (fasta format) along with their gene prediction (gff) and the dS computed from previous step.

```./master.sh -o 4 2>&1 |tee log```

This will perform all steps after the dS.  

Mostly usefull for debugging  

### Only for synteny analysis 

in case you already have a pair of genome (fasta format) along with their gene prediction (gff):

```./master.sh -o 4 2>&1 |tee log```

This will only run GeneSpace and minimap and perform some plots.

Mostly usefull for debugging or if you are not interested in Strata.



### changepoint only 

run simply:  

```./master.sh -o 7 2>&1 |tee log```

this is usefull to explore various parameter settings in the MCP analysis. 
For instance you may want to use prior in the MCP (see below) or tweak the order of the scaffolds. 

### plots only:

run simply:  

```./master.sh -o 8 2>&1 |tee log```

this is usefull to make only the plot after the dS computation, for exemple, you may want to change some options in the R plot after a first pass analysis with default parameters.



### Restart at any step 

the worflow is designed to work at any step in the process. In case of bug you can restart it from whenever it crashes (after fixing the bug) and this should work smoothly.




## Input data 

**/!\ the number of input will vary strongly depending on the analysis you want to perform.**
several files are **Compulsory** 

**1 - genome assemblies:**
Your **input data** may be
\- one genome assembly containing both sex/mating type chromosomes

\- or **ideally:** two separate haplotypes assemblies containing each one of the sex/mating type chromosomes.

**Warning names of fasta and contigs/scaffold/chromosomes :**
we recommend to use short name for each of your genome assemblies name and avoid any special characters appart from underscore.

**For instance:**
species-1.fasta will not be valid in GeneSpace. => Use **species1.fasta** instead.

For the chromosome/contig/scaffold ids we recommand a standard naming including the Species name within it without any thing else than alhpanumeric character.  
\- example: **species1_chr1** or **species1_contigX** or **species1_scaffoldZ**

**2 - other input data** 
additional input data will be highly dependent on the analysis you want to perform (full workflow or not).
These must be specified in the [**config file**](https://github.com/QuentinRougemont/EASYstrata/blob/main/config/config) and will typically include:  

* 2.1 - **ancestral genome**  optional but highly recommended: the link to an ancestral genome to plot the gene order - ideal to infer more accurately single copy orthologs  
* 2.2 - **ancestral gff**  optional but highly recommanded:  the gene prediction associated with the ancestral genome 
* 2.3 - **BUSCO lineage name** : compulsory for genome annotation: the name of the BUSCO lineage from your species (available through busco --list-lineage)) 

**for genome annotation:**  

* 2.4 - **RNAseq** data for each genome (optional): to improve BRAKER annotation  (can also be bam files)  
* 2.5 - **Proteins for annotation (optional)** : if no protein data are available we will use orthoDB12  (downloaded automatically)
* 2.6 - **orthoDB12 species name** : a lineage species name among: 
 "Metazoa" "Vertebrata" "Viridiplantae" "Arthropoda" "Eukaryota" "Fungi" "Alveolata"

**for TE prediction:**  

* 2.7 - **TE database** : compulsory for genome annotation if genome not softmasked already: the name of the TE database (some are available online depending on your taxon) 
* 2.8 - **NCBI taxon** : compulsory: a taxon name for NCBI (used with repeatmasker)  
* 2.9 - **TE bed files** :  optional: a pairs of bed file containing TE for your region of interest if already available (to display on circos plots)

**for evolutionary strata - circos plot - ideogram - etc** 
* 2.10 -  **scaffold.txt**: a tab sepearate file containing the genome name, and scaffold ids of the ancestral genome.  

for mor on the format needed see [example data folder](https://github.com/QuentinRougemont/EASYstrata/blob/main/example_data) 


##  Full config file details:


| option in config | description |
| --- | --- |
| *genome1* | **Compulsory:** Full path to the assembly of the genome or haplotype you wish to analyse. |
| *haplotype1* | **Compulsory:** Name of the genome1 to be used as a /contig/scaffold/chromosome basename |
| \[*genome2*\] | **Optional:** Full path to the assembly of the second haplotype you wish to analyse. Only for the case where you have two haplotypes, each containing one of the sex/mating type chromosomes. |
| \[*haplotype2*\] | **Compulsory:** Name of the second haplotype. This can be a basename of all chromosome if a second assembly is avaiable, or the name of the scaffold/contig/chromosomes corresponding to the sex/MAT chromosome (e.g. species_chrY) |
| annotate | **Compulsory:** a string "YES"/"NO" stating wether genome should be annotated or not|
| *RelatedProt* | **Optional:** Full path to a fasta file containing protein sequences to be used for gene prediction. |
| \[*RNAseqlist*\] | **Optional:**: Full path to a '.txt' file containing the list of RNA-seq data files. |
| \[*bamlist1*\] | **Optional:**. Full path to a .txt file containing the list of bam files for *genome1* (alignment of RNA-seq data onto the fasta of *genome 1*). |
| \[*bamlist2*\] | **Optional:** with option *b* and if *genome2* is given. Full path to a .txt file containing the list of bam files for genome2 (alignment of RNA-seq data onto the fasta of *genome 2*). |
| \[*orthoDBspecies*\] | **Compulsory:** "Metazoa" "Vertebrata" "Viridiplantae" "Arthropoda" "Eukaryota" "Fungi" "Alveolata". Will use a database from **orthoDB** for gene prediction. |
| *fungus* | **Compulsory:** "YES" or "NO" (default), whether your species is a fungus. |
| annotateTE | **Compulsory:** a string "YES"/"NO" stating wether genome should be annotated or not|
| *TEdatabase* | **Compulsory for annotation of TE:** Full path to a database of TE for your species/genus, used in TE prediction, in fasta format. |
| *ncbi_species* | **Compulsory for annotation of TE:** Name of the ncbi species, used in TE prediction. ==list available [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)== |
| \[*gtf1*\] | **Optional**. Full path to a .gtf file for an existing gene prediction on *genome1*. |
| \[*gtf2*\] | **Optional** Full path to a .gtf file for an  existing gene prediction on *genome2*. |
| *busco_lineage* | Lineage used for **busco** analysis. You can access the list of available lineages by typing `busco --list-dataset`. |
| *interpro* | YES or NO (default), whether **interproscan** quality check of the genome annotation should be performed. Warning: this can take up to several days for a big dataset. |
| \[*ancestral_genome*\] | Compulsory if *ancestral* is set as "outgroup". Full path to the genome of the species used as proxy for the ancestral state. |
| \[*ancestral_gff*\] | Compulsory if *ancestral* is set as "outgroup". Full path to the gff (genome annotation) of the species used as proxy for the ancestral state. |
| *scaffolds* | Full path to the list of focal scaffolds (i.e. the scaffolds composing the sex / mating type chromosomes). |


# PART BELOW TO BE UPDATED : 


# Details of the worfklow and results

## I - Perform TE and gene prediction

## RNAseq alignment - TE masking - Gene prediction - Quality assessment

### Parameters set in config

==Give examples of files for RNAseqlist and bamlist HERE==

### Options

**(*a*) - Align RNA & Annotate:**  
Will perform alignment of provided RNA-seq data and use it as additional information for genome annotation.

**(*b*) - Annotate, use BAM of RNA**  
Will use provided BAM of already aligned RNA-seq data and use it as additional information for genome annotation.

**(*c*) - Annotate, no RNA**  
Will perform genome annotation without using RNA information.

**(*d*) - Skip**  
If you have already annotated your genome, will use provided gtf for the following steps, effectively skipping genome annotation.

## Operations of step I

### 1\. Alignment of RNA-seq data (only with option *a*)

Corresponding script: `00_scripts/launch_rnaseq.sh`

- Reads trimming using **trimmomatic**

The script will detect whether the data is Single-End or Paired-End and launch trimmomatic, then count the number of retained reads.

- Creation of database for **gsnap** using **gmap**
- Alignment using **gsnap**

Corresponding scripts: `00_scripts/03_gsnap_PE.sh` for PE ; `00_scripts/03_gsnap_SE.sh` for SE

- Mapping quality assessment

Sequencing depth and MAPQ along the genome will be computed and plotted. The resulting plots can be found in ==XXX/Depth/== and ==XXX/mapq/== .

==Insert example plot here==

### 2\. TE discovery and masking

Corresponding script: `00_scripts/launch_step05_to_08.sh`

- *De novo* repeat identification using **repeatmodeler** on the input genome(s) to be annotated
- Genome masking using **repeatmasker** with known TE libraries

### 3\. Genome annotation, quality assessment and filtering

Corresponding script: ./00_scripts/06_braker.sh

- Five successive rounds of gene prediction based on protein database, using **braker**
    
- One round of gene prediction using RNA-seq data, using **braker** (only with options *a* and *b*)
    
- Quality assessment and reports production for each round of gene prediction
    

Two tools can be used at this stage for quality assessment:  
\- **Busco** (corresponding script: `00_scripts/07_busco_after_braker.sh`)  
\- **Braker** report on the raw hintsfile  
This report includes number of genes, number of introns per gene, gene support, number of complete genes and various histograms useful for error checking.

- Combination of protein-based and RNA-seq-based gene models using **TSEBRA** (only with options *a* and *b*)

Please read the [TSEBRA manual](https://github.com/Gaius-Augustus/TSEBRA) before running the script.  
The best round of protein-based gene prediction and the RNA-seq-based gene prediction are given as input in TSEBRA.  
Warning: TSEBRA parameters *intron_support* and *stasto_support* are set to 0 in this workflow (default in TSEBRA: 1 and 2 respectively). This means that only overlapping genes between the two gene models will be filtered. You can change this parameter and others to adjust to your desired level of stringency in the TSEBRA config file: `config/default.cfg`

- Final genome annotation reshaping

**For options *a* and *b*:** The final genome annotation  is the output from TSEBRA.  
**For option *c*:** The final genome annotation is the best protein-based braker round, as evaluated with busco.

Corresponding script: `00_scripts/08_braker_reshaping.sh`  
The genes will be renamed to insert the scaffold name for clarity in downstream analyses.  
Because the next steps in the workflow involve single copy ortholog identification, only the longest transcript is kept for each gene.

- Final genome annotation quality assessment

Two more in-depth tools can be used at this stage for quality assessment: (==\+ busco on final genome pred ?==)  
\- **Blast** against **Uniprot**  
If you wish to skip this, comment l.295 of the script `00_scripts/08_braker_reshaping.sh`  
\- **InterProScan** (if option interpro is set to "YES" in the config file and Blast against Uniprot successfully ran)  
This tool is more time-consuming.

## II - Identify synteny blocks and rearragements

will enable to infer gene order for dS interpretation 

## Input of step II

### Parameters set in config ==in yellow parameters to be set==

|     |     |
| --- | --- |
| **option in config** | **description** |
| *==ancestral==* | "chromosome" or "outgroup", whether the sequence used as proxy for the ancestral state is one of the sex / mating type chromosomes or an ougroup provided below. |
| ==\[*ancestral_chromosome_scaffolds*\]== | Compulsory if *ancestral* is set as "chromosome". Full path to the list of scaffolds of the chromosome used as proxy for ancestral state. |
| \[*==outgroup_orthofinder==*\] | Advised if *ancestral* is set as "chromosome". Full path to a list of genomes to be used as outgroups in OrthoFinder only. |
| ==\[*ancestral_outgroup_scaffolds*\]== | Compulsory if *ancestral* is set as "outgroup". Full path to the list of focal scaffolds for the outgroup used as proxy for ancestral state. |

==Give examples of files for scaffolds, ancestral_chromosome_scaffolds, outgroup_orthofinder and ancestral_outgroup_scaffolds ?==

### Operations of step II

### 4a. Minimizer alignment and plots of target region

Corresponding script: `00_scripts/11_run_genesSpace_paml_ideogram.sh`

- Alignment between the two haplotypes using **minimap2**

If you provided as input two haplotypes containing each one of the sex/mating type chromosomes, the whole haplotypes will be aligned.  
If you provided one genome containing both sex/mating type chromosomes, only the corresponding focal scaffolds (as indicated with option *scaffold*) will be aligned.

- Alignment between the two haplotypes and an outgroup genome used as proxy for ancestral state if you have one (option B), using **minimap2**
    
- Construction of whole genome dotplot using **pafR** (only with two haplotypes as input)
    
- Construction of synteny plot on the focal scaffolds using **pafR**
    
ex: minimap based divergence along the mating type chromosomes :

![Fig2.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig2.png)

ex: minimap based whole genome alignment : 
	
![Fig3.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig3.png)



### 4b. Ortholog reconstruction

- Launching of **GeneSpace**

In short, this will:  
\- Identify single copy orthologs with **OrthoFinder**  
\- Construct a dotplot and a riparian plot of whole genome (only with two haplotypes as input) \[==GeneSapce==\]  
\- Construct a riparian plot on focal scaffolds \[==GeneSpace==\]  
For more information, consult the [GeneSpace readme](https://github.com/jtlovell/GENESPACE).

ex: Synteny plot from GeneSpace

![Fig4.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig4.png)


# III - Plot dS along the genome

## STEP III - Compute and plot dS - Plot ideogram and rearrangements

## Input of step III

### Parameters set in config

Same as step II (see above).

## Operations of step III

### 5\. Single copy orthologs alignment

Align all coding sequences from the focal scaffolds.

- **TranslatorX**
- **muscle**

### 6\. dS calculation and plotting

- Calculation of dS (& dN) using **PAML**
    
- Plotting dS values using a custom R script
    

Corresponding script: `00_scripts/Rscripts/03_plot_paml.R`  

dS values are plotted along the focal scaffolds, and, if 2 haplotypes were given as input, along the whole genome.  
The gene order will be that of the genome used as proxy for the ancestral state: either one of the two sex/mating type chromosomes, or an outgroup (see option *ancestral*).  
It is possible to modify the R script to adapt the plotting options to your needs (for instance position and direction of scaffolds).

Ex: Ds plot : 

![Fig5.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig5.png)


### === - Plot circos (==step III==)

Corresponding script: `00_scripts/Rscripts/05_plot_circos.R [options]`  

Construction of a circos plot of the focal scaffolds, tracing links between their single copy ortholog genes, using **circlize**.  
* If TE info are available these can also be provided as arguments.

* gene density can be extracted from the bed file in genespace and be plotted if provided as arguments

It is possible to modify the R script to adapt the plotting options to your needs (for instance position and direction of scaffolds).

By default any fused autosome will be plotted but these can be removed from the contig list

See **figure4 panel B** above for example.

# IV - Perform changepoint analysis to identify evolutionary strata

## Step IV

### 1\. Changepoint analyses

Before launching this step, we strongly suggest that you consult the results of the workflow, especially the dS plot. i

Once you have deciphered clear hypotheses as to whether there are strata on your focal scaffolds, and where they occur, you can use the R script.

`00_scripts/Rscripts/06.MCP_model_comp.R` to perform changepoint analyses on the dS, using **mcp**.

To that end, you can automatically launch the code ```master.sh -o7``` and it will launch the MCP, producing several graph as well as colored ideogram according for each model infered by the MCP 

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

**Figure 9:** dS values plotted along the ancestral gene order for all possible models from three to eight changepoints  
each point is a gene dS value colored according to the strata of assignation. 

dS colored by strata along the ancestral genome:

automatically generated for each changepoint values: 
![Fig9.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig9.png)

**Figure 8:** dS values plotted along the ancestral genome for all possible models from three to eight changepoints  
each point is a gene dS value colored according to the strata of assignation

a posterior colored ideogram: 
automatically generated for each changepoint values: 
![Fig10.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig10.png)

**Figure10:**  example ideograms infered for the most likely models here. Links are colored according to their strata of appartenance. 






# Options to run part of the workflow

If you wish to perform only part of the workflow or relaunch it from a given step, use option *\-o*

`bash master.sh -o 2` : perform steps I and II  
now if you have a gtf and and genome assembly (either from running this pipeline or any other annotation tools):
`bash master.sh -o 3` ; perform steps II and III (if step I already ran successfully in a previous run)  
`bash master.sh -o 4` : perform step III only (if steps I and II already ran successfully in a previous run)  
`bash master.sh -o 5` : perform step II only (if step I already ran successfully in a previous run)  
`bash master.sh -o 6` : perform step I only
`bash master.sh -o 7`: perform step IV only




# --------------------------------------------------------------------------

# list of operations and tools


| __Operation__                     |  __Tools__                         |  __data type__  | 
|:---------------------------------:|:------------------------------:|:-----------:| 
| __read trimming__                |  Trimmomatic                   | RNAseq         | 
| __read mapping__                 |  gmap/gsnap                    | RNAseq          | 
| __sorting read__                 |  samtools                      | RNAseq        |
| __mapping quality assement__     |  samtools + R                  | RNAseq        |
| __TE detection and softmasking__ |  RepeatModeler + RepeadMasker  | genome assembly |
| __genome annotation__            |  Braker + tsebra               | genome assembly + protein + database |
| __quality assessment__           |  Busco + Blast + Inter Pro     | genome prediction |
| __Synteny and HOG__              |  GeneSpace (including OrthoFinder/MCScan) | gene prediction and proteins |
| __cds alignement__               |  muscle + translatorX          | gene prediction (single copy orthologs) | 
| __Ds computation__               |  paml                          | CDS alignment |
| __Ds plot/CIRCOS plot__          |  R                             | Ds and genome information |
| __whole genome alignemnt__       |  minimap2                      | genome assemblies |
| __gene microsynteny__            |  R                      | single copy orthologs |
| __changepoint analysis__         |  R                      | Ds values and gene order |



This code has been tested with linux. 


Normally, you should only run the script ```./master.sh```

below we provided a description of what will be done at each steps.


