# genome annotation - Synteny - Ds computation - changepoint analysis - whole genome alignments 
====================================================================================

# TO DO: 
	fill this readme 

# Requirements

This software is suitable only in linux-like systems  (Unfortunately not Windows or MAC) 

# Table of content 

   * [Purpose](#purpose)
   * [Installation](#installation)
   * [Before-launching-the-workflow](#before-launching-the-workflow)
   * [How to use](#how-to-use)
        * [option 1: From genome annotation to strata inference](#option-1-from-genome-annotation-to-strata-inference)
        * [option 2: Anotation and synteny inference only](#option-2-anotation-and-synteny-inference-only)
        * [option 3: Synteny and strata inference from existing data](#option-3-synteny-and-strata-inference-from-existing-data)
        * [option 4: Evolutionary strata inference from existing data](#option-4-evolutionary-strata-inference-from-existing-data)
        * [option 5: Only synteny analysis](#option-5-only-synteny-analysis)
        * [option 6: Genome annotation only](#option-6-genome-annotation-only)
        * [option 7: Only evolutionary strata inference](#option-7-only-evolutionary-strata-inference)
        * [option 8: Synteny plots only](#option-8-synteny-plots-only)  
   * [Input data](#input-data)
   * [Details of the worfklow and results](#details-of-the-worfklow-and-results)
   * [working examples](#working-examples)
        * [1: RNAseq and ancestral genome](#1-rnaseq-and-ancestral-genome)  
        * [2: already mapped RNAseq and ancestral genome](#2-already-mapped-RNAseq-and-ancestral-genome) 
        * [3: no RNAseq nor ancestral genome](#3-no-rnaseq-nor-ancestral-genome) 
        * [4:  already annotated genome](#4-already-annotated-genome)
        * [5: other use cases](#5-other-use-case)


# Purpose:
##  sets of scripts to : 
[I - Perform TE and gene prediction](#I---Perform-TE-and-gene-prediction)

[II - Identify synteny blocks and rearragements](#II---Identify-synteny-blocks-and-rearragements)

[III - Plot dS along the genome](#III---Plot-dS-along-the-genome)

[IV - Perform changepoint analysis to identify evolutionary strata](#IV---Perform-changepoint-analysis-to-identify-evolutionary-strata)

<img src="https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig1.png" width = "490" heigth = "490">


Installation: 
-------

- [Installation instructions](:/f3e842a69b234bb3bd8f17b014c5ec80)



# Before Launching the workflow

After cloning the pipeline, please work from within it to preserve the architecture  

We recommend that you clone the pipeline ***for each of your new project*** and work within it.   

Keep all projects separated otherwise it will be difficult to recover your results.   

All options and full path to input files **must** be set in the **config** file provided in : `config/config` .

PLEASE, carefully read the user guide below before any attempt at running the workflow.


# How to use: 

### In short : 

simply run:
```
./master.sh --help #to see all options 
```

The script provides several options depending on what you want: 

-o 1: to perform all analyses: TE and gene prediction, synteny analysis with GeneSpace including single copy orthologs inference between sex/mating type chromosomes,  synonymous divergence (d<sub>S</sub>) computation, evolutionary strata inference
 
-o 2: to perform TE and gene prediction, as well as synteny analysis with GeneSpace (no d<sub>S</sub> computation or evolutionary strata inference)

-o 3: to perform synteny analysis with GeneSpace and subsequent analyses 

-o 4: to perform d<sub>S</sub> computation and subsequent analysis 

-o 5: to perform only synteny analysis with GeneSpace 

-o 6: to perform only TE and gene prediction

-o 7: to perform only evolutionary strata inference

-o 8: to perform only the plots subsequent to d<sub>S</sub> computation


## Detailed use cases

As described above, the workflow is launched with the command:
```
./master.sh -o X #with X an option from 1 to 8. 
```

Each option and their requirement in the config file are described below :

An example config file is provided [here](https://github.com/QuentinRougemont/EASYstrata/blob/main/example_data/example.config)


### option 1: From genome annotation to strata inference

Basically this means running the whole workflow.  

To do so, run : 

```./master.sh -o1 2>&1 |tee log```

this will:

*	 Perform TE annotation using **repeatmodeller** (de-novo prediction) and **repeatMasker**
*	 Perform gene prediction on the softmasked genome using **BRAKER** (with or without RNAseq). 
*	 Evaluate the quality of the gene prediction (mainly using **BUSCO**)
*	 Run **GeneSpace** between your genomes (and eventual ancestral genome) to infer broad pattern of synteny. This includes the inference of single copy orthologs by **orthofinder**
*	 Run **minimap** between genes to infer gene synteny
*	Run **paml** to estimate synonymous divergence (d<sub>S</sub>) between the sequences/regions of interest
*	Perform various plots: 
	*		 circos plot,
	*		 d<sub>S</sub> along the genome, 
	*		 ideogram, 
	*		 etc. (example plots are given below) 
*	Infer the most likely number of evolutionary strata using a changepoint analysis and produce plots with the results: 
	*		 estimates most likely number of strata (model weight)  ,
	*		 use Bayes-Factors to assess differences among strata ,
	*		 plot changepoint location and uncertainty, convergence of the MCMC, other diagnostic plots
	*		 plot circos colored based on numbers of strata
	*		 plot ideogram colored based on numbers of strata 
	*		 plot strata numbers along the genome and along ancestral order

The following options allow you to run only certain parts of the workflow.

### option 2: Anotation and synteny inference only

If you want to perform synteny inference and are not interested in evolutionary strata inference, run:

```./master.sh -o 2 2>&1 |tee log```

This will run the first part of the workflow, up to the synteny inference step.

### option 3: Synteny and strata inference from existing data

If you already have the gene prediction of your input genomes, run:   

```./master.sh -o 3 2>&1 |tee log```

This will skip TE and gene annotation, start the workflow at synteny inference and perform all following steps.  

### option 4: Evolutionary strata inference from existing data

If you already have the gene prediction, synteny inference and d<sub>S</sub> computation (for instance if you already ran previous step), run:

```./master.sh -o 4 2>&1 |tee log```

This will start the workflow at the plotting step and perform all following steps.  

This is mostly useful for debugging and for customizing plots.

### option 5: Only synteny analysis 

If you already have the gene prediction of your input genomes and are not interested in evolutionary strata inference, run:   

```./master.sh -o 5 2>&1 |tee log```

This will only run synteny inference with **GeneSpace** and **Minimap** and output some plots.

This is mostly useful for debugging.

### option 6: Genome annotation only

To only perform gene and TE annotation of the genome, run:
```./master.sh -o 6 2>&1 |tee log```

NOTE: if you already have a softmasked genome you can directly provide it and set:  

**annotateTE="NO"** in the config file, this will skip TE annotation

### option 7: Only evolutionary strata inference 

If you already ran the workflow up to the d<sub>S</sub> computation, run:

```./master.sh -o 7 2>&1 |tee log```

This will perform the last step of the workflow: the evolutionary strata analysis.

This is useful to explore various parameter settings in the MCP analysis, for instance adding priors (presented below) or tweaking the order of the scaffolds. 

### option 8: Synteny plots only

If you already have the synteny inference and only want to produce synteny plots, run:  

```./master.sh -o 8 2>&1 |tee log```

This is useful for customizing plots.

Table of options:
| Option: | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|:----:| ---| --- | --- | --- | --- | --- | --- | --- |
| TE prediction | X | X |   |   |   | X |   |   |
| Gene prediction | X | X |   |   |   | X |   |   |
| Orthology and synteny (GeneSpace & Minimap2)| X | X | X |   | X |   |   |   |
| Synteny plots | X |   | X |   |   |   |   | X |
| d<sub>S</sub> computation + plots | X |   | X | X |   |   |   |   |
| Evolutionary strata inference | X |   | X | X |   |   | X |   |
| Evolutionary plots | X |   | X | X |   |   | X |X   |

### Restart at any step 

This workflow is designed to work at any step in the process. In case of bug you can restart it from whenever it crashes (after fixing the bug) and this should work smoothly.


## Input data 

**/!\ The input required will vary strongly based on which steps of the workflow you want to perform.**
Several files are **compulsory** 

All input data, including full path to input files, should be provided in the [**config file**](https://github.com/QuentinRougemont/EASYstrata/blob/main/config/config)

### Basic input (all options)
* **Input genome(s)** - compulsory: This may be one genome assembly containing both sex/mating type chromosomes, or **ideally** two separate haplotype assemblies containing each one of the sex/mating-type chromosomes.
* **list of scaffolds** - compulsory: names of the contigs/scaffolds/chromosomes composing the sex/mating-type chromosomes.
* **ancestral genome** - optional but highly recommended: The genome assembly of a species used as a proxy for the ancestral state. This will allow to plot d<sub>S</sub> along 'ancestral' gene order, and to infer more accurately single copy orthologs.
* **ancestral gene prediction** - compulsory with ancestral genome: gene prediction associated with the ancestral genome 

## Warning: names of fasta and contigs/scaffolds/chromosomes**
We recommend you use short names for your genome assemblies and avoid any special characters apart from underscore.
*example:* species-1.fasta will not be valid in GeneSpace. => Use **species1.fasta** instead.

For chromosome/contig/scaffold, you  **MUST** use standardized IDs including the species name, and avoid any special characters apart from underscore.  
example: **species1_chr1** or **species1_contigX** or **species1_scaffoldZ**
**otherwise the code will failed during renaming steps**


### Input for TE prediction (options 1,2,6 if your input genomes are not already softmasked)
* **TE database** - compulsory: the name of the TE database (some are available online depending on your taxon) 
* **NCBI taxon** - compulsory: a taxon name for NCBI (used with repeatmasker)  
* **TE bed files** -  optional: a pair of bed files containing TE for your region of interest if already available (will be displayed on the circos plots)

### Input for gene prediction (options 1,2,6 if your input genomes are not already annotated)
* **BUSCO lineage name** - compulsory: name of the BUSCO lineage corresponding to your species (the list of busco lineages is available with busco --list-lineage)
* **RNAseq** - optional: RNAseq data for each genome, will improve BRAKER annotation  
* **Protein database** - optional: a database of proteins from related species. Alternatively, orthoDB12 can be used (downloaded automatically)
* **orthoDB12 lineage name** - optional: one of "Metazoa" "Vertebrata" "Viridiplantae" "Arthropoda" "Eukaryota" "Fungi" "Alveolata"

For an example of input files, we provide an [example data folder](https://github.com/QuentinRougemont/EASYstrata/blob/main/example_data)

###  Full config file details:

| option in config | description |
| --- | --- |
| *genome1* | **Compulsory:** Full path to the input assembly, either a genome containing both sex/mating-type chromosomes, or a haplotype containing one of the sex/mating-type chromosomes. (ex: /path/to/genome1.hap1.fa.gz) |
| *haplotype1* | **Compulsory:** Name of the genome1 to be used as a /contig/scaffold/chromosome basename. (ex: genome1.hap1) |
| \[*genome2*\] | **Optional:** In the case where a haplotype containing one of the sex/mating-type chromosomes was provided as '*genome1*'. Full path to the input assembly of the haplotype containing the second sex/mating-type chromosome. |
| \[*haplotype2*\] | **Compulsory if genome2 was provided** Name of the second haplotype. This can be a basename of all chromosomes if a second assembly is avaiable, or the name of the scaffold/contig/chromosomes corresponding to the sex/MAT chromosome (e.g. species_chrY). |
| annotate | **Compulsory:** a string "YES"/"NO" stating wether gene prediction should be performed or not. |
| *RelatedProt* | **Optional:** Full path to a protein database from related species, in fasta format. |
| \[*RNAseqlist*\] | **Optional:** Full path to a .txt file containing the list of RNA-seq data files. |
| \[*bamlist1*\] | **Optional, alternative to *RNAseqlist*:** Full path to a .txt file containing the list of bam files for *genome1* (alignment of RNA-seq data onto the fasta of *genome 1*). |
| \[*bamlist2*\] | **Optional, alternative to *RNAseqlist* if genome2 was provided:** Full path to a .txt file containing the list of bam files for *genome2* (alignment of RNA-seq data onto the fasta of *genome 2*). |
| \[*orthoDBspecies*\] | **Compulsory for gene prediction:** One of: "Metazoa" "Vertebrata" "Viridiplantae" "Arthropoda" "Eukaryota" "Fungi" "Alveolata". Will use a database from **orthoDB** for gene prediction. |
| *fungus* | "YES" or "NO" (default), whether your species is a fungus. |
| annotateTE | **Compulsory:** a string "YES"/"NO" stating wether TE prediction should be performed or not|
| *TEdatabase* | **Compulsory for TE prediction:** Full path to a database of TE for your species/genus, in fasta format. (set to NO if you already have a softmasked genome and TE bed files) |
| *ncbi_species* | **Compulsory for TE prediction:** Name of the ncbi species, list available [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) |
| \[*gtf1*\] | **Optional**. Full path to a .gtf file for an existing gene prediction on *genome1*. |
| \[*gtf2*\] | **Optional** Full path to a .gtf file for an  existing gene prediction on *genome2*. |
| *busco_lineage* | **Compulsory for gene prediction:** Lineage used for **busco** analysis. You can access the list of available lineages by typing `busco --list-dataset`. |
| *interpro* | **Optional:** "YES" or "NO" (default), whether **interproscan** quality check of the genome annotation should be performed, additionally to the busco quality check. Warning: this can take up to several days for a big dataset. |
| \[*ancestral_genome*\] | **Optional (but recommended):** Full path to the genome of the species used as proxy for the ancestral state. |
| \[*ancestral_gff*\] | **Compulsory if '*ancestral_genome*' is given** Full path to the gff (genome annotation) of the species used as proxy for the ancestral state. |
| *scaffolds* | **Compulsory** Full path to the list of focal scaffolds (i.e. the scaffolds composing the sex/mating-type chromosomes). |



# PART BELOW TO BE UPDATED : 


# Details of the worfklow and results

## I - Perform TE and gene prediction

## RNAseq alignment - TE masking - Gene prediction - Quality assessment

### Parameters set in config file

**Align RNA & Annotate:**
set: **annotate=YES**

set: **rnaseq="YES"**   #YES/NO a string stating wether rnaseq data is available or not

set: **RNAseqlist="/full/path/to/rnaseq.list.txt"**  

This will perform alignment of provided RNA-seq data and use it for BRAKER.

**Annotate, use BAM of RNA**  

set: **annotate=YES**

set: **rnaseq="YES"** 

set: **bamlist1="/full/path/to/bam.aligned_on_genome1.txt"**
where  *bam.aligned_on_genome1.txt* is a file txt listing all your bam (full path)
set **bamlist2="/full/path/to/bam.aligned_on_genome2.txt"**
	with bam aligned on haplotype2

This will use provided BAM of aligned RNA-seq and use it for BRAKER.

**Annotate, no RNA**  
set **rnaseq="NO"** 
This will perform genome annotation without using RNA information.

**Skip**
set: **annotate=NO**  
set: **gtf1="/full/path/to/haplotype1.gtf"**  
set: **gtf2="/full/path/to/haplotype2.gtf"**  

If you have already annotated your genome, will use provided gtf for the following steps, skipping genome annotation.
set: 

## Operations of step I

### 1\. Alignment of RNA-seq data 

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

**if RNAseq is used:** The final genome annotation  is the output from TSEBRA.  
**without RNAseq*:** The final genome annotation is the best protein-based braker round, as evaluated with busco.

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

- Calculation of d~S~ (& dN) using **PAML**

**NOTE ON GENE NAME**
PAML will fail if special characters occur in the input fasta file, or **if the length of a gene name in the fasta header is above 32 characters.** 
To prevent this, we implemented an automatic renaming procedure to shorten character names and remove special characters.  
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

# working examples 
##  1: RNAseq and ancestral genome 

- [see exemple 1](example_data/example1.md)

## 2: already mapped RNAseq and ancestral genome 

- [see exemple 2](example_data/example2.md)

## 3: no RNAseq nor ancestral genome 

- [see exemple 3](example_data/example3.md)

## 4: already annotated genome  

- [see exemple 4](example_data/example4.md)

## 5: other use cases 

describe here all possibles combinations 

# --------------------------------------------------------------------------

# list of operations and tools


| __Operation__                     |  __Tools__                         |  __data type__  | 
|:---------------------------------:|:------------------------------:|:-----------:| 
| __read trimming__                |  Trimmomatic                   | RNAseq         | 
| __read mapping__                 |  gmap/gsnap                    | RNAseq          | 
| __sorting read__                 |  samtools                      | RNAseq        |
| __mapping quality assement__     |  samtools + R                  | RNAseq        |
| __TE detection and softmasking__ |  RepeatModeler + RepeadMasker  | genome assembly |
| __genome annotation__            |  BRAKER + tsebra               | genome assembly + protein + database |
| __quality assessment__           |  BUSCO + Blast + Inter Pro     | genome prediction |
| __Synteny and HOG__              |  GeneSpace (including OrthoFinder/MCScan) | gene prediction and proteins |
| __cds alignement__               |  muscle + translatorX          | gene prediction (single copy orthologs) | 
| __Ds computation__               |  paml                          | CDS alignment |
| __Ds plot/CIRCOS plot__          |  R                             | Ds and genome information |
| __whole genome alignement__       |  minimap2                      | genome assemblies |
| __gene microsynteny__            |  R                      | single copy orthologs |
| __changepoint analysis__         |  R                      | Ds values and gene order |
