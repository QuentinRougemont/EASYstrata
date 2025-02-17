# EASYstrata - Genome annotation - Synteny - d<sub>S</sub> computation - Changepoint analysis 
====================================================================================

# Requirements

This software is suitable only in linux-like systems (Unfortunately not Windows or MAC) 

# Table of content 

   * [Purpose](#purpose)
   * [Installation](#installation)
   * [Before-launching-the-workflow](#before-launching-the-workflow)
   * [How to use](#how-to-use)
        * [Summary table of options](#summary-table-of-options)
   * [Input data](#input-data)
        * [Basic input](#basic-input)
        * [Input for TE prediction](#input-for-te-prediction)
        * [Input for gene prediction](#input-for-gene-prediction)
   * [Details of the worfklow and outputs](#details-of-the-worfklow-and-outputs)
        * [Operations of step I: TE and gene prediction](#operations-of-step-i-te-and-gene-prediction)
        * [Operations of step II: Identify synteny blocks and rearragements](#operations-of-step-ii-identify-synteny-blocks-and-rearragements)
        * [Operations of step III: Plot d<sub>S</sub> along the genome](#operations-of-step-iii-plot-ds-along-the-genome)
        * [Operations of step IV: Perform changepoint analysis to identify evolutionary strata](#operations-of-step-iv-perform-changepoint-analysis-to-identify-evolutionary-strata)
   * [Working examples](#working-examples)


# Purpose:
##  sets of scripts to : 
[I - Perform TE and gene prediction](##operations-of-step-i:-te-and-gene-prediction)

[II - Identify synteny blocks and rearragements](##operations-of-step-ii:-identify-synteny-blocks-and-rearragements)

[III - Plot d<sub>S</sub> along the genome](##operations-of-step-iii:-Plot-ds-along-the-genomee)

[IV - Perform changepoint analysis to identify evolutionary strata](##operations-of-step-iv:-perform-changepoint-analysis-to-identify-evolutionary-strata)

<img src="https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig1.png" width = "490" heigth = "490">


# Installation: 

- [Installation instructions](INSTALL/INSTALL.md)


# Before launching the workflow

Clone the workflow, then please work from within it to preserve the architecture.

We recommend that you clone the pipeline ***for each of your new project*** and work within it, and to keep all projects separated otherwise it will be difficult to recover your results.   

All options and full paths to input files **must** be set in the **config** file provided in : `config/config` .

PLEASE, carefully read the user guide below before any attempt at running the workflow.

Note that this workflow uses several softwares and packages, notably **BREAKER**, **TSEBRA**, **GeneSpace**, **PAML**, **R mcp**. We recommend that you read their corresponding manuals before launching the workflow.


# How to use: 

The first step is to provide the path to your input files and choose settings in the [config file](config/config). An example config file is provided [here](https://github.com/QuentinRougemont/EASYstrata/blob/main/example_data/example.config)

To launch the workflow, simply run:
```
./master.sh -o X #with X an option from 1 to 8. 
```

There are several options which allow you to choose which steps of the workflow you wish to run. This allows the workflow run from any step in the process. In case of bug you may restart it from whenever it crashes (after fixing the bug) and it should work smoothly.

```
./master.sh --help #to see all options 
```

```./master.sh -o 1 2>&1 |tee log```

**All steps:** performs all steps of the workflow, i.e. gene prediction, synteny analysis with GeneSpace including single copy orthologs inference between sex/mating type chromosomes,  synonymous divergence (d<sub>S</sub>) computation, evolutionary strata inference and production of various plots


The following options allow you to run only certain parts of the workflow.

```./master.sh -o 2 2>&1 |tee log```

**Steps I and II:** performs only gene prediction and synteny analysis with GeneSpace (no d<sub>S</sub> computation or evolutionary strata inference)

```./master.sh -o 3 2>&1 |tee log```

**Steps II to IV:** performs synteny analysis with GeneSpace and subsequent analyses : useful if you already have annotated your genome (either from running this pipeline or any other annotation tools)

```./master.sh -o 4 2>&1 |tee log```

**Steps III to IV:** performs d<sub>S</sub> computation and subsequent analysis : useful if you already ran the synteny analysis with GeneSpace, and for customizing the plots produced at step III

```./master.sh -o 5 2>&1 |tee log```

**Step II:** performs only the synteny analysis with GeneSpace 

```./master.sh -o 6 2>&1 |tee log```

**Step I:** performs only gene prediction

```./master.sh -o 7 2>&1 |tee log```

**Step IV(G and H):** performs only evolutionary strata inference and the production of various plots: useful if you already ran the synteny analysis with GeneSpace and the d<sub>S</sub> computation with PAML. This option is useful and recommanded to explore various parameter settings in the MCP analysis, for instance adding priors or tweaking the order of the scaffolds. 

```./master.sh -o 8 2>&1 |tee log```

**Step IV(H):** performs only the plots subsequent to d<sub>S</sub> computation: useful if you already ran the rest of the workflow and want to customize your plots

### Summary table of options
| Option: | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|:----:| ---| --- | --- | --- | --- | --- | --- | --- |
| I.Gene prediction | X | X |   |   |   | X |   |   |
| II.Orthology and synteny (GeneSpace & Minimap2)| X | X | X |   | X |   |   |   |
| II.Synteny plots | X |   | X |   |   |   |   | X |
| III.d<sub>S</sub> computation + plots | X |   | X | X |   |   |   |   |
| IV.Evolutionary strata inference | X |   | X | X |   |   | X |   |
| IV.Evolutionary strata plots | X |   | X | X |   |   | X |X   |

The different steps of the worflow are detailed [below](#details-of-the-worfklow-and-outputs)

# Input data

**/!\ The input required will vary strongly based on which steps of the workflow you want to perform.**
Several files are **compulsory** 

Again, all input data, including full path to input files, should be provided in the [**config file**](https://github.com/QuentinRougemont/EASYstrata/blob/main/config/config)

### Basic input
all options
* **Input genome(s)** - compulsory: This may be one genome assembly containing both sex/mating type chromosomes, or **ideally** two separate haplotype assemblies containing each one of the sex/mating-type chromosomes.
* **list of scaffolds** - compulsory: names of the contigs/scaffolds/chromosomes composing the sex/mating-type chromosomes.
* **ancestral genome** - optional but highly recommended: The genome assembly of a species used as a proxy for the ancestral state. This will allow to plot d<sub>S</sub> along 'ancestral' gene order, and to infer more accurately single copy orthologs.
* **ancestral gene prediction** - compulsory with ancestral genome: gene prediction associated with the ancestral genome 

#### Warning: names of fasta and contigs/scaffolds/chromosomes
We recommend you use short names for your genome assemblies and avoid any special characters apart from underscore.
*example:* species-1.fasta will not be valid in GeneSpace. => Use **species1.fasta** instead.

For chromosome/contig/scaffold, you  **MUST** use standardized IDs including the species name, and avoid any special characters apart from underscore.  
example: **species1_chr1** or **species1_contigX** or **species1_scaffoldZ**
**otherwise the code will failed during renaming steps**

### Input for TE prediction
options 1,2,6 if your input genomes are not already softmasked
* **TE database** - compulsory: the name of the TE database (some are available online depending on your taxon) 
* **NCBI taxon** - compulsory: a taxon name for NCBI (used with repeatmasker)  
* **TE bed files** -  optional: a pair of bed files containing TE for your region of interest if already available (will be displayed on the circos plots)

### Input for gene prediction
options 1,2,6 if your input genomes are not already annotated
* **BUSCO lineage name** - compulsory: name of the BUSCO lineage corresponding to your species (the list of busco lineages is available with busco --list-lineage)
* **RNAseq** - optional: RNAseq data for each genome, will improve BRAKER annotation  
* **Protein database** - optional: a database of proteins from related species. Alternatively, orthoDB12 can be used (downloaded automatically)
* **orthoDB12 lineage name** - optional: one of "Metazoa" "Vertebrata" "Viridiplantae" "Arthropoda" "Eukaryota" "Fungi" "Alveolata"

Full details on the options in config file are listed in [config folder](https://github.com/QuentinRougemont/EASYstrata/blob/main/config/config_details.md)
For an example of input files, we provide an [example data folder](https://github.com/QuentinRougemont/EASYstrata/blob/main/example_data)

# Details of the worfklow and outputs

## list of operations and tools
| __Operation__                     |  __Tools__                         |  __data type__  | 
|:---------------------------------:|:------------------------------:|:-----------:| 
| __I.A.read trimming__                |  Trimmomatic                   | RNAseq         | 
| __I.A.read mapping__                 |  gmap/gsnap                    | RNAseq          | 
| __I.A.sorting read__                 |  samtools                      | RNAseq        |
| __I.A.mapping quality assement__     |  samtools + R                  | RNAseq        |
| __I.B.TE detection and softmasking__ |  RepeatModeler + RepeadMasker  | genome assembly |
| __I.C.genome annotation__            |  BRAKER + tsebra               | genome assembly + protein + database |
| __I.C.quality assessment__           |  BUSCO + Blast + Inter Pro     | genome prediction |
| __II.D1.whole genome alignement__       |  minimap2                      | genome assemblies |
| __II.D1.gene microsynteny__            |  R                      | single copy orthologs |
| __II.D2.Synteny and orthogroups__              |  GeneSpace (including OrthoFinder/MCScan) | gene prediction and proteins |
| __III.E.cds alignement__               |  muscle + translatorX          | gene prediction (single copy orthologs) | 
| __III.F.d<sub>S</sub> computation__               |  paml                          | CDS alignment |
| __III.F. d<sub>S</sub> plot/CIRCOS plot__          |  R                             | Ds and genome information |
| __IV.G.changepoint analysis__         |  R                      | d<sub>S</sub> values and gene order |

## Operations of step I: TE and gene prediction

### A\. Alignment of RNA-seq data (optional)

:pencil: Corresponding script: `00_scripts/launch_rnaseq.sh`

- Reads trimming using **trimmomatic**, automatic detecting whether the data is Single-End or Paired-End
- Creation of database for **gsnap** using **gmap**
- Alignment using **gsnap**
  :pencil: Corresponding scripts: `00_scripts/03_gsnap_SE.sh` for Single-End or `00_scripts/03_gsnap_SE.sh` for Paired-End
- Mapping quality assessment and plot (sequencing depth and MAPQ)
  plots and output: `haplo1/04_mapped/Depth/` and `haplo2/04_mapped/Depth`

![depth.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/depth.png)
**Figure 1:** example RNAseq depth

### B\. TE discovery and masking

:pencil: Corresponding script: `00_scripts/launch_step05_to_08.sh`

- *De novo* repeat identification using **repeatmodeler** on the input genome(s)
- Genome masking using **repeatmasker** with ncbi TE Library, and custom dataset if provided

### C\. Genome annotation, quality assessment and filtering

:pencil: Corresponding script: `00_scripts/06_braker.sh`

- Five successive rounds of gene prediction based on OrthoDB (and custom protein database if provided) with **BRAKER**
- **If RNA-seq data was provided:** one round of gene prediction based on RNA-seq data with **BRAKER** 
- Quality assessment and reports production for each round of gene prediction
  This uses the raw **BRAKER** hintsfile, and **BUSCO** (:pencil: corresponding script: `00_scripts/07_busco_after_braker.sh`)  
  The report includes number of genes, number of introns per gene, gene support, number of complete genes and various histograms useful for error checking.
- **If RNA-seq data was provided**: combination of the best protein-based gene model and the RNA-seq-based gene model using **TSEBRA**
  **Warning:** TSEBRA parameters *intron_support* and *stasto_support* are set to 0 in this workflow (default in TSEBRA: 1 and 2 respectively). This means that only overlapping genes between the two gene models will be filtered. You can change this parameter and others to adjust to your desired level of stringency in the TSEBRA config file: `config/default.cfg` Please read the [TSEBRA manual](https://github.com/Gaius-Augustus/TSEBRA) before running the script.
- **If no RNA-seq data was provided:** genome annotation is set as the best protein-based gene model, as evaluated with **BUSCO**
- Reshaping of genoma annotation: scaffold name is inserted in gene names (facilitating dowonstream analyses) and the longest transcript of each gene is kept (necessary for single copy ortholog identification)
  :pencil: Corresponding script: `00_scripts/08_braker_reshaping.sh`
- Final genome annotation quality assessment using **Blast** against **Uniprot** database, **BUSCO**, and optionally **InterProScan**
  **Warning** running **InterProScan** is time-consuming and can take up to several days. By default the option is turned off in this workflow. It can be turned on in the config file.


## Operations of step II: Identify synteny blocks and rearragements

### D-1\. Minimizer alignment and plots of target region

:pencil: Corresponding script: `00_scripts/11_run_genesSpace_paml_ideogram.sh`

- Alignment between the two haplotypes using **minimap2**
  If the scaffolds of interest were provided in the *scaffold* table, only those will be aligned. Otherwise, the whole haplotypes will be aligned.  
- **If ancestral genome was provided: alignment between the two haplotypes and ancestral genome using **minimap2**
- Construction of synteny plot on the focal scaffolds using **pafR**
- **If no scaffold table was provided:** construction of whole genome dotplot using **pafR**

![Fig2.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig2.png)
**Figure 2: pafr** synteny plot based on **minimap** between the focal scaffolds

![Fig3.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig3.png)
**Figure 3:** dotplot between the two haplotypes

### D-2\. Ortholog reconstruction

- Launching of **GeneSpace**
\- Identification of single copy orthologs with **OrthoFinder**  
\- **If no scaffold table was provided:** Production of dotplot and a riparian plot for the whole genome 
\- Production of a riparian plot on focal scaffolds  
For more information, consult the [GeneSpace readme](https://github.com/jtlovell/GENESPACE).
- Circos plots production using the R package **circlize**, between the two haplotypes, and the ancestral genome (if provided).
- Micro-synteny plot production using the R package **Rideogram**, examples [below](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig10.png)
:pencil: Corresponding script: `00_scripts/Rscripts/04.ideogram.R`

![Fig4.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig4.png)
**Figure 4:** Synteny plot by **GeneSpace** showing gene synteny between the ancestral genome (ancestral_sp) and the two mating types of *Microbotryum lychnidis-dioiciae*


## Operations of step III: Plot dS along the genome

### E\. Single copy orthologs alignment

- Align all coding sequences from the focal scaffolds using **muscle** from within **TranslatorX**

### F\. d<sub>S</sub> calculation and plotting

- Calculation of d<sub>S</sub> (and d<sub>N</sub>) using **PAML**
**Note: PAML fails if special characters occur in the input fasta file, or if the length of a gene name is above 32 characters. To prevent this, we implemented an automatic renaming procedure to shorten character names and remove special characters. After PAML, genes are converted back to their original names, so this should be transparent to users.**
- Plotting of d<sub>S</sub> values on the focal scaffolds (and on the whole genome if no scaffold table was provided), using a custom R script.
:pencil: Corresponding script: `00_scripts/Rscripts/03_plot_paml.R`
If an ancestral genome was provided, d<sub>S</sub> is plotted along that genome, otherwise it will be plotted along one of the two sex/mating type chromosomes.  
**Note: It is possible to modify this R script to adapt the plotting options to your needs (for instance position and direction of scaffolds).**

![Fig5.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig5.png)
**Figure 5:** A) d<sub>S</sub> values along the position in the ancestral chromosomes. B) d<sub>S</sub> values along the ancestral gene order after inverting the chromosomes and removing the large autosomal part on contig 8.
C) and D) gene rank as a function of ancestral gene order, in genome1 and genome2 respectively. 

- Plotting of circos plots of the focal scaffolds, with links between single copy orthologs genes, using the R package **circlize**
d<sub>S</sub> Corresponding script: `00_scripts/Rscripts/05_plot_circos.R` 
TE and gene density bed files can be provided as arguments. By default fused autosome will be plotted but these can be removed from the contig list.
**Note: It is possible to modify this R script to adapt the plotting options to your needs (for instance position and direction of scaffolds).**

put **figure4 panel B** here
**Figure 6:** Circos plot between the ancestral genome and one of the focal haplotype (left part) and circos plot between the two focal haplotypes (mating types). The most external track displays the position of gene of interest (red, green and light blue) and the centromeres (in purple). The middle track in lightblue displays gene density. The most interior track in green displays TE density. Red and darkblue interior links display single copy orthologs links. Note that the same figure is produced, with inner links colored according to discrete quantile values of d<sub>S</sub> [example here](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Figure12.svg)

## Operations of step IV: Perform changepoint analysis to identify evolutionary strata

### G\. Changepoint analyses

:pencil: Corresponding script: '00_scripts/Rscripts/06.MCP_model_comp.R`

** Warning: This step is automatically performed with default gene order. However, it is unlikely that the default gene order is correct, so we strongly suggest that you consult the results of the workflow, especially the d<sub>S</sub> plot. Once you have deciphered clear hypotheses regarding scaffold order and orientation, you may change them manually in the "scaffold" table, then relaunch the workflow at step IV:
```
bash ./master.sh -o 7
```

- Launch of multiple change point (MCP) analysis based on the d<sub>S</sub> values, with the R package **mcp**
This allows the automatic inference of evolutionary strata through the identification of regions with distinctively different levels of d<sub>S</sub>. Note that regions with low or zero d<sub>S</sub> on the leftmost and rigthmost extremities of the chromosomes of interest usuall correspond to the pseudo-autosomal regions and are therefore not true evolutionary strata.
\- Running MCP analyses from 1 to 8 changepoints
\- Plotting changepoints posterior distribution and models convergence

![Fig6.A.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig6.A.png)
**Figure 7:** Results from the changepoint analysis for 3 (panel A) to 8 changepoints (panel F). Each changepoint panel displays the distribution of d<sub>S</sub> values (black dots) along with 25 draws from the joint posterior distribution (grey lines) and 95% highest density interval (red lines). Posterior distributions of the changepoints are shown in blue with one line for each chain.

![Fig6.B.svg](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig6.B.svg)
**Figure 8:** posterior fit of the model parameter and mixing of the chains. Here only 7-changepoint model posterior fit is shown.

- Export of summary tables
The MCP analysis produces many useful information which are extracted and automatically exported in *.txt* tables:

\- `modelX.chpt.txt` (with X being the number of changepoints tested, from 1 to 8): Contains the output of the summary function from MCP. 
    The columns correspond to:
    1. name: name of the parameter (changepont and interval) 
    2. mean: mean value of d<sub>S</sub> and interval (based on gene order) 
    3. lower/upper: lower and upper boundaries
    4. Rhat: Gelman-Rubin convergence diagnostic which is often taken to be acceptable if <1.1
    5. n.eff: effective sample size computed using effectiveSize. Low effective sample sizes can also be revealed by poor mixing in trace plots.

\- `modelchoice.txt` : Contains information from the loo model choice operation.
    The columns correspond to:
    1. elpd_diff 
    2. se_diff 
    3. elpd_loo 
    4. se_elpd_loo 
    5. p_loo 
    6. se_p_loo looic 
    7. se_looic

\- `weights.txt` : Contains the weights of each tested models. Higher weights indicates higher supports.


\- `HypothesisXstrata.txt` (with X being the number of changepoint tested, from 1 to 8): Contains results from hypothesis testing (BayesFactor and posterior probabilities aiming at testing differences among strata).
    Differences in d<sub>S</sub> values among adjacent strata are tested, from left to right in the gene order. Both directionalities are tested, i.e. both the hypothesis that stratum 1 has higher d<sub>S</sub> than stratum 2, and the hypothesis that stratum 2 has higher d<sub>S</sub> than stratum 1 are tested.

\- `classif.sX.$haplo1.$haplo2` (with X being the number of changepoint tested, from 1 to 8): Contains the assignment of single copy orthologs to strata
    1. column1: gene name in the first haplotype 
    2. column2: gene name in the second haplotype
    3. column3: stratum of appartenance 

\- `df.txt` : Contains a summary of all informations 

### H\. Production of figures based on inferred evolutionary strata
These figures are produced for each number of changepoints tested
- Violin plot of the d<sub>S</sub> for each evolutionary strata inferred by the model with results of statistical tests
- d<sub>S</sub> plot along the ancestral gene order, colored according to strata as inferred by the model
- d<sub>S</sub> plot along the ancestral gene position, colored according to strata as inferred by the model
- Ideograms between the ancestral genome and the haplotypes, colored according to strata as inferred by the model
- Circos plots adding to the exterior track the genes colored according to the strata as inferred by the model
 
![Fig7.svg](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig7.svg)

**Figure 9:** example violinboxplot for the two "most likely" models inferred by the loo analysis.
By default plots are constructed for all models (from 1 to 8 changepoints). Default statiscal test from the ggstats plot package are used
assuming parametric tests. 

![Fig8.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig8.png)

**Figure 10:** d<sub>S</sub> values plotted along the ancestral gene order for all possible models from three to eight changepoints  
each point is a gene d<sub>S</sub> value colored according to the strata of assignation. 

![Fig9.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig9.png)

**Figure 11:** d<sub>S</sub> values plotted along the ancestral genome for all possible models from three to eight changepoints  
each point is a gene d<sub>S</sub> value colored according to the strata of assignation

![Fig10.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig10.png)

**Figure 12:**  example ideograms infered for the most likely models here. Links are colored according to their strata of appartenance. 

![Figure11.svg](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Figure11.svg)

**Figure 13:** circos plot with external links based on inferred strata


# Working examples 
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

## 6: changepoint with priors

- [see example 6](example_data/example6.md)
