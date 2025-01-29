# small workflow to reproduce the main figures from our paper using directly provided gtf files

# including MCP with priors 


# 1 - cloning the repo

```sh
mkdir example6 ; cd example6
git clone https://github.com/QuentinRougemont/EASYstrata/ .

path1=$(pwd) #store variable to edit in config file

#decompress some data:
gunzip example_data/Mlag129.A1.fa.gz
gunzip example_data/Mlag129.A1.gff.gz


```

# 2 - get the data: 

target species: M. violaceum caroliniana
genomes: 1250a1 (MAT 1) and 1250a2 (MAT 2)
genome and gtf files are deposited on zenodo.  

intermediary file are also provided and include mainly the d<sub>S</sub> files

obtaining the data:

```sh
wget wget https://zenodo.org/records/XXXXXXXX/files/data.tar.gz?download=1

tar zxvf data.tar.gz

cd Mvcaliforniana1250
path2=$(pwd)
cd ../


```

# 3 - setting config files:

the following config file is used:

| option in config | description |
| --- | --- |
| *genome1* | ="/path2/to/EASYstrata/data/Mcal1250A1.fa" |
| *haplotype1* | ="Mcal1250A1" |
| \[*genome2*\] | ="/path2/to/EASYstrata/Mcal1250A2.fa" |
| \[*haplotype2*\] | ="Mcal1250A2" |
| annotate | ="NO" |
| *RelatedProt* | ="" |
| \[*RNAseqlist*\] | ="" |
| \[*bamlist1*\] | ="" |
| \[*bamlist2*\] | ="" |
| \[*orthoDBspecies*\] | ="Fungi" |
| *fungus* | ="" |
| *TEdatabase* | =""|
| *ncbi_species* | ="" |
| \[*gtf1*\] | ="path2/Mcal1250A1.gtf"  |
| \[*gtf2*\] | ="path2/Mcal1250A2.gtf" |
| *busco_lineage* | ="" |
| *interpro* | ="" |
| \[*ancestral_genome*\] |  ="path1/example_data/Mlag129.A1.fa" |
| \[*ancestral_gff*\] | ="path1/example_data/Mlag129.A1.gff" |
| *TEancestral* | ="path1/example_data/Mlag129.A1.TE.bed" |
| *scaffolds* | ="path1/example_data/scaffold.txt" |


click [here](example6.config) to see the file 

copy the file in the config folder:

```sh
cp example_data/example6.config config/config
#set the path:
sed -i 's#path1#$path1#g' config/config
sed -i 's#path2#$path2#g' config/config


```

# 4 - first pass analysis : inferring synteny, dS and strata 

to that aim we will run the option 3: 

```sh
./master.sh -o3 2>&1 |tee log
```

this should provide the GeneSpace results and all other important results in `02_results` folder  


# 5 - refining strata with priors : 

to obtain better estimate of strata location we use the mcp with priors based on previous work as well as observed gap in dS distribution along the ancestral genome:

to do so the following script is used:

```R
example_data/MCP_with_priors_example6.R
```
