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
| \[*TEgenome1=*\] | **Optional**. Full path to a bed file of TE on *genome1*. |
| \[*TEgenome2=*\] | **Optional** Full path to a bed file of TE on *genome2*. |
| \[*TEancestral=*\] | **Optional**. Full path to a bed file of TE for the ancestral genome. |
