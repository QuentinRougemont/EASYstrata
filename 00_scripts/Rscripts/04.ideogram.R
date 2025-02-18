#!/usr/bin/env Rscript

#Purpose:
#script to run rideogram plot 
#input required:  
   # 1 - a list of single copy orthologs for the target species pairs
   # 2 - the two bed files for the target species 
   # 3 - the genomes indexes
   # 4 - optionally: a link file of the type "gene1\tgene2\tstatus" with the gene being the single copy orthologs and a status that will be used for colors


#--------------- check if library are installed -------------------------------#
packages <- c('optparse')
install.packages(setdiff(packages, rownames(installed.packages())), repos="https://cloud.r-project.org" )
    
#---------------- load libraries ---------------------------------------------#
invisible(lapply(packages, suppressWarnings(suppressMessages(suppressPackageStartupMessages(library))), character.only = TRUE))


####-------------------------- INITIALISATION ------------------------------####
#------------- read input from the command line -------------------------------#
option_list <- list(
  make_option(c("-c","--single_copy_orthologs_file"), type="character", default=NULL,
              help="txt file of single copy orthologs  [default %default]",
              dest="sco"),
  make_option(c("-f","--fai_haplo1"), type="character", default=NULL,
              help="samtools index file from haplotype1 (generated from previous steps) [default %default]",
              dest="fai1"),
  make_option(c("-g","--fai_haplo2"), type="character", default=NULL,
              help="samtools index file from haplotype2 (generated from previous steps) [default %default]",
              dest="fai2"),
  make_option(c("-i","--bed_haplo1"), type="character", default=NULL,
              help="OPTIONAL: bed file of genes for haplotype1 (generated from previous steps) [default %default]",
              dest="bed1"),
  make_option(c("-j","--bed_haplo2"), type="character", default=NULL,
              help="OPTIONAL: bed file of genes for haplotype2 (generated from previous steps) [default %default]",
              dest="bed2"),
  make_option(c("-l","--links"), type="character", default=NULL,
              help="OPTIONAL: bed file of regions to highlight (gene/centromere/etc) [default %default]",
              ),
  make_option(c("-d","--ds"), type="character", default=NULL,
              help="OPTIONAL: path to a ds file computed from previous steps [default %default]",
              ),
  make_option(c("-s","--scaffold_orientation"), type="character", default=NULL,
              help="OPTIONAL: a table of scaffold orientation: one scaffold per line and a string N/R for Normal or Reverse [default %default]",
              ),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)

#options(error=traceback)

parser <- OptionParser(usage = "%prog -c single-copy-orthologs -i bed_file_species1 
                       -j bed_file species2 -f fai_haplo1 -g fai_haplo2
                       [options (links/orientation/ds)]", 
                       option_list = option_list)
opt = parse_args(parser)

if(opt$v){
  cat(paste0("script to construct ideogram plot\n\n"))
  cat(paste0("\033[0;41m","COMPULSORY PARAMETERS:","\033[0m","\n")) 
  cat(paste0("\t--single_copy_orthologs (-c): ", opt$sco,"\n"))
  cat(paste0("\t--fai_haplo1: index file of haplotype1 (-f): ", opt$fai1,"\n"))
  cat(paste0("\t--fai_haplo2: index file of haplotype2 (-g): ", opt$fai2,"\n"))
  cat(paste0("\t--bed_haplo1 (-i bed file of gene for haplo1): ", opt$bed1,"\n"))
  cat(paste0("\t--bed_haplo2 (-j bed file of gene for haplo2): ", opt$bed2,"\n"))
  cat(paste0("\033[0;42m","optional parameters:\n","\033[0m"))
  cat(paste0("\t--scaffold_orientation (-s): ", opt$scaffold_orientation,"\n"))
  cat(paste0("\t--links (-l): bed file of links  to highlight", opt$links,"\n"))
  cat(paste0("\t--ds (-d): path to the ds file obtained previously use to color the links", opt$ds,"\n\n"))
}


writeLines("\n~~~~~~ loading data ~~~~~~~\n")
if(!is.null(opt$sco)) {
    #syn <- read.table(opt$sco, header=T, as.is=T, sep='\t')
    single_copy_ortho <- opt$sco
} else {
    stop("synteny table is a compulsory parameter. See script usage (--help)")
}

# import contig informations from .fai index
if(!is.null(opt$fai1)) {
    indexA <- opt$fai1
} else {
    stop("index file is a compulsory parameter. See script usage (--help)")
}
if(!is.null(opt$fai2)) {
    indexB <- opt$fai2
} else {
    stop("index file is a compulsory parameter. See script usage (--help)")
}
if(!is.null(opt$bed1)) {
    bedA = opt$bed1
} else {
    stop("index file is a compulsory parameter. See script usage (--help)")
}
if(!is.null(opt$bed2)) {
    bedB = opt$bed2
} else {
    stop("index file is a compulsory parameter. See script usage (--help)")
}

writeLines("\n~~~~~~ compulsory data loaded ~~~~~~~\n")

#    cat(paste0("\033[0;42m","optional: links_file: tab separated file for coloring the genes by groups\n","\033[0m"))
#    cat("\t(automatically generated by modelcomp to color strata for instance)\n
#        column1 and column2: pairs of single copy ortholog genes
#        column3: a grouping structure to be use for color (could be strata id or anything group of gene of interest)\n
#        example:
#        Mlag129A1_contig_8_g7435.t1     Mcal1250A1_MC13_g7409.t1        strata1
#        Mlag129A1_contig_8_g7436.t1     Mcal1250A1_MC13_g7410.t1        strata1
#        Mlag129A1_contig_8_g7437.t1     Mcal1250A1_MC13_g7411.t1        strata1
#        Mlag129A1_contig_8_g7438.t1     Mcal1250A1_MC13_g7412.t1        strata1
#        \n")

#--------------- check if library are installed -------------------------------#
packages <- c('RIdeogram','dplyr','magrittr','data.table','magrittr')
install.packages(setdiff(packages, rownames(installed.packages())), repos="https://cloud.r-project.org" )

#---------------- load libraries ---------------------------------------------#
invisible(lapply(packages, suppressWarnings(suppressMessages(suppressPackageStartupMessages(library))), character.only = TRUE))

sco <- read.table(single_copy_ortho, sep = "\t") %>% select(-V1) 
writeLines("sco ok")

if(!is.null(opt$links)) {
	writeLines("links provided in the links file will be displayed in colors\n")
	writeLines("link file must contain name of gene for species1, name of ortholog for species2 and a status that will be used for coloring the gene\n")
	#to do: add option to provide a coordinate file with status instead of gene file 

    link <- opt$links  

    baselink <-basename(link)
	links <- read.table(link, stringsAsFactors = T) %>%
        set_colnames(.,c("gene1", "gene2","status"))	
	#we will create a vector of color according to the number of status
}
if(!is.null(opt$scaffold_orientation)) {
    scaforder <- read.table(opt$scaffold_orientation, sep="\t") %>% set_colnames(., c("haplo","chr","order"))
}
if(!is.null(opt$ds)) {
    dsfile <- read.table(opt$ds, h=T) %>% na.omit()
}

sp1 <- basename(gsub(".bed", "",bedA)) 
sp2 <- basename(gsub(".bed", "",bedB)) 

#read bed files
#they will be use to create the jointed file:

if(exists("scaforder")){
    bed1 <- read.table(bedA) %>% 
        merge(.,scaforder, by.x="V1", by.y="chr", sort = F) %>%
        select(-haplo) %>%
        merge(sco, ., by.x="V2", by.y = "V4", sort = F ) %>% 
        select(-V3.x) %>%
        set_colnames(., c("gene1", "contig1", "Start_1", "End_1","order") ) %>%
        mutate(species1 = sp1 ) %>%
        group_by(contig1) %>%
        mutate(Start_1 = ifelse(order == "N", Start_1, rev(Start_1) )) %>%
        mutate(End_1 = ifelse(order == "N", End_1, rev(End_1) )) %>%
        select(species1, gene1, contig1, Start_1, End_1) 

     bed2 <- read.table(bedB) %>%
        merge(.,scaforder, by.x="V1", by.y="chr", sort = F) %>%
        select(-haplo) %>%
        merge(sco, ., by.x="V3", by.y = "V4", sort = F ) %>% 
        select(-V2.x) %>%
        set_colnames(., c("gene2", "contig2", "Start_2", "End_2","order") ) %>%
        mutate(species2  = sp2 ) %>%
        group_by(contig2) %>%
        mutate(Start_2 = ifelse(order == "N", Start_2, rev(Start_2) )) %>%
        mutate(End_2 = ifelse(order == "N", End_2, rev(End_2) )) %>%
        select(species2, gene2, contig2, Start_2, End_2) 
} else {
    bed1 <- read.table(bedA) %>% 
        merge(sco, ., by.x="V2", by.y = "V4", sort = F ) %>% 
        select(-V3.x) %>%
        set_colnames(., c("gene1", "contig1", "Start_1", "End_1") ) %>%
        mutate(species1 = sp1 ) %>%
        select(species1, gene1, contig1, Start_1, End_1) 

    bed2 <- read.table(bedB) %>%
        merge(sco, ., by.x="V3", by.y = "V4", sort = F ) %>% 
        select(-V2.x) %>%
        set_colnames(., c("gene2", "contig2", "Start_2", "End_2") ) %>%
        mutate(species2  = sp2 ) %>%
        select(species2, gene2, contig2, Start_2, End_2) 
}

n1 <- nrow(bed1)
n2 <- nrow(bed2)
writeLines(paste0("there is ", n1, " single copy gene in bed1\n"))
writeLines(paste0("there is ", n2, " single copy gene in bed2\n"))


if(identical(n1,n2)=="FALSE"){
    writeLines("error different number of single copy gene in bed1 and bed2!");
    quit("no")
}

#create a pseudo-link file based on dS values:
if (exists("dsfile")){
print("dsfile loaded")
if ("geneX" %in%  colnames(dsfile)) {
  print("subsetting columns")
  #will work for the ancestral only:
  dsfile0 <- select(dsfile, gene, geneX, Ds) %>%
    set_colnames(., c("gene1", "gene2","Ds"))
  dsfile1 <- select(dsfile, geneX, geneY.x, Ds) %>%
    set_colnames(., c("gene1", "gene2","Ds"))

  
} else {
  dsfile0 <- select(dsfile, gene, geneY.x, Ds) %>% 
    set_colnames(., c("gene1", "gene2","Ds"))
}
    #print(head(dsfile0))
    writeLines(paste0("ds file size is :", nrow(dsfile0)))
    
    #create quantile for link:
    dsfile0$quantile <- factor(findInterval(
      dsfile0$Ds, 
      quantile(dsfile0$Ds, na.rm = T, prob=c(0.33, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95))))
    # quantile(syn_ds$Ds, na.rm = T, prob=c(0.33, 0.5, 0.66, 0.75, 0.9,0.95))))

     dsfile1$quantile <- factor(findInterval(
      dsfile1$Ds, 
      quantile(dsfile1$Ds, na.rm = T, prob=c(0.33, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95))))
   
    colnames(sco) <- c("gene1","gene2")
    links <-  merge(dsfile0, sco) %>% 
      select(gene1, gene2, quantile) %>% 
      set_colnames(.,c("gene1","gene2","status"))
    if (nrow(links) == 0){
      links <-  merge(dsfile1, sco) %>% 
      select(gene1, gene2, quantile) %>% 
      set_colnames(.,c("gene1","gene2","status"))
    
    }
}
    
#-------------- merging bed1 and bed2 - fill colors - create rank to match RIdeogram weird requirement 
#---- rename the rank as species
#---- select wanted columns to match RIdeogram requirements: 
#we will merge the bed1 and bed2 
    
    
if(!exists("links")) {
    print("no color")
    all <- cbind(bed1, bed2) %>% 
    group_by(contig1) %>% 
    filter(n()>4) %>% 
    group_by(contig2) %>% 
    filter(n()>4) %>%
    mutate(fill = 'cccccc') %>%
    as.data.frame() %>%  mutate(Species_1 = dense_rank(contig1)) %>%
    mutate(Species_2 = dense_rank(contig2)) %>% 
    #select(Species_1,Start_1,End_1,Species_2,Start_2,End_2,fill) %>%
    as.data.frame(.)
} else {
    #assuming we have a link file that is provided
    print("creating cols")
    #some cols:
    colS <- c("f1bb7b", "fd6467","5b1a18","5b1a88","4575b4",
              "d67236", "fee0d2" , "edf8b1" ,"636363","fc9272", "d73027")
    links$fill <- rep(colS[1:length(levels(links$status))], c(data.frame(table(links$status))[,2]))
    print("cols ok")

    #now merge:
    all <- cbind(bed1, bed2) %>% 
        group_by(contig1) %>% 
        left_join(links, .,  by = join_by(gene1 == gene1, gene2 == gene2) ) %>%
        filter(n()>4) %>% 
        group_by(contig2) %>% 
        filter(n()>4) %>%
        as.data.frame() %>%  
        mutate(Species_1 = dense_rank(contig1)) %>%
        mutate(Species_2 = dense_rank(contig2)) %>% 
        as.data.frame(.) 
   
    all <- na.omit(all)

}

#create dir if not present:
if (!dir.exists("02_results/ideogram")){
  dir.create("02_results/ideogram")
}


#export the joint bed:
write.table(all, paste0("02_results/ideogram/joint_", sp1, "_" , sp2, ".bed" ) , sep="\t", row.names =F, col.names=T, quote = F)
writeLines("joint bed file succesffuly exported\n")
writeLines(paste0("number of lines in bed:" ,nrow(all)) )

#here it would be important to check that the order of the genes in one or the two species is identical
#to the order of the genes in the sco! 
    
#read index to filter chromosome and create a pseudobed file :
#structure of the pseudobed:
#Chr     Start   End      fill    species size    color
#Chr01   0       101369167 969696  Sconica 12      252525
#create pseudo-index1/
index1 <- read.table(indexA)  %>% 
    select(V1,V2) %>%
    filter(V1 %in% all$contig1) %>% 
    mutate(Chr = V1, Start = 0, End = V2, fill = 969696,species = sp1, size = 12, color = 252525) %>%
    select(-V1, -V2)
  
#create pseudo-index2
index2 <- read.table(indexB)  %>% 
    select(V1,V2) %>%
    filter(V1 %in% all$contig2) %>% 
    mutate(Chr = V1, Start = 0,  End = V2, fill = 969696, species = sp2, size = 12, color = 252525) %>%
    select(-V1, -V2)

writeLines("index created")
#combine contig1 and contig 2
karyo <- bind_rows(index1, index2)
#karyo

#RIdeogram plot genomes of the same size
#this is a probelm if genome size if different. 
gapsize <- karyo %>% group_by(species) %>% summarise(size = sum(End) ) %>% summarise(diff = max(size) - min(size))
sp <- karyo %>% group_by(species) %>% summarise(size = sum(End) ) %>%  filter(size == min(size)) %>% select(species)

small  <- data.frame(Chr = "none",  
                     Start = 0, 
                     End = gapsize$diff, 
                     fill = "#0000FF00",  #"#FF0000CC", # "fffffff", #empty fill
                     species = sp$species, #name of 
                     size = 12, 
                     color = 25252525) 

if(small$species==index1$species[1]) {
karyo <- rbind(index1, small, index2)
} else { karyo <- rbind(index1, index2, small) }


#/!\/!\
#bits of code to rework depending on wether species 1 or species 2 is the one with the smallest genome size 
#if not ordered properly the script will fail
all %<>% select(Species_1,Start_1,End_1,Species_2,Start_2,End_2,fill) 


#if (is.null(opt$links)) {
#    if (!exists("dsfile")) {
    if (!exists("links")) {
        writeLines("creating black and white ideogram")
        ideogram(karyotype = karyo, synteny = all, 
            output=paste0('02_results/ideogram/', sp1,sp2,'.svg'))

        convertSVG(paste0('02_results/ideogram/', sp1,sp2,'.svg', sep=''), 
            file = paste0('02_results/ideogram/', sp1,sp2,'.pdf'), device = "pdf")
#    }
} else if(!is.null(opt$links)) {  
   writeLines("exporting ideogram with other colors")
    #assumming links were provided
    ideogram(karyotype = karyo, 
         synteny = all, 
         output=paste0('02_results/ideogram/', sp1,sp2,'_',baselink, '.svg'))

    convertSVG(paste0('02_results/ideogram/', sp1,sp2,'_',baselink, '.svg', sep=''), 
           file = paste0('02_results/ideogram/', sp1,sp2,'_', baselink,'.pdf'), device = "pdf")
} else if(exists("dsfile")) {  
   writeLines("exporting ideogram with quantile colors")
    #assumming links were provided
    ideogram(karyotype = karyo, 
         synteny = all, 
         output=paste0('02_results/ideogram/', sp1,sp2,'_quantile.svg'))

    convertSVG(paste0('02_results/ideogram/', sp1,sp2,'_quantile.svg', sep=''), 
           file = paste0('02_results/ideogram/', sp1,sp2,'_quantile.pdf'), device = "pdf")
}
