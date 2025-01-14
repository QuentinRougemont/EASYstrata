#!/usr/bin/env Rscript

#Purpose:  Script to plot and filter Ds values based on various metadata ##########
#Author: QR
#Date: 12-05-23

#INPUTs -- 3 arguments needed:  
# 1 - bed files for the haplotype1
# 2 - bed files for the haplotype2
# 3 - a 3 column file containing ancestral proxy information for plotting:
#     column1: the haplotype name of the proxy to plot along (ancestral proxy or other)
#     column2: the_chromosome_name to use as proxy (X or ancestral species chromosome id) 
#     column3: a string about chromosome orientation: "N" or "R" for "Normal" or Reversed (R): 

# optional:
# 4 - bed file for the ancestral species 

# - paml results files will be read automatically if previous steps were sucessfull

argv <- commandArgs(T)
if (argv[1]=="-h" || length(argv)==0){
    cat(paste0("\033[0;41m","run script as:\n\tRscript ./03.plot_paml.R haplotype1 haplotype2 chromosome [optional: ancestral_sp]\n\n","\033[0m","\n"))
    cat(paste0("\033[0;42m","COMPULSORY input files should be provided in this exact order\n","\033[0m"))
    cat("\t* haplotype1 basename\n")
    cat("\t* haplotype2 basename\n")
    cat("\t* chromosome: a 3 columns tab separated file containing\n
        \t\tcolumn1: the haplotype name of the proxy to plot along
        \t\tcolumn2: the_chromosome_name to use as proxy (X or ancestral species)
        \t\tcolumn3: a string about chromosome orientation: 'N' or 'R' for 'Normal' or 'Reversed' (R)\n
        example:
        Mlag129A1\tMalg129A1_contig8\tN
        Mlag129A1\tMalg129A1_contig11\tR\n\n")
    cat(paste0("\033[0;44m","optionally:\n","\033[0m"))
    cat("\t*ancestral_sp : basename of the ancestral species\n")
}else{

    #--------------- check if library are installed -------------------------------#
    libs <- c('dplyr','ggplot2','magrittr','cowplot','wesanderson', 'viridis','ggrepel',
              'ggbreak','tidyr','patchwork')
    install.packages(setdiff(libs, rownames(installed.packages())), repos="https://cloud.r-project.org" )
    
    #---------------- load libraries ---------------------------------------------#
    invisible(lapply(libs, suppressWarnings(suppressMessages(suppressPackageStartupMessages(library))), character.only = TRUE))
    
    #--------------------- generic function --------------------------------------#
    
    `%nin%` = Negate(`%in%`) #to negate 
    
    #--------------------- fixed parameters --------------------------------------#
    
    ## ------------------ GGPLOT  CUSTOMISATION ----------------------------------#
    th_plot <- theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
      axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
      axis.title.y=element_text(size=18, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
      axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
      strip.text.x = element_text(size=18),
      panel.grid.major = element_blank(),
      plot.title=element_text(family='Helvetica', face='bold', size=22))
    
    
    mycolor2 <-c("#E69F00",  "#0072B2" ,"#5B1A79",  "#CC79A7", "#D55E00")
    
    #--------------- load the data ----------------------------------------------#
    
    #- common results
    # yn00 results:
    dat <- read.table("02_results/paml/results_YN.txt") %>% 
    	set_colnames(., c("Ds","SEDs","Dn","SEDn", "geneX", "geneY"))
    
    if (length(argv)<3) {
    	  stop("At least the name of 2 species to compare and a txt file containing the name and order of scaffold must be supplied.n", call.=FALSE)
    } else if (length(argv)==3) {
    	writeLines("assuming no ancestral species was used\n")
    	sp1 <- argv[1]     # only the basename is needed !
    	sp2 <- argv[2]     # only the basename is needed !
    	chr <- argv[3]     # table with chr\tstatus [Reversed or Not]
    
    	scaf <- read.table(chr, sep="\t") %>% set_colnames(., c("haplo","chr","order"))
    
    	#orthofinder single copy orthologs:
    	single_cp <- read.table("02_results/paml/single.copy.orthologs", sep = "\t") %>% 
    		     set_colnames(., c("ortho","geneX","geneY" ))
    
    } else {
    	writeLines("assuming an ancestral species exist\n")
    	sp1 <- argv[1]     #only the basename is needed !
    	sp2 <- argv[2]     #only the basename is needed !
    	chr <- argv[3]     # table with chr\tstatus [Reversed or Not]
    	#optional 
    	sp3 <- argv[4]     #the basename of the ancestral species !
    	
    	writeLines("load scaffold info\n")
    	scaf <- read.table(chr, sep ="\t") %>% set_colnames(., c("haplo","chr","order"))
    
    	#orthofinder single copy orthologs:
    	writeLines("load single copy info\n")
    	single_cp <- read.table("02_results/paml/single.copy.orthologs", sep = "\t") %>% 
    		     set_colnames(., c("ortho","gene","geneX","geneY" ))
    
    	#link <- argv[6] 
    	#links <- read.table(link, stringsAsFactors = T) 
        #     %>% set_colnames(.,c("gene1", "gene2","status"))	
    	#we will create a vector of color according to the number of status
    	
    	## read Ancestral species :
    	writeLines("load ancestral species info\n")
    	bedAnc <- read.table(paste0("genespace/bed/",sp3, ".bed", sep = "" )) %>% 
    		set_colnames(., c("scaff","start","end","gene"))
    
    }
    
    ## sp1 + sp2
    bedSp1 <- read.table(paste0("genespace/bed/",sp1, ".bed", sep = "" )) %>% 
    	set_colnames(., c("scaff","start","end","gene" ))
    bedSp2 <- read.table(paste0("genespace/bed/",sp2, ".bed", sep = "" )) %>% 
    	set_colnames(., c("scaffSp2","startSp2","endSp2","geneY"))
    
    ## ------------- arrange the data as needed ----------------------------------------------- ##
    Ds_table <- merge(dat, single_cp, by.x = "geneX", by.y = "geneX")
    
    
    #now we must: 
        #1 - reorder according to the scaffold orientation
        #2 - create an incremenantial gene order accordingly:
    
    if (exists("sp3")) {
        #assuming ancestral species was provided
        writeLines("merging all data\n\n")
        all <- merge(bedAnc, scaf, by.x = "scaff", by.y = "chr", sort =F) %>%
            left_join(., Ds_table, by=join_by(gene == gene)  ) %>%
            arrange(scaff, start) %>%
            group_by(desc(scaff)) %>%
            mutate(St = ifelse(order == "N", start, rev(start) )) %>% 
            arrange(St, .by_group = TRUE) %>%
            ungroup() %>%
            mutate(orderchp = seq(1:nrow(.)))
    } else {
        #assuming non ancestral species 
        #plotting along the X:
        writeLines("merging all data\n\n")
        all <- merge(bedSp1, scaf, by.x = "scaff", by.y = "chr") %>%
            left_join(., Ds_table, by=join_by(geneX == gene) ) %>%
            arrange(scaff, start, sort =F) %>%
            group_by(scaff) %>%
            mutate(St = ifelse(order == "N", start, rev(start) )) %>% 
            arrange(St, .by_group = TRUE) %>%
            ungroup() %>%
            mutate(orderchp = seq(1:nrow(.)))
    }
    
    
    #Ds values above 0.3 will be considered as pseudo-genes for the changepoint analyses. 
    allgood <- all %>% filter((Ds < 0.20) %>% replace_na(TRUE))
    
    #export the df for model comparison on the cluster:
    #write.table(df, "02_results/dS.values.forchangepoint.txt", quote =F, row.names = F, col.names = T, sep = "\t")
    write.table(allgood, "02_results/dS.values.forchangepoint.txt", quote =F, row.names = F, col.names = T, sep = "\t")
    #write.table(all, "02_results/dS.values.metadata.txt", quote =F, row.names = F, col.names = T, sep = "\t")
    allgood2 <- na.omit(allgood) %>% mutate(orderchp = seq(1:nrow(.)))
    
    write.table(allgood2, "02_results/dS.values.forchangepoint_noNA.txt",
                quote =F, row.names = F, col.names = T, sep = "\t")
    
    ########################## make plot now using ALL GENES #######################################
    
    writeLines("making some plots.....\n")
    all$scaff<- factor(all$scaff, levels = c(sort(unique(all$scaff), decreasing=T)))

    Fig1A <- all  %>%   #we plot the D dataframe to obtain the Ds along the order
      ggplot(., aes(x = St, y = Ds )) +
      geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
      facet_wrap(~scaff, scale="free_x") +
      geom_point( size = 1) + 
      theme_classic() +
      scale_y_break(c(1,max(all$Ds, na.rm = T)-0.1)) +
      ylim(c(0,max(all$Ds, na.rm=T)+.1)) +
      xlab("position along chr") +
      ylab( expression(italic("Ds"))) +
      th_plot + theme(legend.position = "none") +
      ggtitle("A") 
    
    Fig1B <- all %>%   #we plot the D dataframe to obtain the Ds along the order
      filter(Ds < 1) %>%
      ggplot(., aes(x = orderchp, y = Ds, colour = scaff)) +
      geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
      geom_point( size = 1) + 
      theme_classic() +
      ylim(c(0,0.4)) +
      xlab("order along reference") +
      ylab( expression(italic("Ds"))) +
      th_plot + theme(legend.position = "none") +
      scale_color_manual(values=wes_palette(n=2, name="GrandBudapest1")) +
      ggtitle("B") 
    
    #create dir if not present:
    if (!dir.exists("02_results/dsplots")){
      dir.create("02_results/dsplots")
    }
    
    
    patch <- Fig1A / Fig1B 
    
    pdf(file = "02_results/dsplots/Ds.pdf",14,8)
    patch 
    dev.off()
    
    writeLines("-------------------------------------------------------")
    writeLines("------- constructing graph with gene order-------------\n")
    writeLines("-------------------------------------------------------")
    
    #only if we have an ancestral reference, otherwise it is a bit meaningless
    
    if(length(argv)==4){
    	# --- now add the order
    	colnames(bedSp1) <- c("scaffSp1","startSp1","endSp1","geneX") 
    	ordSp1<- merge(all, bedSp1, by.x = "geneX", by.y = "geneX", sort = F) %>%
      		group_by(scaff) %>%
      		filter(n()>2) %>%
    	  ungroup() %>%
    	  mutate(rankA1 = dense_rank(startSp1))
    	
    	ordSp2 <-  merge(all, bedSp2, by.x = "geneY.y", by.y = "geneY")%>%
    	  group_by(scaff) %>%
    	  filter(n()>2) %>%
    	  ungroup() %>%
    	  mutate(rankA1 = dense_rank(startSp2))
    	
    	
    	pordSp1 <- ggplot(ordSp1, aes(x = orderchp, y = rankA1, colour = scaffSp1 )) +
    	  geom_point( size = 2) + 
    	  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
    	  theme_classic() +
    	  #ylim(c(0,0.5)) +
    	  xlab("order along reference") +
    	  ylab( expression(italic("gene rank in A1"))) +
    	  th_plot + theme(legend.position = "none") +
    	  theme(axis.text.y=element_blank(),
    	        axis.ticks.y=element_blank() 
    	  ) +
    	  scale_color_viridis(discrete=TRUE) +
          ggtitle("C") 
    	  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1"))   
    	
    	pordSp2 <- ggplot(ordSp2, aes(x = orderchp, y = rankA1, colour = scaffSp2 )) +
    	  geom_point( size = 2) + 
    	  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
    	  theme_classic() +
    	  #ylim(c(0,0.5)) +
    	  xlab("order along reference") +
    	  ylab( expression(italic("gene rank in A2"))) +
    	    th_plot + theme(legend.position = "none") +
    	  theme(axis.text.y=element_blank(),
    	        axis.ticks.y=element_blank() 
    	  ) +
    	  scale_color_viridis(discrete=TRUE) +
          ggtitle("D")
    	  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1"))   
    	
        patch <- Fig1A / Fig1B / pordSp1 / pordSp2 + 
              plot_layout(heights = c(4,4,3,3))
            
    	pdf(file = "02_results/dsplots/Ds_and_arrangements.pdf",18,20)
    	#print(plot_grid(Fig1A, Fig1B, pordSp1, pordSp2, 
            #labels="AUTO", ncol = 1, rel_heights = c(1,1,0.9,0.9)) )
        print(patch)
    	dev.off()
    }
}
