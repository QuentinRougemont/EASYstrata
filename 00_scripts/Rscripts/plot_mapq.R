#!/usr/bin/env Rscript

#date: 01-2024
#Plot samtools depth from mapped read 
#Author: QR

#----- usual checks :
if('dplyr' %in% rownames(installed.packages()) == FALSE)
{install.packages('dplyr', repos='https://cloud.r-project.org') }
if('ggplot2' %in% rownames(installed.packages()) == FALSE)
{install.packages('ggplot2', repos='https://cloud.r-project.org') }
if('data.table' %in% rownames(installed.packages()) == FALSE)
{install.packages('data.table', repos='https://cloud.r-project.org') }
if('ggforce' %in% rownames(installed.packages()) == FALSE)
{install.packages('ggforce', repos='https://cloud.r-project.org') }

#---- load libs: 
library(ggplot2)
library(dplyr)
library(data.table)
library(ggforce)

#---- take input file :
argv <- commandArgs(T)

input <- argv[1]
mapq <- paste0('zcat ', input) 

a <- fread(mapq)

#uncomment this if you want:
m  <- a %>% group_by(V2) %>% summarise(mq = mean(V4))
summary(m$mq)

input2 <- gsub('.txt.gz','', input)

write.table(m, paste0('mean_mapq',input2,'.txt'), quote = F, row.names = F)

p <- a %>%
  ggplot(.) +
  geom_line(aes(x=V3,y=V4)) +
  facet_wrap_paginate(~ V2, ncol = 1, nrow = 8, scales = 'free') +
  ggtitle('mapq') +
  xlab('Genomic position (bins 1 Mb)') +
  ylab('mapq') # + th_plot


#create dir :
if (!dir.exists('02_results/RNA_mapq')){
  dir.create('02_results/RNA_mapq')
}

#export :
for(i in 1:n_pages(p)){
  p_save <-  p +
    facet_wrap_paginate(~ V2, ncol = 1, nrow = 8, scales = 'free', page = i)
  ggsave(plot = p_save, width = 16, height =16, filename = paste0('mapq/',input2, '_', i, '.pdf'))
}
