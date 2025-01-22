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
dp <- paste0('zcat ', input) 

#a <- fread('zcat aln.dp.gz')
a <- fread(dp)

#uncomment this if you want:
m  <- a %>% group_by(V1) %>% summarise(dp = mean(V3))
summary(m$dp)

input2 <- gsub('.dp.gz','', input)

write.table(m, paste0('mean_dp',input2,'.txt'), quote = F, row.names = F)
#filter(m, dp > 10) %>% summarise( meanDP = mean(dp), medianDP = median(dp))

p <- a %>%
  filter(V3 <= quantile(a$V3,0.99) %>%
  ggplot(.) +
  geom_line(aes(x=V2,y=V3)) +
  facet_wrap_paginate(~ V1, ncol = 1, nrow = 8, scales = 'free') +
  ggtitle('depth ') +
  xlab('Genomic position (bins 1 Mb)') +
  ylab('mean DP') # + th_plot


#create dir if not present:
if (!dir.exists("Depth")){
  dir.create("Depth")
}

#export :
for(i in 1:n_pages(p)){
  p_save <-  p +
    facet_wrap_paginate(~ V1, ncol = 1, nrow = 8, scales = 'free', page = i)
  ggsave(plot = p_save, width = 16, height =16, filename = paste0('Depth/',input2, '_', i, '.pdf'))
}

#to do: copy Depth in general results folder
