#!/usr/bin/env Rscript

#nano script to simply run genespace:
#library(devtools)
#devtools::install_github("jtlovell/GENESPACE")
library(GENESPACE)

#some parameters:
blkSize=5   #2 #10
blkRadius=25 #2 #10
UltraSens=TRUE


gpar <- init_genespace(path2mcscanx="mcpath",  wd = "./",
                       blkSize=blkSize, 
                       blkRadius=blkRadius, 
                       diamondUltraSens=UltraSens
                       )

gpar <- run_genespace(gsParam = gpar) 


#€xtra param
#maxOgPlaces = 8,
#nGaps = 5,
#synBuff = 100,
#arrayJump = ceiling(synBuff/2),
#onlyOgAnchors = TRUE,
#nSecondaryHits = 0,
#nGapsSecond = nGaps * 2,
#blkSizeSecond = blkSize,
#blkRadiusSecond = blkRadius,
#onlyOgAnchorsSelf = TRUE,
#onlyOgAnchorsSecond = FALSE,
#maskBuffer = 500,
#onlySameChrs = FALSE,
#dotplots = "check",
#outgroup = ignoreTheseGenomes,
#nSecondHits = nSecondaryHits,
#synBuffSecond = NULL,
#orthofinderMethod = NULL,

