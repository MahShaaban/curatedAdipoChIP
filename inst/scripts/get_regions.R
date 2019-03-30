# load required libraries
library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)

# set options
options(scipen=999)

# read annotations into txdb
# get genes from txdb
txdb <- makeTxDbFromGFF('mm10/annotation.gtf')
grg <- genes(txdb)
grg_standard <- keepStandardChromosomes(grg, pruning.mode = 'coarse')

# set lenght of regions around tss
regions <- c(1,2,3,5,10,100) * 1000

# get regeons
map(regions, function(x) {
  fl <- paste0('mm10/', 'regions_', as.character(x), '.gtf')
  if(!file.exists(fl)) {
    grp <- promoters(grg_standard, upstream = x, downstream = x)
    export(grp, con = fl)
  }
})
