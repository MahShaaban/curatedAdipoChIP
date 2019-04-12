# loading required libraries
library(ChIPseeker)
library(SummarizedExperiment)
library(GenomicFeatures)
library(AnnotationDbi)
library(tidyverse)
library(devtools)
library(S4Vectors)

# loading data
bam_filenames <- read_lines(system.file('extdata', 'bam_filenames.txt', package = 'curatedAdipoChIP'))
count <- read_tsv(system.file('extdata', 'peakset.txt', package = 'curatedAdipoChIP'),
                  col_names = c('seqname', 'start', 'end', 'n', 'peak', bam_filenames))
samples <- read_csv(system.file('extdata', 'samples.csv', package = 'curatedAdipoChIP'))
runs_modified <- read_csv(system.file('extdata', 'runs_modified.csv', package = 'curatedAdipoChIP'))
#runs <- read_csv(system.file('extdata', 'runs.csv', package = 'curatedAdipoChIP'))
txdb <- loadDb(system.file('extdata', 'txdb.db', package = 'curatedAdipoChIP'))
studies <- read_csv(system.file('extdata', 'studies.csv', package = 'curatedAdipoChIP'))
qc <- read_rds(system.file('extdata', 'qc.rds', package = 'curatedAdipoChIP'))

## make counts matrix
mat <- count %>%
  dplyr::select(6:ncol(count)) %>%
  as.matrix()
colnames(mat) <- str_split(colnames(mat), '/|\\.', simplify = TRUE)[,2]

## make features data
peaks <- count[, 1:3]
#peaks$peak <- List(ll)
peaks$peak <- count$peak %>%
  strsplit(',') %>%
  map(function(x) str_replace_all(x, '_peak_\\d+', '') %>% unique) %>%
  List()

peaks$name <- paste('peak', 1:nrow(peaks), sep = '_')

row_ranges <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)

row_ranges <- annotatePeak(row_ranges,
                           TxDb = txdb,
                           level = 'gene') %>%
  as.GRanges()

ind <- peaks$name %in% row_ranges$name

# makeing colData
samples <- runs_modified %>%
  filter(run_average == 1) %>%
  dplyr::select(id, run, library_layout, instrument_model, run_average) %>%
  group_by(id) %>%
  nest(.key = 'runs') %>%
  left_join(samples)

pd <- tibble(id = colnames(mat)) %>%
  inner_join(samples)

rownames(pd) <- pd$id

# add qc to colData
pd$qc <- pd$runs %>%
  map(function(x) {
    if(all(x$library_layout == 'SINGLE -')) {
      ind <- names(qc) %in% x$run
      List(qc[ind])
    } else {
      ind <- names(qc) %in% paste(x$run, 1:2, sep = '_')
      List(qc[ind])
    }
  }) %>%
  List()

# subset matrix to row_ranges
mat2 <- mat[ind,]
rownames(mat2) <- row_ranges$name

nrow(mat2) == length(row_ranges)
ncol(mat2) == nrow(pd)
all(colnames(mat2) == pd$id)

## make a SummarizedExperiment object
peak_counts <- SummarizedExperiment(assays = list(peak_counts = mat2),
                                    colData = pd,
                                    rowRanges = row_ranges,
                                    metadata = list(studies = studies))

# save object to data/
use_data(peak_counts, overwrite = TRUE)
