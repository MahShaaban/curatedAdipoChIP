#!/bin/bash

SCRIPTS=scripts
LOG=log

all: dir_tree \
	make_index \
	get_annotation \
	align_reads \
	get_regions \
	count_features \
	peak_calling \
	peak_annotating \
	peak_object \
	fastqc

# make directory tree
dir_tree:
	test ! -d $(LOG) && mkdir $(LOG) || echo '$(LOG)/ is already there.'

# get fastq files
download_fastq: $(SCRIPTS)/download_fastq.sh runs.csv
	sh $(SCRIPTS)/download_fastq.sh > $(LOG)/download_fastq.out

# quality chekcs
fastqc: $(SCRIPTS)/fastqc.sh $(SCRIPTS)/fastq_filter fastq/*
	sh $(SCRIPTS)/fastqc.sh > $(LOG)/fastqc.out
	R CMD BATCH --vanilla $(SCRIPTS)/fastq_filter.R $(LOG)/fastq_filter.R.out

# make genome index
make_index: $(SCRIPTS)/make_index.sh
	sh $(SCRIPTS)/make_index.sh > $(LOG)/make_index.out

# get genome annotation
get_annotation: $(SCRIPTS)/get_annotation.sh
	sh $(SCRIPTS)/get_annotation.sh > $(LOG)/get_annotation.out

# align reads
align_reads: $(SCRIPTS)/align_reads.sh fastq/* mm10/* samples_runs_singles.tsv samples_runs_paired.csv
	sh $(SCRIPTS)/align_reads.sh > $(LOG)/align_reads.out

# peak calling
peak_calling: $(SCRIPTS)/peak_calling.sh mm10/chromInfo.txt bam/*
	sh $(SCRIPTS)/peak_calling.sh > $(LOG)/peak_calling.out

# get annotation regions
get_regions: $(SCRIPTS)/get_regions.R mm10/annotation.gtf
	R CMD BATCH --vanilla $(SCRIPTS)/get_regions.R $(LOG)/get_regions.R.out

# get peakset
get_peakset: $(SCRIPTS)/get_peakset.sh
	sh $(SCRIPTS)/get_peakset.sh > $(LOG)/get_peakset.out

# count featuresc
count_features: $(SCRIPTS)/count_features.sh bam/* mm10/*.gtf
	sh $(SCRIPTS)/count_features.sh > $(LOG)/count_features.out

# peak object
peak_object: $(SCRIPTS)/peak_object.R counts/peakset.txt mm10/annotation.gtf samples.csv
	R CMD BATCH --vanilla $(SCRIPTS)/peak_object.R $(LOG)/peak_object.R.out
