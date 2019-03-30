#!/bin/bash

# count reads in regions of different length around tss
# count reads in peakset using bedtools multicov
# get the names of the bam files

# the output 'peaks/peakset.txt' goes into `inst/extdata/peakset.txt`
# the output 'bam_filenames.txt' goes into `inst/extdata/bam_filenames.txt`

# count reads in peaks
if [ ! -f "counts/peakset.txt" ]; then
    bedtools multicov \
    -bams ./bam/*.bam \
    -bed mm10/peakset.bed > counts/peakset.txt
    echo "counts/peakset.txt was created."
fi

# make list of file names counted in peaks
ls bam/*.bam > bam_filenames.txt

# define variables
REG=$(ls mm10/ | grep regions | cut -d '.' -f1)

# make directory of counts
test ! -d counts && mkdir counts || echo 'Already exists'

# count reads in gene body
if [ ! -f "counts/gene_body.txt" ]; then
    featureCounts \
    -T 8 \
    -F GTF \
    -t exon \
    -g gene_id \
    -a ./mm10/annotation.gtf \
    -o ./counts/gene_body.txt \
    ./bam/*.bam
    echo "counts/gene_body.txt was created."
fi

# count regions around tss
for i in $REG; do
  GTF=$(printf "./mm10/%s.gtf" "$i")
  OUT=$(printf "./counts/%s.txt" "$i")
  if [ ! -f $OUT ]; then
    featureCounts \
    -T 8 \
    -F GTF \
    -t sequence_feature \
    -g gene_id \
    -a $GTF \
    -o $OUT \
    ./bam/*.bam
    echo "$OUT was created."
  fi
done
