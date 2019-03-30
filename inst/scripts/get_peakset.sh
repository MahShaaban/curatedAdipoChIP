#!/bin/bash

# merge overlapping peaks into a peakset using bedtools

# merge peaks
if [ ! -f "mm10/peakset.bed" ]; then
cat peaks/*.narrowPeak \
  | bedtools sort -i - \
  | bedtools merge -d -1 -c 4,4 -o count_distinct,distinct -i - > mm10/peakset.bed
fi
