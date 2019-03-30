#!/bin/bash

# download fastq files from sra ftp

# make directory of the downloaded fastq files
test ! -d fastq && mkdir fastq || echo 'fastq/ is already there.'

# download all files/ ftp links using wget
cat runs.csv \
  | cut -d ',' -f10 \
  | tail -n +2 \
  | xargs -P100 -n 1 wget -q -c -P fastq/
