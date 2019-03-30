#!/bin/bash

# align single and paired end reads to mm10 using bowtie2
# transform the output .sam files to .bam
# sort and index the output .bam files

# define variables
INDEX='mm10/mm10'

# make directory of the alignment output
test ! -d bam && mkdir bam || echo 'Already exists'

# make sample runs file for single end
tail -n+2 runs_modified.csv \
  | grep SINGLE \
  | grep TRUE \
  | awk -F, '{print "bam/"$1".bam" "\t" "fastq/"$3".fastq.gz"}' \
  | awk '{if(a[$1])a[$1]=a[$1]","$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS='\t' \
  > samples_runs_singles.tsv

# processing single end reads
while read l; do
  OUT=`echo $l | cut -f1 -d ' '`
  IN=`echo $l | cut -f2 -d ' '`
  if [ ! -f $OUT ]; then
    bowtie2 \
      --no-unal \
      -p 8 \
      -x $INDEX \
      -U $IN \
      | samtools view -Sb - \
      | samtools sort -o - > $OUT
    samtools index $OUT
    echo $out was created.
  fi
done < samples_runs_singles.tsv

# make sample runs file for paired-end
tail -n+2 runs.csv \
  | grep PAIRED \
  | awk -F, '{print "bam/"$1".bam" "," "fastq/"$2"_1.fastq.gz" " " "fastq/"$2"_2.fastq.gz"}' \
  | awk -F, '{if(a[$1])a[$1]=a[$1]" "$2; else a[$1]=$2;}END{for (i in a)print i, a[i];}' OFS=, \
  > samples_runs_paired.csv

# processing paired end reads
while read l; do
  OUT=`echo $l | cut -d ',' -f1`
  IN=`echo $l | cut -d ',' -f2`
  IN1=`echo $IN | cut -b1- | tr ' ' '\n' | grep '_1' | paste -sd ','`
  IN2=`echo $IN | cut -b1- | tr ' ' '\n' | grep '_2' | paste -sd ','`
  if [ ! -f $OUT ]; then
    bowtie2 \
      --no-unal \
      -p 8 \
      -x $INDEX \
      -1 $IN1 -2 $IN2 \
      | samtools view -Sb - \
      | samtools sort -o - > $OUT
    samtools index $OUT
    echo $out was created.
  fi
done < samples_runs_paired.csv
