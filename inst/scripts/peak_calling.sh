#!/bin/bash

# call peaks using macs2
# build signal track files

# make directory of the peaks and tracks
test ! -d peaks && mkdir peaks || echo 'Already exists'

# make control sample table
cat samples.csv | cut -d ',' -f1,8 | tail -n+2 > sample_control.csv

# calling peaks
while read l;
do
  BAM=`echo $l | cut -d ',' -f 1`
  CTL=`echo $l | cut -d ',' -f 2`
  EF=$(printf "peaks/%s_EF.bdg" "$BAM")
  if [ ! -f $EF ]; then
    macs2 callpeak \
    -t bam/$BAM.bam \
    -c bam/$CTL.bam \
    -B --nomodel --SPMR -g mm \
    --outdir peaks/ \
    -n $BAM
    PILEUP=$(printf "peaks/%s_treat_pileup.bdg" "$BAM")
    LAMBDA=$(printf "peaks/%s_control_lambda.bdg" "$BAM")
    macs2 bdgcmp \
    -t $PILEUP \
    -c $LAMBDA \
    -o $EF \
    -m FE
  fi
done < samples_control.csv
