[![Travis build status](https://travis-ci.org/MahShaaban/curatedAdipoChIP.svg?branch=master)](https://travis-ci.org/MahShaaban/curatedAdipoChIP)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/MahShaaban/curatedAdipoChIP?branch=master&svg=true)](https://ci.appveyor.com/project/MahShaaban/curatedAdipoChIP)

# curatedAdipoChIP

A Curated ChIP-Seq Dataset of MDI-induced Differentiated Adipocytes (3T3-L1)

# Overview

In this documnet, we introduce the purpose of `curatedAdipoChIP` package,
its contents and its potential use cases. This package is a curated dataset of
ChIP-Seq samples. The samples are MDI-induced pre-adipocytes (3T3-L1) at
different time points/stage of differentiation. The ChIP-antibodies used in 
this dataset are of transcription factors, chromatin remodelers and histone 
modifications. The package document the data collection, pre-processing and 
processing. In addition to the documentation the package contains the scripts
that were used to generated the data in `inst/scripts` and the final
`RangedSummarizedExperiment` object in can be accessed through the
`ExperimentHub`.

# Introduction

## What is `curatedAdipoChIP`?

It is an R package for documenting and distributing a curated dataset. The 
package doesn't contain any R functions.

## What is contained in `curatedAdipoChIP`?

The package contains two different things:

1. Scripts for documenting/reproducing the data in `inst/scripts`
2. Access to a `RangedSummarizedExperiment` object in through the
`ExperimentHub`.

## What is `curatedAdipoChIP` for?

The `RangedSummarizedExperiment` object contains the `peak_counts`, `colData`,
`rowRanges` and `metadata` which can be used for the purposes of conducting 
differential peak binding or gene set enrichment analysis on the cell line 
model.

# Installation

The package can be installed using `BiocManager`.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("curatedAdipoChIP")
```

# Docker image

The pre-processing and processing of the data setup environment is available as
a `docker` image. This image is also suitable for reproducing this document. 
The `docker` image can be obtained using the `docker` CLI client.

```
$ docker pull bcmslab/adiporeg_chip:latest
```

# Citing `curatedAdipoChIP`

For citing the package use:

```r
# citing the package
citation("curatedAdipoChIP")
```
