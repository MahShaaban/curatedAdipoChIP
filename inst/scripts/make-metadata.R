meta <- data.frame(
  Title = "A Curated ChIP-Seq Dataset of MDI-induced Differentiated Adipocytes (3T3-L1)",
  Description = "A curated dataset of publicly available ChIP-sequencing of transcription factors, chromatin remodelers and histone modifications in the 3T3-L1 pre-adipocyte cell line",
  BiocVersion = "3.6",
  Genome = "mm10",
  SourceType = "FASTA",
  SourceUrl = "https://github.com/MahShaaban/curatedAdipoChIP",
  SourceVersion = "March 29 2019",
  Species = "Mus musculus",
  TaxonomyId = "10090",
  Coordinate_1_based = TRUE,
  DataProvider = "SRA",
  Maintainer = "Bioconductor Package Maintainer <mahmoud.s.fahmy@students.kasralainy.edu.eg>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "RDA",
  RDataPath = "curatedAdipoChIP/data/peak_counts.rda",
  Tags = "",
  Notes = ""
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
