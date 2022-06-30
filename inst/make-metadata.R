# Create metadata.csv
metadata <- data.frame(
  Title = c("Visium mouse skeletal muscle 2 days post notexin injury",
            "Small demo subset of mouse skeletal muscle Visium data",
            "Small demo subset of mouse skeletal muscle Visium data 2"),
  Description = c("Space Ranger processed Visium data from 'Large-scale integration of single-cell transcriptomic data captures transitional progenitor states in mouse skeletal muscle regeneration' was downloaded from GEO. The gene count matrix, Visium spot polygons, tissue boundary polygons, and myofiber and nuclei segmentation are stored in a SpatialFeatureExperiment (SFE) object. Morphological metrics of the segmented myofibers and nuclei are also included.",
                  "A small subset of the first dataset in this package for quick demos.",
                  "A second small subset of the first dataset to demonstrate working with multiple samples in the same SFE object."),
  BiocVersion = "3.16",
  Genome = "mm10",
  SourceType = "MTX",
  SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4904759",
  SourceVersion = "",
  Species = "Mus musculus",
  TaxonomyId = "10090",
  Coordinate_1_based = NA,
  DataProvider = "Cornell University",
  Maintainer = "Lambda Moses <dlu2@caltech.edu>",
  RDataClass = "SpatialFeatureExperiment",
  DispatchClass = "Rds",
  RDataPath = "placeholder"
)
write.csv(metadata, "inst/extdata/metadata.csv")
