# Create metadata.csv
metadatas <- list()
metadatas$mckellar <- data.frame(
  Title = c("Visium mouse skeletal muscle 2 days post notexin injury",
            "Small demo subset of mouse skeletal muscle Visium data",
            "Small demo subset of mouse skeletal muscle Visium data 2"),
  Description = c("Space Ranger processed Visium data from 'Large-scale integration of single-cell transcriptomic data captures transitional progenitor states in mouse skeletal muscle regeneration' was downloaded from GEO. The gene count matrix, Visium spot polygons, tissue boundary polygons, and myofiber and nuclei segmentation are stored in a SpatialFeatureExperiment (SFE) object. Morphological metrics of the segmented myofibers and nuclei are also included.",
                  "A small subset of the first dataset in this package for quick demos.",
                  "A second small subset of the first dataset to demonstrate working with multiple samples in the same SFE object."),
  BiocVersion = "3.16",
  Genome = "GRCm38",
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
  RDataPath = file.path("SFEData", c("sfe_vis5a.rds", "sfe_vis5a_small.rds", "sfe_small2.rds"))
)

metadatas$biermann <- data.frame(
    Title = c("Human melanoma brain metastasis slide-seq2 data",
              "Human melanoma extracranial metastasis slide-seq2 data"),
    Description = c("Gene count matrix and bead locations from one melanoma brain metastasis sample as processed by the authors as in the paper 'Dissecting the treatment-naive ecosystem of human melanoma brain metastasis'",
                    "Gene count matrix and bead locations from one melanoma extracranial metastasis sample as processed by the authors as in the paper 'Dissecting the treatment-naive ecosystem of human melanoma brain metastasis'"),
    BiocVersion = "3.16",
    Genome = "GRCh38",
    SourceType = "CSV",
    SourceUrl = c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6025935",
                  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6025946"),
    SourceVersion = "",
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = NA,
    DataProvider = "Columbia University Irving Medical Center",
    Maintainer = "Lambda Moses <dlu2@caltech.edu>",
    RDataClass = "SpatialFeatureExperiment",
    DispatchClass = "Rds",
    RDataPath = file.path("SFEData", c("mbm_slide_seq.rds", "ecm_slide_seq.rds"))
)

metadatas$he <- data.frame(
    Title = "Nanostring FFPE CosMX human NSCLC data",
    Description = "One of the CosMX formalin fixed paraffin embedded (FFPE) example datasets for human non small cell lung cancer (NSCLC, Lung5_Rep1) from the Nanostring website, described in the paper 'High-plex Multiomic Analysis in FFPE at Subcellular Level by Spatial Molecular Imaging'",
    BiocVersion = "3.16",
    Genome = "GRCh38",
    SourceType = "CSV",
    SourceUrl = "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/",
    SourceVersion = "",
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = NA,
    DataProvider = "Nanostring Technologies Inc.",
    Maintainer = "Lambda Moses <dlu2@caltech.edu>",
    RDataClass = "SpatialFeatureExperiment",
    DispatchClass = "Rds",
    RDataPath = file.path("SFEData", "cosmx1.rds")
)

metadatas$janesick <- data.frame(
    Title = paste0("Xenium FFPE human breast cancer data (rep", 1:2, ")"),
    Description = paste0("Example Xenium dataset of formalin fixed paraffin embedded (FFPE) human breast cancer from 10X Genomics, rep", 1:2, ", described in the paper 'High resolution mapping of the breast cancer tumor microenvironment using integrated single cell'"),
    BiocVersion = "3.16",
    Genome = "GRCh38",
    SourceType = "HDF5",
    SourceUrl = "https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast",
    SourceVersion = "",
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = NA,
    DataProvider = "10X Genomics",
    Maintainer = "Lambda Moses <dlu2@caltech.edu>",
    RDataClass = "SpatialFeatureExperiment",
    DispatchClass = "Rds",
    RDataPath = file.path("SFEData", c("xenium1.rds", "xenium2.rds"))
)

metadatas$vizgen <- data.frame(
    Title = "Vizgen MERFISH mouse liver data",
    Description = "This is one of the example datasets from Vizgen's website",
    BiocVersion = "3.16",
    Genome = "GRCm38",
    SourceType = "CSV",
    SourceUrl = "https://console.cloud.google.com/storage/browser/vz-liver-showcase/Liver1Slice1;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false&pli=1",
    SourceVersion = "",
    Species = "Mus musculus",
    TaxonomyId = "10090",
    Coordinate_1_based = NA,
    DataProvider = "Vizgen Inc.",
    Maintainer = "Lambda Moses <dlu2@caltech.edu>",
    RDataClass = "SpatialFeatureExperiment",
    DispatchClass = "Rds",
    RDataPath = file.path("SFEData", "merfish_liver1.rds")
)

metadatas$lohoff <- data.frame(
  Title = paste0("seqFISH mouse gastrulation data (rep", 1:3, ")"),
  Description = paste0(
    "Example seqFISH dataset of tissue sections of mouse embryos at the 8–12 somite stage by the Marioni lab at the Cancer Research UK Cambridge Institute, rep", 1:3, ", described in the paper 'Integration of spatial and single-cell transcriptomic data elucidates mouse organogenesis'"),
  BiocVersion = "3.17",
  Genome = "GRCm38",
  SourceType = "RDS",
  SourceUrl = "https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/",
  SourceVersion = "",
  Species = "Mus musculus",
  TaxonomyId = "10090",
  Coordinate_1_based = NA,
  DataProvider = "Cancer Research UK Cambridge Institute",
  Maintainer = "Kayla Jackson <kaylajac@caltech.edu>",
  RDataClass = "SpatialFeatureExperiment",
  DispatchClass = "Rds",
  RDataPath = file.path("SFEData", c("seqfish_em1.rds", "seqfish_em2.rds", "seqfish_em3.rds"))
)


metadata <- do.call(rbind, metadatas)

write.csv(metadata, "inst/extdata/metadata.csv", row.names = FALSE)
