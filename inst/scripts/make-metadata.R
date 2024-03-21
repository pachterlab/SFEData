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
  Title = paste0("seqFISH mouse gastrulation data (rep ", 1:3, ")"),
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

metadatas$vizgen_out <- data.frame(
  Title = c("Small subset of Vizgen output from human brain cancer data",
            "Small subset of Vizgen output from human brain cancer data, Cellpose output"),
  Description = c("Small subset of Vizgen output from in house unpublished brain cancer data, used to demonstrate and test readVizgen()",
                  "Small subset of Vizgen output from in house unpublished brain cancer data, with CellPose cell segmentation directory, used to demonstrate and test readVizgen()"),
  BiocVersion = "3.19",
  Genome = "GRCh38",
  SourceType = "tar.gz",
  SourceUrl = "https://www.dkfz.de/en/single-cell-sequencing/open-lab.html",
  SourceVersion = "",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = NA,
  DataProvider = "scOpenLab DKFZ",
  Maintainer = "Lambda Moses <dlu2@caltech.edu>",
  RDataClass = "SpatialFeatureExperiment",
  DispatchClass = "FilePath",
  RDataPath = file.path("SFEData", c("vizgen.tar.gz", "vizgen_cellpose.tar.gz"))
)

metadatas$cosmx <- data.frame(
  Title = "Small subset of CosMX output from mouse quarter brain data",
  Description = "Small subset of CosMX output from mouse quarter brain data including the hippocampus downloaded from the Nanostring website, used to demonstrate and test readCosMX()",
  BiocVersion = "3.19",
  Genome = "GRCm38",
  SourceType = "tar.gz",
  SourceUrl = "https://nanostring.com/resources/coronal-hippocampus-and-cortex-basic-data-files/",
  SourceVersion = "",
  Species = "Mus musculus",
  TaxonomyId = "10090",
  Coordinate_1_based = NA,
  DataProvider = "Nanostring",
  Maintainer = "Lambda Moses <dlu2@caltech.edu>",
  RDataClass = "SpatialFeatureExperiment",
  DispatchClass = "FilePath",
  RDataPath = file.path("SFEData", "cosmx.tar.gz")
)

metadatas$xenium_out <- data.frame(
  Title = c("Xenium Onboarding Analysis v1 output from mouse brain",
            "Small subset of Xenium Onboarding Analysis v2 output from human pancreas"),
  Description = c("'Tiny subset' of Xenium mouse brain data from 10X website, generated with Xenium Onboarding Analysis v1, without the zarr files and only with lower resolution images, to demonstrate and test readXenium()",
                  "Small subset of Xenium human pancreas data from 10X website, generated with Xenium Onboarding Analysis v2, without the zarr files and only with lower resolution images, to demonstrate and test readXenium()"),
  BiocVersion = "3.19",
  Genome = c("GRCm38", "GRCh38"),
  SourceType = "Zip",
  SourceUrl = c("https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip",
                "https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Pancreas_FFPE/Xenium_V1_human_Pancreas_FFPE_outs.zip"),
  SourceVersion = "",
  Species = c("Mus musculus", "Homo sapiens"),
  TaxonomyId = c("10090", "9606"),
  Coordinate_1_based = NA,
  DataProvider = "10X Genomics",
  Maintainer = "Lambda Moses <dlu2@caltech.edu>",
  RDataClass = "SpatialFeatureExperiment",
  DispatchClass = "FilePath",
  RDataPath = file.path("SFEData", c("xenium.tar.gz", "xenium2.tar.gz"))
)

metadata <- do.call(rbind, metadatas)

write.csv(metadata, "inst/extdata/metadata.csv", row.names = FALSE)
