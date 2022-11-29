library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(vroom)
library(sf)
library(terra)
library(EBImage)
library(rmapshaper)
library(lwgeom)
library(tidyverse)
library(jsonlite)
library(Matrix)
library(Voyager)
library(R.utils)
library(stringr)
library(glue)
library(rlang)
library(arrow)
library(DropletUtils)
library(BiocParallel)
library(DelayedArray)
# Mouse skeletal muscle Visium dataset==========================

# The data comes from this paper:
# McKellar DW, Walter LD, Song LT, Mantri M et
# al. Large-scale integration of single-cell transcriptomic data captures
# transitional progenitor states in mouse skeletal muscle regeneration. Commun
# Biol 2021 Nov 12;4(1):1280.

# The dataset included in the first version of this version comes from the first
# time point (2 days) after notexin injury. The gene count matrix of all Visium
# spots, the metadata, and the H&E image were downloaded from here:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4904759

# Myofiber segmentation-------------
# Myofibers were manually segmented in LabKit on a 4x downsized image, downsized
# so the image can be loaded into LabKit. The results are originally in a TIFF.
# Here I convert the myofiber masks in the TIFF into vector format.

raster2polygon <- function(seg, keep = 0.1) {
  r <- rast(t(as.array(seg)), extent = ext(0, nrow(seg), 0, ncol(seg)))
  r <- terra::flip(r)
  r[r < 1] <- NA
  contours <- st_as_sf(as.polygons(r, dissolve = TRUE))
  # Use mapshaper to keep polygon contiguity.
  #browser()
  simplified <- ms_simplify(contours, keep = keep, keep_shapes = TRUE, no_repair = FALSE)
  list(full = contours,
       simplified = simplified)
}

myofiber_tiff <- readImage("data/tissue_downsized_5a_myofibers.tiff", as.is = TRUE)
#myofiber_tiff <- fillHull(myofiber_tiff)
# Get rid of holes, somehow fillHull of the whole thing doesn't work
myofiber_poly <- raster2polygon(myofiber_tiff)
invalid_inds <- which(!st_is_valid(myofiber_poly$full))
for (i in invalid_inds) {
  val <- myofiber_poly$full[i,]$lyr.1
  ind <- myofiber_tiff == val
  ind <- fillHull(ind)
  ind2 <- ind > 0
  myofiber_tiff[ind2] <- val
}
# I also realized that sometimes two myofibers have the same label. That needs to be fixed.
# Sometimes the MULTIPOLYGON is from a dangling pixel, which should be removed.
kern <- makeBrush(5, shape = "disc")
multi_inds <- which(st_geometry_type(myofiber_poly$full) == "MULTIPOLYGON")
for (i in multi_inds) {
  val <- myofiber_poly$full[i,]$lyr.1
  ind <- myofiber_tiff == val
  ind2 <- opening(ind, kern)
  ind2 <- bwlabel(ind2)
  myofiber_tiff[ind] <- 0
  for (j in seq_len(max(ind2))) {
    ind3 <- ind2 == j
    if (j == 1L) myofiber_tiff[ind3] <- val
    else {
      myofiber_tiff[ind3] <- max(myofiber_tiff) + j - 1L
    }
  }
}

myofiber_poly2 <- raster2polygon(myofiber_tiff)
all(st_is_valid(myofiber_poly2$simplified))
# Some simplified polygons are invalid due to a protruding line.
myofiber_poly2$simplified <- st_buffer(myofiber_poly2$simplified, 0)
unique(st_geometry_type(myofiber_poly2$full))
# Now all polygons, full and simplified, are valid. No more MULTIPOLYGONs.
# Compute morphological metrics
myofiber_morphology <- computeFeatures.moment(myofiber_tiff)
# Compute Haralick metrics, needs the image
myofiber_he <- readImage("data/tissue_downsized.tif")
myofiber_haralick <- computeFeatures.haralick(myofiber_tiff, myofiber_he)

# Prepare the annotGeometries
# Use the simplified polygons for analyses
myofiber_simplified <- myofiber_poly2$simplified %>%
  mutate(area = st_area(geometry),
         perimeter = st_perimeter(geometry),
         eccentricity = myofiber_morphology[,"m.eccentricity"],
         theta = myofiber_morphology[,"m.theta"], # Already in radians
         sine_theta = sin(theta),
         convexity = area / st_area(st_convex_hull(geometry)))
# There're empty labels, remove them from Haralick results
myofiber_haralick <- myofiber_haralick[rowSums(myofiber_haralick) > 0,]
myofiber_simplified <- cbind(myofiber_simplified, as.data.frame(myofiber_haralick))

# To scale the coordinates up to match coordinates of full image
myofiber_simplified$geometry <- myofiber_simplified$geometry * 4
myofiber_full <- myofiber_poly2$full
myofiber_full$geometry <- myofiber_full$geometry * 4
myofiber_simplified$sample_id <- "Vis5A"
myofiber_full$sample_id <- "Vis5A"

# Nuclei segmentation----------------
# I manually annotated about 3000 nuclei in the H&E image in LabKit as training
# data for StarDist. Then the trained model is used to segment all the nuclei.
# StarDist returns a TIFF mask, which I have converted to polygons with OpenCV.
# I have also computed morphological metrics with OpenCV. The reason why I did
# not use the OpenCV polygons for myofibers is that OpenCV does not preserve
# contiguity, which is less relevant to nuclei.

# For nuclei
nuc_poly <- read_csv("data/cell_contours.csv")
nuc_poly <- nuc_poly[,c("index", "x", "y")]
names(nuc_poly)[1] <- "ID"
nuc_poly <- df2sf(nuc_poly, geometryType = "POLYGON")

nuc_centroids <- read_csv("data/cell_centroids.csv")
nuc_centroids <- df2sf(nuc_centroids, geometryType = "POINT")
# bbox as nuclei was performed on a cropped image of the bouding box of the tissue,
# while coordinates of spots are of the full image.
nuc_bbox <- read_lines("data/tissue_bbox.txt")
nuc_bbox <- str_extract_all(nuc_bbox, "\\d+", simplify = TRUE)
nuc_bbox <- as.integer(nuc_bbox)
nuc_poly$geometry <- nuc_poly$geometry + nuc_bbox[1:2]
nuc_centroids$geometry <- nuc_centroids$geometry + nuc_bbox[1:2]

# Get the morphological metrics
nuc_metrics <- read_csv("data/shape_metrics.csv")
nuc_poly <- cbind(nuc_poly, nuc_metrics)
nuc_centroids <- cbind(nuc_centroids, nuc_metrics)
nuc_poly$sample_id <- "Vis5A"
nuc_centroids$sample_id <- "Vis5A"

# Tissue boundary--------------
# The tissue boundary was found with thresholding in OpenCV.
tb <- read_csv("data/contours_tissue_vis5a.csv")
names(tb)[1] <- "ID"
tb <- df2sf(tb, geometryType = "POLYGON")
tb$geometry <- tb$geometry + nuc_bbox[1:2]
tb$sample_id <- "Vis5A"

# Construct SFE object--------------
barcodes <- readLines("data/barcodes.tsv")
barcodes <- str_remove(barcodes, "-\\d+$")
genes <- read_tsv("data/features.tsv", col_names = c("Ensembl", "symbol", "type"))
mat <- readMM("data/matrix.mtx")
mat <- as(mat, "dgCMatrix")
rownames(mat) <- genes$Ensembl
colnames(mat) <- barcodes
# The authors only provided the Loupe manual alignment JSON
spot_info <- read_json("data/spot_info.json", simplifyVector = TRUE)
spots <- spot_info$oligo
data("visium_row_col")
col_data <- tibble(barcode = barcodes) %>%
  left_join(visium_row_col, by = "barcode") %>%
  mutate(row = row - 1L, # 0 based indexing
         col = col - 1L) %>%
  left_join(spots, by = c("row", "col"))

sfe <- SpatialFeatureExperiment(assays = list(counts = mat),
                                colData = col_data, rowData = genes,
                                sample_id = "Vis5A",
                                spatialCoordsNames = c("imageX", "imageY"),
                                spotDiameter = spots$dia[1])
# Add the geometries
tissueBoundary(sfe) <- tb
annotGeometry(sfe, "myofiber_full") <- myofiber_full
annotGeometry(sfe, "myofiber_simplified") <- myofiber_simplified
annotGeometry(sfe, "nuclei") <- nuc_poly
annotGeometry(sfe, "nuclei_centroid") <- nuc_centroids

# Some basic QC
colData(sfe)$nCounts <- colSums(mat)
colData(sfe)$nGenes <- colSums(mat > 0)
colData(sfe)$prop_mito <- colSums(mat[str_detect(genes$symbol, "^mt-"),]) / colData(sfe)$nCounts
rowData(sfe)$means <- rowMeans(mat)
rowData(sfe)$vars <- rowVars(mat)

sfe <- sfe[rowData(sfe)$means > 0,]
rowData(sfe)$cv2 <- rowData(sfe)$vars/rowData(sfe)$means^2

# Whether spots are on tissue
colData(sfe)$in_tissue <- annotPred(sfe, "spotPoly", "tissueBoundary")

saveRDS(sfe, "data/sfe_vis5a.rds")

# Small subset for function examples
sfe_small <- crop(sfe, xmin = 5000, xmax = 7000, ymin = 13000, ymax = 15000)
saveRDS(sfe_small, "data/sfe_vis5a_small.rds")

# I have only manually annotated myofibers for Vis5A.
# For demonstration purposes, I'll use two different portions of this one
# dataset and pretend that they're different samples.
sfe_small2 <- crop(sfe, xmin = 6000, xmax = 8000, ymin = 7000, ymax = 9000)
colData(sfe_small2)$sample_id <- "sample02"
for (n in names(int_metadata(sfe_small2)$annotGeometries)) {
  int_metadata(sfe_small2)$annotGeometries[[n]]$sample_id <- "sample02"
}
names(int_metadata(sfe_small2)$spatialGraphs) <- "sample02"
# Maybe I need a convenience function to change sample_id.
saveRDS(sfe_small2, "data/sfe_small2.rds")

# Melanoma brain metastasis slide-seq2 data====================

# The data comes from
# Jana Biermann, Johannes C. Melms, Amit Dipak Amin, Yiping Wang, Lindsay A.
# Caprio, Alcida Karz, Somnath Tagore, Irving Barrera, Miguel A.
# Ibarra-Arellano, Massimo Andreatta, Benjamin T. Fullerton, Kristjan H.
# Gretarsson, Varun Sahu, Vaibhav S. Mangipudy, Trang T.T. Nguyen, Ajay Nair,
# Meri Rogava, Patricia Ho, Peter D. Koch, Matei Banu, Nelson Humala, Aayushi
# Mahajan, Zachary H. Walsh, Shivem B. Shah, Daniel H. Vaccaro, Blake Caldwell,
# Michael Mu, Florian Wünnemann, Margot Chazotte, Simon Berhe, Adrienne M.
# Luoma, Joseph Driver, Matthew Ingham, Shaheer A. Khan, Suthee Rapisuwon, Craig
# L. Slingluff, Thomas Eigentler, Martin Röcken, Richard Carvajal, Michael B.
# Atkins, Michael A. Davies, Albert Agustinus, Samuel F. Bakhoum, Elham Azizi,
# Markus Siegelin, Chao Lu, Santiago J. Carmona, Hanina Hibshoosh, Antoni Ribas,
# Peter Canoll, Jeffrey N. Bruce, Wenya Linda Bi, Praveen Agrawal, Denis
# Schapiro, Eva Hernando, Evan Z. Macosko, Fei Chen, Gary K. Schwartz, Benjamin
# Izar,
# Dissecting the treatment-naive ecosystem of human melanoma brain metastasis,
# Cell,
# Volume 185, Issue 14,
# 2022,
# Pages 2591-2608.e30,
# ISSN 0092-8674,
# https://doi.org/10.1016/j.cell.2022.06.007.
# (https://www.sciencedirect.com/science/article/pii/S0092867422007127)

# The gene count matrix and spatial coordinates were downloaded from GEO
mbm <- vroom("GSM6025935_MBM05_rep1_slide_raw_counts.csv")
mbm_coords <- vroom("GSM6025935_MBM05_rep1_slide_spatial_coordinates.csv")
genes <- mbm[,1]
mbm <- as.matrix(mbm[,-1])
mbm <- as(mbm, "dgCMatrix")
rownames(mbm) <- genes$...1
mbm <- mbm[rowSums(mbm) > 0,]
# The genes are already in symbols rather than Ensembl IDs
sfe_mbm <- SpatialFeatureExperiment(assays = list(counts = mbm),
                                    spatialCoords = as.matrix(mbm_coords[,c("xcoord", "ycoord")]))
addQC <- function(sfe, species = "human") {
  mat <- counts(sfe)
  colData(sfe)$nCounts <- colSums(mat)
  colData(sfe)$nGenes <- colSums(mat > 0)
  mt_regex <- if (species == "human") "^MT-" else "^Mt-"
  if (any(str_detect(rownames(mat), mt_regex))) {
      mito_genes <- str_detect(rownames(mat), "^MT-")
      colData(sfe)$prop_mito <- colSums(mat[mito_genes,])/colData(sfe)$nCounts
  }
  rowData(sfe)$means <- rowMeans(mat)
  rowData(sfe)$vars <- rowVars(mat)
  rowData(sfe)$cv2 <- rowData(sfe)$vars/rowData(sfe)$means^2
  sfe
}
sfe_mbm <- addQC(sfe_mbm)
saveRDS(sfe_mbm, "mbm_slide_seq.rds")

# Extracraneal metastasis
ecm <- vroom("GSM6025946_ECM01_rep1_slide_raw_counts.csv")
ecm_coords <- vroom("GSM6025946_ECM01_rep1_slide_spatial_coordinates.csv")
genes <- ecm[,1]
ecm <- as.matrix(ecm[,-1])
ecm <- as(ecm, "dgCMatrix")
rownames(ecm) <- genes$...1
ecm <- ecm[rowSums(ecm) > 0,]

sfe_ecm <- SpatialFeatureExperiment(assays = list(counts = ecm),
                                    spatialCoords = as.matrix(ecm_coords[,c("xcoord", "ycoord")]))
sfe_ecm <- addQC(sfe_ecm)
saveRDS(sfe_ecm, "ecm_slide_seq.rds")

# Xenium FFPE human breast cancer data=========================
# Downloaded from 10x website
# https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
# Rep 1
sce <- read10xCounts("Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5")
counts(sce) <- as(realize(counts(sce)), "dgCMatrix")
cell_info <- vroom("Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz")
cell_poly <- read_parquet("Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_boundaries.parquet")
nuc_poly <- read_parquet("Xenium_FFPE_Human_Breast_Cancer_Rep1_nucleus_boundaries.parquet")
names(cell_poly)[1] <- "ID"
names(nuc_poly)[1] <- "ID"
cells_sf <- df2sf(cell_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
nuc_sf <- df2sf(nuc_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
all(st_is_valid(cells_sf))
# there're some cases of self-intersection
all(st_is_valid(nuc_sf))
# In that case, I'll do buffer 0 and then remove the holes
# Since nuclei in one z plane shouldn't have holes
ind_invalid <- !st_is_valid(nuc_sf)
nuc_sf[ind_invalid,] <- nngeo::st_remove_holes(st_buffer(nuc_sf[ind_invalid,], 0))
colData(sce) <- cbind(colData(sce), cell_info[,-1])
spe <- toSpatialExperiment(sce, spatialCoordsNames = c("x_centroid", "y_centroid"))
sfe <- toSpatialFeatureExperiment(spe)
cellSeg(sfe, withDimnames = FALSE) <- cells_sf
nucSeg(sfe, withDimnames = FALSE) <- nuc_sf
sfe <- addQC(sfe)
colData(sfe)$total_counts <- NULL
saveRDS(sfe, file = "xenium1.rds")

# Rep 2
sce <- read10xCounts("Xenium_FFPE_Human_Breast_Cancer_Rep2_cell_feature_matrix.h5")
counts(sce) <- as(realize(counts(sce)), "dgCMatrix")
cell_info <- read_parquet("Xenium_FFPE_Human_Breast_Cancer_Rep2_cells.parquet")
cell_poly <- read_parquet("Xenium_FFPE_Human_Breast_Cancer_Rep2_cell_boundaries.parquet")
nuc_poly <- read_parquet("Xenium_FFPE_Human_Breast_Cancer_Rep2_nucleus_boundaries.parquet")
names(cell_poly)[1] <- "ID"
names(nuc_poly)[1] <- "ID"
cells_sf <- df2sf(cell_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
nuc_sf <- df2sf(nuc_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
all(st_is_valid(cells_sf))
# there're some cases of self-intersection
all(st_is_valid(nuc_sf))

ind_invalid <- !st_is_valid(nuc_sf)
nuc_sf[ind_invalid,] <- nngeo::st_remove_holes(st_buffer(nuc_sf[ind_invalid,], 0))
colData(sce) <- cbind(colData(sce), cell_info[,-1])
spe <- toSpatialExperiment(sce, spatialCoordsNames = c("x_centroid", "y_centroid"))
sfe <- toSpatialFeatureExperiment(spe)
cellSeg(sfe, withDimnames = FALSE) <- cells_sf
nucSeg(sfe, withDimnames = FALSE) <- nuc_sf
sfe <- addQC(sfe)
saveRDS(sfe, file = "xenium2.rds")

# CosMX lung sample=====================
# Data downloaded from https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/
# Lung5_Rep1 is used here. 
# Cell polygons were downloaded from the link on Seurat's vignette: https://www.dropbox.com/s/hl3peavrx92bluy/Lung5_Rep1-polygons.csv?dl=0
mat <- vroom("Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_exprMat_file.csv")
metadata <- vroom("Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_metadata_file.csv")
cell_poly <- vroom("Lung5_Rep1-polygons.csv")
cell_poly <- cell_poly %>% 
    unite("ID", fov:cellID)
mat <- mat %>% unite("ID", fov:cell_ID)
metadata <- metadata %>% unite("ID", fov:cell_ID)
mat <- mat[match(metadata$ID, mat$ID),]
cell_poly <- cell_poly %>% 
    filter(ID %in% metadata$ID)
cell_sf <- df2sf(cell_poly, spatialCoordsNames = c("x_global_px", "y_global_px"),
                 geometryType = "POLYGON")
all(st_is_valid(cell_sf))
# 2 cells have fewer than 3 vertices. Remove those
mat <- mat[mat$ID %in% cell_sf$ID,]
metadata <- metadata[metadata$ID %in% cell_sf$ID,]
m <- t(as.matrix(mat[,-1]))
colnames(m) <- mat$ID
m <- as(m, "dgCMatrix")
cell_sf <- cell_sf[match(mat$ID, cell_sf$ID),]
sfe <- SpatialFeatureExperiment(assays = list(counts = m),
                                colData = metadata[, c(2:3, 6:19)],
                                spatialCoordsNames = c("CenterX_global_px",
                                                       "CenterY_global_px"))
cellSeg(sfe) <- cell_sf
sfe <- addQC(sfe)
saveRDS(sfe, "cosmx1.rds")

# Vizgen MERFISH liver data ======================
# Downloaded from https://console.cloud.google.com/storage/browser/vz-liver-showcase/Liver1Slice1;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false&pli=1
library(rhdf5)
# The cell boundaries are in hdf5 files
# Each one is for an FOV, and there are 7 z-planes.
# However, from the few I checked, all the 7 are the same
# So I wonder why they even bothered to save 7 copies. I'll just use the first one.
# I think I'll save each z-plane as a separate sf data frame, and keep the sf
# data frames separate from the main SFE object, which has the centroids.
# I'm not going through df2sf here.
h52poly_fov <- function(fn, i) {
    l <- rhdf5::h5dump(fn)[[1]]
    cell_ids <- names(l)
    geometries <- lapply(l, function(m) sf::st_polygon(list(t(m[["zIndex_0"]]$p_0$coordinates[,,1]))))
    df <- data.frame(geometry = sf::st_sfc(geometries),
                     ID = cell_ids,
                     fov = i)
    sf::st_sf(df)
}
fns <- list.files("cell_boundaries", "*.hdf5", full.names = TRUE)
polys <- bpmapply(h52poly_fov, fn = fns, i = seq_along(fns), SIMPLIFY = FALSE, 
                  BPPARAM = SnowParam(20, progressbar = TRUE))
polys <- do.call(bind_rows, polys)

mat <- vroom("Liver1Slice1_cell_by_gene.csv", col_types = cols(...1 = "c"))
polys <- polys[match(mat$...1, polys$ID),]
metadata <- vroom("Liver1Slice1_cell_metadata.csv", col_types = cols(...1 = "c"))
metadata <- metadata[match(mat$...1, metadata$...1),]

m <- as.matrix(mat[,-1])
m <- as(m, "dgCMatrix")
rownames(m) <- mat$...1
m <- t(m)

polys$fov <- NULL
# Takes a while to make the POINT geometry for the centroids, not too bad
sfe <- SpatialFeatureExperiment(assays = list(counts = m),
                                colData = metadata[,-1], 
                                spatialCoordsNames = c("center_x", "center_y"))
rownames(polys) <- polys$ID
polys$ID <- NULL
cellSeg(sfe) <- polys
sfe <- addQC(sfe)
saveRDS(sfe, "merfish_liver1.rds")


# Mouse gastrulation embryo seqFISH data====================
# Data downloaded from hhttps://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/

library(magrittr)

counts <- readRDS('counts.Rds')
meta <- readRDS('metadata.Rds')
seg <- readRDS('segmentation_vertices.Rds')

# create sf dataframe for all cells
polygos <- seg %>% 
  select(ID = uniqueID, contains("segmentation"))

names(polygos)[-1] <- c('x', 'y')

# 6 low qual cells are removed because they lack segmentation data
polygos <- df2sf(polygos, geometryType = "POLYGON")


# make valid polygons for objs with self-intersections
all(st_is_valid(polygos))

inv_inds <- which(!st_is_valid(polygos))
fixed_poly <- st_make_valid(polygos[inv_inds,], geos_method = 'valid_structure')

polygos[inv_inds,] <- fixed_poly

# add back cells lacking segmentation data
# add empty geom for dropped cells
mask <- !(meta$uniqueID %in% polygos$ID)
if (sum(mask)) {
  df <- data.frame(ID = meta$uniqueID[mask])
  nr <- nrow(df)
  
  expr_ <- glue::glue("st_sfc({str_flatten(rep('st_polygon()', nr), ',')})")
  st_geometry(df) <- eval(parse_expr(expr_))
  
  polygos <- rbind(polygos, df)
}

polygos %<>%  mutate(embryo = str_extract(ID, "(embryo[0-9]+)"))

# create sfe object for each embryo
bio_reps <- unique(meta$embryo)

for (em in bio_reps){
  mask_em <- meta$embryo %in% em
  
  # subset meta, counts, segmentation for embryo
  meta_em <- meta %>% filter(embryo %in% em) %>% select(-contains('segmentation'))
  counts_em <- counts[, mask_em]
  
  poly_em <- polygos %>% filter(embryo %in% em)  %>% select(-embryo)
  
  sfe <- SpatialFeatureExperiment(
    assays = list(counts = counts_em),
    colData = meta_em,
    colGeometries = list(seg_coords = poly_em)
  )
  
  fn = glue::glue("seqfish_em{str_extract(em, '[0-9]')}.Rds")

  saveRDS(sfe, file = fn)
}


