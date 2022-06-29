library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(sf)
library(terra)
library(EBImage)
library(rmapshaper)
library(lwgeom)
library(tidyverse)
library(jsonlite)
library(Matrix)
library(Voyager)
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
colData(sfe)$in_tissue <- annotOp(sfe, "spotPoly", "tissueBoundary")

saveRDS(sfe, "data/sfe_vis5a.rds")
