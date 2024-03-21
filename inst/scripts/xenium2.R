# Make xenium test data
# Test readXenium w/o the image part first
# Then do the subsetting. The images are cropped separately
# TBH I can very well use the "tiny" subset from the 10X website, but it's nearly 300 MB
# even after compression. My concern is download time.
# The purpose of this even tinier subset is fast download for testing and examples

library(SingleCellExperiment)
library(sf)
library(Voyager)
library(ggplot2)
options( java.parameters = "-Xmx4g" )
library(RBioFormats)
library(EBImage)
library(terra)
library(xml2)
devtools::load_all() # version 1.5.2 of SFE

# Downloaded from https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
# https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard
sfe <- readXenium("inst/extdata/xenium_toy", sample_id = "test_xenium",
                  add_molecules = TRUE)
bbox(sfe)
plotSpatialFeature(sfe, "total_counts", colGeometryName = "cellSeg") + theme_bw()
img <- toEBImage(getImg(sfe), resolution = 4L)
display(imgRaster(img) |> normalize(), method = "raster")
img@image <- normalize(img@image, ft = c(0,255))
spi <- toSpatRasterImage(img, save_geotiff = FALSE)
imgData(sfe)$data[[1]] <- spi
bbox_use <- c(xmin = 1933, xmax = 2433, ymin = -2190, ymax = -1690)
plotSpatialFeature(sfe, "total_counts", bbox = bbox_use,
                   colGeometryName = "nucSeg",
                   image_id = "morphology_focus.ome", dark = TRUE,
                   aes_use = "color", fill = NA) + theme_dark()
# Bug: With lower resolution images, the cells are shifted to the bottom right
# compared to nuclei boundaries. Could be because the lower resolution images
# actually have smaller extent compared to the full resolution one.
# Possibly this is because when 10X made this subset, they only included pixels
# that are entirely contained within the bounding box. Nothing to do with cropping

# The challenge: crop the OME-TIFF file at different resolutions.
# I can edit the XML metadata in R, write it out, and replace the original
# as in https://bio-formats.readthedocs.io/en/latest/users/comlinetools/edit.html
# How am I going to do it?
# 1. Select an arbitrary small bbox I want, crop the SFE object
# 2. Read each resolution of the image with the toEBImage function, can read the
# highest resolution first, get actual extent. Convert to SpatRasterImage
# 3. Read the subsequent resolutions with the extent in toEBImage.
# I need to keep the full extent of all the resolutions the same. I can create
# an empty rast with just the extent from the full resolution and the same number
# of pixels in each dimension as the lower resolution. Then I'll resample the
# full resolution image with this empty rast to make lower resolution images that
# would actually correctly align, as if it were full extent.

# Also need to flip the geometries back before writing
sfe <- readXenium("inst/extdata/xenium_toy", sample_id = "test_xenium",
                  add_molecules = TRUE)
bbox_use <- c(xmin = 1000, xmax = 1500, ymin = -2000, ymax = -1500)
sfe_sub <- crop(sfe, bbox_use, colGeometryName = "cellSeg", keep_whole = "col",
                cover = TRUE)

# Just to see what AnnotatedImage for multiple series should look like
foo <- read.image("inst/extdata/xenium_toy/morphology_focus.ome.tif", resolution = c(4,6))
# AnnotatedImageList whose elements are the different resolutions
# Each element is of class AnnotatedImage. Use constructor.
# coreMetadata: list of length 18. I only need to change sizeX and sizeY
metas <- metadata(foo) |> coreMetadata()
rm(foo)
gc()
xml_meta <- read.omexml("img45.ome.tiff")

# Maybe I can try simply using lower resolution images rather than cropping
xml_full <- read.omexml("inst/extdata/xenium_toy/morphology_focus.ome.tif")
xml_full <- read_xml(xml_full)
xml_attr(xml_child(xml_full, 3) |> xml_child(2), "PhysicalSizeX")
# I think I might not need to crop after all. Just keep low resolutions, reset PhysicalSize
# Also downsample the transcript spots
img1 <- getImg(sfe_sub, image_id = "morphology_focus.ome")
ebi1 <- toEBImage(img1, resolution = 1L)
display(imgRaster(ebi1) |> normalize(), method = "raster")
ebi1@image <- normalize(ebi1@image, ft = c(0,255))
spi <- toSpatRasterImage(ebi1, save_geotiff = FALSE)
imgData(sfe_sub)$data[[1]] <- spi
plotSpatialFeature(sfe_sub, "total_counts", colGeometryName = "nucSeg",
                   image_id = "morphology_focus.ome", aes_use = "color", fill = NA)

ebis <- lapply(2:6, toEBImage, x = img1)

img_f <- getImg(sfe)
bbox_use2 <- c(xmin = 1000, xmax = 1500, ymin = -1500, ymax = -1000)
img_c <- cropImg(img_f, bbox_use2)
ebi_c <- toEBImage(img_c, resolution = 4L)
display(imgRaster(ebi_c) |> normalize(), method = "raster")
