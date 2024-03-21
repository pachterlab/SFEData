options( java.parameters = "-Xmx4g" )
library(RBioFormats)
library(xml2)
library(data.table)
library(arrow)
library(Voyager)
library(terra)
library(EBImage)
library(ggplot2)
library(BiocParallel)
devtools::load_all() # version 1.5.2 of SFE

# Downloaded from https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
# https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard
file <- "inst/extdata/xenium_toy/morphology_mip.ome.tif"
imgs <- read.image(file, resolution = c(4,6), normalize = FALSE)
cm1 <- metadata(imgs) |> coreMetadata(1)
cm1$resolutionLevel <- 1L
cm2 <- metadata(imgs) |> coreMetadata(2)
cm2$resolutionLevel <- 2L
coreMetadata(imgs[[1]]) <- cm1
coreMetadata(imgs[[2]]) <- cm2
metadata(imgs) |> coreMetadata()
write.image(imgs, "xenium_lr/morphology_mip.ome.tif", force = TRUE,
            littleEndian = TRUE, pixelType = "uint16")

m <- RBioFormats::read.metadata(file)
cm <- coreMetadata(m)
# Get the downsizing factor
sizeXs <- vapply(cm, function(x) x$sizeX, FUN.VALUE = numeric(1))
sizeYs <- vapply(cm, function(x) x$sizeY, FUN.VALUE = numeric(1))

fcts_x <- sizeXs[1:5] / sizeXs[2:6]
fcts_y <- sizeYs[1:5] / sizeYs[2:6]
# It's usually a little larger than 2 which is the intended factor, but easy to round

sfs <- .get_fullres_scale_factor(file)
sfx <- sfs[1]; sfy <- sfs[2] # From microns to full res pixel
size_full <- .get_fullres_size(file)
sizeX_full <- size_full[1]; sizeY_full <- size_full[2]

meta <- coreMetadata(m, series = 4L)
fct_x <- sizeX_full/meta$sizeX
fct_y <- sizeY_full/meta$sizeY
fct_round <- round(fct_x) # Should be the same for x and y
fctx2 <- fct_x/fct_round
fcty2 <- fct_y/fct_round

sfx2 <- meta$sizeX*fctx2/sizeX_full # Multiply to get from full res pixel to low res pixel
sfy2 <- meta$sizeY*fcty2/sizeY_full
psx4 <- 1/(sfx*sfx2)
psy4 <- 1/(sfy*sfy2)

# Doesn't work==================
# Create new xml
xml_orig <- read.omexml(file) |> read_xml()
xml_meta <- read.omexml("xenium_lr/morphology_mip.ome.tif") |> read_xml()
ps_child <- xml_child(xml_meta) |> xml_child()
# Actually this modifies the original
xml_set_attr(ps_child, "PhysicalSizeX", value = psx4)
xml_set_attr(ps_child, "PhysicalSizeY", value = psy4)
#xml_set_attr(ps_child, "SizeX", value = dim(imgs[[1]])[1])
#xml_set_attr(ps_child, "SizeY", value = dim(imgs[[1]])[2])

sa_child <- xml_child(xml_orig, 4) |> xml_child() |> xml_child()
xml_remove(xml_children(sa_child)[1:4])
xml_set_attr(xml_child(sa_child), "K", "1")

sa_top <- xml_child(xml_orig, 4)
xml_add_child(xml_meta, sa_top)

xml2::write_xml(xml_meta, "meta4.xml")
# https://bio-formats.readthedocs.io/en/latest/users/comlinetools/edit.html
# tiffcomment -set meta4.xml img46.ome.tiff
# See if it worked
xml_new <- read.omexml("xenium_lr/morphology_mip.ome.tif") |> read_xml()
xml_structure(xml_new)
bar <- read.image("xenium_lr/morphology_mip.ome.tif")
m <- metadata(bar)
cm2 <- coreMetadata(m, 2)
# Doesn't work. It turns the two resolutions into two series instead no matter
# how I hack the xml metadata.================= End of what doesn't work

# Trying bfconvert
write.image(imgs[[1]], "morphology_mip_res4.tiff", force = TRUE)
xml_meta <- read.omexml("morphology_mip_res4.tiff") |> read_xml()
xml_structure(xml_meta)
ps_child <- xml_child(xml_meta) |> xml_child(2)
# Actually this modifies the original
xml_set_attr(ps_child, "PhysicalSizeX", value = psx4)
xml_set_attr(ps_child, "PhysicalSizeY", value = psy4)
write_xml(xml_meta, "meta4.xml")
# tiffcomment -set meta4.xml morphology_mip_res4.tiff
# Download bfconvert: https://downloads.openmicroscopy.org/bio-formats/6.7.0/artifacts/bftools.zip
# bfconvert -noflat -pyramid-resolutions 2 -pyramid-scale 3 morphology_mip_res4.tiff morphology_mip.ome.tif
foo <- read.image("morphology_mip.ome.tiff")
metadata(foo) |> coreMetadata() # Worked!
# Do the same to the focus image
write.image(read.image("inst/extdata/xenium_toy/morphology_focus.ome.tif",
                       resolution = 4, normalize = FALSE),
            file = "morphology_focus_res4.tiff", force = TRUE)

# Downsample the transcript spots
#spots_parq <- read_parquet("xenium_lr/transcripts.parquet")
spots_dt <- fread("xenium_lr/transcripts.csv.gz")
set.seed(29)
n <- seq_len(nrow(spots_dt))
inds <- sample(n, size = round(nrow(spots_dt)/10), replace = FALSE)

#sub_parq <- spots_parq[inds,]
sub_dt <- spots_dt[inds,]

write_parquet(sub_dt, "xenium_lr/transcripts.parquet")
fwrite(sub_dt, "xenium_lr/transcripts_sub.csv.gz")

# Convert the arrow binary in cell and nuclei segmentation
# The actual conversion is super slow, taking too long for testing purposes
# I'll create another parquet file from the csv without the binaries
# Then only keep a small subset of cells in the parquet files that do have binaries
# for testing purposes only
# parquet file is faster to read
cell_parq <- read_parquet("xenium_lr/cell_boundaries.parquet")
cell <- fread("xenium_lr/cell_boundaries.csv.gz")
#cells_nobinary <- .rawToChar_df(cells, BPPARAM = MulticoreParam(2, progressbar = TRUE))
nuc_parq <- read_parquet("xenium_lr/nucleus_boundaries.parquet")
nuc <- fread("xenium_lr/nucleus_boundaries.csv.gz")
#nuc_nobinary <- .rawToChar_df(nuc, BPPARAM = MulticoreParam(2, progressbar = TRUE))

write_parquet(cell, "xenium_lr/cell_boundaries_nobinary.parquet")
write_parquet(nuc, "xenium_lr/nucleus_boundaries_nobinary.parquet")
cell_ids <- unique(cell$cell_id)
cells_sub <- sample(cell_ids, round(length(cell_ids)/100), replace = FALSE)
inds <- which(cell$cell_id %in% cells_sub)
cell_parq_sub <- cell_parq[inds,]
inds <- which(nuc$cell_id %in% cells_sub)
nuc_parq_sub <- nuc_parq[inds,]

write_parquet(cell_parq_sub, "xenium_lr/cell_boundaries_binary.parquet")
write_parquet(nuc_parq_sub, "xenium_lr/nucleus_boundaries_binary.parquet")

# cell2 <- .rawToChar_df(cell_parq_sub, SerialParam(progressbar = TRUE))
# Not taking too long
# Will rename the appropriate files when testing
file.remove("xenium_lr/nucleus_boundaries.parquet")
file.remove("xenium_lr/cell_boundaries.parquet")

# See if it works
sfe <- readXenium("xenium_lr", add_molecules = TRUE)
bbox(sfe)
ext(getImg(sfe))
bfi <- getImg(sfe)
img <- toEBImage(bfi, resolution = 2L)
display(imgRaster(img) |> normalize(), method = "raster")
img@image <- normalize(img@image, ft = c(0,255))
spi <- toSpatRasterImage(img, save_geotiff = FALSE)

imgData(sfe)$data[[1]] <- spi
bbox_use <- c(xmin = 0, xmax = 500, ymin = -500, ymax = 0)
bbox_lr <- c(xmin = 1932, xmax = 2432, ymin = -2189, ymax = -1689)
plotSpatialFeature(sfe, "total_counts", bbox = bbox_lr,
                   colGeometryName = "nucSeg",
                   image_id = "morphology_focus.ome", dark = TRUE,
                   aes_use = "color", fill = NA) + theme_dark()

# See the size of the tar.gz and if it's acceptable
tar("xenium.tar.gz", "xenium_lr", compression = "gzip", compression_level = 9)
# 24 MB, which I consider acceptable, comparable to the Visium datasets I
# downloaded from 10X website for testing purposes
