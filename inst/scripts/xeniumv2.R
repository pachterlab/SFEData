# Create Xenium XOA v2 test data
options( java.parameters = "-Xmx4g" )
devtools::load_all("~/SpatialFeatureExperiment")
devtools::load_all("~/Voyager")
library(sf)
library(RBioFormats)
library(DropletUtils)
library(xml2)
library(jsonlite)
library(arrow)
library(data.table)
library(sfheaders)
library(tidyverse)
sfe <- readXenium("/Volumes/Archives/xenium2_pancreas/", add_molecules = TRUE)
bbox(sfe, include_images = TRUE)
plotSpatialFeature(sfe, "total_counts", colGeometryName = "cellSeg", dark = TRUE,
                   show_axes = TRUE)
bbox_use <- c(xmin = 1000, ymin = -1000, xmax = 2000, ymax = 0)
# From profiling, most of the time was spent cropping rowGeometries
sfe_sub <- crop(sfe, bbox_use, keep_whole = "col")
plotSpatialFeature(sfe_sub, "total_counts", colGeometryName = "cellSeg", dark = TRUE,
                   show_axes = TRUE)
bbox(sfe_sub, include_images = TRUE)
dim(sfe_sub)

# Write to disk in XOA output format---------------
dir.create("xenium2")
# Write the gene count matrix h5 file
write10xCounts("xenium2/cell_feature_matrix.h5", counts(sfe_sub), gene.id = rowData(sfe_sub)$ID,
               gene.symbol = rowData(sfe_sub)$Symbol, gene.type = rowData(sfe_sub)$Type,
               version = "3", type = "HDF5", overwrite = TRUE)

# Write image, I'll get the pixel ranges and then use them to crop with bfconvert
dir.create("xenium2/morphology_focus")
# Stop it half way just to get the pixel ranges
pss <- .get_pixel_size(imgSource(getImg(sfe_sub)), resolution = 3L)
debug(.toEBImage)
img <- toEBImage(getImg(sfe_sub), resolution = 3L)
# so should be bfconvert -compression LZW -crop 1163,3,1208,1187 -series 2 -overwrite \
# /Volumes/Archives/xenium2_pancreas/morphology_focus/morphology_focus_0000.ome.tif \
# xenium2/morphology_focus/morphology_focus_000%c.ome.tif
# Won't make pyramid since there's the buffer too small error and I have no idea
# what to do

# Still need to add pixel size
xml_meta <- read.omexml("xenium2/morphology_focus/morphology_focus_0000.ome.tif") |> read_xml()
xml_structure(xml_meta)
ps_child <- xml_child(xml_meta, 2) |> xml_child()
# Actually this modifies the original
xml_set_attr(ps_child, "PhysicalSizeX", value = pss[1])
xml_set_attr(ps_child, "PhysicalSizeY", value = pss[2])
xml_set_attr(ps_child, "PhysicalSizeXUnit", value = "µm")
xml_set_attr(ps_child, "PhysicalSizeYUnit", value = "µm")
write_xml(xml_meta, "meta_new.xml")
# tiffcomment -set meta_new.xml xenium2/morphology_focus/morphology_focus_0000.ome.tif
# Do the same for the other files

# Write experiment.xenium
experiment <- fromJSON(file = "/Volumes/Archives/xenium2_pancreas/experiment.xenium",
                       simplify = TRUE)
experiment$run_start_time <- format(Sys.time(), format = "%Y-%m-%dT%H:%M:%S%Z")
experiment$num_cells <- ncol(sfe_sub)
experiment$pixel_size <- pss[1]
write_json(experiment, "xenium2/experiment.xenium", auto_unbox = TRUE, pretty = TRUE)

# Write cell vertices parquet file (not GeoParquet)
cells <- read_parquet("/Volumes/Archives/xenium2_pancreas/cell_boundaries.parquet")
cells <- cells[cells$cell_id %in% colnames(sfe_sub),]
# Translate the extent to pretend that it's full extent,
bbox_use <- bbox(sfe_sub, include_images = TRUE)
v <- c(bbox_use["xmin"], -bbox_use["ymax"])
cells <- cells |>
    mutate(vertex_x = vertex_x - v[1],
           vertex_y = vertex_y - v[2])

# Make sure there are polygons that are invalid
cells_sf <- df2sf(cells, spatialCoordsNames = c("vertex_x", "vertex_y"),
                  geometryType = "POLYGON", id_col = "cell_id")
st_bbox(cells_sf)
all(st_is_valid(cells_sf)) # There are invalid ones
write_parquet(cells, "xenium2/cell_boundaries.parquet")
fwrite(cells, "xenium2/cell_boundaries.csv.gz")

nuclei <- read_parquet("/Volumes/Archives/xenium2_pancreas/nucleus_boundaries.parquet")
nuclei <- nuclei[nuclei$cell_id %in% colnames(sfe_sub),]
nuclei <- nuclei |>
    mutate(vertex_x = vertex_x - v[1],
           vertex_y = vertex_y - v[2])
write_parquet(nuclei, "xenium2/nucleus_boundaries.parquet")
fwrite(nuclei, "xenium2/nucleus_boundaries.csv.gz")

# Write transcript spot parquet file (not GeoParquet)
tx_orig <- read_parquet("/Volumes/Archives/xenium2_pancreas/transcripts.parquet")
tx_orig <- tx_orig |>
    filter(between(x_location, bbox_use["xmin"], bbox_use["xmax"]),
           between(y_location, -bbox_use["ymax"], -bbox_use["ymin"])) |>
    mutate(x_location = x_location - v[1],
           y_location = y_location - v[2])
write_parquet(tx_orig, "xenium2/transcripts.parquet")
fwrite(tx_orig, "xenium2/transcripts.csv.gz")

# Write cell metadata
cell_metadata <- read_parquet("/Volumes/Archives/xenium2_pancreas/cells.parquet")
cell_metadata <- cell_metadata |>
    filter(cell_id %in% colnames(sfe_sub)) |>
    mutate(x_centroid = x_centroid - v[1],
           y_centroid = y_centroid - v[2])
write_parquet(cell_metadata, "xenium2/cells.parquet")
fwrite(cell_metadata, "xenium2/cells.csv.gz")

# See if it worked
sfe2 <- readXenium("xenium2", add_molecules = TRUE)
getImg(sfe2)
plotSpatialFeature(sfe2, "total_counts", colGeometryName = "cellSeg", dark = TRUE,
                   show_axes = TRUE)
ext(getImg(sfe2))
st_bbox(cellSeg(sfe2))
st_bbox(nucSeg(sfe2))
st_bbox(txSpots(sfe2))

# Seems to work, but I'll further check after I fix the Voyager plotting function
# for images that are not SpatRasterImage

# Check tar.gz file size
# tar -czvf xenium2.tar.gz xenium2
# 35.7 MB, not too bad
