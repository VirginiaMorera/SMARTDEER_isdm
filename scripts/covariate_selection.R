library(rasterVis)

covars <- stack("large_env_data/all_covars1km.grd")
plot(covars)

x <- covars$landCover
plot(x)
r <- ratify(x)
rat <- levels(r)[[1]]
rat$landcover <- c("Builtup", "Saltwater_related", "No_vegetation", "Freshwater_related", 
                   "Other_vegetation", "Agricultural", "Pasture", "Grassland", "Transitional", 
                   "Coniferous_forest", "Mixed_forest", "Broadleaf_forest")
levels(r) <- rat
levelplot(r)

subs <- subset(covars, c("landCover", "tree_cover_density", "elevation", "slope", "human_footprint_index"))

plot(subs)

corsub <-layerStats(subs,'pearson', na.rm = T)

corrSub <- ggcorrplot(corsub$`pearson correlation coefficient`, method = "circle", type = "upper", 
                      legend.title = "Pear Cor 1km", show.diag = T, lab = T)

# project each layer separately to use different method for landCover
env_data_ITM <- stack() 
env_data_ITM$tree_cover_density <- projectRaster(subs$tree_cover_density, crs = CRS("+init=epsg:2157"), method = "bilinear")
env_data_ITM$elevation <- projectRaster(subs$elevation, crs = CRS("+init=epsg:2157"), method = "bilinear")
env_data_ITM$slope <- projectRaster(subs$slope, crs = CRS("+init=epsg:2157"), method = "bilinear")
env_data_ITM$human_footprint_index <- projectRaster(subs$human_footprint_index, crs = CRS("+init=epsg:2157"), method = "bilinear")
env_data_ITM$landCover <- projectRaster(subs$landCover, crs = CRS("+init=epsg:2157"), method = "ngb")

writeRaster(env_data_ITM, filename = "large_env_data/covar_subset_ITM.grd", format = "raster")
