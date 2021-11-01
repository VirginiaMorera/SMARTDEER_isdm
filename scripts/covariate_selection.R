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
levelplot(r, maxpixels = 10e9)

subs <- subset(covars, c("landCover", "tree_cover_density", "elevation", "slope", "human_footprint_index"))

plot(subs)

corsub <-layerStats(subs,'pearson', na.rm = T)

corrSub <- ggcorrplot(corsub$`pearson correlation coefficient`, method = "circle", type = "upper", 
                      legend.title = "Pear Cor 1km", show.diag = T, lab = T)

writeRaster(subs, filename = "large_env_data/covar_subset.grd", format = "raster")
