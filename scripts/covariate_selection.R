covars <- stack("large_env_data/all_covars1km.grd")
plot(covars)

ireland <- st_read("data/ireland_ITM.shp")
ireland_proj <- st_transform(ireland, covars@crs)
ireland_sp <- as_Spatial(ireland_proj)


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


levelplot(r, col.regions = viridis(12), 
          xlab = "Longitude", ylab = "Latitude", 
          # xlim = c(-22, -7), ylim = c(15, 33), 
          margin = FALSE, zscaleLog = FALSE, main = "Land use") + 
  latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "black"))


levelplot(covars$tree_cover_density, col.regions = viridis(100), 
          xlab = "Longitude", ylab = "Latitude", 
          # xlim = c(-22, -7), ylim = c(15, 33), 
          margin = FALSE, zscaleLog = FALSE, main = "Tree cover density") + 
  latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "white"))
