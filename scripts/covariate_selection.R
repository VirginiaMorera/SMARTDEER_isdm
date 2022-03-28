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


# merge dist to forest and dist to forest edge

dist_to_fEdge <- (-1)*covars$dist_to_fEdge
dist_to_fEdge[is.na(dist_to_fEdge)] <- 0

dist_to_forest <- covars$dist_to_forest
dist_to_forest[is.na(dist_to_forest)] <- 0

forest_distances <- sum(dist_to_forest, dist_to_fEdge)

forest_distances <- mask(forest_distances, ireland_sp)
plot(forest_distances)
covars$forest_distances <- forest_distances
plot(covars)

subs <- subset(covars, c("tree_cover_density", "elevation", "slope", "human_footprint_index", "forest_distances", 
                         "small_woody_features"))

plot(subs)

corsub <-layerStats(subs,'pearson', na.rm = T)

colnames(corsub$`pearson correlation coefficient`) <- c("Tree cover density", "Elevation", "Slope", "Human footprint index", "Distance to forest edge", 
                                                           "Density of small woody features")
rownames(corsub$`pearson correlation coefficient`) <- c("Tree cover density", "Elevation", "Slope", "Human footprint index", "Distance to forest edge", 
                                                                "Density of small woody features")
                  
corrSub <- ggcorrplot(corsub$`pearson correlation coefficient`, method = "circle", type = "upper", 
                      legend.title = "Pearson Correlation", show.diag = F, lab = T)

Cairo::CairoPDF(file = "Fig5Sup.pdf", width = 8, height = 6)
corrSub
dev.off()


# project each layer separately to use different method for landCover
env_data_ITM <- stack() 
env_data_ITM$tree_cover_density <- projectRaster(subs$tree_cover_density, crs = CRS("+init=epsg:2157"), method = "bilinear")
env_data_ITM$elevation <- projectRaster(subs$elevation, crs = CRS("+init=epsg:2157"), method = "bilinear")
env_data_ITM$slope <- projectRaster(subs$slope, crs = CRS("+init=epsg:2157"), method = "bilinear")
env_data_ITM$human_footprint_index <- projectRaster(subs$human_footprint_index, crs = CRS("+init=epsg:2157"), method = "bilinear")
env_data_ITM$landCover <- projectRaster(subs$landCover, crs = CRS("+init=epsg:2157"), method = "ngb")
env_data_ITM$forest_distances <- projectRaster(subs$forest_distances, crs = CRS("+init=epsg:2157"), method = "bilinear")
env_data_ITM$small_woody_features <- projectRaster(subs$small_woody_features, crs = CRS("+init=epsg:2157"), method = "bilinear")

writeRaster(env_data_ITM, filename = "large_env_data/covar_subset_ITM.grd", format = "raster", overwrite = T)


levelplot(r, col.regions = viridis(12), 
          xlab = "Longitude", ylab = "Latitude", 
          # xlim = c(-22, -7), ylim = c(15, 33), 
          margin = FALSE, zscaleLog = FALSE, main = "Land use") + 
  latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "black"))


levelplot(covars$forest_distances, col.regions = viridis(100), 
          xlab = "Longitude", ylab = "Latitude", 
          # xlim = c(-22, -7), ylim = c(15, 33), 
          margin = FALSE, zscaleLog = FALSE, main = "Tree cover density") + 
  latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "white"))


covar_stack <- stack("large_env_data/covar_subset_ITM.gri")
covar_stackKM <- projectRaster(covar_stack, crs = projKM)

writeRaster(covar_stackKM, filename = "large_env_data/covar_subset_KM.grd", format = "raster", overwrite = T)

covar_stackKM <- stack("large_env_data/covar_subset_KM.grd")



covar_stackKM <- projectRaster(covar_stackKM, crs = CRS("EPSG:4326"))
ireland_sp <- spTransform(ireland_sp, CRSobj = CRS("EPSG:4326"))

covar_stackKM <- mask(covar_stackKM, ireland_sp)

levelplot(covar_stackKM$forest_distances, col.regions = viridis(100), 
          xlab = "Longitude", ylab = "Latitude", 
          xlim = c(-11, -5), ylim = c(51.4, 55.5),
          margin = FALSE, zscaleLog = FALSE, main = "Distance to forest") + 
  latticeExtra::layer(sp.polygons(ireland_sp, col = "white", alpha = 0.5))
