# load data
ireland <- st_read("data/ireland_ITM.shp")

spde <- readRDS("cluster_outputs/spde_range100_sd2.RDS")
ireland_proj <- st_transform(ireland, sp_pred@proj4string)
ireland_sp <- as_Spatial(ireland_proj)

sp_pred <- readRDS("cluster_outputs/sppred_range100_sd2.RDS")


pred <- stack(raster(sp_pred['mean']), raster(sp_pred['sd']), raster(sp_pred['median']), 
              raster(sp_pred['q0.025']), raster(sp_pred['q0.975']))

spderas <- stack(raster(spde['mean']), raster(spde['sd']), raster(spde['median']), 
              raster(spde['q0.025']), raster(spde['q0.975']))

spderas$IQR <- spderas$q0.975 - spderas$q0.025

# (plt <- levelplot(pred$median, col.regions = viridis(16), 
#                   xlab = "Longitude", ylab = "Latitude", 
#                   # xlim = c(-22, -7), ylim = c(15, 33), 
#                   margin = TRUE, zscaleLog = FALSE) +
#     latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray")) +
#     latticeExtra::layer(sp.polygons(nha_sp, lwd = 0.1, col = "pink")))

spderas2 <- subset(spderas, c(1,2,3,6))

nl <- nlayers(spderas2)
m <- matrix(1:nl, nrow=2)
for (i in 1:nl){
  p <- levelplot(spderas2, layers=i, col.regions = viridis(16), 
                 xlab = "Longitude", ylab = "Latitude", 
                 # xlim = c(-22, -7), ylim = c(15, 33), 
                 margin = FALSE, zscaleLog = FALSE, main = names(spderas2)[i]) + 
    latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
  print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
}
  
