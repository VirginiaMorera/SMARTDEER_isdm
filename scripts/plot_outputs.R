
ireland <- st_read("C:/Users/Virginia/Dropbox/Virginia_post_UB/01_Smartdeer/Counties.shp")
mask <- readRDS("data/barrier.RDS")
sp_pred <- readRDS("cluster_outputs/sppredresp_range200_sd0.1.RDS")
spde <- readRDS("cluster_outputs/spde_range100_sd1.RDS")
ireland_proj <- st_transform(ireland, sp_pred@proj4string)
ireland_sp <- as_Spatial(ireland_proj)
mask_sf <- st_as_sf(mask, crs = sp_pred@proj4string)

plot(sp_pred['mean'])
plot(ireland_sp, add = T)

plot(sp_pred['sd'], )
plot(ireland_sp, add = T)



plot(spde['mean'])
plot(ireland_sp, add = T)

pred_mean <- raster(sp_pred['mean'])
pred_sd <- raster(sp_pred['sd'])

pred <- stack(raster(sp_pred['mean']), raster(sp_pred['sd']), raster(sp_pred['median']), 
              raster(sp_pred['q0.025']), raster(sp_pred['q0.975']))


pred$IQR <- pred$q0.975 - pred$q0.025

(plt <- levelplot(pred$median, col.regions = viridis(16), 
                  xlab = "Longitude", ylab = "Latitude", 
                  # xlim = c(-22, -7), ylim = c(15, 33), 
                  margin = TRUE, zscaleLog = FALSE) +
    latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray")) )


nl <- nlayers(pred)
m <- matrix(1:nl, nrow=2)
for (i in 1:nl){
  p <- levelplot(pred, layers=i, col.regions = viridis(16), 
                 xlab = "Longitude", ylab = "Latitude", 
                 # xlim = c(-22, -7), ylim = c(15, 33), 
                 margin = FALSE, zscaleLog = FALSE, main = names(pred)[i]) + 
    latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
  print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
}
  