# load ireland data
ireland <- st_read("data/ireland_ITM.shp")


## Red deer ##

# load model
mdl_RD <- readRDS("server_outputs/mdl_RD_def.RDS")

fixed.effects <- mdl_RD$summary.fixed
fixed.effects$variable <- row.names(fixed.effects)

fixed.effects <- fixed.effects %>% 
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_sel_intercept"))

names(fixed.effects) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")

# fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
#                             "Human footprint index", "Land cover")

fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
                            "Human footprint index", "Forest distance", "Small Woody Features")


fixedRD <- fixed.effects %>% 
  # filter(variable != "Human footprint index") %>%
  ggplot(aes(x = variable, 
             # y = median,
             y = mean,
             colour = variable)) + 
  # geom_errorbar(aes(ymin = lower, ymax = higher), width =.1) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width =.1) +
  geom_line() +
  geom_point() + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  coord_flip() + 
  labs(title = "Covariate effects Red deer", y = "Median", x = "", colour = "Covariate") +
  NULL


# map outputs
RD_pred_lin <- readRDS("server_outputs/RD_pred_lin_def.RDS")
RD_pred_resp <- readRDS("server_outputs/RD_pred_resp_def.RDS")
RD_spde_PO <- readRDS("server_outputs/RD_spde_PO_def.RDS")
RD_spde_PA <- readRDS("server_outputs/RD_spde_PA_def.RDS")

ireland_proj <- st_transform(ireland, RD_pred_lin@proj4string)
ireland_sp <- as_Spatial(ireland_proj)


pred_lin <- stack(raster(RD_pred_lin['mean']), raster(RD_pred_lin['sd']))
pred_resp <- stack(raster(RD_pred_resp['mean']), raster(RD_pred_resp['sd']))
spde_PO <- stack(raster(RD_spde_PO['mean']), raster(RD_spde_PO['sd']))
spde_PA <- stack(raster(RD_spde_PA['mean']), raster(RD_spde_PA['sd']))


# Plot 

plotlayer <- pred_resp

nl <- nlayers(plotlayer)
m <- matrix(1:nl, ncol = 2)
for (i in 1:nl){
  p <- levelplot(plotlayer, layers=i, col.regions = viridis(16), 
                 xlab = "Longitude", ylab = "Latitude", 
                 zscaleLog = FALSE,
                 # xlim = c(-22, -7), ylim = c(15, 33), 
                 margin = FALSE, main = names(plotlayer)[i]) + 
    latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
  print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
}
  
spdes <- stack(spde_PO$mean, spde_PA$mean)
names(spdes) <- c("PO", "PA")

levelplot(spdes, col.regions = viridis(16), 
          xlab = "Longitude", ylab = "Latitude", 
          zscaleLog = FALSE,
          # xlim = c(-22, -7), ylim = c(15, 33), 
          margin = FALSE, main = names(spdes)) + 
  latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))


## Sika deer ##

# load model
mdl_SD <- readRDS("server_outputs/mdl_SD_def.RDS")

fixed.effects <- mdl_SD$summary.fixed
fixed.effects$variable <- row.names(fixed.effects)

fixed.effects <- fixed.effects %>% 
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_sel_intercept"))

names(fixed.effects) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")

# fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
#                             "Human footprint index", "Land cover")

fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
                            "Human footprint index", "Forest distance", "Small Woody Features")


fixedSD <- fixed.effects %>% 
  # filter(variable != "Human footprint index") %>%
  ggplot(aes(x = variable, 
             # y = median,
             y = mean,
             colour = variable)) + 
  # geom_errorbar(aes(ymin = lower, ymax = higher), width =.1) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width =.1) +
  geom_line() +
  geom_point() + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  coord_flip() + 
  labs(title = "Covariate effects Sika deer", y = "Median", x = "", colour = "Covariate") +
  NULL


# map outputs
SD_pred_lin <- readRDS("server_outputs/SD_pred_lin_def.RDS")
SD_pred_resp <- readRDS("server_outputs/SD_pred_resp_def.RDS")
SD_spde_PO <- readRDS("server_outputs/SD_spde_PO_def.RDS")
SD_spde_PA <- readRDS("server_outputs/SD_spde_PA_def.RDS")

ireland_proj <- st_transform(ireland, SD_pred_lin@proj4string)
ireland_sp <- as_Spatial(ireland_proj)


pred_lin <- stack(raster(SD_pred_lin['mean']), raster(SD_pred_lin['sd']))
pred_resp <- stack(raster(SD_pred_resp['mean']), raster(SD_pred_resp['sd']))
spde_PO <- stack(raster(SD_spde_PO['mean']), raster(SD_spde_PO['sd']))
spde_PA <- stack(raster(SD_spde_PA['mean']), raster(SD_spde_PA['sd']))


# Plot 

plotlayer <- pred_resp

nl <- nlayers(plotlayer)
m <- matrix(1:nl, ncol = 2)
for (i in 1:nl){
  p <- levelplot(plotlayer, layers=i, col.regions = viridis(16), 
                 xlab = "Longitude", ylab = "Latitude", 
                 zscaleLog = FALSE,
                 # xlim = c(-22, -7), ylim = c(15, 33), 
                 margin = FALSE, main = names(plotlayer)[i]) + 
    latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
  print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
}

spdes <- stack(spde_PO$mean, spde_PA$mean)
names(spdes) <- c("PO", "PA")

levelplot(spdes, col.regions = viridis(16), 
          xlab = "Longitude", ylab = "Latitude", 
          zscaleLog = FALSE,
          # xlim = c(-22, -7), ylim = c(15, 33), 
          margin = FALSE, main = names(spdes)) + 
  latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))

## Fallow deer ##

# load model
mdl_FD <- readRDS("server_outputs/mdl_FD_def.RDS")

fixed.effects <- mdl_FD$summary.fixed
fixed.effects$variable <- row.names(fixed.effects)

fixed.effects <- fixed.effects %>% 
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_sel_intercept"))

names(fixed.effects) <- c("mean", "FD", "lower", "median", "higher", "mode", "kld", "variable")

# fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
#                             "Human footprint index", "Land cover")

fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
                            "Human footprint index", "Forest distance", "Small Woody Features")


fixedFD <- fixed.effects %>% 
  # filter(variable != "Human footprint index") %>%
  ggplot(aes(x = variable, 
             # y = median,
             y = mean,
             colour = variable)) + 
  # geom_errorbar(aes(ymin = lower, ymax = higher), width =.1) +
  geom_errorbar(aes(ymin = mean - FD, ymax = mean + FD), width =.1) +
  geom_line() +
  geom_point() + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  coord_flip() + 
  labs(title = "Covariate effects Fallow deer", y = "Median", x = "", colour = "Covariate") +
  NULL

cowplot::plot_grid(plotlist = list(fixedRD, fixedSD, fixedFD), ncol = 2)

# map outputs
FD_pred_lin <- readRDS("server_outputs/FD_pred_lin_def.RDS")
FD_pred_resp <- readRDS("server_outputs/FD_pred_resp_def.RDS")
FD_spde_PO <- readRDS("server_outputs/FD_spde_PO_def.RDS")
FD_spde_PA <- readRDS("server_outputs/FD_spde_PA_def.RDS")

ireland_proj <- st_transform(ireland, FD_pred_lin@proj4string)
ireland_sp <- as_Spatial(ireland_proj)


pred_lin <- stack(raster(FD_pred_lin['mean']), raster(FD_pred_lin['sd']))
pred_resp <- stack(raster(FD_pred_resp['mean']), raster(FD_pred_resp['sd']))
spde_PO <- stack(raster(FD_spde_PO['mean']), raster(FD_spde_PO['sd']))
spde_PA <- stack(raster(FD_spde_PA['mean']), raster(FD_spde_PA['sd']))


# Plot 

plotlayer <- pred_resp

nl <- nlayers(plotlayer)
m <- matrix(1:nl, ncol = 2)
for (i in 1:nl){
  p <- levelplot(plotlayer, layers=i, col.regions = viridis(16), 
                 xlab = "Longitude", ylab = "Latitude", 
                 zscaleLog = FALSE,
                 # xlim = c(-22, -7), ylim = c(15, 33), 
                 margin = FALSE, main = names(plotlayer)[i]) + 
    latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
  print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
}

spdes <- stack(spde_PO$mean, spde_PA$mean)
names(spdes) <- c("PO", "PA")

levelplot(spdes, col.regions = viridis(16), 
          xlab = "Longitude", ylab = "Latitude", 
          zscaleLog = FALSE,
          # xlim = c(-22, -7), ylim = c(15, 33), 
          margin = FALSE, main = names(spdes)) + 
  latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))



### Plot predictions in the log scale

preds <- stack(raster(RD_pred_resp["mean"]), raster(SD_pred_resp["mean"]), 
              raster(FD_pred_resp["mean"]))

names(preds) <- c("Red deer", "Sika deer", "Fallow deer")

plotlayer <- preds

nl <- nlayers(plotlayer)
m <- matrix(1:nl, ncol = 3)
for (i in 1:nl){
  p <- levelplot(plotlayer, layers=i, col.regions = viridis(16), 
                 xlab = "Longitude", ylab = "Latitude", 
                 zscaleLog = FALSE,
                 # xlim = c(-22, -7), ylim = c(15, 33), 
                 margin = FALSE, main = names(plotlayer)[i]) + 
    latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
  print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
}

