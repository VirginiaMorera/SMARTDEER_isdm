
#### load all data ####
ireland <- st_read("data/ireland_ITM.shp")

mdl_RD <- readRDS("server_outputs/mdl_RD_def.RDS")
RD_pred_lin <- readRDS("server_outputs/RD_pred_lin_def.RDS")
RD_pred_resp <- readRDS("server_outputs/RD_pred_resp_def.RDS")
RD_spde_PO <- readRDS("server_outputs/RD_spde_PO_def.RDS")
RD_spde_PA <- readRDS("server_outputs/RD_spde_PA_def.RDS")

mdl_SD <- readRDS("server_outputs/mdl_SD_def.RDS")
SD_pred_lin <- readRDS("server_outputs/SD_pred_lin_def.RDS")
SD_pred_resp <- readRDS("server_outputs/SD_pred_resp_def.RDS")
SD_spde_PO <- readRDS("server_outputs/SD_spde_PO_def.RDS")
SD_spde_PA <- readRDS("server_outputs/SD_spde_PA_def.RDS")

mdl_FD <- readRDS("server_outputs/mdl_FD_def.RDS")
FD_pred_lin <- readRDS("server_outputs/FD_pred_lin_def.RDS")
FD_pred_resp <- readRDS("server_outputs/FD_pred_resp_def.RDS")
FD_spde_PO <- readRDS("server_outputs/FD_spde_PO_def.RDS")
FD_spde_PA <- readRDS("server_outputs/FD_spde_PA_def.RDS")

ireland_proj <- st_transform(ireland, RD_pred_lin@proj4string)
ireland_sp <- as_Spatial(ireland_proj)
ireland_latlon <- spTransform(ireland_sp, CRSobj = CRS("EPSG:4326"))


#### Plot fixed effects ####

## Red deer 

fixed.effectsRD <- mdl_RD$summary.fixed
fixed.effectsRD$variable <- row.names(fixed.effectsRD)

fixed.effectsRD <- fixed.effectsRD %>% 
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_sel_intercept"))

names(fixed.effectsRD) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")

fixed.effectsRD$variable <- c("Tree cover density", "Elevation", "Slope",
                            "Human footprint index", "Forest distance", "Small Woody Features")

fixed.effectsRD %>% 
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

## Sika deer 

fixed.effectsSD <- mdl_SD$summary.fixed
fixed.effectsSD$variable <- row.names(fixed.effectsSD)

fixed.effectsSD <- fixed.effectsSD %>% 
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_sel_intercept"))

names(fixed.effectsSD) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")

# fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
#                             "Human footprint index", "Land cover")

fixed.effectsSD$variable <- c("Tree cover density", "Elevation", "Slope",
                              "Human footprint index", "Forest distance", "Small Woody Features")


fixed.effectsSD %>% 
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

## Fallow deer

fixed.effectsFD <- mdl_FD$summary.fixed
fixed.effectsFD$variable <- row.names(fixed.effectsFD)

fixed.effectsFD <- fixed.effectsFD %>% 
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_sel_intercept"))

names(fixed.effectsFD) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")

# fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
#                             "Human footprint index", "Land cover")

fixed.effectsFD$variable <- c("Tree cover density", "Elevation", "Slope",
                              "Human footprint index", "Forest distance", "Small Woody Features")


fixed.effectsFD %>% 
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

fixed.effectsRD$Species <- "Red deer"
fixed.effectsSD$Species <- "Sika deer"
fixed.effectsFD$Species <- "Fallow deer"

fixed.effects_all <- bind_rows(fixed.effectsFD, fixed.effectsRD, fixed.effectsSD)

fixed.effects_all <- fixed.effects_all %>% 
  mutate(Species = factor(Species, levels = c("Red deer", "Sika deer", "Fallow deer")))

covarplot <- fixed.effects_all %>% 
  ggplot(aes(x = variable, 
             y = mean,
             colour = Species)) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width =.1, position = position_dodge(width = -0.3)) +
  geom_line(position = position_dodge(width = -0.3)) +
  geom_point(position = position_dodge(width = -0.3)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  scale_colour_colorblind() + 
  coord_flip() + 
  labs(title = "Covariate effects", y = "Effect size", x = "", colour = "Species") +
  # facet_wrap(~Species) +
  NULL

Cairo::CairoPDF(file = "Fig2.pdf", width = 8, height = 6)
covarplot
dev.off()


#### Get predictions ####

## Red deer ##

pred_linRD <- stack(raster(RD_pred_lin['mean']), raster(RD_pred_lin['sd']))
pred_respRD <- stack(raster(RD_pred_resp['mean']), raster(RD_pred_resp['sd']))
spde_PORD <- stack(raster(RD_spde_PO['mean']), raster(RD_spde_PO['sd']))
spde_PARD <- stack(raster(RD_spde_PA['mean']), raster(RD_spde_PA['sd']))

## Sika deer ##

pred_linSD <- stack(raster(SD_pred_lin['mean']), raster(SD_pred_lin['sd']))
pred_respSD <- stack(raster(SD_pred_resp['mean']), raster(SD_pred_resp['sd']))
spde_POSD <- stack(raster(SD_spde_PO['mean']), raster(SD_spde_PO['sd']))
spde_PASD <- stack(raster(SD_spde_PA['mean']), raster(SD_spde_PA['sd']))

## Fallow deer ##

pred_linFD <- stack(raster(FD_pred_lin['mean']), raster(FD_pred_lin['sd']))
pred_respFD <- stack(raster(FD_pred_resp['mean']), raster(FD_pred_resp['sd']))
spde_POFD <- stack(raster(FD_spde_PO['mean']), raster(FD_spde_PO['sd']))
spde_PAFD <- stack(raster(FD_spde_PA['mean']), raster(FD_spde_PA['sd']))


#### Rescale lin predictions from 0 to 1 ####

rescaled_lins <- stack(rescale0to1(pred_linRD$mean), 
                       rescale0to1(pred_linSD$mean), 
                       rescale0to1(pred_linFD$mean))

rescaled_lins_latlon <- projectRaster(rescaled_lins, crs = CRS("EPSG:4326"))
rescaled_lins_latlon <- mask(rescaled_lins_latlon, ireland_latlon)


#### Plot ####

## prediction means 

predplot <- levelplot(rescaled_lins_latlon, col.regions = viridis(100), 
                      ylab = "Latitude", 
                      zscaleLog = FALSE,
                      names.attr = c("", "", ""),
                      xlab=NULL,
                      margin = FALSE) + 
  latticeExtra::layer(sp.polygons(ireland_latlon, lwd = 1, col = "darkgray"))

Cairo::CairoPDF(file = "Fig3a.pdf", height = 5, width = 10)
predplot
dev.off()

## sd 
sds <- stack(pred_linRD$sd, pred_linSD$sd, pred_linFD$sd)

sds_latlon <- projectRaster(sds, crs = CRS("EPSG:4326"))
sds_latlon <- mask(sds_latlon, ireland_latlon)

sdplot <- levelplot(sds_latlon, col.regions = viridis(100), 
                      ylab = "Latitude", 
                      zscaleLog = FALSE,
                      names.attr = c("", "", ""),
                      xlab = "Longitude",
                      margin = FALSE) + 
  latticeExtra::layer(sp.polygons(ireland_latlon, lwd = 1, col = "darkgray"))

Cairo::CairoPDF(file = "Fig3b.pdf", height = 5, width = 10)
sdplot
dev.off()

Cairo::CairoPDF(file = "Fig3.pdf", height = 10, width = 12)
cowplot::plot_grid(predplot, sdplot, nrow = 2)
dev.off()

#### Plot spatial fields####

## spde means

spdes <- stack(raster(RD_spde_PO['mean']), raster(SD_spde_PO['mean']), raster(FD_spde_PO['mean']), 
               raster(RD_spde_PA['mean']), raster(SD_spde_PA['mean']), raster(FD_spde_PA['mean']))

spdes_latlon <- projectRaster(spdes, crs = CRS("EPSG:4326"))
spdes_latlon <- mask(spdes_latlon, ireland_latlon)


spdesplot <- levelplot(spdes_latlon, col.regions = viridis(100), 
                      ylab = "Latitude", 
                      zscaleLog = FALSE,
                      names.attr = c("", "", "", "", "", ""),
                      xlab = "Longitude",
                      layout=c(3, 2),
                      margin = FALSE) + 
  latticeExtra::layer(sp.polygons(ireland_latlon, lwd = 1, col = "darkgray"))

Cairo::CairoPDF(file = "FigS1.pdf", height = 10, width = 10)
predplot
dev.off()

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

lins <- stack(raster(RD_pred_lin["mean"]), raster(SD_pred_lin["mean"]), 
               raster(FD_pred_lin["mean"]))

names(preds) <- c("Red deer", "Sika deer", "Fallow deer")
names(lins) <- c("Red deer", "Sika deer", "Fallow deer")

plotlayer <- lins

nl <- nlayers(plotlayer)
m <- matrix(1:nl, ncol = 3)
for (i in 1:nl){
  p <- levelplot(plotlayer, layers=i, col.regions = viridis(16), 
                 xlab = "Longitude", ylab = "Latitude", 
                 zscaleLog = F,
                 # xlim = c(-22, -7), ylim = c(15, 33), 
                 margin = FALSE, main = names(plotlayer)[i]) + 
    latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
  print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
}

Cairo::CairoPDF(file = "fallowdeer.pdf", width = 8, height = 10)
p
dev.off()




### Binary maps and multispecies map 

quants <- c(0.25, 0.5, 0.75, 0.95)

sumlist <- stack()
for (i in seq_along(quants)) {
  # i = 1
  binaries <- rescaled_lins
  binaries[binaries >= quants[i]] <- 1
  binaries[binaries < quants[i]] <- 0
  binsum <- sum(binaries)
  sumlist <- stack(sumlist, binsum)
}

names(sumlist) <- paste0("Quantile ", quants)
plot(sumlist)

levelplot(sumlist,
          col.regions = viridis(6),
          layout=c(2, 2),
          xlab = "Longitude", ylab = "Latitude", 
          zscaleLog = FALSE,
          # xlim = c(-22, -7), ylim = c(15, 33), 
          margin = FALSE) + 
  latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))

