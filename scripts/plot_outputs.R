source("scripts/setups.R")

#### load all data ####
ireland <- st_read("data/ireland_ITM.shp")

mdl_RD <- readRDS("server_outputs/RedDeer_mdl_copy.RDS")
RD_pred <- readRDS("server_outputs/RedDeer_prediction.RDS")
RD_spde <- readRDS("server_outputs/RedDeer_spatial_field.RDS")
RD_covar <- readRDS("server_outputs/RedDeer_covariate_effect.RDS")

mdl_SD <- readRDS("server_outputs/SikaDeer_mdl_copy.RDS")
SD_pred <- readRDS("server_outputs/SikaDeer_prediction.RDS")
SD_spde <- readRDS("server_outputs/SikaDeer_spatial_field.RDS")
SD_covar <- readRDS("server_outputs/SikaDeer_covariate_effect.RDS")

mdl_FD <- readRDS("server_outputs/FallowDeer_mdl_copy.RDS")
FD_pred <- readRDS("server_outputs/FallowDeer_prediction.RDS")
FD_spde <- readRDS("server_outputs/FallowDeer_spatial_field.RDS")
FD_covar <- readRDS("server_outputs/FallowDeer_covariate_effect.RDS")

ireland_proj <- st_transform(ireland, RD_pred$predictions@proj4string)
ireland_sp <- as_Spatial(ireland_proj)
ireland_latlon <- spTransform(ireland_sp, CRSobj = CRS("EPSG:4326"))


#### Plot fixed effects ####

## Red deer 

fixed.effectsRD <- mdl_RD$summary.fixed
fixed.effectsRD$variable <- row.names(fixed.effectsRD)

fixed.effectsRD <- fixed.effectsRD %>% 
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_NI_intercept", "PA_data_ROI_intercept"))

names(fixed.effectsRD) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")

fixed.effectsRD$variable <- c("Tree cover density", "Elevation", "Slope",
                            "Human footprint index", "Distance to forest edge", "Small Woody Features")

fixed.effectsRD %>% 
  # filter(variable != "Human footprint index") %>%
  ggplot(aes(x = variable, 
             y = median,
             # y = mean,
             colour = variable)) + 
  geom_errorbar(aes(ymin = lower, ymax = higher), width =.1) +
  # geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width =.1) +
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
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_NI_intercept", "PA_data_ROI_intercept"))

names(fixed.effectsSD) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")

# fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
#                             "Human footprint index", "Land cover")

fixed.effectsSD$variable <- c("Tree cover density", "Elevation", "Slope",
                              "Human footprint index", "Distance to forest edge", "Small Woody Features")


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
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_NI_intercept", "PA_data_ROI_intercept"))

names(fixed.effectsFD) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")

# fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
#                             "Human footprint index", "Land cover")

fixed.effectsFD$variable <- c("Tree cover density", "Elevation", "Slope",
                              "Human footprint index", "Distance to forest edge", "Small Woody Features")


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
  mutate(Species = factor(Species, levels = c("Red deer", "Sika deer", "Fallow deer")), 
         variable = factor(variable, levels = c("Tree cover density", "Distance to forest edge", 
                                                "Small Woody Features", "Elevation", "Slope", 
                                                "Human footprint index")))

covarplot <- fixed.effects_all %>% 
  ggplot(aes(x = variable, 
             # y = mean,
             y = median, 
             colour = Species)) + 
  geom_errorbar(aes(ymin = lower, ymax = higher), width =.1, position = position_dodge(width = -0.3)) +
  geom_errorbar(aes(ymin = lower, ymax = higher), width =.1, position = position_dodge(width = -0.3)) +
  geom_line(position = position_dodge(width = -0.3)) +
  geom_point(position = position_dodge(width = -0.3)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  scale_colour_colorblind() + 
  coord_flip() + 
  labs(title = "Covariate effects", y = "Effect size", x = "", colour = "Species") +
  # facet_wrap(~Species) +
  NULL

Cairo::CairoPDF(file = "server_outputs/Fig2_rev.pdf", width = 8, height = 6)
covarplot
dev.off()


#### Get predictions ####

## Red deer ##

pred_respRD <- stack(raster(RD_pred$predictions['mean']), raster(RD_pred$predictions['sd']))
spde_RD <- stack(raster(RD_spde$predictions['mean']), raster(RD_spde$predictions['sd']))
covars_RD <- stack(raster(RD_covar$predictions['mean']), raster(RD_covar$predictions['sd']))

## Sika deer ##

pred_respSD <- stack(raster(SD_pred$predictions['mean']), raster(SD_pred$predictions['sd']))
spde_SD <- stack(raster(SD_spde$predictions['mean']), raster(SD_spde$predictions['sd']))
covars_SD <- stack(raster(SD_covar$predictions['mean']), raster(SD_covar$predictions['sd']))

## Fallow deer ##

pred_respFD <- stack(raster(FD_pred$predictions['mean']), raster(FD_pred$predictions['sd']))
spde_FD <- stack(raster(FD_spde$predictions['mean']), raster(FD_spde$predictions['sd']))
covars_FD <- stack(raster(FD_covar$predictions['mean']), raster(FD_covar$predictions['sd']))


#### Rescale predictions from 0 to 1 ####

rescaled_preds <- stack(rescale0to1(pred_respRD$mean), 
                       rescale0to1(pred_respSD$mean), 
                       rescale0to1(pred_respFD$mean))

rescaled_preds_latlon <- projectRaster(rescaled_preds, crs = CRS("EPSG:4326"))
rescaled_preds_latlon <- mask(rescaled_preds_latlon, ireland_latlon)


#### Plot ####

## prediction means 

predplot <- levelplot(rescaled_preds_latlon, col.regions = viridis(100), 
                      ylab = "Latitude", 
                      zscaleLog = FALSE,
                      names.attr = c("Red deer", "Sika deer", "Fallow deer"),
                      xlab=NULL,
                      margin = FALSE) + 
  latticeExtra::layer(sp.polygons(ireland_latlon, lwd = 1, col = "darkgray"))

Cairo::CairoPDF(file = "outputs/Fig3a_rev.pdf", height = 5, width = 10)
predplot
dev.off()

## sd 
sds <- stack(pred_respRD$sd, pred_respSD$sd, pred_respFD$sd)

sds_latlon <- projectRaster(sds, crs = CRS("EPSG:4326"))
sds_latlon <- mask(sds_latlon, ireland_latlon)

sdplot <- levelplot(sds_latlon, col.regions = viridis(100), 
                      ylab = "Latitude", 
                      zscaleLog = FALSE,
                      names.attr = c("", "", ""),
                      xlab = "Longitude",
                      margin = FALSE) + 
  latticeExtra::layer(sp.polygons(ireland_latlon, lwd = 1, col = "darkgray"))

Cairo::CairoPDF(file = "outputs/Fig3b_rev.pdf", height = 5, width = 10)
sdplot
dev.off()

Cairo::CairoPDF(file = "server_outputs/Fig3_rev.pdf", height = 10, width = 12)
cowplot::plot_grid(predplot, sdplot, nrow = 2, labels = c('A', 'B'))
dev.off()

#### Plot spatial fields####

rescaled_spde_means <- stack(rescale0to1(spde_RD$mean), 
                             rescale0to1(spde_SD$mean), 
                             rescale0to1(spde_FD$mean))

rescaled_spde_means_latlon <- projectRaster(rescaled_spde_means, crs = CRS("EPSG:4326"))
rescaled_spde_means_latlon <- mask(rescaled_spde_means_latlon, ireland_latlon)


## spde means

spdeplot <- levelplot(rescaled_spde_means_latlon, col.regions = viridis(100), 
                      ylab = "Latitude", 
                      zscaleLog = FALSE,
                      names.attr = c("", "", ""),
                      xlab=NULL,
                      margin = FALSE) + 
  latticeExtra::layer(sp.polygons(ireland_latlon, lwd = 1, col = "darkgray"))

spdeplot

## sd 

spde_sds <- stack(spde_RD$sd, spde_SD$sd, spde_FD$sd)

spde_sds_latlon <- projectRaster(spde_sds, crs = CRS("EPSG:4326"))
spde_sds_latlon <- mask(spde_sds_latlon, ireland_latlon)

sd_spde_plot <- levelplot(spde_sds_latlon, col.regions = viridis(100), 
                    ylab = "Latitude", 
                    zscaleLog = FALSE,
                    names.attr = c("", "", ""),
                    xlab = "Longitude",
                    margin = FALSE) + 
  latticeExtra::layer(sp.polygons(ireland_latlon, lwd = 1, col = "darkgray"))

sd_spde_plot

Cairo::CairoPDF(file = "outputs/FigS1_rev.pdf", height = 10, width = 12)
cowplot::plot_grid(spdeplot, sd_spde_plot, nrow = 2)
dev.off()


# # Plot 
# 
# plotlayer <- pred_resp
# 
# nl <- nlayers(plotlayer)
# m <- matrix(1:nl, ncol = 2)
# for (i in 1:nl){
#   p <- levelplot(plotlayer, layers=i, col.regions = viridis(16), 
#                  xlab = "Longitude", ylab = "Latitude", 
#                  zscaleLog = FALSE,
#                  # xlim = c(-22, -7), ylim = c(15, 33), 
#                  margin = FALSE, main = names(plotlayer)[i]) + 
#     latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
#   print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
# }
#   
# spdes <- stack(spde_PO$mean, spde_PA$mean)
# names(spdes) <- c("PO", "PA")
# 
# levelplot(spdes, col.regions = viridis(16), 
#           xlab = "Longitude", ylab = "Latitude", 
#           zscaleLog = FALSE,
#           # xlim = c(-22, -7), ylim = c(15, 33), 
#           margin = FALSE, main = names(spdes)) + 
#   latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
# 
# 
# ## Sika deer ##
# 
# plotlayer <- pred_resp
# 
# nl <- nlayers(plotlayer)
# m <- matrix(1:nl, ncol = 2)
# for (i in 1:nl){
#   p <- levelplot(plotlayer, layers=i, col.regions = viridis(16), 
#                  xlab = "Longitude", ylab = "Latitude", 
#                  zscaleLog = FALSE,
#                  # xlim = c(-22, -7), ylim = c(15, 33), 
#                  margin = FALSE, main = names(plotlayer)[i]) + 
#     latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
#   print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
# }
# 
# spdes <- stack(spde_PO$mean, spde_PA$mean)
# names(spdes) <- c("PO", "PA")
# 
# levelplot(spdes, col.regions = viridis(16), 
#           xlab = "Longitude", ylab = "Latitude", 
#           zscaleLog = FALSE,
#           # xlim = c(-22, -7), ylim = c(15, 33), 
#           margin = FALSE, main = names(spdes)) + 
#   latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
# 
# ## Fallow deer ##
# 
# plotlayer <- pred_resp
# 
# nl <- nlayers(plotlayer)
# m <- matrix(1:nl, ncol = 2)
# for (i in 1:nl){
#   p <- levelplot(plotlayer, layers=i, col.regions = viridis(16), 
#                  xlab = "Longitude", ylab = "Latitude", 
#                  zscaleLog = FALSE,
#                  # xlim = c(-22, -7), ylim = c(15, 33), 
#                  margin = FALSE, main = names(plotlayer)[i]) + 
#     latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
#   print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
# }
# 
# spdes <- stack(spde_PO$mean, spde_PA$mean)
# names(spdes) <- c("PO", "PA")
# 
# levelplot(spdes, col.regions = viridis(16), 
#           xlab = "Longitude", ylab = "Latitude", 
#           zscaleLog = FALSE,
#           # xlim = c(-22, -7), ylim = c(15, 33), 
#           margin = FALSE, main = names(spdes)) + 
#   latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
# 
# 
# 
# ### Plot predictions in the log scale
# 
# preds <- stack(raster(RD_pred_resp["mean"]), raster(SD_pred_resp["mean"]), 
#               raster(FD_pred_resp["mean"]))
# 
# lins <- stack(raster(RD_pred_lin["mean"]), raster(SD_pred_lin["mean"]), 
#                raster(FD_pred_lin["mean"]))
# 
# names(preds) <- c("Red deer", "Sika deer", "Fallow deer")
# names(lins) <- c("Red deer", "Sika deer", "Fallow deer")
# 
# plotlayer <- lins
# 
# nl <- nlayers(plotlayer)
# m <- matrix(1:nl, ncol = 3)
# for (i in 1:nl){
#   p <- levelplot(plotlayer, layers=i, col.regions = viridis(16), 
#                  xlab = "Longitude", ylab = "Latitude", 
#                  zscaleLog = F,
#                  # xlim = c(-22, -7), ylim = c(15, 33), 
#                  margin = FALSE, main = names(plotlayer)[i]) + 
#     latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
#   print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
# }
# 
# Cairo::CairoPDF(file = "fallowdeer.pdf", width = 8, height = 10)
# p
# dev.off()
# 
# 
# 
# 
# ### Binary maps and multispecies map 
# 
# quants <- c(0.25, 0.5, 0.75, 0.95)
# 
# sumlist <- stack()
# for (i in seq_along(quants)) {
#   # i = 1
#   binaries <- rescaled_lins_latlon
#   binaries[binaries >= quants[i]] <- 1
#   binaries[binaries < quants[i]] <- 0
#   binsum <- sum(binaries)
#   sumlist <- stack(sumlist, binsum)
# }
# 
# names(sumlist) <- paste0("Quantile ", quants)
# plot(sumlist)
# 
# levelplot(sumlist$Quantile.0.5,
#           # col.regions = viridis(3),
#           # layout=c(2, 2),
#           xlab = "Longitude", ylab = "Latitude", 
#           zscaleLog = FALSE,
#           # xlim = c(-22, -7), ylim = c(15, 33), 
#           margin = FALSE) + 
#   latticeExtra::layer(sp.polygons(ireland_sp, lwd = 1, col = "darkgray"))
# 
# 
# x <- sumlist$Quantile.0.5
# # x <- projectRaster(x, crs = ireland_latlon@proj4string)
# plot(x)
# r <- ratify(x)
# rat <- levels(r)[[1]]
# rat$landcover <- c("0", "1", "2", "3")
# levels(r) <- rat
# levelplot(r)
# 
# 
# levelplot(r,
#           col.regions = viridis(4),
#           # layout=c(2, 2),
#           xlab = "Longitude", ylab = "Latitude", 
#           zscaleLog = FALSE,
#           # xlim = c(-22, -7), ylim = c(15, 33), 
#           margin = FALSE) + 
#   latticeExtra::layer(sp.polygons(ireland_latlon, lwd = 1, col = "darkgray"))
