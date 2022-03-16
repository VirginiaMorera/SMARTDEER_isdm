rm(list = ls())
devtools::install_github('PhilipMostert/inlabruSDMs', ref = "main")
library(inlabruSDMs)
library(raster)
library(tidyverse)
library(ggthemes)
`%!in%` = Negate(`%in%`)
# library(sf)
setwd("~/scratch/ISDM")
in_bound <- readRDS("data/inner_boundary.RDS")
mesh0 <- readRDS("data/mesh.RDS")
mesh1 <- readRDS("data/meshLarge.RDS")
mesh2 <- readRDS("data/meshVLarge.RDS")
covar_stack <- stack("data/covar_subset_KM.gri")
covar_scaled <- scale(covar_stack)
spatialcovs <- as(covar_scaled, "SpatialPixelsDataFrame")


PO_data <- read.csv("data/PO_data_all.csv", row.names = NULL)
PO_data_sel <- PO_data %>% 
        filter(Source %!in% c("CameraTraps")) %>% 
        filter(Species %in% c("SikaDeer"))
PO_data_sp <- SpatialPointsDataFrame(coords = cbind(PO_data_sel$X, PO_data_sel$Y), 
                                     data = data.frame(Count = PO_data_sel$Count, 
                                                       PO = PO_data_sel$Count),
                                     proj4string = in_bound@proj4string)

PA_data <- read.csv("data/PA_data_all.csv", row.names = NULL)
PA_data_sel <- PA_data %>%
        filter(Species %in% c("SikaDeer"))

PA_data_sp <- SpatialPointsDataFrame(coords = cbind(PA_data_sel$X, PA_data_sel$Y), 
                                     data = data.frame(PA = PA_data_sel$PA),
                                     proj4string = in_bound@proj4string)


# select the mesh we're going to use
sel_mesh <- mesh1

pix <- pixels(mesh = sel_mesh, mask = in_bound)
pix@proj4string <- sel_mesh$crs

pix_large <- pixels(mesh = sel_mesh, mask = in_bound, nx = 603, ny = 659)
pix_large@proj4string <- sel_mesh$crs

SD_data_for_model <- organize_data(PA_data_sel, PO_data_sel, 
                                   countresp = "PO", paresp = "PA", 
                                   coords = c("X", "Y"), 
                                   proj = sel_mesh$crs,
                                   marks = F,  
                                   mesh = sel_mesh, 
                                   speciesname = "Species",
                                   boundary = in_bound)


SD_spdemodel_PO = INLA::inla.spde2.pcmatern(mesh = sel_mesh, alpha = 3/2,  ### mesh and smoothness parameter
                                            prior.range = c(40, 0.05), ### P(range < 40km) = 0.05 Very small probability the range is smaller than 40 (triangle edge size) (with high likelihood the range is between 20 (res of the model) and infinity i.e. flat prior)
                                            prior.sigma = c(10, 0.05)) ### P(sigma > 2) = 0.05 Very small probability that sigma is larger than 0.2. Sigma is most likely between 0 and 2

SD_spdemodel_PA = INLA::inla.spde2.pcmatern(mesh = sel_mesh, alpha = 3/2,  ### mesh and smoothness parameter
                                            prior.range = c(100, 0.05), ### P(range < 40km) = 0.05 Very small probability the range is smaller than 40 (triangle edge size) (with high likelihood the range is between 20 (res of the model) and infinity i.e. flat prior)
                                            prior.sigma = c(6, 0.05)) ### P(sigma > 2) = 0.05 Very small probability that sigma is larger than 0.2. Sigma is most likely between 0 and 2


mdl_SD <- bru_sdm(data = SD_data_for_model,
                  spatialcovariates = covar_scaled,
                  covariatestoinclude = c("elevation", "tree_cover_density", "human_footprint_index", #"landCover",
                                          "forest_distances", "small_woody_features", "slope"),
                  # covariatesbydataset = list(PO_data_sel = c("elevation", "tree_cover_density", "human_footprint_index",
                  #                                            "forest_distances", "small_woody_features", "slope"),
                  #                            PA_data_sel = c("elevation", "tree_cover_density",
                  #                                            "forest_distances", "small_woody_features", "slope")),
                  specieseffects = TRUE,
                  pointsintercept = TRUE,
                  marksintercept = FALSE,
                  sharedspatial = F,
                  spdemodel =  list(PO_data_sel = SD_spdemodel_PO, 
                                    PA_data_sel = SD_spdemodel_PA),
                  # spdemodel = RD_spdemodel_shared, 
                  pointsspatial = TRUE,
                  marksspatial = FALSE,
                  spatialdatasets = NULL,
                  tolerance = NULL)

# saveRDS(mdl_SD, file = "mdl_SD_def.RDS")
# mdl_SD <- readRDS("mdl_SD_def.RDS")

summary(mdl_SD)


fixed.effects <- mdl_SD$summary.fixed
fixed.effects$variable <- row.names(fixed.effects)

fixed.effects <- fixed.effects %>% 
        filter(variable %!in% c("PO_data_sel_intercept", "PA_data_sel_intercept"))

names(fixed.effects) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")
fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
                            "Human footprint index", "Forest distances", "Small woody features")

ggplot(fixed.effects, aes(x = variable, y = median, colour = variable)) + 
        geom_errorbar(aes(ymin = lower, ymax = higher), width =.1) +
        geom_line() +
        geom_point() + 
        geom_hline(yintercept = 0) + 
        theme_bw() + 
        coord_flip() + 
        labs(title = "Sika deer", y = "Median", x = "", colour = "Covariate") +
        NULL

SD_pred_lin <- predict(mdl_SD, 
                       data = pix, 
                       formula = ~(tree_cover_density + elevation + slope + human_footprint_index +
                                           forest_distances + small_woody_features +
                                           PA_data_sel_spde +
                                           PO_data_sel_spde +
                                           PO_data_sel_intercept + PA_data_sel_intercept))

plot(SD_pred_lin['mean'], main = "Mean of prediction (linear scale) Sika deer")
plot(SD_pred_lin['sd'], main = "SD of prediction (linear scale) Sika deer")
saveRDS(SD_pred_lin, file = "SD_pred_lin_def.RDS")

plot(SD_pred_lin ['mean'])

SD_pred_resp <- predict(mdl_SD, 
                        data = pix,
                        formula = ~exp(tree_cover_density + elevation + slope + human_footprint_index + 
                                               forest_distances + small_woody_features +
                                               # shared_spatial +
                                               PO_data_sel_intercept + PA_data_sel_intercept +
                                               PA_data_sel_spde +
                                               PO_data_sel_spde))

plot(SD_pred_resp['mean'], main = "Mean of prediction (response scale) Sika deer")
plot(SD_pred_resp['sd'], main = "SD of prediction (response scale) Sika deer")
saveRDS(SD_pred_resp, file = "SD_pred_resp_def.RDS")

SD_covars <- predict(mdl_SD, 
                     data = pix, 
                     formula = ~exp(tree_cover_density + elevation + slope + human_footprint_index + 
                                            forest_distances + small_woody_features +
                                            # shared_spatial 
                                            PO_data_sel_intercept + PA_data_sel_intercept))

plot(SD_covars['mean'], main = "Mean of the covariate effect (response scale)")
plot(SD_covars['sd'], main = "SD of the covariate effect (response scale)")

SD_spde_PO <- predict(mdl_SD, 
                      data = pix, 
                      formula = ~(PO_data_sel_spde))

plot(SD_spde_PO["mean"], main = "Mean of spatial effect PO")
plot(SD_spde_PO["sd"], main = "Sd of spatial effect PO")
saveRDS(SD_spde_PO, file = "SD_spde_PO_def.FDS")


SD_spde_PA <- predict(mdl_SD, 
                      data = pix, 
                      formula = ~exp(PA_data_sel_spde))

plot(SD_spde_PA["mean"], main = "Mean of spatial effect PA")
plot(SD_spde_PA["sd"], main = "Sd of spatial effect PA")
saveRDS(SD_spde_PA, file = "SD_spde_PA_def.FDS")


ips <- ipoints(in_bound, mesh1)
mdl1_SD_Ab <- predict(mdl_SD, 
                   ips,
                   formula = ~sum(weight*exp((tree_cover_density + elevation + slope + human_footprint_index + 
                                                      forest_distances + small_woody_features +
                                                      PA_data_sel_spde +
                                                      PO_data_sel_spde))))


loo_PO <- leave_one_out(mdl_SD, dataset = "PO_data_sel")
loo_PA <- leave_one_out(mdl_SD, dataset = "PA_data_sel")

loo_PO$Leaving_out_PO_data_sel
loo_PA$Leaving_out_PA_data_sel


PA.fixed <- loo_PO$Leaving_out_PO_data_sel$summary.fixed
PO.fixed <- loo_PA$Leaving_out_PA_data_sel$summary.fixed


PO.fixed$variable <- row.names(PO.fixed)
PA.fixed$variable <- row.names(PA.fixed)


PO.fixed <- PO.fixed %>% 
        filter(variable %!in% c("PO_data_sel_intercept"))

PA.fixed <- PA.fixed %>% 
        filter(variable %!in% c("PA_data_sel_intercept"))

names(PO.fixed) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")
names(PA.fixed) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")

PO.fixed$variable <- c("Tree cover density", "Elevation", "Slope",
                       "Human footprint index", "Forest distance", "Small Woody Features")

PA.fixed$variable <- c("Tree cover density", "Elevation", "Slope",
                       "Human footprint index", "Forest distance", "Small Woody Features")

fixed.effects$model <- "Joint"
PO.fixed$model <- "Count"
PA.fixed$model <- "PA"


all_fixed <- bind_rows(fixed.effects, PO.fixed, PA.fixed)

all_fixed %>% 
        # filter(model != "Count") %>% 
        ggplot(aes(x = variable, y = median, colour = model)) + 
        geom_errorbar(aes(ymin = lower, ymax = higher), width =.1, position = position_dodge(width=0.3)) +
        geom_point(position = position_dodge(width=0.3)) + 
        geom_hline(yintercept = 0) + 
        theme_bw() + 
        coord_flip() +
        scale_colour_colorblind() +
        labs(y = "Effect size", x = "", colour = "Model") +
        NULL

# CairoPNG(filename = "coefficients.png", bg = "transparent", dpi = 800, width = 15, height = 8, units = "cm")
# print(coeffs)
# dev.off()

