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
covar_stack <- stack("large_env_data/covar_subset_KM.gri")
covar_scaled <- scale(covar_stack)
spatialcovs <- as(covar_scaled, "SpatialPixelsDataFrame")


PO_data <- read.csv("data/PO_data_all.csv", row.names = NULL)
PO_data_sel <- PO_data %>% 
        filter(Source %in% c("CullReturnsNI", "MammalSurvey", "Bycatch", "NBDC", "CitizenScience", "webSurvey", "Other")) %>% 
        filter(Species %in% c("RedDeer"))
PO_data_sp <- SpatialPointsDataFrame(coords = cbind(PO_data_sel$X, PO_data_sel$Y), 
                                     data = data.frame(PO = PO_data_sel$PO, 
                                                       Count = PO_data_sel$Count),
                                     proj4string = in_bound@proj4string)

PA_data <- read.csv("data/PA_data_all.csv", row.names = NULL)
PA_data_sel <- PA_data %>%
        filter(Species %in% c("RedDeer"))

PA_data_sp <- SpatialPointsDataFrame(coords = cbind(PA_data_sel$X, PA_data_sel$Y), 
                                     data = data.frame(PA = PA_data_sel$PA),
                                     proj4string = in_bound@proj4string)


# select the mesh we're going to use
sel_mesh <- mesh1

pix <- pixels(mesh = sel_mesh, mask = in_bound)
pix@proj4string <- sel_mesh$crs

pix_large <- pixels(mesh = sel_mesh, mask = in_bound, nx = 603, ny = 659)
pix_large@proj4string <- sel_mesh$crs

RD_data_for_model <- organize_data(PA_data_sel, PO_data_sel, 
                                   countresp = "PO", paresp = "PA", 
                                   coords = c("X", "Y"), 
                                   proj = sel_mesh$crs,
                                   marks = F,  
                                   mesh = sel_mesh, 
                                   speciesname = "Species",
                                   boundary = in_bound)

RD_spdemodel_PO = INLA::inla.spde2.pcmatern(mesh = sel_mesh, alpha = 3/2,  ### mesh and smoothness parameter
                                            prior.range = c(40, 0.05), ### P(range < 40km) = 0.05 Very small probability the range is smaller than 40 (triangle edge size) (with high likelihood the range is between 20 (res of the model) and infinity i.e. flat prior)
                                            prior.sigma = c(1, 0.05)) ### P(sigma > 2) = 0.05 Very small probability that sigma is larger than 0.2. Sigma is most likely between 0 and 2

RD_spdemodel_PA = INLA::inla.spde2.pcmatern(mesh = sel_mesh, alpha = 3/2,  ### mesh and smoothness parameter
                                            prior.range = c(40, 0.05), ### P(range < 40km) = 0.05 Very small probability the range is smaller than 40 (triangle edge size) (with high likelihood the range is between 20 (res of the model) and infinity i.e. flat prior)
                                            prior.sigma = c(1, 0.05)) ### P(sigma > 2) = 0.05 Very small probability that sigma is larger than 0.2. Sigma is most likely between 0 and 2

# RD_spdemodel_shared = INLA::inla.spde2.pcmatern(mesh = mesh1, alpha = 3/2,  ### mesh and smoothness parameter
#                                                 prior.range = c(100, 0.5), ### P(range < 40km) = 0.05 Very small probability the range is smaller than 40 (triangle edge size) (with high likelihood the range is between 20 (res of the model) and infinity i.e. flat prior) 
#                                                 prior.sigma = c(5, 0.5)) ### P(sigma > 2) = 0.05 Very small probability that sigma is larger than 0.2. Sigma is most likely between 0 and 2


mdl_RD <- bru_sdm(data = RD_data_for_model,
                  spatialcovariates = covar_scaled,
                  covariatestoinclude = c("elevation", "tree_cover_density", "human_footprint_index", #"landCover",
                                          "forest_distances", "small_woody_features", "slope"),
                  specieseffects = TRUE,
                  pointsintercept = TRUE,
                  marksintercept = FALSE,
                  sharedspatial = F,
                  spdemodel =  list(PO_data_sel = RD_spdemodel_PO, 
                                    PA_data_sel = RD_spdemodel_PA),
                  pointsspatial = TRUE,
                  marksspatial = FALSE,
                  spatialdatasets = NULL,
                  tolerance = NULL)

summary(mdl_RD)

# saveRDS(mdl_RD, file = "mdl1_RD_def.RDS")
# mdl_RD <- readRDS("mdl1_RD_def.RDS")

fixed.effects <- mdl_RD$summary.fixed
fixed.effects$variable <- row.names(fixed.effects)

fixed.effects <- fixed.effects %>% 
        filter(variable %!in% c("PO_data_sel_intercept", "PA_data_sel_intercept"))

names(fixed.effects) <- c("mean", "sd", "lower", "median", "higher", "mode", "kld", "variable")

# fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
#                             "Human footprint index", "Land cover")

fixed.effects$variable <- c("Tree cover density", "Elevation", "Slope",
                            "Human footprint index", "Forest distance", "Small Woody Features")


fixed.effects %>% 
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


RD_pred_lin <- predict(mdl_RD, 
                       data = pix, 
                       formula = ~(tree_cover_density + elevation + slope + human_footprint_index +
                                           forest_distances + small_woody_features +
                                           PA_data_sel_spde +
                                           PO_data_sel_spde +
                                           PO_data_sel_intercept + PA_data_sel_intercept))


plot(RD_pred_lin['mean'], main = "Mean of prediction (linear scale) Red Deer")
plot(RD_pred_lin['sd'], main = "SD of prediction (linear scale) Red Deer")
saveRDS(RD_pred_lin, file = "RD_pred_lin_def.RDS")

RD_pred_resp <- predict(mdl_RD, 
                        data = pix,
                        formula = ~exp(tree_cover_density + elevation + slope + human_footprint_index + 
                                               forest_distances + small_woody_features +
                                               PO_data_sel_intercept + PA_data_sel_intercept +
                                               PA_data_sel_spde +
                                               PO_data_sel_spde))

plot(RD_pred_resp['mean'], main = "Mean of prediction (response scale) Red Deer")
plot(RD_pred_resp['sd'], main = "SD of prediction (response scale) Red Deer")
saveRDS(RD_pred_resp, file = "RD_pred_resp_def.RDS")

RD_covars <- predict(mdl_RD, 
                     data = pix, 
                     formula = ~exp(tree_cover_density + elevation + slope + human_footprint_index + 
                                         forest_distances + small_woody_features +
                                         # shared_spatial 
                                         PO_data_sel_intercept + PA_data_sel_intercept))

plot(RD_covars['mean'], main = "Mean of the covariate effect (response scale)")
plot(RD_covars['sd'], main = "SD of the covariate effect (response scale)")

RD_hfi <- predict(mdl_RD, 
                     data = pix, 
                     formula = ~exp(human_footprint_index))

plot(RD_hfi['mean'], main = "Mean of the human footprint index effect (response scale)")
plot(RD_hfi['sd'], main = "SD of the human footprint index effect (response scale)")

# RD_spde_shared <- predict(mdl_RD, 
#                           data = pix, 
#                           formula = ~(shared_spatial))
# 
# plot(RD_spde_shared["mean"], main = "Shared spatial effect")


RD_spde_PO <- predict(mdl_RD, 
                   data = pix, 
                   formula = ~(PO_data_sel_spde))

plot(RD_spde_PO["mean"], main = "Mean of spatial effect PO")
plot(RD_spde_PO["sd"], main = "Sd of spatial effect PO")
saveRDS(RD_spde_PO, file = "RD_spde_PO_def.RDS")


RD_spde_PA <- predict(mdl_RD, 
                      data = pix, 
                      formula = ~(PA_data_sel_spde))

plot(RD_spde_PA["mean"], main = "Mean of spatial effect PA")
plot(RD_spde_PA["sd"], main = "Sd of spatial effect PA")
saveRDS(RD_spde_PA, file = "RD_spde_PA_def.RDS")

ips <- ipoints(in_bound, mesh1)
mdl1_Ab <- predict(mdl_RD, 
                   ips,
                   formula = ~sum(weight*exp((tree_cover_density + elevation + slope + human_footprint_index + 
                                                      forest_distances + small_woody_features +
                                                      PA_data_sel_spde +
                                                      PO_data_sel_spde))))


loo_PO <- leave_one_out(mdl_RD, dataset = "PO_data_sel")
loo_PA <- leave_one_out(mdl_RD, dataset = "PA_data_sel")

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
        filter(model != "Count") %>% 
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

