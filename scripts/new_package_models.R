##Script using PointedSDMs
# install.packages("PointedSDMs")
library(PointedSDMs)
library(raster)
library(tidyverse)
library(ggthemes)

`%!in%` = Negate(`%in%`)

in_bound <- readRDS("data/inner_boundary.RDS")
mesh1 <- readRDS("data/meshLarge.RDS")
covar_stack <- stack("large_env_data/covar_subset_KM.gri")
covar_stack <- subset(covar_stack, c(1:4, 6:7))
covar_scaled <- scale(covar_stack)
ireland <- sf::st_read("data/ireland_ITM.shp")

ireland_km  <-   sf::st_transform(ireland, crs = in_bound@proj4string)
ireland_km_sp <- sf::as_Spatial(ireland_km)

PO_data <- read.csv("data/PO_data_all.csv", row.names = NULL)

PO_data_sel <- PO_data %>% 
  filter(Source %in% c("CullReturnsNI", "MammalSurvey", "Bycatch", "NBDC", "CitizenScience", "webSurvey", "Other")) %>% 
  filter(Species %in% c("RedDeer")) %>% 
  filter(Y < 965) %>% 
  select(PO, X, Y) %>% 
  drop_na()

PO_data_sel_sp <- SpatialPointsDataFrame(coords = PO_data_sel[,c(2:3)], 
                                         data = data.frame(PO = PO_data_sel[,1]),
                                         proj4string = mesh1$crs)


PA_data <- read.csv("data/PA_data_all.csv", row.names = NULL)

PA_data_sel <- PA_data %>%
  filter(Species %in% c("RedDeer")) %>% 
  select(PA = Deer.Presence, X, Y) %>% 
  drop_na()

PA_data_sel_sp <- SpatialPointsDataFrame(coords = PA_data_sel[,c(2:3)], 
                                         data = data.frame(PA = PA_data_sel[,1]),
                                         proj4string = mesh1$crs)


PA_data_NI <- PA_data %>%
  filter(Species %in% c("RedDeer")) %>% 
  filter(Source == "BDS_survey") %>% 
  select(PA = Deer.Presence, X, Y) %>% 
  drop_na()

PA_data_ROI <- PA_data %>%
  filter(Species %in% c("RedDeer")) %>% 
  filter(Source != "BDS_survey") %>% 
  select(PA = Deer.Presence, X, Y) %>% 
  drop_na()


pix <- pixels(mesh = mesh1, mask = in_bound)
pix@proj4string <- mesh1$crs

#### Shared field model ####

## set up the model ####
shared_field_model <- PointedSDMs::intModel(PA_data_NI, PA_data_ROI, PO_data_sel,
                                            spatialCovariates = covar_scaled,
                                            Coordinates = c('X', 'Y'),
                                            Mesh = mesh1,
                                            pointsSpatial = 'shared',
                                            responsePA = 'PA',
                                            #responseCounts = 'PO',
                                            speciesName = NULL,
                                            Projection = mesh1$crs)

## specifiy the shared spatial field ####

shared_field_model$specifySpatial(sharedSpatial = TRUE,
                                  alpha = 3/2,
                                  prior.range = c(200, 0.01), 
                                  prior.sigma = c(0.1, 0.01))

## run the model ####

mdl_RD_shared <- runModel(shared_field_model)
summary(mdl_RD_shared)

fixed.effectsRD <- mdl_RD_shared$summary.fixed
fixed.effectsRD$variable <- row.names(fixed.effectsRD)

fixed.effectsRD <- fixed.effectsRD %>% 
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_ROI_intercept", "PA_data_NI_intercept"))

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
#saveRDS(object = mdl_RD_shared, file = 'mdl_shared_field.RDS')


## predict ####
mdl_shared_pred = predict(mdl_RD_shared, mesh = mesh1, mask = in_bound, predictor = TRUE, fun = 'linear')
mdl_shared_spde = predict(mdl_RD_shared, mesh = mesh1, mask = in_bound, 
                          formula = ~shared_spatial,
                          # predictor = TRUE, 
                          fun = 'linear')
#saveRDS(object = mdl_sep_pred, 'mdl_sep_pred.RDS')

par(mfrow = c(1,2))
plot(raster(mdl_shared_pred$predictions["mean"]), main = "Mean, shared spde model")
plot(raster(mdl_shared_pred$predictions["sd"]), main = "SD, shared spde model")
par(mfrow = c(1,1))

plot(mdl_shared_spde$predictions["sd"])


## Separate field model ####

## set up the model 

sep_field_model <- PointedSDMs::intModel(PA_data_NI, PA_data_ROI, PO_data_sel,
                                         spatialCovariates = covar_scaled,
                                         Coordinates = c('X', 'Y'),
                                         Mesh = mesh1,
                                         pointsSpatial = "individual",
                                         responsePA = 'PA',
                                         responseCounts = 'PO',
                                         speciesName = NULL,
                                         Projection = mesh1$crs)

## specify the separate spatial effects ####

sep_field_model$specifySpatial(datasetName = 'PO_data_sel',
                               alpha = 3/2, 
                               prior.range = c(40, 0.05),
                               prior.sigma = c(1, 0.05))

sep_field_model$specifySpatial(datasetName = 'PA_data_NI',
                               alpha = 3/2, 
                               prior.range = c(40, 0.05), 
                               prior.sigma = c(1, 0.05))

sep_field_model$specifySpatial(datasetName = 'PA_data_ROI',
                               alpha = 3/2, 
                               prior.range = c(40, 0.05), 
                               prior.sigma = c(1, 0.05))

## run the model ####
mdl_RD_sep <- runModel(sep_field_model)
#saveRDS(object = mdl_RD_sep, file = 'mdl_sep_fields.RDS')
summary(mdl_RD_sep)

fixed.effectsRD <- mdl_RD_sep$summary.fixed
fixed.effectsRD$variable <- row.names(fixed.effectsRD)

fixed.effectsRD <- fixed.effectsRD %>% 
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_ROI_intercept", "PA_data_NI_intercept"))

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



## predict ####
mdl_sep_pred = predict(mdl_RD_sep, mesh = mesh1, mask = in_bound, predictor = TRUE, fun = 'linear')
#saveRDS(object = mdl_sep_pred, 'mdl_sep_pred.RDS')

par(mfrow = c(1,2))
plot(raster(mdl_sep_pred$predictions["mean"]), main = "Mean, separate spde model")
plot(raster(mdl_sep_pred$predictions["sd"]), main = "SD, separate spde model")
par(mfrow = c(1,1))

## Copy field model ####

## set up the model
copy_field_model <- PointedSDMs::intModel(PA_data_NI, PA_data_ROI, PO_data_sel,
                                          spatialCovariates = covar_scaled,
                                          Coordinates = c('X', 'Y'),
                                          Mesh = mesh1,
                                          pointsSpatial = 'individual',
                                          responsePA = 'PA',
                                          # responseCounts = 'PO',
                                          speciesName = NULL,
                                          Projection = mesh1$crs)

## specify the copied spatial fields ####

copy_field_model$specifySpatial(datasetName = 'PO_data_sel',
                                alpha = 3/2, 
                                prior.range = c(100, 0.5),
                                prior.sigma = c(1, 0.5))

copy_field_model$changeComponents(addComponent = 'PA_data_NI_spatial(main = coordinates, copy = "PO_data_sel_spatial", fixed = FALSE)')
copy_field_model$changeComponents(addComponent = 'PA_data_ROI_spatial(main = coordinates, copy = "PO_data_sel_spatial", fixed = FALSE)')

## run model ####
mdl_RD_copy <- runModel(copy_field_model)

summary(mdl_RD_copy)
fixed.effectsRD <- mdl_RD_copy$summary.fixed
fixed.effectsRD$variable <- row.names(fixed.effectsRD)

fixed.effectsRD <- fixed.effectsRD %>% 
  filter(variable %!in% c("PO_data_sel_intercept", "PA_data_ROI_intercept", "PA_data_NI_intercept"))

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


#saveRDS(object = mdl_RD_copy, file = 'mdl_copy_fields.RDS')

## predict ####
mdl_copy_pred = predict(mdl_RD_copy, mesh = mesh1, mask = in_bound, 
                        # formula = ~PO_data_NI_spatial,
                        predictor = TRUE,
                        fun = 'linear')


mdl_copy_pred_PO = predict(mdl_RD_copy, mesh = mesh1, mask = in_bound, 
                           formula = ~PO_data_sel_spatial,
                           # predictor = TRUE, 
                           fun = 'linear')

mdl_copy_pred_ROI = predict(mdl_RD_copy, mesh = mesh1, mask = in_bound, 
                            formula = ~PA_data_ROI_spatial,
                            # predictor = TRUE, 
                            fun = 'linear')

mdl_copy_pred_NI = predict(mdl_RD_copy, mesh = mesh1, mask = in_bound, 
                           formula = ~PA_data_NI_spatial,
                           # predictor = TRUE, 
                           fun = 'linear')

mdl_copy_pred_covars = predict(mdl_RD_copy, mesh = mesh1, mask = in_bound, 
                        formula = ~tree_cover_density + elevation + slope + human_footprint_index + forest_distances + small_woody_features,
                        # predictor = TRUE,
                        fun = 'linear')

par(mfrow = c(1,3))
plot(raster(mdl_copy_pred$predictions["mean"]), main = "Mean, copied spde model")
plot(raster(mdl_copy_pred$predictions["sd"]), main = "SD, copied spde model")
plot(raster(mdl_copy_pred_covars$predictions["mean"]), main = "covar effect")
par(mfrow = c(1,1))


#saveRDS(object = mdl_sep_pred, 'mdl_sep_pred.RDS')

PO <- raster(mdl_copy_pred_PO$predictions["mean"])
ROI <- raster(mdl_copy_pred_ROI$predictions["mean"])
NI <- raster(mdl_copy_pred_NI$predictions["mean"])

par(mfrow = c(1,3))
plot(PO, main = "PO spatial field")
plot(ROI, main = "ROI PA 'copy' spatial field")
plot(NI, main = "NI PA 'copy' spatial field")
par(mfrow = c(1,1))
