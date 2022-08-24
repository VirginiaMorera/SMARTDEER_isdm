source("scripts/setups.R")

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
  dplyr::select(PO, X, Y) %>% 
  drop_na()

PO_data_sel_sp <- SpatialPointsDataFrame(coords = PO_data_sel[,c(2:3)], 
                                         data = data.frame(PO = PO_data_sel[,1]),
                                         proj4string = mesh1$crs)


PA_data <- read.csv("data/PA_data_all.csv", row.names = NULL)

PA_data_sel <- PA_data %>%
  filter(Species %in% c("RedDeer")) %>% 
  dplyr::select(PA = Deer.Presence, X, Y) %>% 
  drop_na()

PA_data_sel_sp <- SpatialPointsDataFrame(coords = PA_data_sel[,c(2:3)], 
                                         data = data.frame(PA = PA_data_sel[,1]),
                                         proj4string = mesh1$crs)


PA_data_NI <- PA_data %>%
  filter(Species %in% c("RedDeer")) %>% 
  filter(Source == "BDS_survey") %>% 
  dplyr::select(PA = Deer.Presence, X, Y) %>% 
  drop_na()

PA_data_ROI <- PA_data %>%
  filter(Species %in% c("RedDeer")) %>% 
  filter(Source != "BDS_survey") %>% 
  dplyr::select(PA = Deer.Presence, X, Y) %>% 
  drop_na()


pix <- pixels(mesh = mesh1, mask = in_bound)
pix@proj4string <- mesh1$crs


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
                                prior.range = c(100, 0.1),
                                prior.sigma = c(0.5, 0.1))

copy_field_model$changeComponents(addComponent = 'PA_data_NI_spatial(main = coordinates, copy = "PO_data_sel_spatial", fixed = FALSE)')
copy_field_model$changeComponents(addComponent = 'PA_data_ROI_spatial(main = coordinates, copy = "PO_data_sel_spatial", fixed = FALSE)')

## run model ####
mdl_RD_copy <- runModel(copy_field_model, options = list(num.threads = 1))
beep(2)


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

saveRDS(object = mdl_RD_copy, file = 'server_outputs/RedDeer_mdl_copy.RDS')
mdl_RD_copy <- readRDS('server_outputs/RedDeer_mdl_copy.RDS')
## predict ####
mdl_copy_pred = predict(mdl_RD_copy, mesh = mesh1, mask = in_bound, 
                        # formula = ~PO_data_NI_spatial,
                        predictor = TRUE,
                        fun = 'linear')
beep(9)

mdl_copy_pred_PO = predict(mdl_RD_copy, mesh = mesh1, mask = in_bound, 
                           formula = ~PO_data_sel_spatial,
                           # predictor = TRUE, 
                           fun = 'linear')

mdl_copy_pred_covars = predict(mdl_RD_copy, mesh = mesh1, mask = in_bound, 
                        formula = ~tree_cover_density + elevation + slope + human_footprint_index + forest_distances + small_woody_features,
                        # predictor = TRUE,
                        fun = 'linear')

par(mfrow = c(2,2))
plot(raster(mdl_copy_pred$predictions["mean"]), main = "Mean, copied spde model")
plot(raster(mdl_copy_pred$predictions["sd"]), main = "SD, copied spde model")
plot(raster(mdl_copy_pred_covars$predictions["mean"]), main = "Covariate effect")
plot(raster(mdl_copy_pred_PO$predictions["mean"]), main = "Spatial field")
par(mfrow = c(1,1))

saveRDS(mdl_copy_pred, file = "server_outputs/RedDeer_prediction.RDS")
saveRDS(mdl_copy_pred_PO, file = "server_outputs/RedDeer_spatial_field.RDS")
saveRDS(mdl_copy_pred_covars, file = "server_outputs/RedDeer_covariate_effect.RDS")

#saveRDS(object = mdl_sep_pred, 'mdl_sep_pred.RDS')

