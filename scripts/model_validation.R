# load data
ireland <- st_read("data/ireland_ITM.shp")
dens_data <- read.csv("data/density_latlon_coillte_dens_surveys.csv", row.names = NULL)

ireland_proj <- st_transform(ireland, projKM)

dens_sp <- dens_data %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(ireland)) %>% 
  st_transform(projKM)

dens_avg <- dens_sp %>% 
  group_by(Site_id, Species) %>% 
  summarise(Dens.avg = mean(Deer.Density, na.rm = T)) %>% 
  ungroup()

## Red deer ##

# load prediction 

RD_pred_resp <- readRDS("server_outputs/RD_pred_resp_def.RDS")
pred_resp <- stack(raster(RD_pred_resp['mean']), raster(RD_pred_resp['sd']))

# extract values
RD_dens_pred <- dens_avg %>% 
  mutate(prediction = raster::extract(pred_resp$mean, dens_avg), 
         predictionLog = log(raster::extract(pred_resp$mean, dens_avg))) %>% 
  filter(Species == "RedDeer") 
  

# plot prediction
my.formula <- y ~ x

RD_plot <- RD_dens_pred %>% 
  # filter(prediction > 0.1) %>%
  # filter(Dens.avg > 0) %>%
  ggplot(aes(x = Dens.avg, y = prediction)) +
  geom_point() + 
  geom_smooth(method='lm', formula= my.formula) + 
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  theme_bw() + 
  labs(x = "Tim Burkitt density", y = "Model prediction", title = "Red deer prediction validation")


## Fallow deer ##

# load prediction 

FD_pred_resp <- readRDS("server_outputs/FD_pred_resp_def.RDS")
pred_resp <- stack(raster(FD_pred_resp['mean']), raster(FD_pred_resp['sd']))

# extract values
FD_dens_pred <- dens_avg %>% 
  mutate(prediction = raster::extract(pred_resp$mean, dens_avg), 
         predictionLog = log(raster::extract(pred_resp$mean, dens_avg))) %>% 
  filter(Species == "FallowDeer") 


# plot prediction
my.formula <- y ~ x

FD_plot <- FD_dens_pred %>% 
  # filter(prediction > 0.1) %>%
  # filter(Dens.avg > 0) %>%
  ggplot(aes(x = Dens.avg, y = prediction)) +
  geom_point() + 
  geom_smooth(method='lm', formula= my.formula) + 
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  theme_bw() + 
  labs(x = "Tim Burkitt density", y = "Model prediction", title = "Fallow deer prediction validation")

## Sika deer ##

# load prediction 

SD_pred_resp <- readRDS("server_outputs/SD_pred_resp_def.RDS")
pred_resp <- stack(raster(SD_pred_resp['mean']), raster(SD_pred_resp['sd']))

# extract values
SD_dens_pred <- dens_avg %>% 
  mutate(prediction = raster::extract(pred_resp$mean, dens_avg), 
         predictionLog = log(raster::extract(pred_resp$mean, dens_avg))) %>% 
  filter(Species == "SikaDeer") 


# plot prediction
my.formula <- y ~ x

SD_plot <- SD_dens_pred %>% 
  # filter(prediction > 0.1) %>%
  # filter(Dens.avg > 0) %>%
  ggplot(aes(x = Dens.avg, y = prediction)) +
  geom_point() + 
  geom_smooth(method='lm', formula= my.formula) + 
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  theme_bw() + 
  labs(x = "Tim Burkitt density", y = "Model prediction", title = "Sika deer prediction validation")



cowplot::plot_grid(RD_plot, FD_plot, SD_plot, nrow = 2)


