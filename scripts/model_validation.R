### VALIDATION WITH TIM BURKITT DENSITY DATA ###

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
RD_dens_pred <- dens_sp %>% 
  mutate(prediction = raster::extract(pred_resp$mean, dens_sp), 
         predictionLog = log(raster::extract(pred_resp$mean, dens_sp))) %>% 
  filter(Species == "RedDeer") %>% 
  group_by(Year) %>% 
  mutate(count = n()) %>% 
  ungroup() %>% 
  filter(count > 10) # only retain years with more than 10 sites sampled
  

# plot prediction
my.formula <- y ~ x

(RD_plot <- RD_dens_pred %>% 
    # filter(prediction > 0.1) %>%
    # filter(Deer.Density > 0) %>%
    ggplot(aes(x = Deer.Density, y = prediction)) +
    geom_point() + 
    geom_smooth(method='lm', formula= my.formula) + 
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE) +
    theme_bw() + 
    facet_wrap(~Year, scales = "free") + 
    labs(x = "Tim Burkitt density", y = "Model prediction", title = "Red deer prediction validation"))
  

## Fallow deer ##

# load prediction 

FD_pred_resp <- readRDS("server_outputs/FD_pred_resp_def.RDS")
pred_resp <- stack(raster(FD_pred_resp['mean']), raster(FD_pred_resp['sd']))

# extract values
FD_dens_pred <- dens_sp %>% 
  mutate(prediction = raster::extract(pred_resp$mean, dens_sp), 
         predictionLog = log(raster::extract(pred_resp$mean, dens_sp))) %>% 
  filter(Species == "FallowDeer") %>% 
  group_by(Year) %>% 
  mutate(count = n()) %>% 
  ungroup() %>% 
  filter(count > 10) # only retain years with more than 10 sites sampled


# plot prediction
my.formula <- y ~ x

(FD_plot <- FD_dens_pred %>% 
  # filter(prediction > 0.1) %>%
  # filter(Dens.avg > 0) %>%
  ggplot(aes(x = Deer.Density, y = prediction)) +
  geom_point() + 
  geom_smooth(method='lm', formula= my.formula) + 
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  theme_bw() + 
  facet_wrap(~Year, scales = "free") + 
  labs(x = "Tim Burkitt density", y = "Model prediction", title = "Fallow deer prediction validation"))

## Sika deer ##

# load prediction 

SD_pred_resp <- readRDS("server_outputs/SD_pred_resp_def.RDS")
pred_resp <- stack(raster(SD_pred_resp['mean']), raster(SD_pred_resp['sd']))

# extract values
SD_dens_pred <- dens_sp %>% 
  mutate(prediction = raster::extract(pred_resp$mean, dens_sp), 
         predictionLog = log(raster::extract(pred_resp$mean, dens_sp))) %>% 
  filter(Species == "SikaDeer") %>% 
  group_by(Year) %>% 
  mutate(count = n()) %>% 
  ungroup() %>% 
  filter(count > 10) # only retain years with more than 10 sites sampled



# plot prediction
my.formula <- y ~ x

(SD_plot <- SD_dens_pred %>% 
    # filter(prediction > 0.1) %>%
    # filter(Dens.avg > 0) %>%
    ggplot(aes(x = Deer.Density, y = prediction)) +
    geom_point() + 
    geom_smooth(method='lm', formula= my.formula) + 
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE) +
    theme_bw() + 
    facet_wrap(~Year, scales = "free") + 
    labs(x = "Tim Burkitt density", y = "Model prediction", title = "Sika deer prediction validation"))



cowplot::plot_grid(RD_plot, FD_plot, SD_plot, nrow = 2)


### VALIDATION WITH COILLTE CATEGORICAL DATA ###

# load Coillte data 

coillte_data <- readRDS("data/PA_data_all.RDS")

coillte_data <- coillte_data %>% 
  select(Red, Fallow, Sika, Year) %>% 
  st_transform(projKM) %>% 
  st_centroid()



## Red deer ##

pred_resp <- stack(raster(RD_pred_resp['mean']), raster(RD_pred_resp['sd']))

# extract values
RD_dens_pred <- coillte_data %>% 
  mutate(prediction = raster::extract(pred_resp$mean, coillte_data))


# plot prediction
my.formula <- y ~ x

(RD_plot <- RD_dens_pred %>% 
    mutate(Red = factor(Red, levels = c("No", "Low", "Moderate", "High")), 
           Red.num = as.numeric(Red)) %>% 
    ggplot(aes(x = Red.num, y = prediction)) +
    geom_boxplot(aes(x = Red, y = prediction)) +
    # geom_point() +
    geom_smooth(method='lm', formula= my.formula) +
    stat_poly_eq(formula = my.formula,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw() + 
    facet_wrap(~Year, scales = "free") + 
    labs(x = "Coillte density estimates", y = "Model prediction", title = "Red deer prediction validation"))



## Fallow deer ##

pred_resp <- stack(raster(FD_pred_resp['mean']), raster(FD_pred_resp['sd']))

# extract values
FD_dens_pred <- coillte_data %>% 
  mutate(prediction = raster::extract(pred_resp$mean, coillte_data))


# plot prediction
my.formula <- y ~ x

(FD_plot <- FD_dens_pred %>% 
    mutate(Fallow = factor(Fallow, levels = c("No", "Low", "Moderate", "High")), 
           Fallow.num = as.numeric(Fallow)) %>% 
    ggplot(aes(x = Fallow.num, y = prediction)) +
    geom_boxplot(aes(x = Fallow, y = prediction)) +
    # geom_point() +
    geom_smooth(method='lm', formula= my.formula) +
    stat_poly_eq(formula = my.formula,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw() + 
    facet_wrap(~Year, scales = "free") + 
    labs(x = "Coillte density estimates", y = "Model prediction", title = "Fallow deer prediction validation"))

## Sika deer ##

pred_resp <- stack(raster(SD_pred_resp['mean']), raster(SD_pred_resp['sd']))

# extract values
SD_dens_pred <- coillte_data %>% 
  mutate(prediction = raster::extract(pred_resp$mean, coillte_data))


# plot prediction
my.formula <- y ~ x

(SD_plot <- SD_dens_pred %>% 
    mutate(Sika = factor(Sika, levels = c("No", "Low", "Moderate", "High")), 
           Sika.num = as.numeric(Sika)) %>% 
    ggplot(aes(x = Sika.num, y = prediction)) +
    geom_boxplot(aes(x = Sika, y = prediction)) +
    # geom_point() +
    geom_smooth(method='lm', formula= my.formula) +
    stat_poly_eq(formula = my.formula,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE) +
    theme_bw() + 
    facet_wrap(~Year, scales = "free") + 
    labs(x = "Coillte density estimates", y = "Model prediction", title = "Sika deer prediction validation"))

