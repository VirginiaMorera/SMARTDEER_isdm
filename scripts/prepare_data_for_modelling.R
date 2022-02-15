##---------------##
#### Deer data ####
##---------------##

all_data <- read.csv("data/all_data.csv", row.names = NULL)

ireland <- st_read("data/ireland_ITM.shp") 

ireland_km <- ireland %>% 
  st_transform(projKM)

sel_data <- all_data %>% 
  dplyr::select(County, Year, Species, Deer.Presence, Deer.Count, Y = Latitude, X = Longitude, Source, Type) %>% 
  filter(X > 421874) %>% 
  st_as_sf(coords = c("X", "Y")) %>% 
  st_set_crs(st_crs(ireland)) %>% 
  st_transform(projKM) # IRENET in km

sel_data %>% 
  filter(Type == "PO") %>%
  filter(Species == "FallowDeer") %>% 
  filter(Source != "CameraTraps") %>%
  ggplot + 
  geom_sf(aes(col = Source), alpha = 0.5) + 
  geom_sf(data = ireland_km, fill = NA, col = "darkgray") + 
  # facet_wrap(~Type) + 
  theme_bw()

PA_data <- sel_data %>% 
  filter(Type  == "PA") %>% 
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                Y = sf::st_coordinates(.)[,2], 
                Deer.Presence = if_else(Deer.Presence == "Yes", 1, 0)) %>% 
  rename(PA = Deer.Presence) %>% 
  dplyr::select(-Deer.Count) %>% 
  st_set_geometry(NULL)

write.csv(PA_data, file = "data/PA_data_all.csv", row.names = F)

PO_data <- sel_data %>% 
  filter(Type == "PO") %>%
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                Y = sf::st_coordinates(.)[,2], 
                Deer.Presence = if_else(Deer.Presence == "Yes", 1, 0)) %>% 
  rename(PO = Deer.Presence, 
         Count = Deer.Count) %>% 
  st_set_geometry(NULL)

write.csv(PO_data, file = "data/PO_data_all.csv", row.names = F)


##------------------------##
#### Environmental data ####
##------------------------##

# We're only going to load the data that, a priori, we want to enter in the model
env_data <- stack("large_env_data/covar_subset_ITM.grd")
plot(env_data)

# Previsualisation with Presence/Absence data
PA_data_sf <- PA_data %>% 
  st_as_sf(coords = c("X", "Y"), crs = st_crs(ireland)) 

PA_covars <- raster::extract(env_data, PA_data_sf)  

PA_covars <- as.data.frame(cbind(PA = PA_data$PA, PA_covars)) 
PA_covars <- PA_covars %>% 
  mutate(landCover2 = recode(landCover, "1" = "Built area", "2" = "Saltwater related", "3" = "No vegetation", 
                             "4" = "Freshwater related", "5" = "Other vegetation", "6" = "Agricultural land", 
                             "7" = "Pastures", "8" = "Grassland", "9" = "Transitional", "10" = "Coniferous forest", 
                             "11" = "Mixed forest", "12" = "Broadleaf forest")) %>%
  mutate(landCover2 = factor(landCover2, levels = c("Built area", "Saltwater related", "No vegetation", 
                                                    "Freshwater related", "Other vegetation", "Agricultural land", 
                                                    "Pastures", "Grassland", "Transitional", "Coniferous forest", 
                                                    "Mixed forest", "Broadleaf forest")))
  

tcd_plot <- ggplot(PA_covars) + 
  geom_jitter(aes(x = factor(PA), y = tree_cover_density, col = factor(PA)), 
              alpha = 0.1, width = 0.1) + 
  geom_boxplot(aes(x = factor(PA), y = tree_cover_density, fill = factor(PA)), 
               alpha = 0.5, outlier.shape = NA) + 
  scale_fill_tableau() +
  scale_color_tableau() +
  labs(x = "", y = "Tree Cover Density", fill = "PA") +
  theme_bw() + 
  guides(color = "none") +
  # facet_wrap(~Classification, scales = "free_x") +
  # theme(legend.position = "None") + 
  NULL

ele_plot <- ggplot(PA_covars) + 
  geom_jitter(aes(x = factor(PA), y = elevation, col = factor(PA)), 
              alpha = 0.1, width = 0.1) + 
  geom_boxplot(aes(x = factor(PA), y = elevation, fill = factor(PA)), 
               alpha = 0.5, outlier.shape = NA) + 
  scale_fill_tableau() +
  scale_color_tableau() +
  labs(x = "", y = "elevation", fill = "PA") +
  theme_bw() + 
  guides(color = "none") +
  # facet_wrap(~Classification, scales = "free_x") +
  # theme(legend.position = "None") + 
  NULL

slo_plot <- ggplot(PA_covars) + 
  geom_jitter(aes(x = factor(PA), y = slope, col = factor(PA)), 
              alpha = 0.1, width = 0.1) + 
  geom_boxplot(aes(x = factor(PA), y = slope, fill = factor(PA)), 
               alpha = 0.5, outlier.shape = NA) + 
  scale_fill_tableau() +
  scale_color_tableau() +
  labs(x = "", y = "Slope", fill = "PA") +
  theme_bw() + 
  guides(color = "none") +
  # facet_wrap(~Classification, scales = "free_x") +
  # theme(legend.position = "None") + 
  NULL

hfi_plot <- ggplot(PA_covars) + 
  geom_jitter(aes(x = factor(PA), y = human_footprint_index, col = factor(PA)), 
              alpha = 0.1, width = 0.1) + 
  geom_boxplot(aes(x = factor(PA), y = human_footprint_index, fill = factor(PA)), 
               alpha = 0.5, outlier.shape = NA) + 
  scale_fill_tableau() +
  scale_color_tableau() +
  labs(x = "", y = "Human footprint index", fill = "PA") +
  theme_bw() + 
  guides(color = "none") +
  # facet_wrap(~Classification, scales = "free_x") +
  # theme(legend.position = "None") + 
  NULL

lcv_plot <- PA_covars %>%
  ggplot + 
  geom_bar(aes(fill = landCover2, x = factor(PA)), position = "fill") + 
  scale_fill_viridis_d() + 
  labs(y = "Proportion", x = "", fill = "Land Cover") +
  theme_bw()


plot_grid(tcd_plot, ele_plot, slo_plot, hfi_plot, lcv_plot)
ggsave("outputs/fallowDeer_covar_eval.png", scale = 2)


# Previsualisation with PO data
PO_data_sf <- PO_data %>% 
  st_as_sf(coords = c("X", "Y"), crs = st_crs(ireland)) 

PO_win <- as.owin(st_bbox(ireland))
PO_ppp <- ppp(x = st_coordinates(PO_data_sf)[,1], y = st_coordinates(PO_data_sf)[,2], 
              window = PO_win)
PO_dens <- density(PO_ppp, sigma = bw.diggle, eps = c(500, 500))
PO_raster <- raster(PO_dens)
PO_raster@crs <- CRS("+init=epsg:2157")
PO_raster <- raster::crop(PO_raster, y = extent(env_data))
PO_raster2 <- raster::mask(PO_raster, ireland)
plot(PO_raster2)
PO_raster2[PO_raster2 < 0.0000001] <- NA

RL <- list(env_data$landCover, env_data$tree_cover_density, env_data$elevation, env_data$slope, env_data$human_footprint_index, PO_raster2)
# devtools::source_url("https://github.com/VirginiaMorera/Useful-little-functions/blob/master/list_to_stack.R?raw=TRUE")

new_stack <- list_to_stack(RL, new_res = raster::res(env_data), dest_crs = env_data@crs, turn_0_to_NA = T)

plot(new_stack)
names(new_stack)[6] <- "Fallow_deer_per_cell"
cor <-layerStats(new_stack, 'pearson', na.rm = T)
ggcorrplot(cor$`pearson correlation coefficient`, method = "circle", show.diag = F, type = "upper", lab = T) 
ggsave("outputs/fallowDeerCorr.png")

