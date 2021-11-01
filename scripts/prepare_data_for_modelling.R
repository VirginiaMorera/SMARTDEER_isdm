##---------------##
#### Deer data ####
##---------------##

pre_data <- read.csv("data/all_data.csv", row.names = NULL)
pre_data_NI <- read.csv("data/all_NI_data.csv", row.names = NULL)

pre_data_NI <- pre_data_NI %>% 
  rename(Longitude = X, Latitude = Y)
ireland <- st_read("data/ireland_ITM.shp") %>% 
  st_transform(st_crs("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=0.99982 +x_0=600000 +y_0=750000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"))  # IRENET in km

# all_data <- pre_data  
all_data <- bind_rows(pre_data, pre_data_NI)

sel_data <- all_data %>% 
  dplyr::select(County, Year, Species, Deer.Presence, Y = Latitude, X = Longitude, Source) %>% 
  st_as_sf(coords = c("X", "Y")) %>% 
  st_set_crs(st_crs("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=0.99982 +x_0=600000 +y_0=750000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")) %>%  # IRENET in m
  st_transform(ITM_km)  # IRENET in km

PA_data <- sel_data %>% 
  filter(Source %in% c("Coillte_density", "Coillte_desk")) %>% 
  filter(Species == "RedDeer") %>% 
  dplyr::select(-Source)  %>% 
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                Y = sf::st_coordinates(.)[,2], 
                Deer.Presence = if_else(Deer.Presence == "Yes", 1, 0)) %>% 
  rename(PA = Deer.Presence) %>% 
  st_set_geometry(NULL)

write.csv(PA_data, file = "data/PA_data_RD.csv", row.names = F)

PO_data <- sel_data %>% 
  # filter(Source %in% c("NBDC", "webSurvey")) %>%
  filter(Source == "webSurvey") %>%
  filter(Species == "RedDeer") %>% 
  dplyr::select(-Source) %>% 
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                Y = sf::st_coordinates(.)[,2], 
                Deer.Presence = if_else(Deer.Presence == "Yes", 1, 0)) %>% 
  rename(PO = Deer.Presence) %>% 
  st_set_geometry(NULL)

write.csv(PO_data, file = "data/PO_webSurvey_data_RD.csv", row.names = F)


##------------------------##
#### Environmental data ####
##------------------------##

# We're only going to load the data that, a priori, we want to enter in the model
env_data <- stack("large_env_data/covar_subset.grd")
plot(env_data)


# Previsualisation with Presence/Absence data
PA_data_sf <- PA_data %>% 
  st_as_sf(coords = c("X", "Y"), crs = ITM_km) %>% 
  st_transform(env_data@crs)

PA_covars <- extract(env_data, PA_data_sf)  

PA_covars <- as.data.frame(cbind(PA = PA_data$PA, PA_covars))

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
  mutate(landCover2 = recode(landCover, "1" = "Built area", "4" = "Freshwater related", "5" = "Other vegetation",
                             "6" = "Agricultura land", "7" = "Pastures", "8" = "Grassland", "9" = "Transitional",
                             "10" = "Coniferous forest", "11" = "Mixed forest", "12" = "Broadleaf forest")) %>%
  ggplot + 
  geom_bar(aes(fill = landCover2, x = factor(PA)), position = "fill") + 
  scale_fill_tableau() + 
  labs(y = "Proportion", x = "", fill = "Land Cover") +
  theme_bw()


plot_grid(tcd_plot, ele_plot, slo_plot, hfi_plot, lcv_plot)
ggsave("covar_eval.png", scale = 2)

PO_data_sf <- PO_data %>% 
  st_as_sf(coords = c("X", "Y"), crs = ITM_km) %>% 
  st_transform(env_data@crs)

PO_data_sp <- SpatialPoints(st_coordinates(PO_data_sf), proj4string = env_data@crs)

PO_raster <- rasterize(PO_data_sp, env_data$landCover, fun='count', na.rm = T)

plot(PO_raster)

env_data2 <- stack(env_data, PO_raster)

names(env_data2)[6] <- "Deer_per_cell"

cor <-layerStats(env_data2, 'pearson', na.rm = T)
ggcorrplot(cor$`pearson correlation coefficient`, method = "circle", show.diag = F, type = "upper", lab = T) 
