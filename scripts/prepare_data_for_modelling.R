##---------------##
#### Deer data ####
##---------------##

pre_data <- read.csv("data/all_data.csv", row.names = NULL)
pre_data_NI <- read.csv("data/all_NI_data.csv", row.names = NULL)

pre_data_NI <- pre_data_NI %>% 
  rename(Longitude = X, Latitude = Y)

ireland <- st_read("data/ireland_ITM.shp") 

# all_data <- pre_data  
all_data <- bind_rows(pre_data, pre_data_NI)

sel_data <- all_data %>% 
  dplyr::select(County, Year, Species, Deer.Presence, Y = Latitude, X = Longitude, Source) %>% 
  st_as_sf(coords = c("X", "Y")) %>% 
  st_set_crs(st_crs(ireland))  # IRENET 

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
env_data <- stack("large_env_data/covar_subset_ITM.grd")
plot(env_data)

# Previsualisation with Presence/Absence data
PA_data_sf <- PA_data %>% 
  st_as_sf(coords = c("X", "Y"), crs = st_crs(ireland)) 

PA_covars <- raster::extract(env_data, PA_data_sf)  

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
  mutate(landCover2 = recode(landCover, "1" = "Built area", "2" = "Saltwater related", "3" = "No vegetation", 
                             "4" = "Freshwater related", "5" = "Other vegetation", "6" = "Agricultura land", 
                             "7" = "Pastures", "8" = "Grassland", "9" = "Transitional", "10" = "Coniferous forest", 
                             "11" = "Mixed forest", "12" = "Broadleaf forest")) %>%
  ggplot + 
  geom_bar(aes(fill = landCover2, x = factor(PA)), position = "fill") + 
  scale_fill_viridis_d() + 
  labs(y = "Proportion", x = "", fill = "Land Cover") +
  theme_bw()


plot_grid(tcd_plot, ele_plot, slo_plot, hfi_plot, lcv_plot)
ggsave("fSikaDeer_covar_eval.png", scale = 2)



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

RL <- list(env_data$landCover, env_data$tree_cover_density, env_data$elevation, env_data$slope, env_data$human_footprint_index, PO_raster2)
# devtools::source_url("https://github.com/VirginiaMorera/Useful-little-functions/blob/master/list_to_stack.R?raw=TRUE")

new_stack <- list_to_stack(RL, new_res = raster::res(env_data), dest_crs = env_data@crs, turn_0_to_NA = F)

plot(new_stack)
names(new_stack)[6] <- "Sika_deer_per_cell"

cor <-layerStats(new_stack, 'pearson', na.rm = T)
ggcorrplot(cor$`pearson correlation coefficient`, method = "circle", show.diag = F, type = "upper", lab = T) 
ggsave("sikaDeerCorr.png")



##--------------------------------##
#### Package data for modelling ####
##--------------------------------##

# Functions from the package inlabruSDMs by the wonderful Philip Mostert
in_bound <- readRDS("data/inner_boundary.RDS")
mesh0 <- readRDS("data/mesh.RDS")
PO_data <- read.csv("data/PO_webSurvey_data_RD.csv", row.names = NULL)
PA_data <- read.csv("data/PA_data_RD.csv", row.names = NULL)

RD_data_for_model <- organize_data(PA_data, PO_data, poresp = "PO", paresp = "PA", coords = c("X", "Y"), proj = CRS("+init=epsg:2157"),
                                marks = F,  mesh = mesh0, boundary = in_bound)

saveRDS(RD_data_for_model, file = "data/RD_data_for_model.RDS")
