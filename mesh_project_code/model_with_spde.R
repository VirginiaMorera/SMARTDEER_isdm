rm(list = ls())
source("scripts/setups.R")

# 1. Load data ####

ireland <- st_read("data/ireland_ITM.shp") %>% 
  st_transform(projKM)

PO_data <- read.csv("data/PO_data_all.csv", row.names = NULL)

PO_data_sel <- PO_data %>% 
  filter(Species %in% c("FallowDeer")) %>% 
  filter(Y < 965) %>% 
  dplyr::select(PO, X, Y) %>% 
  drop_na() %>% 
  st_as_sf(coords = c("X", "Y")) %>% 
  st_set_crs(projKM)

covar_stack <- rast("large_env_data/covar_subset_KM.grd")
covar_stack <- subset(covar_stack, 
                      c("landCover", "forest_distances", "human_footprint_index"),
                      negate = T) # remove factor var
names(covar_stack) <- c("Tree cover density", "Elevation", "Slope", "Small woody features")
covar_scaled <- scale(covar_stack)
covar_scaled <- trim(covar_scaled)
plot(covar_scaled)

covar_scaled <- project(covar_scaled, projKM)

(pp <- ggplot() + 
    geom_spatraster(data = covar_scaled, aes(fill = Elevation)) +
    scale_fill_viridis_c(na.value = NA) + 
    geom_sf(data = ireland, fill = NA) + 
    geom_sf(data = PO_data_sel) + 
    theme_bw())


# 2. Create mesh ####

boundary <- ireland %>% 
  summarise() %>%
  st_simplify(dTolerance = 10) %>% 
  st_buffer(dist = 20)

boundary2 <- boundary %>% 
  st_buffer(dist = 200)

mesh <- fm_mesh_2d_inla(boundary = list(boundary, boundary2),
                        max.edge = c(10,100), cutoff = 10, crs = projKM)


ipoints <-  fm_int(fmesher:::fm_subdivide(mesh, 2))

## 9.1 Set up SPDE ####

pcmatern <- inla.spde2.pcmatern(mesh, 
                                prior.range = c(120, 0.01), #in km                                
                                prior.sigma = c(1, 0.01)) 
## 9.2 Formula ####
form1 <- geometry ~  Intercept(1)  +
  Eff.elevation(covar_scaled$Elevation, model = "linear") +
  # Eff.slope(covar_scaled$Slope, model = "linear",
  #           mean.linear = -1, prec.linear = 1) +
  # Eff.tree_cover(covar_scaled$`Tree cover density`, model = "linear", 
  #                mean.linear = 2, prec.linear = 1) + 
  # Eff.swf(covar_scaled$`Small woody features`, model = "linear", 
          # mean.linear = 1, prec.linear = 1) +
  Eff.spde(geometry, model = pcmatern) +
  NULL

## 9.3 Run model and summary ####

m1 <- lgcp(components = form1,
           data = PO_data_sel,
           # samplers = boundary,
           ips = ipoints)

summary(m1)


df2 <- fm_pixels(mesh,
                 dims = c(72*10, 81*10),
                 mask = boundary,
                 format = "sf")

pred <- predict(
  object = m1, 
  newdata = df2, 
  samples = 100,
  formula = ~ Intercept +
    Eff.spde +
    Eff.elevation)
    # Eff.slope +
    # Eff.tree_cover +
    # Eff.swf +

spde <- predict(
  object = m1, 
  newdata = df2, 
  samples = 100,
  formula = ~ Eff.spde )



ggplot() + 
  gg(data = pred, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = ireland, fill = NA) +
  # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
  labs(x = "", y = "", fill = "Median",
       title = "Red deer distribution (linear scale)") +
  theme_bw() +
  scale_fill_viridis_c() +
  NULL +
  
ggplot() + 
  gg(data = spde, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = ireland, fill = NA) +
  geom_sf(data = PO_data_sel, alpha = 0.5, size = 1) +
  labs(x = "", y = "", fill = "Median",
       title = "Spatial Random Effect") +
  theme_bw() +
  scale_fill_distiller(palette = 'RdBu') + 
  NULL +
  
plot_layout(ncol = 2)




(Lambda <- predict(
  m1,
  fm_int(mesh, boundary),
  ~ sum(weight * exp(Intercept +
                       Eff.spde +
                       # Eff.slope +
                       # Eff.tree_cover +
                       # Eff.swf + 
                       Eff.elevation))
))



sample <- predict(
  object = m1, 
  newdata = df2, 
  samples = 1,
  formula = ~ Intercept +
    Eff.spde +
    # Eff.slope +
    # Eff.tree_cover +
    # Eff.swf + 
    Eff.elevation)

ggplot() + 
  gg(data = sample, aes(fill = mean), geom = "tile") +
  geom_sf(data = ireland, fill = NA) +
  # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
  labs(x = "", y = "", fill = "Median",
       title = "Red deer distribution (sample)") +
  theme_bw() +
  scale_fill_viridis_c() +
  NULL +
  
ggplot() + 
  gg(data = pred, aes(fill = mean), geom = "tile") +
  geom_sf(data = ireland, fill = NA) +
  # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
  labs(x = "", y = "", fill = "Median",
       title = "Red deer distribution (mean)") +
  theme_bw() +
  scale_fill_viridis_c() +
  NULL +
  
plot_layout(ncol = 2)










