rm(list = ls())
source("scripts/setups.R")
library(GeoThinneR)

# 1. Load data ####

ireland <- st_read("data/ireland_ITM.shp") %>% 
  st_transform(projKM)

PO_data <- read.csv("data/PO_data_all.csv", row.names = NULL)

PO_data_ll <- PO_data %>% 
  filter(Species %in% c("RedDeer")) %>% 
  filter(Y < 965) %>% 
  dplyr::select(PO, X, Y) %>% 
  drop_na() %>% 
  st_as_sf(coords = c("X", "Y")) %>% 
  st_set_crs(projKM) %>% 
  st_transform(PO_data_thinned, crs = 4326)

coords <- st_coordinates(PO_data_ll)
thinning <- distance_thinning(thin_dist = 2.5, 
                              coordinates = coords, 
                              trials = 1000)

PO_data_thinned <- PO_data_ll %>% 
  mutate(thin = thinning[[1]]) %>% 
  filter(thin == TRUE) %>% 
  select(-thin) %>% 
  st_transform(projKM)

pp <- ggplot() + 
  geom_sf(data = ireland, fill = NA) + 
  geom_sf(data = PO_data_ll, size = 2) +
  geom_sf(data = PO_data_thinned, col = "red", shape = 1, size = 2) +
  theme_bw()

ggsave(filename = "mesh_outputs/point_data.pdf", 
       height = 12, width = 8.5, units = "in") 


covar_stack <- rast("large_env_data/covar_subset_KM.grd")
covar_stack <- subset(covar_stack, 
                      c("landCover", "forest_distances", "human_footprint_index"),
                      negate = T) # remove factor var
names(covar_stack) <- c("Tree cover density", "Elevation", "Slope", "Small woody features")
covar_scaled <- scale(covar_stack)
covar_scaled <- trim(covar_scaled)
plot(covar_scaled)

covar_scaled <- project(covar_scaled, projKM)


# 2. Create mesh ####

boundary <- ireland %>% 
  summarise() %>%
  st_simplify(dTolerance = 10) %>% 
  st_buffer(dist = 20)

boundary2 <- boundary %>% 
  st_buffer(dist = 200)

library(fmesher)
mesh <- fm_mesh_2d_inla(boundary = list(boundary, boundary2),
                        max.edge = c(20,100), cutoff = 10, crs = projKM)

## 2.1 Create prediction dataset ####

df2 <- fm_pixels(mesh,
                 dims = c(72*10, 81*10),
                 mask = boundary,
                 format = "sf")


# 3. Set up SPDE ####

pcmatern <- inla.spde2.pcmatern(mesh, 
                                prior.range = c(120, 0.01), #in km                                
                                prior.sigma = c(0.1, 0.01)) 
# 4. Formula ####
form1 <- geometry ~  Intercept(1)  +
  Eff.elevation(covar_scaled$Elevation, model = "linear") +
  Eff.slope(covar_scaled$Slope, model = "linear") +
  Eff.tree_cover(covar_scaled$`Tree cover density`, model = "linear") +
  Eff.swf(covar_scaled$`Small woody features`, model = "linear") +
  Eff.spde(geometry, model = pcmatern) +
  NULL

# 5. FM_INT 1 ####

## 5.1 Subdivide mesh ####
load("mesh_outputs/objects/m1.Rdata")

ipoints1 <-  fm_int(fm_subdivide(mesh, 1))

## 5.2 Run model and summary ####

m1 <- lgcp(components = form1,
           data = PO_data_thinned,
           # samplers = boundary,
           ips = ipoints1)

summary(m1)


plot(spde.posterior(m1, "Eff.spde", what = "range")) +
plot(spde.posterior(m1, "Eff.spde", what = "variance"))


## 5.3 Predict #### 

pred1 <- predict(
  object = m1, 
  newdata = df2, 
  samples = 100,
  formula = ~ Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation)

resp1 <- predict(
  object = m1, 
  newdata = df2, 
  samples = 100,
  formula = ~ exp(Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation))

spde1 <- predict(
  object = m1, 
  newdata = df2, 
  samples = 100,
  formula = ~ Eff.spde )

(Lambda1 <- predict(
  m1,
  fm_int(mesh, boundary),
  ~ sum(weight * exp(Intercept +
                       Eff.spde +
                       Eff.slope +
                       Eff.tree_cover +
                       Eff.swf +
                       Eff.elevation))
))

## 5.4. Plot ####
ggplot() + 
  gg(data = pred1, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = ireland, fill = NA) +
  # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
  labs(x = "", y = "", fill = "Median",
       title = "Red deer distribution (linear scale)") +
  theme_bw() +
  scale_fill_viridis_c() +
  NULL +
  
ggplot() + 
  gg(data = spde1, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = ireland, fill = NA) +
  # geom_sf(data = PO_data_thinned, alpha = 0.5, size = 1) +
  labs(x = "", y = "", fill = "Median",
       title = "Spatial Random Effect") +
  theme_bw() +
  scale_fill_distiller(palette = 'RdBu') + 
  NULL +
  
plot_layout(ncol = 2)

fixed1 <-  tibble::rownames_to_column(m1$summary.fixed, "Variable")

fixed1 <- fixed1 %>% 
  select(Variable, Median = '0.5quant', 
         Low = '0.025quant', High = '0.975quant') %>% 
  mutate(Model = "fm_subdivide_1")

ggplot(fixed1) + 
  geom_point(aes(x = Variable, y = Median, col = Model), 
             position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = Low, ymax = High, 
                    x = Variable, col = Model), 
                position=position_dodge(width=0.5), 
                width = 0.1) + 
  scale_color_brewer(type = "qual", palette = "Set1") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() + 
  labs(title = "fm_int_1") 

save(m1, ipoints1, pred1, resp1, spde1, Lambda1, file = "mesh_outputs/objects/m1.Rdata")

# 6. FM_INT 2 ####

## 6.1 Subdivide mesh ####

load("mesh_outputs/objects/m2.Rdata")

ipoints2 <-  fm_int(fm_subdivide(mesh, 2))

## 6.2 Run model and summary ####

m2 <- lgcp(components = form1,
           data = PO_data_thinned,
           # samplers = boundary,
           ips = ipoints2)

summary(m2)


## 6.3 Predict #### 

pred2 <- predict(
  object = m2, 
  newdata = df2, 
  samples = 100,
  formula = ~ Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation)

resp2 <- predict(
  object = m2, 
  newdata = df2, 
  samples = 100,
  formula = ~ exp(Intercept +
                    Eff.spde +
                    Eff.slope +
                    Eff.tree_cover +
                    Eff.swf +
                    Eff.elevation))

spde2 <- predict(
  object = m2, 
  newdata = df2, 
  samples = 100,
  formula = ~ Eff.spde)

Lambda2 <- predict(
  m2,
  fm_int(mesh, boundary),
  ~ sum(weight * exp(Intercept +
                       Eff.spde +
                       Eff.slope +
                       Eff.tree_cover +
                       Eff.swf +
                       Eff.elevation))
)

bind_rows(Lambda1, Lambda2)

## 6.4. Plot ####
ggplot() + 
  gg(data = resp2, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = ireland, fill = NA) +
  # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
  labs(x = "", y = "", fill = "Median",
       title = "Red deer distribution (linear scale)") +
  theme_bw() +
  scale_fill_viridis_c() +
  NULL 

ggplot() + 
  gg(data = pred2, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = ireland, fill = NA) +
  # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
  labs(x = "", y = "", fill = "Median",
       title = "Red deer distribution (linear scale)") +
  theme_bw() +
  scale_fill_viridis_c() +
  NULL +
  
ggplot() + 
  gg(data = spde2, aes(fill = q0.5), geom = "tile") +
  geom_sf(data = ireland, fill = NA) +
  labs(x = "", y = "", fill = "Median",
       title = "Spatial Random Effect") +
  theme_bw() +
  scale_fill_distiller(palette = 'RdBu') + 
  NULL +
  
plot_layout(ncol = 2)

fixed2 <-  tibble::rownames_to_column(m2$summary.fixed, "Variable") %>% 
  select(Variable, Median = '0.5quant', 
         Low = '0.025quant', High = '0.975quant') %>% 
  mutate(Model = "fm_subdivide_2")

all_fixed <- bind_rows(fixed1, fixed2)

ggplot(all_fixed) + 
  geom_point(aes(x = Variable, y = Median, col = Model), 
             position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = Low, ymax = High, 
                    x = Variable, col = Model), 
                position=position_dodge(width=0.5), 
                width = 0.1) + 
  scale_color_brewer(type = "qual", palette = "Set1") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() + 
  labs(title = "all") 

save(m2, ipoints2, pred2, resp2, spde2, Lambda2, file = "mesh_outputs/objects/m2.Rdata")

# 7. FM_INT 3 ####

## 7.1 Subdivide mesh ####
load("mesh_outputs/objects/m3.Rdata")

ipoints3 <-  fm_int(fm_subdivide(mesh, 3))

## 7.2 Run model and summary ####

m3 <- lgcp(components = form1,
           data = PO_data_thinned,
           # samplers = boundary,
           ips = ipoints3)

summary(m3)

# plot(spde.posterior(m1, "Eff.spde", what = "range")) + 
#   plot(spde.posterior(m2, "Eff.spde", what = "range")) + 
#   plot(spde.posterior(m3, "Eff.spde", what = "range")) 
  

# plot(spde.posterior(m3, "Eff.spde", what = "variance"))


## 7.3 Predict #### 

pred3 <- predict(
  object = m3, 
  newdata = df2, 
  samples = 100,
  formula = ~ Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation)

resp3 <- predict(
  object = m3, 
  newdata = df2, 
  samples = 100,
  formula = ~ exp(Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation))

spde3 <- predict(
  object = m3, 
  newdata = df2, 
  samples = 100,
  formula = ~ Eff.spde)

Lambda3 <- predict(
  m3,
  fm_int(mesh, boundary),
  ~ sum(weight * exp(Intercept +
                       Eff.spde +
                       Eff.slope +
                       Eff.tree_cover +
                       Eff.swf +
                       Eff.elevation))
)

bind_rows(Lambda1, Lambda2, Lambda3)

## 7.4. Plot ####
# ggplot() + 
#   gg(data = pred3, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
#   labs(x = "", y = "", fill = "Median",
#        title = "Red deer distribution (linear scale)") +
#   theme_bw() +
#   scale_fill_viridis_c() +
#   NULL +
#   
# ggplot() + 
#   gg(data = spde3, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   geom_sf(data = PO_data_thinned, alpha = 0.5, size = 1) +
#   labs(x = "", y = "", fill = "Median",
#        title = "Spatial Random Effect") +
#   theme_bw() +
#   scale_fill_distiller(palette = 'RdBu') + 
#   NULL +
#   
#   plot_layout(ncol = 2)
# 
fixed3 <-  tibble::rownames_to_column(m3$summary.fixed, "Variable")

fixed3 <- fixed3 %>%
  select(Variable, Median = '0.5quant',
         Low = '0.025quant', High = '0.975quant') %>%
  mutate(Model = "fm_subdivide_3")
# 
# all_fixed <- bind_rows(fixed1, fixed2, fixed3)
# 
# ggplot(all_fixed) + 
#   geom_point(aes(x = Variable, y = Median, col = Model), 
#              position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin = Low, ymax = High, 
#                     x = Variable, col = Model), 
#                 position=position_dodge(width=0.5), 
#                 width = 0.1) + 
#   scale_color_brewer(type = "qual", palette = "Set1") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_bw() + 
#   labs(title = "all") 

save(m3, ipoints3, pred3, resp3, spde3, Lambda3, file = "mesh_outputs/objects/m3.Rdata")

# 8. FM_INT 4 ####

## 8.1 Subdivide mesh ####
load("mesh_outputs/objects/m4.Rdata")

ipoints4 <-  fm_int(fm_subdivide(mesh, 4))

## 8.2 Run model and summary ####

m4 <- lgcp(components = form1,
           data = PO_data_thinned,
           # samplers = boundary,
           ips = ipoints4)

summary(m4)


## 8.3 Predict #### 

pred4 <- predict(
  object = m4, 
  newdata = df2, 
  samples = 100,
  formula = ~ Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation)

resp4 <- predict(
  object = m4, 
  newdata = df2, 
  samples = 100,
  formula = ~ exp(Intercept +
                    Eff.spde +
                    Eff.slope +
                    Eff.tree_cover +
                    Eff.swf +
                    Eff.elevation))

spde4 <- predict(
  object = m4, 
  newdata = df2, 
  samples = 100,
  formula = ~ Eff.spde)

Lambda4 <- predict(
  m4,
  fm_int(mesh, boundary),
  ~ sum(weight * exp(Intercept +
                       Eff.spde +
                       Eff.slope +
                       Eff.tree_cover +
                       Eff.swf +
                       Eff.elevation))
)

bind_rows(Lambda1, Lambda2, Lambda3, Lambda4)

## 8.4. Plot ####
# ggplot() + 
#   gg(data = pred4, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
#   labs(x = "", y = "", fill = "Median",
#        title = "Red deer distribution (linear scale)") +
#   theme_bw() +
#   scale_fill_viridis_c() +
#   NULL +
#   
# ggplot() + 
#   gg(data = spde4, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   geom_sf(data = PO_data_thinned, alpha = 0.5, size = 1) +
#   labs(x = "", y = "", fill = "Median",
#        title = "Spatial Random Effect") +
#   theme_bw() +
#   scale_fill_distiller(palette = 'RdBu') + 
#   NULL +
#   
# plot_layout(ncol = 2)
# 
fixed4 <-  tibble::rownames_to_column(m4$summary.fixed, "Variable")

fixed4 <- fixed4 %>%
  select(Variable, Median = '0.5quant',
         Low = '0.025quant', High = '0.975quant') %>%
  mutate(Model = "fm_subdivide_4")
# 
# all_fixed <- bind_rows(fixed1, fixed2, fixed3, fixed4) 
# 
# ggplot(all_fixed) + 
#   geom_point(aes(x = Variable, y = Median, col = Model), 
#              position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin = Low, ymax = High, 
#                     x = Variable, col = Model), 
#                 position=position_dodge(width=0.5), 
#                 width = 0.1) + 
#   scale_color_brewer(type = "qual", palette = "Set1") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_bw() + 
#   labs(title = "all") 

save(m4, ipoints4, pred4, resp4, spde4, Lambda4, file = "mesh_outputs/objects/m4.Rdata")

# 9. FM_INT 5 ####

## 9.1 Subdivide mesh ####
load("mesh_outputs/objects/m5.Rdata")

ipoints5 <-  fm_int(fm_subdivide(mesh, 5))

## 9.2 Run model and summary ####

m5 <- lgcp(components = form1,
           data = PO_data_thinned,
           # samplers = boundary,
           ips = ipoints5)

summary(m5)

## 9.3 Predict #### 

pred5 <- predict(
  object = m5, 
  newdata = df2, 
  samples = 100,
  formula = ~ Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation)

resp5 <- predict(
  object = m5, 
  newdata = df2, 
  samples = 100,
  formula = ~ exp(Intercept +
                    Eff.spde +
                    Eff.slope +
                    Eff.tree_cover +
                    Eff.swf +
                    Eff.elevation))

spde5 <- predict(
  object = m5, 
  newdata = df2, 
  samples = 100,
  formula = ~ Eff.spde)

Lambda5 <- predict(
  m5,
  fm_int(mesh, boundary),
  ~ sum(weight * exp(Intercept +
                       Eff.spde +
                       Eff.slope +
                       Eff.tree_cover +
                       Eff.swf +
                       Eff.elevation))
)

bind_rows(Lambda1, Lambda2, Lambda3, Lambda4, Lambda5)

## 9.4. Plot ####
# ggplot() + 
#   gg(data = pred5, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
#   labs(x = "", y = "", fill = "Median",
#        title = "Red deer distribution (linear scale)") +
#   theme_bw() +
#   scale_fill_viridis_c() +
#   NULL +
#   
# ggplot() + 
#   gg(data = spde5, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   geom_sf(data = PO_data_thinned, alpha = 0.5, size = 1) +
#   labs(x = "", y = "", fill = "Median",
#        title = "Spatial Random Effect") +
#   theme_bw() +
#   scale_fill_distiller(palette = 'RdBu') + 
#   NULL +
#   
#   plot_layout(ncol = 2)
# 
fixed5 <-  tibble::rownames_to_column(m5$summary.fixed, "Variable")

fixed5 <- fixed5 %>%
  select(Variable, Median = '0.5quant',
         Low = '0.025quant', High = '0.975quant') %>%
  mutate(Model = "fm_subdivide_5")
# 
# all_fixed <- bind_rows(fixed1, fixed2, fixed3, fixed4, fixed5) 
# 
# ggplot(all_fixed) + 
#   geom_point(aes(x = Variable, y = Median, col = Model), 
#              position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin = Low, ymax = High, 
#                     x = Variable, col = Model), 
#                 position=position_dodge(width=0.5), 
#                 width = 0.1) + 
#   scale_color_brewer(type = "qual", palette = "Set1") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_bw() + 
#   labs(title = "all") 

save(m5, ipoints5, pred5, spde5, resp5, Lambda5, file = "mesh_outputs/objects/m5.Rdata")

# 10. FM_INT 6 ####

## 10.1 Subdivide mesh ####
load("mesh_outputs/objects/m6.Rdata")

ipoints6 <-  fm_int(fm_subdivide(mesh, 6))

## 10.2 Run model and summary ####

m6 <- lgcp(components = form1,
           data = PO_data_thinned,
           # samplers = boundary,
           ips = ipoints6)

summary(m6)


## 10.3 Predict #### 

pred6 <- predict(
  object = m6, 
  newdata = df2, 
  samples = 100,
  formula = ~ Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation)

resp6 <- predict(
  object = m6, 
  newdata = df2, 
  samples = 100,
  formula = ~ exp(Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation))

spde6 <- predict(
  object = m6, 
  newdata = df2, 
  samples = 100,
  formula = ~ Eff.spde)

Lambda6 <- predict(
  m6,
  fm_int(mesh, boundary),
  ~ sum(weight * exp(Intercept +
                       Eff.spde +
                       Eff.slope +
                       Eff.tree_cover +
                       Eff.swf +
                       Eff.elevation))
)

bind_rows(Lambda1, Lambda2, Lambda3, Lambda4, Lambda5, Lambda6)

## 10.4. Plot ####
# ggplot() + 
#   gg(data = pred6, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
#   labs(x = "", y = "", fill = "Median",
#        title = "Red deer distribution (linear scale)") +
#   theme_bw() +
#   scale_fill_viridis_c() +
#   NULL +
#   
#   ggplot() + 
#   gg(data = spde6, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   geom_sf(data = PO_data_thinned, alpha = 0.5, size = 1) +
#   labs(x = "", y = "", fill = "Median",
#        title = "Spatial Random Effect") +
#   theme_bw() +
#   scale_fill_distiller(palette = 'RdBu') + 
#   NULL +
#   
#   plot_layout(ncol = 2)
# 
fixed6 <-  tibble::rownames_to_column(m6$summary.fixed, "Variable")

fixed6 <- fixed6 %>%
  select(Variable, Median = '0.5quant',
         Low = '0.025quant', High = '0.975quant') %>%
  mutate(Model = "fm_subdivide_6")
 
# all_fixed <- bind_rows(fixed1, fixed2, fixed3, fixed4, fixed5, fixed6) 
# 
# covars <- all_fixed %>% 
#   filter(Variable != "Intercept")
# 
# intercept <- all_fixed %>% 
#   filter(Variable == "Intercept")
# 
# ggplot(covars) + 
#   geom_point(aes(x = Variable, y = Median, col = Model), 
#              position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin = Low, ymax = High, 
#                     x = Variable, col = Model), 
#                 position=position_dodge(width=0.5), 
#                 width = 0.1) + 
#   scale_color_brewer(type = "qual", palette = "Set1") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_bw() + 
#   labs(title = "Covariates") +
#   
# ggplot(intercept) + 
#   geom_point(aes(x = Variable, y = Median, col = Model), 
#              position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin = Low, ymax = High, 
#                     x = Variable, col = Model), 
#                 position=position_dodge(width=0.5), 
#                 width = 0.1) + 
#   scale_color_brewer(type = "qual", palette = "Set1") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_bw() + 
#   labs(title = "Intercept")

save(m6, ipoints6, pred6, resp6, spde6, Lambda6, file = "mesh_outputs/objects/m6.Rdata")

# 11. FM_INT 7 ####

## 11.1 Subdivide mesh ####
load("mesh_outputs/objects/m7.Rdata")

ipoints7 <-  fm_int(fm_subdivide(mesh, 7))

## 11.2 Run model and summary ####

m7 <- lgcp(components = form1,
           data = PO_data_thinned,
           # samplers = boundary,
           ips = ipoints7)

summary(m7)

## 11.3 Predict #### 

pred7 <- predict(
  object = m7, 
  newdata = df2, 
  samples = 100,
  formula = ~ Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation)

resp7 <- predict(
  object = m7, 
  newdata = df2, 
  samples = 100,
  formula = ~ exp(Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation))

spde7 <- predict(
  object = m7, 
  newdata = df2, 
  samples = 100,
  formula = ~ Eff.spde)

Lambda7 <- predict(
  m7,
  fm_int(mesh, boundary),
  ~ sum(weight * exp(Intercept +
                       Eff.spde +
                       Eff.slope +
                       Eff.tree_cover +
                       Eff.swf +
                       Eff.elevation))
)

# bind_rows(Lambda1, Lambda2, Lambda3, Lambda4, Lambda5, Lambda6, Lambda7)
# 
## 11.4. Plot ####
# ggplot() + 
#   gg(data = resp7, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
#   labs(x = "", y = "", fill = "Median",
#        title = "Red deer distribution (linear scale)") +
#   theme_bw() +
#   scale_fill_viridis_c() +
#   NULL 
# 
# 
# ggplot() + 
#   gg(data = pred7, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
#   labs(x = "", y = "", fill = "Median",
#        title = "Red deer distribution (linear scale)") +
#   theme_bw() +
#   scale_fill_viridis_c() +
#   NULL +
#   
# ggplot() + 
#   gg(data = spde7, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   geom_sf(data = PO_data_thinned, alpha = 0.5, size = 1) +
#   labs(x = "", y = "", fill = "Median",
#        title = "Spatial Random Effect") +
#   theme_bw() +
#   scale_fill_distiller(palette = 'RdBu') + 
#   NULL +
#   
#   plot_layout(ncol = 2)
# 
fixed7 <-  tibble::rownames_to_column(m7$summary.fixed, "Variable")

fixed7 <- fixed7 %>%
  select(Variable, Median = '0.5quant',
         Low = '0.025quant', High = '0.975quant') %>%
  mutate(Model = "fm_subdivide_7")

# all_fixed <- bind_rows(fixed1, fixed2, fixed3, fixed4, fixed5, fixed6, fixed7) 
# 
# covars <- all_fixed %>% 
#   filter(Variable != "Intercept")
# 
# intercept <- all_fixed %>% 
#   filter(Variable == "Intercept")
# 
# ggplot(covars) + 
#   geom_point(aes(x = Variable, y = Median, col = Model), 
#              position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin = Low, ymax = High, 
#                     x = Variable, col = Model), 
#                 position=position_dodge(width=0.5), 
#                 width = 0.1) + 
#   scale_color_brewer(type = "qual", palette = "Set1") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_bw() + 
#   labs(title = "Covariates") +
#   
#   ggplot(intercept) + 
#   geom_point(aes(x = Variable, y = Median, col = Model), 
#              position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin = Low, ymax = High, 
#                     x = Variable, col = Model), 
#                 position=position_dodge(width=0.5), 
#                 width = 0.1) + 
#   scale_color_brewer(type = "qual", palette = "Set1") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_bw() + 
#   labs(title = "Intercept") 

save(m7, ipoints7, pred7, spde7, resp7, Lambda7, file = "mesh_outputs/objects/m7.Rdata")

# 12. FM_INT 8 ####

## 12.1 Subdivide mesh ####
load("mesh_outputs/objects/m8.Rdata")

ipoints8 <-  fm_int(fm_subdivide(mesh, 8))

## 11.2 Run model and summary ####

m8 <- lgcp(components = form1,
           data = PO_data_thinned,
           # samplers = boundary,
           ips = ipoints8)

summary(m8)

## 11.3 Predict #### 

pred8 <- predict(
  object = m8, 
  newdata = df2, 
  samples = 100,
  formula = ~ Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation)

resp8 <- predict(
  object = m8, 
  newdata = df2, 
  samples = 100,
  formula = ~ exp(Intercept +
    Eff.spde +
    Eff.slope +
    Eff.tree_cover +
    Eff.swf +
    Eff.elevation))

spde8 <- predict(
  object = m8, 
  newdata = df2, 
  samples = 100,
  formula = ~ Eff.spde)

Lambda8 <- predict(
  m8,
  fm_int(mesh, boundary),
  ~ sum(weight * exp(Intercept +
                       Eff.spde +
                       Eff.slope +
                       Eff.tree_cover +
                       Eff.swf +
                       Eff.elevation))
)

bind_rows(Lambda1, Lambda2, Lambda3, Lambda4, Lambda5, Lambda6, Lambda7, Lambda8)

## 11.4. Plot ####
# ggplot() + 
#   gg(data = resp8, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
#   labs(x = "", y = "", fill = "Median",
#        title = "Red deer distribution (linear scale)") +
#   theme_bw() +
#   scale_fill_viridis_c() +
#   NULL
# 
# ggplot() + 
#   gg(data = pred8, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   # geom_sf(data = sett_subset, alpha = 0.5, size = 1, col = "white") +
#   labs(x = "", y = "", fill = "Median",
#        title = "Red deer distribution (linear scale)") +
#   theme_bw() +
#   scale_fill_viridis_c() +
#   NULL +
#   
#   ggplot() + 
#   gg(data = spde8, aes(fill = q0.5), geom = "tile") +
#   geom_sf(data = ireland, fill = NA) +
#   geom_sf(data = PO_data_thinned, alpha = 0.5, size = 1) +
#   labs(x = "", y = "", fill = "Median",
#        title = "Spatial Random Effect") +
#   theme_bw() +
#   scale_fill_distiller(palette = 'RdBu') + 
#   NULL +
#   
#   plot_layout(ncol = 2)
# 
fixed8 <-  tibble::rownames_to_column(m8$summary.fixed, "Variable")

fixed8 <- fixed8 %>%
  select(Variable, Median = '0.5quant',
         Low = '0.025quant', High = '0.975quant') %>%
  mutate(Model = "fm_subdivide_8")

all_fixed <- bind_rows(fixed1, fixed2, fixed3, fixed4, fixed5, fixed6, fixed7, fixed8)

# covars <- all_fixed %>% 
#   filter(Variable != "Intercept")
# 
# intercept <- all_fixed %>% 
#   filter(Variable == "Intercept")
# 
# ggplot(covars) + 
#   geom_point(aes(x = Variable, y = Median, col = Model), 
#              position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin = Low, ymax = High, 
#                     x = Variable, col = Model), 
#                 position=position_dodge(width=0.5), 
#                 width = 0.1) + 
#   # scale_color_brewer(type = "qual", palette = "Set3") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_bw() + 
#   labs(title = "Covariates") +
#   
#   ggplot(intercept) + 
#   geom_point(aes(x = Variable, y = Median, col = Model), 
#              position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin = Low, ymax = High, 
#                     x = Variable, col = Model), 
#                 position=position_dodge(width=0.5), 
#                 width = 0.1) + 
#   # scale_color_brewer(type = "qual", palette = "Set3") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   theme_bw() + 
#   labs(title = "Intercept") 

save(m8, ipoints8, pred8, spde8, resp8, Lambda8, file = "mesh_outputs/objects/m8.Rdata")


# 13. final plots ####

## 13.1 Predicted means ####

pred1_raster <- data.frame(st_coordinates(pred1), z = pred1$q0.5) %>% 
  rast(.)

pred2_raster <- data.frame(st_coordinates(pred2), z=pred2$q0.5) %>% 
  rast(.)

pred3_raster <- data.frame(st_coordinates(pred3), z=pred3$q0.5) %>% 
  rast(.)

pred4_raster <- data.frame(st_coordinates(pred4), z=pred4$q0.5) %>% 
  rast(.)

pred5_raster <- data.frame(st_coordinates(pred5), z=pred5$q0.5) %>% 
  rast(.)

pred6_raster <- data.frame(st_coordinates(pred6), z=pred6$q0.5) %>% 
  rast(.)

pred7_raster <- data.frame(st_coordinates(pred7), z=pred7$q0.5) %>% 
  rast(.)

pred8_raster <- data.frame(st_coordinates(pred8), z=pred8$q0.5) %>% 
  rast(.)


all_preds <- c(pred1_raster, pred2_raster, pred3_raster, pred4_raster, 
               pred5_raster, pred6_raster, pred7_raster, pred8_raster)
crs(all_preds) <- crs(ireland)
all_preds <- mask(all_preds, ireland)
names(all_preds) <- c("1 Subdivision", "2 subdivisions","3 subdivisions",
                      "4 subdivisions", "5 subdivisions", "6 subdivisions", 
                      "7 subdivisions", "8 subdivisions") 

ggplot() + 
  geom_spatraster(data = all_preds) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA) + 
  theme_bw() +
  facet_wrap(~lyr, ncol = 4)

ggsave(filename = "mesh_outputs/predictions.pdf", 
       width = 12, height = 8.5, units = "in") 

## 13.1bis Predicted in response scale ####

resp1_raster <- data.frame(st_coordinates(resp1), z = resp1$q0.5) %>% 
  rast(.)

resp2_raster <- data.frame(st_coordinates(resp2), z=resp2$q0.5) %>% 
  rast(.)

resp3_raster <- data.frame(st_coordinates(resp3), z=resp3$q0.5) %>% 
  rast(.)

resp4_raster <- data.frame(st_coordinates(resp4), z=resp4$q0.5) %>% 
  rast(.)

resp5_raster <- data.frame(st_coordinates(resp5), z=resp5$q0.5) %>% 
  rast(.)

resp6_raster <- data.frame(st_coordinates(resp6), z=resp6$q0.5) %>% 
  rast(.)

resp7_raster <- data.frame(st_coordinates(resp7), z=resp7$q0.5) %>% 
  rast(.)

resp8_raster <- data.frame(st_coordinates(resp8), z=resp8$q0.5) %>% 
  rast(.)


all_resps <- c(resp1_raster, resp2_raster, resp3_raster, resp4_raster, 
               resp5_raster, resp6_raster, resp7_raster, resp8_raster)
crs(all_resps) <- crs(ireland)
all_resps <- mask(all_resps, ireland)


names(all_resps) <- c("1 Subdivision", "2 subdivisions","3 subdivisions",
                      "4 subdivisions", "5 subdivisions", "6 subdivisions", 
                      "7 subdivisions", "8 subdivisions") 

ggplot() + 
  geom_spatraster(data = all_resps$`1 Subdivision`) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA, option = "H") + 
  ggtitle("1 subdivision") + 
  theme_bw() +

ggplot() + 
  geom_spatraster(data = all_resps$`2 subdivisions`) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA, option = "H") + 
  ggtitle("2 subdivisions") +
  theme_bw() +
  
ggplot() + 
  geom_spatraster(data = all_resps$`3 subdivisions`) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA, option = "H") + 
  ggtitle("3 subdivisions") +
  theme_bw() +
  
ggplot() + 
  geom_spatraster(data = all_resps$`4 subdivisions`) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA, option = "H") + 
  ggtitle("4 subdivisions") +
  theme_bw() + 
  
ggplot() + 
  geom_spatraster(data = all_resps$`5 subdivisions`) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA, option = "H") + 
  ggtitle("5 subdivisions") +
  theme_bw() +
  
ggplot() + 
  geom_spatraster(data = all_resps$`6 subdivisions`) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA, option = "H") + 
  ggtitle("6 subdivisions") +
  theme_bw() +
  
ggplot() + 
  geom_spatraster(data = all_resps$`7 subdivisions`) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA, option = "H") + 
  ggtitle("7 subdivisions") +
  theme_bw() +
  
ggplot() + 
  geom_spatraster(data = all_resps$`8 subdivisions`) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA, option = "H") + 
  ggtitle("8 subdivisions") +
  theme_bw() +
  
plot_layout(nrow = 3)

ggsave(filename = "mesh_outputs/response_pred.pdf", 
       width = 12, height = 10, units = "in") 

## 13.1.3 Prediction aggregated by county ####   

pred_by_county <- st_as_sf(extract(all_resps, ireland, fun = sum, bind = TRUE))

pred_by_county_agg <- pred_by_county %>% 
  select(NAME_TAG, X1.Subdivision:X8.subdivisions) %>% 
  pivot_longer(X1.Subdivision:X8.subdivisions, names_to = "Subdivision",
               values_to = "Prediction") %>% 
  mutate(Subdivision = str_remove(Subdivision, "X")) %>% 
  mutate(Subdivision = str_replace(Subdivision, "\\.", " "))

ggplot() + 
  geom_sf(data = pred_by_county_agg, 
          aes(fill = Prediction)) + 
  scale_fill_viridis() + 
  facet_wrap(~Subdivision) + 
  theme_bw()

ggsave(filename = "mesh_outputs/response_pred_agg.pdf", 
       width = 8.5, height = 12, units = "in") 


## 13.1.4 Prediction uncertainty ####

CI1_raster <- data.frame(st_coordinates(pred1), z = pred1$q0.975-pred1$q0.025) %>% 
  rast(.)

CI2_raster <- data.frame(st_coordinates(pred2), z = pred2$q0.975-pred1$q0.025) %>% 
  rast(.)

CI3_raster <- data.frame(st_coordinates(pred3), z = pred3$q0.975-pred1$q0.025) %>% 
  rast(.)

CI4_raster <- data.frame(st_coordinates(pred4), z = pred4$q0.975-pred1$q0.025) %>% 
  rast(.)

CI5_raster <- data.frame(st_coordinates(pred5), z = pred5$q0.975-pred1$q0.025) %>% 
  rast(.)

CI6_raster <- data.frame(st_coordinates(pred6), z = pred6$q0.975-pred1$q0.025) %>% 
  rast(.)

CI7_raster <- data.frame(st_coordinates(pred7), z = pred7$q0.975-pred1$q0.025) %>% 
  rast(.)

CI8_raster <- data.frame(st_coordinates(pred8), z = pred8$q0.975-pred1$q0.025) %>% 
  rast(.)

all_CIs <- c(CI1_raster, CI2_raster, CI3_raster, CI4_raster, 
             CI5_raster, CI6_raster, CI7_raster, CI8_raster)

crs(all_CIs) <- crs(ireland)

all_CIs <- mask(all_CIs, ireland)

names(all_CIs) <- c("1 Subdivision", "2 subdivisions","3 subdivisions",
                    "4 subdivisions", "5 subdivisions", "6 subdivisions", 
                    "7 subdivisions", "8 subdivisions") 

ggplot() + 
  geom_spatraster(data = all_CIs) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA) + 
  theme_bw() +
  facet_wrap(~lyr, ncol = 4)

ggsave(filename = "mesh_outputs/prediction_CIs.pdf", 
       width = 12, height = 8.5, units = "in") 


## 13.2 Spatial fields ####

spde1_raster <- data.frame(st_coordinates(spde1), z = spde1$q0.5) %>% 
  rast(.)

spde2_raster <- data.frame(st_coordinates(spde2), z=spde2$q0.5) %>% 
  rast(.)

spde3_raster <- data.frame(st_coordinates(spde3), z=spde3$q0.5) %>% 
  rast(.)

spde4_raster <- data.frame(st_coordinates(spde4), z=spde4$q0.5) %>% 
  rast(.)

spde5_raster <- data.frame(st_coordinates(spde5), z=spde5$q0.5) %>% 
  rast(.)

spde6_raster <- data.frame(st_coordinates(spde6), z=spde6$q0.5) %>% 
  rast(.)

spde7_raster <- data.frame(st_coordinates(spde7), z=spde7$q0.5) %>% 
  rast(.)

spde8_raster <- data.frame(st_coordinates(spde8), z=spde8$q0.5) %>% 
  rast(.)


all_spdes <- c(spde1_raster, spde2_raster, spde3_raster, spde4_raster, 
               spde5_raster, spde6_raster, spde7_raster, spde8_raster)

crs(all_spdes) <- crs(ireland)

all_spdes <- mask(all_spdes, ireland)

names(all_spdes) <- c("1 Subdivision", "2 subdivisions","3 subdivisions",
                      "4 subdivisions", "5 subdivisions", "6 subdivisions", 
                      "7 subdivisions", "8 subdivisions") 

ggplot() + 
  geom_spatraster(data = all_spdes ) +
  geom_sf(data = ireland, fill = NA) + 
  # scale_fill_distiller(
  #   na.value = NA,
  #   palette = "RdBu",
  #   limits = c(-1, 1) * 11.5, # Symmetrical limits
  #   direction = -1) +
  scale_fill_viridis_c(na.value = NA, option = "A") +
  theme_bw() +
  facet_wrap(~lyr, ncol = 4)

ggsave(filename = "mesh_outputs/spdes.pdf", 
       width = 12, height = 8.5, units = "in") 



ggplot() + 
  geom_spatraster(data = all_spdes[[1]]) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_viridis_c(na.value = NA, option = "A") + 
  labs(title = "1 subdivision") +
  theme_bw() + 
  
  
ggplot() + 
  geom_spatraster(data = all_spdes[[2]]-all_spdes[[1]]) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_distiller(
    na.value = NA,
    palette = "RdBu",
    limits = c(-1, 1) *1, # Symmetrical limits
    direction = -1) +
  labs(title = "Difference 2 subdivisions - 1 subdivision") + 
  theme_bw() +
  
ggplot() + 
  geom_spatraster(data = all_spdes[[3]]-all_spdes[[2]]) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_distiller(
    na.value = NA,
    palette = "RdBu",
    limits = c(-1, 1) *1 , # Symmetrical limits
    direction = -1) +
  labs(title = "Difference 3 subdivisions - 2 subdivisions") + 
  theme_bw() +

ggplot() + 
  geom_spatraster(data = all_spdes[[4]]-all_spdes[[3]]) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_distiller(
    na.value = NA,
    palette = "RdBu",
    limits = c(-1, 1) * 1, # Symmetrical limits
    direction = -1) +
  labs(title = "Difference 4 subdivisions - 3 subdivisions") + 
  theme_bw() +
  
ggplot() + 
  geom_spatraster(data = all_spdes[[5]]-all_spdes[[4]]) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_distiller(
    na.value = NA,
    palette = "RdBu",
    limits = c(-1, 1) * 1, # Symmetrical limits
    direction = -1) +
  labs(title = "Difference 5 subdivisions - 4 subdivisions") + 
  theme_bw() +
  
ggplot() + 
  geom_spatraster(data = all_spdes[[6]]-all_spdes[[5]]) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_distiller(
    na.value = NA,
    palette = "RdBu",
    limits = c(-1, 1) * 1, # Symmetrical limits
    direction = -1) +
  labs(title = "Difference 6 subdivisions - 5 subdivisions") + 
  theme_bw() +
  
ggplot() + 
  geom_spatraster(data = all_spdes[[7]]-all_spdes[[6]]) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_distiller(
    na.value = NA,
    palette = "RdBu",
    limits = c(-1, 1) * 1, # Symmetrical limits
    direction = -1) +
  labs(title = "Difference 7 subdivisions - 6 subdivisions") + 
  theme_bw() +
  
ggplot() + 
  geom_spatraster(data = all_spdes[[8]]-all_spdes[[7]]) +
  geom_sf(data = ireland, fill = NA) + 
  scale_fill_distiller(
    na.value = NA,
    palette = "RdBu",
    limits = c(-1, 1) * 1, # Symmetrical limits
    direction = -1) +
  labs(title = "Difference 8 subdivisions - 7 subdivisions") + 
  theme_bw() +

plot_layout(ncol = 3)
  


ggsave(filename = "mesh_outputs/spdes_dif.pdf", 
       width = 12, height = 8.5, units = "in") 


## 12.3 Predicted abundances plot  ####

pred_points <- bind_rows(Lambda1, Lambda2, Lambda3, Lambda4, 
                         Lambda5, Lambda6, Lambda7, Lambda8) %>% 
  mutate(model = paste(1:8, "subdivisions", sep = " "))


ggplot(pred_points) + 
  geom_point(aes(x = model, y = q0.5, col = model), 
             position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = q0.025, ymax = q0.975, 
                    x = model, col = model), 
                position=position_dodge(width=0.5), 
                width = 0.1) + 
  # scale_color_brewer(type = "qual", palette = "Set3") +
  geom_hline(yintercept = 574, linetype = 2) +
  theme_bw() + 
  labs(title = "Predicted total abundance", x = "Model", y = "Abundance") 

ggsave(filename = "mesh_outputs/TotAb.pdf", 
       width = 12, height = 8.5, units = "in") 

## 12.4 Fixed effects plot ####

all_fixed <- bind_rows(fixed1, fixed2, fixed3, fixed4, fixed5, fixed6, fixed7, fixed8) 

all_fixed <- all_fixed %>% 
  mutate(Model = recode(Model, 
                        fm_subdivide_1 = "1 subdivision", 
                        fm_subdivide_2 = "2 subdivisions",
                        fm_subdivide_3 = "3 subdivisions",
                        fm_subdivide_4 = "4 subdivisions",
                        fm_subdivide_5 = "5 subdivisions",
                        fm_subdivide_6 = "6 subdivisions",
                        fm_subdivide_7 = "7 subdivisions",
                        fm_subdivide_8 = "8 subdivisions"), 
         Variable = recode(Variable, 
                           Eff.elevation = "Elevation", 
                           Eff.slope = "Slope", 
                           Eff.tree_cover = "Tree cover", 
                           Eff.swf = "Small woody \n features"))

covars <- all_fixed %>% 
  filter(Variable != "Intercept") 

intercept <- all_fixed %>% 
  filter(Variable == "Intercept")

ggplot(covars) + 
  geom_point(aes(x = Variable, y = Median, col = Model), 
             position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = Low, ymax = High, 
                    x = Variable, col = Model), 
                position=position_dodge(width=0.5), 
                width = 0.1) + 
  # scale_color_brewer(type = "qual", palette = "Set3") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() + 
  labs(title = "Covariates") +
  
ggplot(intercept) + 
  geom_point(aes(x = Variable, y = Median, col = Model), 
             position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = Low, ymax = High, 
                    x = Variable, col = Model), 
                position=position_dodge(width=0.5), 
                width = 0.1) + 
  # scale_color_brewer(type = "qual", palette = "Set3") +
  # geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() + 
  labs(title = "Intercept") 

ggsave(filename = "mesh_outputs/FixedEffects.pdf", 
       width = 12, height = 6, units = "in") 


## 12.5 Integration points ####

ipoints1f <- st_filter(ipoints1, ireland)
ipoints2f <- st_filter(ipoints2, ireland)
ipoints3f <- st_filter(ipoints3, ireland)
ipoints4f <- st_filter(ipoints4, ireland)
ipoints5f <- st_filter(ipoints5, ireland)
ipoints6f <- st_filter(ipoints6, ireland)
ipoints7f <- st_filter(ipoints7, ireland)
ipoints8f <- st_filter(ipoints8, ireland)


ggplot() + 
  geom_sf(data = ireland, fill = NA, col = "darkred") + 
  geom_sf(data = ipoints1f, size = 0.5) + 
  theme_bw() + 

ggplot() + 
  geom_sf(data = ireland, fill = NA, col = "darkred") + 
  geom_sf(data = ipoints2f, size = 0.5) + 
  theme_bw() + 

ggplot() + 
  geom_sf(data = ireland, fill = NA, col = "darkred") + 
  geom_sf(data = ipoints3f, size = 0.5) + 
  theme_bw() + 
  
ggplot() + 
  geom_sf(data = ireland, fill = NA, col = "darkred") + 
  geom_sf(data = ipoints4f, size = 0.5) + 
  theme_bw() + 
  
ggplot() + 
  geom_sf(data = ireland, fill = NA, col = "darkred") + 
  geom_sf(data = ipoints5f, size = 0.5) + 
  theme_bw() + 
  
ggplot() + 
  geom_sf(data = ireland, fill = NA, col = "darkred") + 
  geom_sf(data = ipoints6f, size = 0.5) + 
  theme_bw() + 
  
ggplot() + 
  geom_sf(data = ireland, fill = NA, col = "darkred") + 
  geom_sf(data = ipoints7f, size = 0.5) + 
  theme_bw() + 
  
ggplot() + 
  geom_sf(data = ireland, fill = NA, col = "darkred") + 
  geom_sf(data = ipoints8f, size = 0.5) + 
  theme_bw() + 
  
plot_layout(ncol = 4)



ggsave(filename = "mesh_outputs/ipoints.pdf", 
       width = 12, height = 6, units = "in") 








