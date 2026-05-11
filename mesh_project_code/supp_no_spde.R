rm(list = ls())
source("scripts/setups.R")

# 1. Load data ####

ireland <- st_read("data/ireland_ITM.shp") %>% 
  st_transform(projKM)

PO_data <- read.csv("data/PO_data_all.csv", row.names = NULL)

PO_data_sel <- PO_data %>% 
  filter(Species %in% c("RedDeer")) %>% 
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

pp <- ggplot() + 
  geom_spatraster(data = covar_scaled, aes(fill = Elevation)) +
  scale_fill_viridis_c(na.value = NA) + 
  geom_sf(data = ireland, fill = NA) + 
  geom_sf(data = PO_data_sel) + 
  theme_bw()

Cairo::CairoPDF(file = "mesh_outputs/point_data.pdf", width = 6, height = 9)
pp
dev.off()

Cairo::CairoPDF(file = "mesh_outputs/covar_data.pdf", width = 9, height = 6)
plot(covar_stack)
dev.off()

# 2. Create large meshes ####

## 2.1 Inner boundaries ####

boundary <- ireland %>% 
  summarise() %>%
  st_simplify(dTolerance = 10) %>% 
  st_buffer(dist = 20)

set.seed(123)
shift2 <- data.frame(x = rnorm(1, mean = 5, sd = 5), y = rnorm(1, mean = 5, sd = 5)) %>% 
  st_as_sf(coords = c("x", "y"))

shift3 <- data.frame(x = rnorm(1, mean = 5, sd = 5), y = rnorm(1, mean = 5, sd = 5)) %>% 
  st_as_sf(coords = c("x", "y"))

shift4 <- data.frame(x = rnorm(1, mean = 5, sd = 5), y = rnorm(1, mean = 5, sd = 5)) %>% 
  st_as_sf(coords = c("x", "y"))

boundary2 <- boundary$geometry + shift2$geometry
boundary3 <- boundary$geometry + shift3$geometry
boundary4 <- boundary$geometry + shift4$geometry

boundary2 <- boundary2 %>% st_set_crs(st_crs(boundary))
boundary3 <- boundary3 %>% st_set_crs(st_crs(boundary))
boundary4 <- boundary4 %>% st_set_crs(st_crs(boundary))

Cairo::CairoPDF(file = "mesh_outputs/boundaries.pdf", width = 6, height = 9)
ggplot() + 
  geom_sf(data = ireland, col = "black", fill = "gray") + 
  geom_sf(data = boundary, fill = NA) + 
  geom_sf(data = boundary2, col = "red", fill = NA) + 
  geom_sf(data = boundary3, col = "orange", fill = NA) + 
  geom_sf(data = boundary4, col = "purple", fill = NA) + 
  theme_bw()
dev.off()

bound_list <- list(boundary, boundary2, boundary3, boundary4)

## 2.2 Outer boundaries ####

bound2_list <- list()

for(i in seq_along(bound_list)) {
  boundary <- bound_list[[i]]
  boundary2 <- boundary %>% 
    st_buffer(dist = 200)

  bound2_list[i] <- list(boundary2)  
}

## 2.3 Create meshes ####
meshes <-  list()
ipoints <- list()

for(i in 1:4) {
  meshes[[i]] =  fm_mesh_2d_inla(boundary = list(bound_list[[i]], bound2_list[[i]]),
                                 max.edge = c(20,100), cutoff = 20, crs = projKM)
  ipoints[[i]] <- fm_int(meshes[[i]])
}

ipb <- ggplot() + 
  geom_spatraster(data = covar_scaled, aes(fill = Elevation), alpha = 0.5) + 
  scale_fill_viridis_c(na.value = NA) + 
  gg(meshes[[1]]) + 
  geom_sf(data = ipoints[[1]], col = "red") + 
  geom_sf(data = ipoints[[2]], col = "darkgreen") +
  geom_sf(data = ipoints[[3]], col = "orange") + 
  geom_sf(data = ipoints[[4]], col = "purple") +
  coord_sf(xlim = st_bbox(covar_scaled)[c(1,3)], 
           ylim = st_bbox(covar_scaled)[c(2,4)], 
           expand = F) + 
  labs(title = "30 km edge") +
  theme_bw() 


# 3. Model with large mesh ####

model_list_b <- list()
fixed_list_b <- list()

for(i in 1:4) {
  print(i)
  mesh <- meshes[[i]]
  boundary <- bound_list[[i]] 

  ## 3.1 Formula ####
  form1 <- geometry ~  Intercept(1)  +
    Eff.elevation(covar_scaled$Elevation, model = "linear") +
    Eff.slope(covar_scaled$Slope, model = "linear") +
    Eff.tree_cover(covar_scaled$`Tree cover density`, model = "linear") + 
    Eff.swf(covar_scaled$`Small woody features`, model = "linear") +
    NULL
  
  ## 3.2 Run model and summary ####
  m1 <- lgcp(components = form1,
             data = PO_data_sel,
             # samplers = boundary,
             domain =  list(geometry = mesh))
  
  model_list_b[i] <- list(m1)
  
  fixed <-  tibble::rownames_to_column(m1$summary.fixed, "Variable")
  fixed_list_b[i] <- list(fixed)
}  

## 3.4 Plot fixed effects ####

big_x <- data.table::rbindlist(fixed_list_b, idcol = "Model")

big_x2 <-  big_x %>% 
  select(Variable, Model, Median = '0.5quant', 
         Low = '0.025quant', High = '0.975quant') %>% 
  mutate(Model = paste("Mesh", Model, sep = "_")) 

(fixed_big <- ggplot(big_x2) + 
  geom_point(aes(x = Variable, y = Median, col = Model), 
             position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = Low, ymax = High, 
                    x = Variable, col = Model), 
                position=position_dodge(width=0.5), 
                width = 0.2) + 
  scale_color_brewer(type = "qual", palette = "Set1") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() + 
  labs(title = "20km edge") +
  # facet_wrap(~Variable, scales = "free") + 
  NULL)


# 4. Create medium meshes ####

med_meshes <- list()
med_ipoints <- list()

for(i in seq_along(meshes)) {
  med_meshes[[i]] <-  fmesher:::fm_subdivide(meshes[[i]])
  med_ipoints[[i]] <- fm_int(med_meshes[[i]])
} 

ipm <- ggplot() + 
  geom_spatraster(data = covar_scaled, aes(fill = Elevation), alpha = 0.5) + 
  gg(med_meshes[[1]]) + 
  scale_fill_viridis_c(na.value = NA) + 
  geom_sf(data = med_ipoints[[1]], col = "red")  +
  geom_sf(data = med_ipoints[[2]], col = "darkgreen") +
  geom_sf(data = med_ipoints[[3]], col = "orange") + 
  geom_sf(data = med_ipoints[[4]], col = "purple") +
  coord_sf(xlim = st_bbox(covar_scaled)[c(1,3)], 
           ylim = st_bbox(covar_scaled)[c(2,4)], 
           expand = F) + 
  labs(title = "20/2 km edge") +
  theme_bw()


# 5. Model with medium mesh ####

med_model_list <- list()
med_fixed_list <- list()

for(i in 1:4) {
  print(i)
  mesh <- med_meshes[[i]]
  boundary <- bound_list[[i]] 
  
  ## 5.1 Formula ####
  form1 <- geometry ~  Intercept(1)  +
    Eff.elevation(covar_scaled$Elevation, model = "linear") +
    Eff.slope(covar_scaled$Slope, model = "linear") +
    Eff.tree_cover(covar_scaled$`Tree cover density`, model = "linear") + 
    Eff.swf(covar_scaled$`Small woody features`, model = "linear") +
    NULL
  ## 5.2 Run model and summary ####
  
  m1 <- lgcp(components = form1,
             data = PO_data_sel,
             # samplers = boundary,
             domain =  list(geometry = mesh))
  
  med_model_list[i] <- list(m1)
  
  fixed <-  tibble::rownames_to_column(m1$summary.fixed, "Variable")
  med_fixed_list[i] <- list(fixed)

}  


## 5.4 Plot fixed effects ####

med_x <- data.table::rbindlist(med_fixed_list, idcol = "Model")
med_x2 <-  med_x %>% 
  select(Variable, Model, Median = '0.5quant', 
         Low = '0.025quant', High = '0.975quant') %>% 
  mutate(Model = paste("Mesh", Model, sep = "_")) 

(fixed_med <- ggplot(med_x2) + 
  geom_point(aes(x = Variable, y = Median, col = Model), 
             position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = Low, ymax = High, 
                    x = Variable, col = Model), 
                position=position_dodge(width=0.5), 
                width = 0.2) + 
  scale_color_brewer(type = "qual", palette = "Set1") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() + 
  labs(title = "20/2 km edge") +
  NULL)

# 6. Create small meshes ####

small_meshes <- list()
small_ipoints <- list()

for(i in seq_along(meshes)) {
  small_meshes[[i]] <-  fmesher:::fm_subdivide(med_meshes[[i]])
  small_ipoints[[i]] <- fm_int(small_meshes[[i]])
} 

ips <- ggplot() + 
  geom_spatraster(data = covar_scaled, aes(fill = Elevation), alpha = 0.5) + 
  gg(small_meshes[[1]]) + 
  scale_fill_viridis_c(na.value = NA) + 
  geom_sf(data = small_ipoints[[1]], col = "red")  +
  geom_sf(data = small_ipoints[[2]], col = "darkgreen") +
  geom_sf(data = small_ipoints[[3]], col = "orange") + 
  geom_sf(data = small_ipoints[[4]], col = "purple") +
  coord_sf(xlim = st_bbox(covar_scaled)[c(1,3)], 
           ylim = st_bbox(covar_scaled)[c(2,4)], 
           expand = F) + 
  labs(title = "20/4 km edge") +
  theme_bw()

Cairo::CairoPDF(file = "mesh_outputs/ipoints.pdf", width = 18, height = 10)
ipb + ipm + ips + plot_layout(ncol = 3)
dev.off()

# 7. Model with small mesh ####

small_model_list <- list()
small_fixed_list <- list()

for(i in 1:4) {
  print(i)
  mesh <- small_meshes[[i]]
  boundary <- bound_list[[i]] 
  
  ## 7.1 Formula ####
  form1 <- geometry ~  Intercept(1)  +
    Eff.elevation(covar_scaled$Elevation, model = "linear") +
    Eff.slope(covar_scaled$Slope, model = "linear") +
    Eff.tree_cover(covar_scaled$`Tree cover density`, model = "linear") + 
    Eff.swf(covar_scaled$`Small woody features`, model = "linear") +
    NULL
  
  ## 7.2 Run model and summary ####
  
  m1 <- lgcp(components = form1,
             data = PO_data_sel,
             # samplers = boundary,
             domain =  list(geometry = mesh))
  
  small_model_list[i] <- list(m1)
  
  fixed <-  tibble::rownames_to_column(m1$summary.fixed, "Variable")
  small_fixed_list[i] <- list(fixed)

}  


## 7.4 Plot fixed effects ####

small_x <- data.table::rbindlist(small_fixed_list, idcol = "Model")
small_x2 <-  small_x %>% 
  select(Variable, Model, Median = '0.5quant', 
         Low = '0.025quant', High = '0.975quant') %>% 
  mutate(Model = paste("Mesh", Model, sep = "_")) 

(fixed_small <- ggplot(small_x2) + 
  geom_point(aes(x = Variable, y = Median, col = Model), 
             position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = Low, ymax = High, 
                    x = Variable, col = Model), 
                position=position_dodge(width=0.5), 
                width = 0.2) + 
  scale_color_brewer(type = "qual", palette = "Set1") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() + 
  labs(title = "20/4 km edge") +
  # facet_wrap(~Variable, scales = "free") + 
  NULL)

# 8. Create tiny meshes ####

tiny_meshes <- list()
tiny_ipoints <- list()

for(i in seq_along(meshes)) {
  tiny_meshes[[i]] <-  fmesher:::fm_subdivide(small_meshes[[i]])
  tiny_ipoints[[i]] <- fm_int(tiny_meshes[[i]])
} 

ipt <- ggplot() + 
  geom_spatraster(data = covar_scaled, aes(fill = Elevation), alpha = 0.5) + 
  gg(tiny_meshes[[1]]) + 
  scale_fill_viridis_c(na.value = NA) + 
  geom_sf(data = tiny_ipoints[[1]], col = "red")  +
  geom_sf(data = tiny_ipoints[[2]], col = "darkgreen") +
  geom_sf(data = tiny_ipoints[[3]], col = "orange") + 
  geom_sf(data = tiny_ipoints[[4]], col = "purple") +
  coord_sf(xlim = st_bbox(covar_scaled)[c(1,3)], 
           ylim = st_bbox(covar_scaled)[c(2,4)], 
           expand = F) + 
  labs(title = "20/4 km edge") +
  theme_bw()

Cairo::CairoPDF(file = "mesh_outputs/ipoints.pdf", width = 18, height = 16)
ipb + ipm + ips + ipt + plot_layout(ncol = 2)
dev.off()


# 9. Model with tiny mesh ####

tiny_model_list <- list()
tiny_fixed_list <- list()

for(i in 1:4) {
  print(i)
  mesh <- tiny_meshes[[i]]
  boundary <- bound_list[[i]] 
  
  ## 9.1 Formula ####
  form1 <- geometry ~  Intercept(1)  +
    Eff.elevation(covar_scaled$Elevation, model = "linear") +
    Eff.slope(covar_scaled$Slope, model = "linear") +
    Eff.tree_cover(covar_scaled$`Tree cover density`, model = "linear") + 
    Eff.swf(covar_scaled$`Small woody features`, model = "linear") +
    NULL
  
  ## 7.2 Run model and summary ####
  
  m1 <- lgcp(components = form1,
             data = PO_data_sel,
             # samplers = boundary,
             domain =  list(geometry = mesh))
  
  tiny_model_list[i] <- list(m1)
  
  fixed <-  tibble::rownames_to_column(m1$summary.fixed, "Variable")
  tiny_fixed_list[i] <- list(fixed)
  
}  


## 7.4 Plot fixed effects ####

tiny_x <- data.table::rbindlist(tiny_fixed_list, idcol = "Model")
tiny_x2 <-  tiny_x %>% 
  select(Variable, Model, Median = '0.5quant', 
         Low = '0.025quant', High = '0.975quant') %>% 
  mutate(Model = paste("Mesh", Model, sep = "_")) 

(fixed_tiny <- ggplot(tiny_x2) + 
    geom_point(aes(x = Variable, y = Median, col = Model), 
               position=position_dodge(width=0.5)) + 
    geom_errorbar(aes(ymin = Low, ymax = High, 
                      x = Variable, col = Model), 
                  position=position_dodge(width=0.5), 
                  width = 0.2) + 
    scale_color_brewer(type = "qual", palette = "Set1") +
    geom_hline(yintercept = 0, linetype = 2) +
    theme_bw() + 
    labs(title = "20/8 km edge") +
    # facet_wrap(~Variable, scales = "free") + 
    NULL)



# 8. plot all fixed ####

Cairo::CairoPDF(file = "mesh_outputs/linear_effects.pdf", width = 16, height = 8)
fixed_big + fixed_med + fixed_small + fixed_tiny + plot_layout(ncol = 2)  
dev.off()

