source("scripts/setups.R")

ded <- st_read("kk_ty_ls.shp")
st_crs(ded) <- st_crs(29902)



counties <- ded %>% 
  group_by(COUNTY) %>% 
  summarise()
st_crs(counties) <- st_crs(29902)

ggplot() + 
  geom_sf(data = ded, col = "darkgray") +
  geom_sf(data = counties, col = "black", fill = NA) + 
  theme_bw()


# SMARTDEER rasters ####
fallow <- terra::rast("fallow_deer.grd")
sika <- terra::rast("sika_deer.grd")
red <- terra::rast("red_deer.grd")

fallow_IG <- project(fallow, crs("EPSG:29902"))
fallow_mean <- rescale0to1(fallow_IG$mean)

fallow_mask <- mask(crop(fallow_mean, counties), counties)

ggplot() +
  geom_spatraster(data = fallow_mask, aes(fill = mean)) +
  geom_sf(data = ded, fill = NA, col = "darkgray") + 
  geom_sf(data = counties, col = "black", fill = NA) +
  scale_fill_viridis_c(na.value = NA) + 
  theme_bw() + 
  ggtitle("Fallow deer")

sika_IG <- project(sika, crs("EPSG:29902"))
sika_mean <- rescale0to1(sika_IG$mean)

sika_mask <- mask(crop(sika_mean, counties), counties)

ggplot() +
  geom_spatraster(data = sika_mask, aes(fill = mean)) +
  geom_sf(data = ded, fill = NA, col = "darkgray") + 
  geom_sf(data = counties, col = "black", fill = NA) +
  scale_fill_viridis_c(na.value = NA) + 
  theme_bw() + 
  ggtitle("Sika deer")

red_IG <- project(red, crs("EPSG:29902"))
red_mean <- rescale0to1(red_IG$mean)

red_mask <- mask(crop(red_mean, counties), counties)

ggplot() +
  geom_spatraster(data = red_mask, aes(fill = mean)) +
  geom_sf(data = ded, fill = NA, col = "darkgray") + 
  geom_sf(data = counties, col = "black", fill = NA) +
  scale_fill_viridis_c(na.value = NA) + 
  theme_bw() + 
  ggtitle("Red deer")

# disag rasters ####

fallow18 <- terra::rast("Fallow18_prediction.tif")
crs(fallow18) <- "EPSG:2157"
sika18 <- terra::rast("Sika18_prediction.tif")
crs(sika18) <- "EPSG:2157"
red18 <- terra::rast("Red18_prediction.tif")
crs(red18) <- "EPSG:2157"

fallow_IG <- project(fallow18, crs("EPSG:29902"))

fallow_mask <- mask(crop(fallow_IG, counties), counties)

ggplot() +
  geom_spatraster(data = fallow_mask, aes(fill = Fallow18_Prediction)) +
  geom_sf(data = ded, fill = NA, col = "darkgray") + 
  geom_sf(data = counties, col = "white", fill = NA) +
  scale_fill_viridis_c(na.value = NA) + 
  theme_bw() + 
  ggtitle("Culled fallow deer in 2018")

sika_IG <- project(sika18, crs("EPSG:29902"))

sika_mask <- mask(crop(sika_IG, counties), counties)

ggplot() +
  geom_spatraster(data = sika_mask, aes(fill = Sika18_Prediction)) +
  geom_sf(data = ded, fill = NA, col = "darkgray") + 
  geom_sf(data = counties, col = "white", fill = NA) +
  scale_fill_viridis_c(na.value = NA) + 
  theme_bw() + 
  ggtitle("Culled sika deer in 2018")

red_IG <- project(red18, crs("EPSG:29902"))

red_mask <- mask(crop(red_IG, counties), counties)

ggplot() +
  geom_spatraster(data = red_mask, aes(fill = Red18_Prediction )) +
  geom_sf(data = ded, fill = NA, col = "darkgray") + 
  geom_sf(data = counties, col = "white", fill = NA) +
  scale_fill_viridis_c(na.value = NA) + 
  theme_bw() + 
  ggtitle("Culled red deer in 2018")



# calculate values per county and per total area

library(exactextractr)

counties <- counties %>% 
  mutate(
    sum_red = exact_extract(x = red_mask, counties, 'sum'),
    mean_red = exact_extract(x = red_mask, counties, 'mean'),
    min_red = exact_extract(x = red_mask, counties, 'min'),
    max_red = exact_extract(x = red_mask, counties, 'max'),
    sum_sika = exact_extract(x = sika_mask, counties, 'sum'),
    mean_sika = exact_extract(x = sika_mask, counties, 'mean'),
    min_sika = exact_extract(x = sika_mask, counties, 'min'),
    max_sika = exact_extract(x = sika_mask, counties, 'max'),
    sum_fallow = exact_extract(x = fallow_mask, counties, 'sum'),
    mean_fallow = exact_extract(x = fallow_mask, counties, 'mean'),
    min_fallow = exact_extract(x = fallow_mask, counties, 'min'),
    max_fallow = exact_extract(x = fallow_mask, counties, 'max')
  )

library(patchwork)

# Fallow 
ggplot(counties) +
  geom_sf(aes(fill = sum_fallow)) +
  scale_fill_viridis_c() + 
  labs(title = "Total culled animals", fill = "Fallow deer") +
  theme_bw() + 
  
  ggplot(counties) +
  geom_sf(aes(fill = mean_fallow)) +
  scale_fill_viridis_c() + 
  labs(title = "Average per 25 sqkm", fill = "Fallow deer") +
  theme_bw() + 
  
  ggplot(counties) +
  geom_sf(aes(fill = min_fallow)) +
  scale_fill_viridis_c() + 
  labs(title = "Minimum per 25 sqkm", fill = "Fallow of deer") +
  theme_bw() + 
  
  ggplot(counties) +
  geom_sf(aes(fill = max_fallow)) +
  scale_fill_viridis_c() + 
  labs(title = "Maximum per 25 sqkm", fill = "Fallow deer") +
  theme_bw() + 
  
  plot_annotation(title = 'Fallow deer culling numbers', 
                  subtitle = "2894 deer culled in total in the area in 2018")

# Sika 
ggplot(counties) +
  geom_sf(aes(fill = sum_sika)) +
  scale_fill_viridis_c() + 
  labs(title = "Total culled animals", fill = "sika deer") +
  theme_bw() + 
  
  ggplot(counties) +
  geom_sf(aes(fill = mean_sika)) +
  scale_fill_viridis_c() + 
  labs(title = "Average per 25 sqkm", fill = "sika deer") +
  theme_bw() + 
  
  ggplot(counties) +
  geom_sf(aes(fill = min_sika)) +
  scale_fill_viridis_c() + 
  labs(title = "Minimum per 25 sqkm", fill = "sika of deer") +
  theme_bw() + 
  
  ggplot(counties) +
  geom_sf(aes(fill = max_sika)) +
  scale_fill_viridis_c() + 
  labs(title = "Maximum per 25 sqkm", fill = "sika deer") +
  theme_bw() + 
  
  plot_annotation(title = 'Sika deer culling numbers', 
                  subtitle = "507 deer culled in total in the area in 2018")



# Red 
ggplot(counties) +
  geom_sf(aes(fill = sum_red)) +
  scale_fill_viridis_c() + 
  labs(title = "Total culled animals", fill = "Red deer") +
  theme_bw() + 
  
ggplot(counties) +
  geom_sf(aes(fill = mean_red)) +
  scale_fill_viridis_c() + 
  labs(title = "Average per 25 sqkm", fill = "Red deer") +
  theme_bw() + 
  
ggplot(counties) +
  geom_sf(aes(fill = min_red)) +
  scale_fill_viridis_c() + 
  labs(title = "Minimum per 25 sqkm", fill = "Red deer") +
  theme_bw() + 
  
ggplot(counties) +
  geom_sf(aes(fill = max_red)) +
  scale_fill_viridis_c() + 
  labs(title = "Maximum per 25 sqkm", fill = "Red deer") +
  theme_bw() + 

plot_annotation(title = 'Red deer culling numbers', 
                subtitle = "429 deer culled in total in the area in 2018")



