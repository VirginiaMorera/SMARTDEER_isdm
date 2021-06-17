##---------------##
#### Deer data ####
##---------------##

pre_data <- read.csv("data/all_data.csv", row.names = NULL)

library(tidyverse)
library(sf)
library(sp)
library(raster)

sel_data <- pre_data %>% 
  dplyr::select(Year, Species, Deer.Presence, X = Longitude, Y = Latitude, Source)

PA_data <- sel_data %>% 
  filter(Source %in% c("Coillte_density", "Coillte_desk")) %>% 
  dplyr::select(-Source)  %>% 
  st_as_sf(coords = c("X", "Y")) %>% 
  st_set_crs(st_crs("+proj=tmerc +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")) %>%  # IRENET in m
  st_transform(st_crs("+proj=tmerc +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")) %>%  # IRENET in km
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                Y = sf::st_coordinates(.)[,2], 
                Deer.Presence = if_else(Deer.Presence == "Yes", 1, 0)) %>% 
  rename(PA = Deer.Presence) %>% 
  st_set_geometry(NULL)

write.csv(PA_data, file = "data/PA_data.csv", row.names = F)

PO_data <- sel_data %>% 
  filter(Source %in% c("NBDC", "webSurvey")) %>% 
  dplyr::select(-Source) %>% 
  st_as_sf(coords = c("X", "Y")) %>% 
  st_set_crs(st_crs("+proj=tmerc +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")) %>%  # IRENET in m
  st_transform(st_crs("+proj=tmerc +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")) %>%  # IRENET in km
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                Y = sf::st_coordinates(.)[,2], 
                Deer.Presence = if_else(Deer.Presence == "Yes", 1, 0)) %>% 
  rename(PO = Deer.Presence) %>% 
  st_set_geometry(NULL)

write.csv(PO_data, file = "data/PO_data.csv", row.names = F)

##------------------------##
#### Environmental data ####
##------------------------##

# Obtain integration points (mesh nodes) to extract covariate values at those points
mesh <- readRDS("data/mesh.RDS") 

ipoints <- data.frame(X = mesh$loc[,1], 
                      Y = mesh$loc[,2])

ipoints_sp <- SpatialPoints(coords = mesh$loc, proj4string = CRS("+proj=tmerc +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"))

# Load outer boundary to crop the rasters
bound <- readRDS("data/outer_boundary.RDS")

# Human footprint index 2009

HFI <- raster("large_env_data/wildareas-v3-2009-human-footprint.tif")
bound_temp <- spTransform(bound, CRSobj = HFI@crs)
HFI_crop <- crop(HFI, bound_temp)
HFI_crop <- mask(HFI_crop, bound_temp)
HFI_crop[HFI_crop > 100] <- NA
plot(HFI_crop)
writeRaster(HFI_crop, "large_env_data/HFI_crop.tif", format = "GTiff")

# Digital elevation model

DEM_files <- list.files(path = "large_env_data/DEM", pattern = ".TIF$", full.names = T)
DEM <- mosaic(raster(DEM_files[1]), raster(DEM_files[2]), fun = "mean")

bound_temp <- spTransform(bound, CRSobj = DEM@crs)
DEM_crop <- crop(DEM, bound_temp)
DEM_crop <- mask(DEM_crop, bound_temp)
plot(DEM_crop)

writeRaster(DEM_crop, "large_env_data/DEM_crop.tif", format = "GTiff")

##--------------------------##
#### Extract covar values ####
##--------------------------##

# put all rasters in a stack using custom function 
# (this takes ages and it's better done in the RStudio server or Sonic cluster)
cropped_rasters <- list.files(path = "large_env_data", pattern = "crop.tif$", full.names = T)
covar_list <- map(cropped_rasters, raster)

source("scripts/aux_functions.R")

covar_stack <- list_to_stack(covar_list, new_res = c(1, 1), 
                             dest_crs = CRS("+proj=tmerc +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"))

saveRDS(covar_stack, file = "data/covar_stack.RDS") # we've actually cheated and done this in the RStudio server, saved it, and we load it here

covar_stack <- readRDS("data/covar_stack.RDS")

# extract values and turn into SpatialPointsDataFrame
covar_values <- extract(covar_stack, ipoints_sp, sp = TRUE)
