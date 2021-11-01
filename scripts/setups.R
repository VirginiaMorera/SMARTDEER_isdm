# install all necessary packages # 
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(inlabru)
library(INLA)
library(rasterVis)
library(ggcorrplot)
'%!in%' <- function(x,y)!('%in%'(x,y))

ITM_km <- CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=0.99982 +x_0=600000 +y_0=750000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")

