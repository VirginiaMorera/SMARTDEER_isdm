# install all necessary packages # 
# devtools::install_github('PhilipMostert/inlabruSDMs', force = T)
# install.packages("INLA",repos=c(getOption("repos"),
#                                 INLA="https://inla.r-inla-download.org/R/stable"), 
#                  dep=TRUE)
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(inlabru)
library(INLA)
library(rasterVis)
library(ggcorrplot)
library(ggthemes)
library(cowplot)
library(spatstat)
library(PointedSDMs)
library(rasterVis)
library(viridis)
library(INLA)
library(ggpmisc)
library(climateStability)
library(ggrepel)
library(ggthemes)
library(beepr)
<<<<<<< HEAD
library(rstatix)
library(magrittr)
=======
>>>>>>> bb9d3ef22acdf933b04a11b776b83bf2a7a5fb15


projKM <- "+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=0.99982 +x_0=600000 +y_0=750000
           +ellps=GRS80 +units=km +no_defs"

`%!in%` = Negate(`%in%`)


