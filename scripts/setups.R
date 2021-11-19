# install all necessary packages # 
devtools::install_github('PhilipMostert/inlabruSDMs', force = T)
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
library(inlabruSDMs)
'%!in%' <- function(x,y)!('%in%'(x,y))


