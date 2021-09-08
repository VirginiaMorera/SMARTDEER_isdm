# install all necessary packages # 
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(inlabru)
library(INLA)

'%!in%' <- function(x,y)!('%in%'(x,y))
