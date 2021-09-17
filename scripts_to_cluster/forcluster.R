
components_joint <- readRDS("components_joint.RDS")
likelihoods <- readRDS("likelihoods.RDS")
options <- readRDS("options.RDS")

install.packages("inlabru")
library("inlabru")

##---------------------------------------##
#### Load objects needed for modelling ####
##---------------------------------------##

in_bound <- readRDS("data/inner_boundary.RDS")
mesh0 <- readRDS("data/mesh.RDS")
covar_stack <- readRDS("data/covar_stack.RDS")
PO_data <- read.csv("data/PO_data_RD.csv", row.names = NULL)
PA_data <- read.csv("data/PA_data_RD.csv", row.names = NULL)

library(inlabru)
library(INLA)

##------------------------------##
#### source necessary scripts ####
##------------------------------##
source("scripts/params_bru_sdm.R")

source("scripts/bru_sdm.R")

##------------------##
#### run function ####
##------------------##
model_out <-  bru_sdm(spatialcovariates = spatialcovs, marks = FALSE, markfamily = 'gaussian',
                      inclmarks = NULL, coords = c('X','Y'), poresp = 'PO', paresp = 'PA',
                      trialname= NULL, inclcoords = FALSE, mesh = mesh0, 
                      meshpars = NULL, ips = NULL, bdry = bdry, 
                      proj = CRS("+proj=tmerc +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"),
                      predictions = TRUE, residuals = 'response', intercept = FALSE,
                      indivintercepts = TRUE, pointsspatial = TRUE, marksspatial = FALSE, 
                      options = list(), poformula = NULL, paformula = NULL, tol = 0, 
                      PO_data, PA_data) # other specifications

