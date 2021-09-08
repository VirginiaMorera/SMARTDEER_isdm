in_bound <- readRDS("data/inner_boundary.RDS")
mesh0 <- readRDS("data/mesh.RDS")
covar_stack <- readRDS("data/covar_stack.RDS")
spatialcovariates <- as(covar_stack$HFI_crop, "SpatialPixelsDataFrame")
source("scripts/params_bru_sdm.R")
source("scripts/bru_sdm.R")

PO_data <- read.csv("data/PO_data_RD.csv", row.names = NULL)
PA_data <- read.csv("data/PA_data_RD.csv", row.names = NULL)

red_deer_data <- organize_data(PO_data, PA_data, poresp = "PO", paresp = "PA",
                               trialname = NULL, coords = c("X", "Y"), proj = ITM,
                               marks = FALSE, inclmarks = NULL,
                               markfamily = 'gaussian', timevariable = NULL,
                               ips = NULL, mesh = mesh0,
                               meshpars = NULL, boundary = in_bound)



mdl1 <- bru_sdm(PO_data, PA_data, spatialcovariates = NULL, marks = FALSE, markfamily = 'gaussian',
                inclmarks = NULL, coords = c('X','Y'), poresp = "PO", paresp = "PA",
                trialname = NULL, inclcoords = FALSE, mesh = mesh, meshpars = NULL,
                spdemodel = NULL, ips = NULL, bdry = bdry,
                proj = CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=0.99982 +x_0=600000 +y_0=750000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"),
                predictions = FALSE,
                residuals = NULL, intercept = FALSE, indivintercepts = TRUE,
                pointsspatial = TRUE, marksspatial = FALSE, options = list(),
                poformula = NULL, paformula = NULL, tol = 0)


mdl2 <- bru_sdm(PO_data, PA_data, spatialcovariates = NULL, marks = FALSE, markfamily = 'gaussian',
                inclmarks = NULL, coords = c('X','Y'), poresp = "PO", paresp = "PA",
                trialname = NULL, inclcoords = FALSE, mesh = mesh, meshpars = NULL,
                spdemodel = spdemodel_pc, ips = NULL, bdry = bdry,
                proj = CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=0.99982 +x_0=600000 +y_0=750000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"),
                predictions = FALSE,
                residuals = NULL, intercept = FALSE, indivintercepts = TRUE,
                pointsspatial = TRUE, marksspatial = FALSE, options = list(),
                poformula = NULL, paformula = NULL, tol = 0)



source("scripts/outputs.R")
summary.bru_sdm(mdl1)
summary.bru_sdm(mdl2)


mdl3 <- bru_sdm(PO_data_sub, PA_data_sub, spatialcovariates = spatialcovariates, marks = FALSE, markfamily = 'gaussian',
                inclmarks = NULL, coords = c('X','Y'), poresp = "PO", paresp = "PA",
                trialname = NULL, inclcoords = FALSE, mesh = mesh, meshpars = NULL,
                spdemodel = spdemodel_pc, ips = NULL, bdry = bdry,
                proj = CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=0.99982 +x_0=600000 +y_0=750000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"),
                predictions = FALSE,
                residuals = NULL, intercept = FALSE, indivintercepts = TRUE,
                pointsspatial = TRUE, marksspatial = FALSE, options = list(),
                poformula = NULL, paformula = NULL, tol = 0.9)
