in_bound <- readRDS("data/inner_boundary.RDS")
mesh0 <- readRDS("data/mesh.RDS")
covar_stack <- stack("large_env_data/covar_subset_ITM.gri")
spatialcovariates <- as(covar_stack, "SpatialPixelsDataFrame")
red_deer_data <- readRDS("data/data_for_model.RDS")

spdemodel_pc = INLA::inla.spde2.pcmatern(mesh = mesh0, alpha = 3/2,  ### mesh and smoothness parameter
                                         prior.range = c(20, 0.01), ### P(practic.range < 20) Very small probability the range is smaller than 20 (with high likelihood the range is between 20 (res of the model) and infinity i.e. flat prior) 
                                         prior.sigma = c(0.2, 0.01)) ### P(sigma > 0.2) Very small probability that sigma is larger than 0.2. Sigma is most likely between 0 and 0.2



mdl1 <- bru_sdm(data = red_deer_data, 
                spatialcovariates, 
                covariatestoinclude = c("tree_cover_density", "elevation", "slope", "human_footprint_index", "landCover"),
                spdemodel = spdemodel_pc,
                sharedspatial = T)

predict(mdl1)
