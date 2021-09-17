# aim of this script is to extract environmental covariates at the study area level

# libraries ----
install.packages("pacman")
source("functions_covar_extraction.R")
install.packages("pacman")
pacman::p_load(sp, sdmvspecies, tidyverse, rgdal, maptools, rgeos, rJava, dismo, 
               sf, raster, mapview, virtualspecies, spdplyr, ENMeval, lubridate, 
               amt, osmextract, egg, spatialEco, viridis, pool, dbplyr, DBI, 
               RPostgreSQL, dbplyr, dbplot, rpostgis)

# data import ---
# sas = study areas
#sas <- readOGR("D:/EUROBOAR/PAGERS/PAGER_pig_farm_interaction/GIS",
#              layer="study_areas_210212", stringsAsFactors = F)

ireland <- st_read("~/ISDM/data/ireland_ITM.shp") %>% 
  st_transform(st_crs("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=0.99982 +x_0=600000 +y_0=750000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"))  # IRENET in km

ireland_utm <- ireland %>% 
  st_transform("EPSG:3035")

ireland_bbox <- st_sf(st_envelope(ireland))
class(ireland_bbox)

sas_sf <- ireland_bbox
sas52 <- st_transform(sas_sf, "EPSG:3035")

# Landcover ----
# res 100m 
landcov <- raster("copernicus/corine/DATA/U2018_CLC2018_V2020_20u1.tif")
names(landcov) <- "landcover"
# list of habitat class
df_cov <- data.frame(ID=c(1:44,48), landcover=c("Continuous	urban	fabric",
                                                "Discontinuous	urban	fabric",
                                                "Industrial	or	commercial	units",
                                                "Road	and	rail	networks	and	associated	land",
                                                "Port	areas",
                                                "Airports",#
                                                "Mineral	extraction	sites",
                                                "Dump	sites",
                                                "Construction	sites",
                                                "Green	urban	areas",
                                                "Sport	and	leisure	facilities",
                                                "Non-irrigated	arable	land",
                                                "Permanently	irrigated	land",
                                                "Rice	fields",
                                                "Vineyards",
                                                "Fruit	trees	and	berry	plantations",
                                                "Olive	groves",
                                                "Pastures",
                                                "Annual	crops	associated	with	permanent	crops",
                                                "Complex	cultivation	patterns",
                                                "Land	principally	occupied	by	agriculture	with	significant	areas	of	natural	vegetation",
                                                "Agro-forestry	areas",
                                                "Broad-leaved	forest",
                                                "Coniferous	forest",
                                                "Mixed	forest",
                                                "Natural	grasslands",
                                                "Moors	and	heathland",
                                                "Sclerophyllous	vegetation",
                                                "Transitional	woodland-shrub",
                                                "Beaches	-	dunes	-	sands",
                                                "Bare	rocks",
                                                "Sparsely vegetated	areas",
                                                "Burnt	areas",
                                                "Glaciers	and	perpetual	snow",
                                                "Inland	marshes",
                                                "Peat	bogs",
                                                "Salt	marshes",
                                                "Salines",
                                                "Intertidal	flats",
                                                "Water	courses",
                                                "Water	bodies",
                                                "Coastal	lagoons",
                                                "Estuaries",
                                                "Sea	and	ocean",
                                                "NODATA"))
library(janitor)
df_cov$landcover <- make_clean_names(df_cov$landcover)
write.csv(df_cov, "copernicus/corine/legend_corine.csv")


# Tree cover density ----
# res 100m: to be treated as CONTINUOUS, or eventually classes TCD_2015_020m_eu_03035_d04_E40N20/TCD_2015_020m_eu_03035_d04_E40N20
tcd <- raster("copernicus/tree_cover_density/TCD_2015_100m_eu_03035_d04_full.tif")
names(tcd) <- "tcd"
#tcd <- raster::crop(tcd, spdf52)
#plot(tcd)

# Forest type ----
# non-forest (0), broadleaf (1), coniferous (2), unclass (>2)
fortype <- raster("copernicus/forest_type/FTY_2015_100m_eu_03035_d02_full.tif")
names(fortype) <- "forest_type"
# plot(fortype)

# Dominant leaf type ----
dlt <- raster("copernicus/dominant_leaf_type/DLT_2015_020m_eu_03035_d04_full.tif")
names(dlt) <- "dominant_leaf"

# Small woody cover ----
swf <- raster("copernicus/small_woody_forest/Data/swfawf_2015_100m_EU_3035_V013.tif")
names(swf) <- "small_woody_forest"
# plot(swf)

# DEM ----
# res 100m from https://zenodo.org/record/4057883#.YMH4fkxCRPY
dem <- raster("Europe_Digital_Terrain_Model/dtm_elev.lowestmode_gedi.eml_m_100m_0..0cm_2000..2018_eumap_epsg3035_v0.1.tif")

# human foot print ----
hfp <- raster("human_foot_print/HFP2009.tif")


# loop ----
# out <- list()
# for (i in  unique(sas52$study_name)) { 
#  
#   sas_i <- sas52 %>% filter(study_name==i)
sas_i <- sas52
dForest <- NULL; dEdges<- NULL; dUrban<- NULL;tcd_j<- NULL; dlt_j<- NULL;fortype_j<- NULL; swf_j<- NULL; dRoads<- NULL; dPaths<- NULL; dem_j<- NULL; slope<- NULL; hli<- NULL;
twi<- NULL; hfp_j<- NULL

# landcover
landcov_i <- raster::crop(landcov, sas_i)
# landcov_i <- as.factor(landcov_i)
# rat <- levels(landcov_i)[[1]]
# rat <- left_join(rat, df_cov, by="ID")
# levels(landcov_i) <- rat

forest <- landcov_i==23|landcov_i==24|landcov_i==25
forest[forest == 0] <- NA
names(forest) <- "forest"

edges <- forest
edges[is.na(edges)] <- 0
edges[edges > 0] <- NA
edges[edges == 0] <- 1
names(edges) <- "edges"

urban <- landcov_i==1|landcov_i==2|landcov_i==3
urban[urban == 0] <- NA
names(urban) <- "urban"
# distance forest
dForest <- distance(forest)
names(dForest) <- "distance_forest"
#distance_edge
dEdges <- distance(edges)
names(dEdges) <- "distance_edges"
#distance_urban
true_false <- data.frame(table(!is.na(values(urban))))
if (length(true_false$Var1) > 1) {
  dUrban <- distance(urban)
  names(dUrban) <- "distance_urban"
}

# tree cover density
tcd_i <- raster::crop(tcd, sas_i)

# dominant leaf type
dlt_i <- raster::crop(dlt, sas_i)
dlt_i <- resample(dlt_i, tcd_i, method='ngb')
# forest type
fortype_i <- raster::crop(fortype, sas_i)
# small woody
swf_i <- raster::crop(swf, sas_i)

# Linear infrastructure
osm_roads <- NULL
osm_paths <- NULL

# sas84 <- sas_sf %>% filter(study_name==i)
sas84 <- ireland_bbox
centroid <- st_centroid(sas84)
its_osm_fr <- oe_match(centroid, provider = "openstreetmap_fr")
its_geof <- oe_match(centroid, provider = "geofabrik")
its_details <- ifelse(its_osm_fr$file_size > its_geof$file_size, its_geof[[1]], its_osm_fr[[1]])

its_pbf <- oe_download(
  #file_url= "https://download.geofabrik.de/europe/spain-latest.osm.pbf",
  file_url = its_details, 
  #file_size = its_details$file_size,
  provider = ifelse(its_osm_fr$file_size > its_geof$file_size, "geofabrik", "openstreetmap_fr"),
  download_directory = "osm",
  force_download = F
)
its_gpkg = oe_vectortranslate(its_pbf)
#+ mapview(osm_roads) + mapview(sas84)

osm_roads = oe_get(centroid, stringsAsFactors = FALSE, quiet = TRUE,
                   provider = ifelse(its_osm_fr$file_size > its_geof$file_size, "geofabrik", "openstreetmap_fr"),
                   download_directory="osm",
                   query = "SELECT * FROM 'lines' WHERE highway IN ('motorway','trunk','primary',
                   'secondary','tertiary','unclassified','residential','service')")
osm_roads <- st_crop(osm_roads, sas84)

if (nrow(osm_roads) > 0) {
  osm_roads <-  st_transform(osm_roads,3035)
  
  roadsr <- rasterize(osm_roads, dForest)
  dRoads <- distance(roadsr)
  names(dRoads) <- "distance_roads"
}

osm_paths = oe_get(centroid, stringsAsFactors = FALSE, quiet = TRUE,
                   provider = ifelse(its_osm_fr$file_size > its_geof$file_size, "geofabrik", "openstreetmap_fr"),
                   download_directory="osm",
                   query = "SELECT * FROM 'lines' WHERE highway IN ('path','track','footway')")
osm_paths <- st_crop(osm_paths, sas84)
if (nrow(osm_paths) > 0) {
  osm_paths <-  st_transform(osm_paths,3035)
  
  pathsr <- rasterize(osm_paths, dForest)
  dPaths <- distance(pathsr)
  names(dPaths) <- "distance_paths"
}

# DEM
dem_i <- raster::crop(dem, sas_i)
dem_i <- resample(dem_i, dForest, method='ngb')
names(dem_i) <- "elevation"
# Slope 
slope <- terrain(dem_i, "slope")
names(slope) <- "slope"
# Aspect
aspect <- terrain(dem_i, "aspect")
names(aspect) <- "aspect"
# heat load index
hli <- spatialEco::hli(dem_i)
names(hli) <- "hli"
# Topographic Wetness Index
topo <- create_layers(dem_i, fill.sinks = TRUE, deg = .1)
twi <- subset(topo, "atb")
names(twi) <- "twi"

# human foot print index
hfp_i <- raster::crop(hfp, sas_i)
hfp_i <- resample(hfp_i, dForest, method='ngb')
names(hfp_i) <- "human_print"

stack_i <- stack(landcov_i,dForest, dEdges, dUrban,tcd_i, dlt_i,fortype_i, swf_i, dRoads, dPaths, dem_i, slope, hli, twi, hfp_i)
writeRaster(stack_i, paste("derived_data",i, sep=""), overwrite=T)





# road density (if needed) ----
library(spatstat)
library(maptools)
pspSl <- as.psp(osm_roads$geometry)
# Pixellate with resolution of 0.5, i.e. 2x2 pixels
W <- as(pext, "owin")
px <- pixellate(pspSl, W=W, eps=1000)
# This can be converted to raster as desired
road_dens <- raster(px)
road_dens <- resample(road_dens, r, method='ngb')
#road.dens <- mask(road.dens, gcp)
names(road_dens) <- "density_roads"
plot(road_dens)

# Path density
pspSl <- as.psp(osm_paths$geometry)
# Pixellate with resolution of 0.5, i.e. 2x2 pixels
px <- pixellate(pspSl, W=W, eps=1000)
# This can be converted to raster as desired
path_dens <- raster(px)
path_dens <- resample(path_dens, r, method='ngb')
#path.dens <- mask(path.dens, gcp)
names(path_dens) <- "density_paths"
plot(path_dens)








