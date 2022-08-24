rm(list = ls())
source("scripts/setups.R")

# Model validation
## To try and validate the results of the integrated species distribution models, 
## we're going to use the density data based on fecal sampling by Tim Burkitt, and the estimates of low, moderate and high density collected during Coillte's desk surveys. We'll extract the value of the prediction at the centroid of each sampled property and, using linear models, we'll see how well the ISDM predictions can predict the observed values. 

## Load and modify data ####

# load data

## map
ireland <- st_read("data/ireland_ITM.shp")

ireland_proj <- st_transform(ireland, projKM)

## tim burkitt densities
dens_data <- readRDS("data/coillte_dens_surveys_area_and_point.RDS")

dens_latlon <- dens_data %>%
  st_set_geometry(NULL) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(ireland)) %>% 
  st_transform(projKM)

dens_area <- dens_data %>% 
  dplyr::select(-Longitude, -Latitude) %>% 
  st_transform(projKM)

## coillte surveys
coillte_data <- readRDS("data/PA_data_all.RDS")

coillte_data <- coillte_data %>% 
  dplyr::select(Red, Fallow, Sika, Year) %>% 
  st_transform(projKM) %>% 
  st_centroid()

## modelled densities
RD_pred_resp <- readRDS("server_outputs/RedDeer_prediction.RDS")
SD_pred_resp <- readRDS("server_outputs/SikaDeer_prediction.RDS")
FD_pred_resp <- readRDS("server_outputs/FallowDeer_prediction.RDS")

RD_pred_scaled <- rescale0to1(raster(RD_pred_resp$predictions["mean"]))
SD_pred_scaled <- rescale0to1(raster(SD_pred_resp$predictions["mean"]))
FD_pred_scaled <- rescale0to1(raster(FD_pred_resp$predictions["mean"]))

preds <- stack(RD_pred_scaled,
               SD_pred_scaled,
               FD_pred_scaled)

names(preds) <- c("Red", "Sika", "Fallow")

## culling returns
culling_data <- read.csv("data/culling_data.csv", row.names = NULL)

culling_data <- culling_data %>% 
  dplyr::select(Year = SeasonStart, County, Species, Deer_corrected, Deer_killed) %>% 
  filter(Year >= 2010)

## -----------------------------------------------------------------------------

## Validate models with Tim Burkitt data, extracting densities at centroid ####

# Average
dens_latlon <- dens_latlon %>% 
  group_by(Site_id, Species) %>% 
  summarise(Dens.avg = mean(Deer_density, na.rm = T)) %>% 
  ungroup()

# Extract and merge
dens_latlon2 <- dens_latlon %>% 
  distinct(Site_id, .keep_all = T) %>% 
  st_as_sf(sf_column_name = "geometry") 

preds_at_points <- raster::extract(preds, dens_latlon2)

extracted_dens <- data.frame(
  Site_id = rep(dens_latlon2$Site_id, times = 3),
  Species = c(rep("RedDeer", 418), rep("SikaDeer", 418), rep("FallowDeer", 418)),
  Model_density = c(preds_at_points[,1], preds_at_points[,2], preds_at_points[,3]))

dens_latlon %<>%
  left_join(extracted_dens)
         

# Calculate Kendall correlation
corr_data <- dens_latlon %>% 
  st_drop_geometry(NULL) %>% 
  group_by(Species) %>% 
  # filter(Species == "RedDeer") %>% 
  cor_test(Dens.avg, Model_density, method = "kendall") %>% 
  mutate(x = 0.30,
         y = 50)


# And finally we generate the plots
dens_latlon %>% 
  ggplot(aes(y = Dens.avg, x = Model_density)) +
  geom_point() +
  # geom_smooth(method='lm', formula= my.formula) + 
  # stat_poly_eq(formula = my.formula, 
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #              parse = TRUE) +
  geom_text(data = corr_data, aes(x = x, y = y, 
                                  label = paste("Kendall corr. =", cor, 
                                  "\n p-value = ", p, sep = " "))) +
  theme_bw() + 
  facet_wrap(~Species, nrow = 2) + 
  labs(y = "Tim Burkitt density", x = "Model prediction")


## -----------------------------------------------------------------------------

## Validate models by county with total culling returns ####

# We first have to add up the predicted intensity by county
 
bycounty_abundances <- data.frame(
  County = ireland$NAME_TAG, 
  Red_model = NA, 
  Sika_model = NA, 
  Fallow_model = NA
)

for (i in seq_along(bycounty_abundances$County)) {
  # i = 1
  county <- ireland_proj %>% filter(NAME_TAG == bycounty_abundances$County[i])
  cropped <- mask(crop(preds, county), county)
  # print(levelplot(cropped, main = bycounty_abundances$County[i]) +
  #         latticeExtra::layer(sp.polygons(as_Spatial(county))))

  sum_raster_cells <-cellStats(cropped, 'sum')  
  bycounty_abundances[i, 2:4] <- sum_raster_cells
}

bycounty_abundances %>% 
  kableExtra::kable(digits = 2, align = "rrrr", caption = "Added abundances by county") %>% 
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")


 
# We now merge that dataset with the culling returns (averaging the last 10 years)
# (Remember that there is no culling data for the NI counties in this dataset)

cull_validating_avg <- culling_data %>% 
  group_by(County, Species) %>% 
  summarise(Deer_total = mean(Deer_killed, na.rm = T)) %>% 
  ungroup() %>% 
  full_join(bycounty_abundances)


# And finally plot 
my.formula <-  y~x

(RD_cull <- cull_validating_avg %>% 
  filter(Species == "Red") %>% 
  ggplot(aes(y = Deer_total, x = Red_model)) + 
  geom_point() + 
  geom_smooth(method='lm', formula= my.formula) +
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  geom_text(aes(label = County), check_overlap = TRUE, nudge_y = 1.5) +
  theme_bw() + 
  labs(y = "Average culling returns (total)", x = "Model prediction", 
       title = "Red deer prediction validation by county"))

cull_validating_avg %>% 
  filter(Species == "Red") %>% 
  cor_test(Deer_total, Red_model, method = "kendall") 
  

(FD_cull <- cull_validating_avg %>% 
  filter(Species == "Fallow") %>% 
  ggplot(aes(y = Deer_total, x = Fallow_model)) + 
  geom_point() + 
  geom_smooth(method='lm', formula= my.formula) +
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  geom_text(aes(label = County), check_overlap = TRUE, nudge_y = 1.5) +
  theme_bw() + 
  labs(y = "Average culling returns (total)", x = "Model prediction", 
       title = "Fallow deer prediction validation by county"))

cull_validating_avg %>% 
  filter(Species == "Fallow") %>% 
  cor_test(Deer_total, Fallow_model, method = "kendall") 


(SD_cull <- cull_validating_avg %>% 
  filter(Species == "Sika") %>% 
  # filter(County != "Wicklow") %>%
  ggplot(aes(y = Deer_total, x = Sika_model)) + 
  geom_point() + 
  geom_smooth(method='lm', formula= my.formula) +
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  geom_text(aes(label = County), check_overlap = TRUE, nudge_y = 1.5) +
  theme_bw() + 
  labs(y = "Average culling returns (total)", x = "Model prediction", 
       title = "Sika deer prediction validation by county"))

cull_validating_avg %>% 
  filter(Species == "Sika") %>% 
  cor_test(Deer_total, Sika_model, method = "kendall") 


(SD2_cull <- cull_validating_avg %>% 
  filter(Species == "Sika") %>% 
  filter(County != "Wicklow") %>%
  ggplot(aes(y = Deer_total, x = Sika_model)) + 
  geom_point() + 
  geom_smooth(method='lm', formula= my.formula) +
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  geom_text(aes(label = County), check_overlap = TRUE, nudge_y = 1.5) +
  theme_bw() + 
  labs(y = "Average culling returns (total)", x = "Model prediction", 
       title = "Sika deer prediction validation by county"))

cull_validating_avg %>% 
  filter(Species == "Red") %>% 
  filter(County != "Wicklow") %>%
  cor_test(Deer_total, Red_model, method = "kendall") 

 
# Finally plot them all together
cowplot::plot_grid(RD_cull, FD_cull, SD_cull, SD2_cull, ncol = 2)
