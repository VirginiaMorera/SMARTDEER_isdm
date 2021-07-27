####--------------------------------####
#### source output summary function ####
####--------------------------------####

source("scripts/outputs.R")


#### model output summary ####
summary.bru_sdm(model_joint)


df <-  pixels(mesh0)
pr <- predict(model_joint, df, ~ exp(Intercept + spde2))
ggplot() + gg(pr)