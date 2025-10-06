rm(list = ls(envir = globalenv()), envir = globalenv())
gc()

list_of_packages <- c("picante", "here", "phytools", "ape", "readxl", "dplyr",
                      "ade4", "FD", "PVR", "broom")
lapply(list_of_packages, library, character.only = TRUE)

# LOAD DATA
path_datafile <- here('data/Seger_etal_JVS_data.RData')
load(path_datafile)

# ORGANIZING DATA ####

tree
is.ultrametric(tree)
is.binary(tree)

peak_modal_angle_year_flower
peak_modal_angle_year_fruit

flor1 <- peak_modal_angle_year_flower %>% 
  select(year1 = modal_angle_flower1) %>% 
  filter(!is.na(as.numeric(year1))) %>%
  mutate(year1 = round(as.numeric(year1), digits = 2))
dim(flor1)
flor1

flor2 <- peak_modal_angle_year_flower %>% 
  select(year2 = modal_angle_flower2) %>% 
  filter(!is.na(as.numeric(year2))) %>%
  mutate(year2 = round(as.numeric(year2), digits = 2))
dim(flor2)
flor2

fruto1 <- peak_modal_angle_year_fruit %>% 
  select(year1 = modal_angle_fruit1) %>% 
  filter(!is.na(as.numeric(year1))) %>%
  mutate(year1 = round(as.numeric(year1), digits = 2))
dim(fruto1)
fruto1

fruto2 <- peak_modal_angle_year_fruit %>% 
  select(year2 = modal_angle_fruit2) %>% 
  filter(!is.na(as.numeric(year2))) %>%
  mutate(year2 = round(as.numeric(year2), digits = 2))
dim(fruto2)
fruto2

env2
env <- sqrt(env2)
env$month1 <- rownames(env)
env$month2 <- rownames(env)
env


#
##### PHENOLOGICAL AXES ####

# FLOWERING - YEAR 1
## species ordering matching phylotree
phylo_traits <- match.phylo.data(tree, flor1)
phy <- phylo_traits$phy
angles <- phylo_traits$data
angles<- as.matrix(angles)

# calculate matrix of angular distances
matriz.Y <- matrix(angles, nrow = length(angles),
                   ncol = length(angles),
                   byrow = FALSE)
rownames(matriz.Y) <- rownames(angles)
colnames(matriz.Y) <- rownames(angles)
matriz.Y[1:6, 1:6]

step_1 <- matriz.Y-t(matriz.Y) #difference between each pair of species
step_2 <- abs(step_1) # transforming differences to absolute values
step_3 <- ifelse(step_2>180, step_2-360, step_2)
step_4 <- abs(step_3) # transforming to absolute values
dist_angular <- as.dist(step_4) # matrix showing the angular distances
as.matrix(dist_angular)
dim(step_4)
step_4[1:6, 1:6]

# PCoA to extract the scores that represent the angles on a linear scale
pcoa <- dudi.pco(dist_angular, scannf = FALSE, nf = 2)  
explic.ax1 <- pcoa$eig[1]/sum(pcoa$eig)
explic.ax2 <- pcoa$eig[2]/sum(pcoa$eig)
explic.ax1
explic.ax2
scores_PCoA <- round(pcoa$li, 3)
rownames(scores_PCoA) <- rownames(angles)

# Matching Phylo and Pheno(scores)
phylo_traits <- match.phylo.data(phy, scores_PCoA)
phy_flor1 <- phylo_traits$phy
scores_flower1 <- phylo_traits$data
colnames(scores_flower1) <- colnames(scores_PCoA)
head(scores_flower1)


# FLOWERING - YEAR 2
## species ordering matching phylotree
phylo_traits <- match.phylo.data(tree, flor2)
phy <- phylo_traits$phy
angles <- phylo_traits$data
angles<- as.matrix(angles)

# calculate matrix of angular distances
matriz.Y <- matrix(angles, nrow = length(angles),
                   ncol = length(angles),
                   byrow = FALSE)
rownames(matriz.Y) <- rownames(angles)
colnames(matriz.Y) <- rownames(angles)
matriz.Y[1:6, 1:6]

step_1 <- matriz.Y-t(matriz.Y) #difference between each pair of species
step_2 <- abs(step_1) # transforming differences to absolute values
step_3 <- ifelse(step_2>180, step_2-360, step_2)
step_4 <- abs(step_3) # transforming to absolute values
dist_angular <- as.dist(step_4) # matrix showing the angular distances
as.matrix(dist_angular)
dim(step_4)
step_4[1:6, 1:6]

# PCoA to extract the scores that represent the angles on a linear scale
pcoa <- dudi.pco(dist_angular, scannf = FALSE, nf = 2)  
explic.ax1 <- pcoa$eig[1]/sum(pcoa$eig)
explic.ax2 <- pcoa$eig[2]/sum(pcoa$eig)
explic.ax1
explic.ax2
scores_PCoA <- round(pcoa$li, 3)
rownames(scores_PCoA) <- rownames(angles)

# Matching Phylo and Pheno(scores)
phylo_traits <- match.phylo.data(phy, scores_PCoA)
phy_flor2 <- phylo_traits$phy
scores_flower2 <- phylo_traits$data
colnames(scores_flower2) <- colnames(scores_PCoA)
head(scores_flower2)


# FRUITING - YEAR 1
## species ordering matching phylotree
phylo_traits <- match.phylo.data(tree, fruto1)
phy <- phylo_traits$phy
angles <- phylo_traits$data
angles<- as.matrix(angles)

# calculate matrix of angular distances
matriz.Y <- matrix(angles, nrow = length(angles),
                   ncol = length(angles),
                   byrow = FALSE)
rownames(matriz.Y) <- rownames(angles)
colnames(matriz.Y) <- rownames(angles)
matriz.Y[1:6, 1:6]

step_1 <- matriz.Y-t(matriz.Y) #difference between each pair of species
step_2 <- abs(step_1) # transforming differences to absolute values
step_3 <- ifelse(step_2>180, step_2-360, step_2)
step_4 <- abs(step_3) # transforming to absolute values
dist_angular <- as.dist(step_4) # matrix showing the angular distances
as.matrix(dist_angular)
dim(step_4)
step_4[1:6, 1:6]

# PCoA to extract the scores that represent the angles on a linear scale
pcoa <- dudi.pco(dist_angular, scannf = FALSE, nf = 2)  
explic.ax1 <- pcoa$eig[1]/sum(pcoa$eig)
explic.ax2 <- pcoa$eig[2]/sum(pcoa$eig)
explic.ax1
explic.ax2
scores_PCoA <- round(pcoa$li, 3)
rownames(scores_PCoA) <- rownames(angles)

# Matching Phylo and Pheno(scores)
phylo_traits <- match.phylo.data(phy, scores_PCoA)
phy_fruit1 <- phylo_traits$phy
scores_fruit1 <- phylo_traits$data
colnames(scores_fruit1) <- colnames(scores_PCoA)
head(scores_fruit1)


# FRUITING - YEAR 2
## species ordering matching phylotree
phylo_traits <- match.phylo.data(tree, fruto2)
phy <- phylo_traits$phy
angles <- phylo_traits$data
angles<- as.matrix(angles)

# calculate matrix of angular distances
matriz.Y <- matrix(angles, nrow = length(angles),
                   ncol = length(angles),
                   byrow = FALSE)
rownames(matriz.Y) <- rownames(angles)
colnames(matriz.Y) <- rownames(angles)
matriz.Y[1:6, 1:6]

step_1 <- matriz.Y-t(matriz.Y) #difference between each pair of species
step_2 <- abs(step_1) # transforming differences to absolute values
step_3 <- ifelse(step_2>180, step_2-360, step_2)
step_4 <- abs(step_3) # transforming to absolute values
dist_angular <- as.dist(step_4) # matrix showing the angular distances
as.matrix(dist_angular)
dim(step_4)
step_4[1:6, 1:6]

# PCoA to extract the scores that represent the angles on a linear scale
pcoa <- dudi.pco(dist_angular, scannf = FALSE, nf = 2)  
explic.ax1 <- pcoa$eig[1]/sum(pcoa$eig)
explic.ax2 <- pcoa$eig[2]/sum(pcoa$eig)
explic.ax1
explic.ax2
scores_PCoA <- round(pcoa$li, 3)
rownames(scores_PCoA) <- rownames(angles)

# Matching Phylo and Pheno(scores)
phylo_traits <- match.phylo.data(phy, scores_PCoA)
phy_fruit2 <- phylo_traits$phy
scores_fruit2 <- phylo_traits$data
colnames(scores_fruit2) <- colnames(scores_PCoA)
head(scores_fruit2)


#
# ENVIRONMENT PER YEAR ####

# ENVIRONMENT FLOWERING - YEAR 1

flor1.1 <- peak_modal_angle_year_flower %>% 
  rownames_to_column(var = "species") %>% 
  select(species, angle1 = modal_angle_flower1, month1) %>% 
  filter(!is.na(as.numeric(angle1))) %>% 
  mutate(angle1 = round(as.numeric(angle1), digits = 2))
dim(flor1.1)
flor1.1

env1 <- env %>% 
  select(D, month1)

flor_env1 <- left_join(flor1.1, env1, by = "month1") %>% 
  tibble::column_to_rownames(var = "species") %>% 
  select(D)
flor_env1

# Matching Phylo and environment
phylo_traits <- match.phylo.data(tree, flor_env1)
env_flor1 <- phylo_traits$data
colnames(env_flor1) <- colnames(flor_env1)
head(env_flor1)


# ENVIRONMENT FLOWERING - YEAR 2

flor1.2 <- peak_modal_angle_year_flower %>% 
  rownames_to_column(var = "species") %>% 
  select(species, angle2 = modal_angle_flower2, month2) %>% 
  filter(!is.na(as.numeric(angle2))) %>% 
  mutate(angle2 = round(as.numeric(angle2), digits = 2))
dim(flor1.2)
flor1.2

env1 <- env %>% 
  select(D, month2)

flor_env2 <- left_join(flor1.2, env1, by = "month2") %>% 
  tibble::column_to_rownames(var = "species") %>% 
  select(D)
flor_env2

# Matching Phylo and environment
phylo_traits <- match.phylo.data(tree, flor_env2)
env_flor2 <- phylo_traits$data
colnames(env_flor2) <- colnames(flor_env2)
head(env_flor2)


# ENVIRONMENT FRUITING - YEAR 1

fruto1.1 <- peak_modal_angle_year_fruit %>% 
  rownames_to_column(var = "species") %>% 
  select(species, angle1 = modal_angle_fruit1, month1) %>% 
  filter(!is.na(as.numeric(angle1))) %>% 
  mutate(angle1 = round(as.numeric(angle1), digits = 2))
dim(fruto1.1)
fruto1.1

env1 <- env %>% 
  select(D_lag3, CO, D_lag1, month1)

fruit_env1 <- left_join(fruto1.1, env1, by = "month1") %>% 
  tibble::column_to_rownames(var = "species") %>% 
  select(D_lag3, CO, D_lag1)
fruit_env1

# Matching Phylo and environment
phylo_traits <- match.phylo.data(tree, fruit_env1)
env_fruit1 <- phylo_traits$data
colnames(env_fruit1) <- colnames(fruit_env1)
head(env_fruit1)


# ENVIRONMENT FRUITING - YEAR 2

fruto1.2 <- peak_modal_angle_year_fruit %>% 
  rownames_to_column(var = "species") %>% 
  select(species, angle2 = modal_angle_fruit2, month2) %>% 
  filter(!is.na(as.numeric(angle2))) %>% 
  mutate(angle2 = round(as.numeric(angle2), digits = 2))
dim(fruto1.2)
fruto1.2

env1 <- env %>% 
  select(D_lag3, CO, D_lag1, month2)

fruit_env2 <- left_join(fruto1.2, env1, by = "month2") %>% 
  tibble::column_to_rownames(var = "species") %>% 
  select(D_lag3, CO, D_lag1)
fruit_env2

# Matching Phylo and environment
phylo_traits <- match.phylo.data(tree, fruit_env2)
env_fruit2 <- phylo_traits$data
colnames(env_fruit2) <- colnames(fruit_env2)
head(env_fruit2)


#
# PVR ####

############# FLOWERING - YEAR 1 - DAYLENGTH
env_flor1
phy_flor1
scores_flower1


# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_flor1)

pvr_flor1 <- PVR::PVR(pvr_class,
                      trait = scores_flower1$A1,
                      envVar = env_flor1$D,
                      method = "moran")
pvr_flor1@PVR$R2
pvr_flor1@PVR$p

test_p_value <- cbind(scores_flower1$A1, pvr_flor1@Selection$Vectors, env_flor1$D)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

flor_pc1y1_daylength <- c(list(r.squared = summary(mod)$r.squared,
                         p.value = unname(glance(mod)$p.value)),
                         pvr_flor1@VarPart)
flor_pc1y1_daylength 


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_flor1)

pvr_flor1.2 <- PVR::PVR(pvr_class, 
                        trait = scores_flower1$A2,
                        envVar = env_flor1$D,
                        method = "moran")
pvr_flor1.2@PVR$R2
pvr_flor1.2@PVR$p

test_p_value <- cbind(scores_flower1$A2, pvr_flor1.2@Selection$Vectors, env_flor1$D)
colnames(test_p_value) <- c("T", "Phy1", "Phy2", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+Phy2+E,data=test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

flor_pc2y1_daylength <- c(list(r.squared = summary(mod)$r.squared,
                               p.value = unname(glance(mod)$p.value)),
                          pvr_flor1.2@VarPart)
flor_pc2y1_daylength 


############# FLOWERING YEAR 2 - DAYLENGTH
env_flor2
phy_flor2
scores_flower2

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_flor2)

pvr_flor2 <- PVR::PVR(pvr_class, 
                      trait = scores_flower2$A1,
                      envVar = env_flor2$D,
                      method = "moran")

pvr_flor2@PVR$R2
pvr_flor2@PVR$p

test_p_value <- cbind(scores_flower2$A1, pvr_flor2@Selection$Vectors, env_flor2$D)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

flor_pc1y2_daylength <- c(list(r.squared = summary(mod)$r.squared,
                               p.value = unname(glance(mod)$p.value)),
                          pvr_flor2@VarPart)
flor_pc1y2_daylength 


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_flor2)

pvr_flor2.2 <- PVR::PVR(pvr_class, 
                        trait = scores_flower2$A2,
                        envVar = env_flor2$D,
                        method = "moran")
pvr_flor2.2@PVR
pvr_flor2.2@PVR$R2
pvr_flor2.2@PVR$p

test_p_value <- cbind(scores_flower2$A2, pvr_flor2.2@Selection$Vectors, env_flor2$D)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

flor_pc2y2_daylength <- c(list(r.squared = summary(mod)$r.squared,
                               p.value = unname(glance(mod)$p.value)),
                          pvr_flor2.2@VarPart)
flor_pc2y2_daylength 


########## FRUTING - YEAR 1 - daylength lag3
env_fruit1
phy_fruit1
scores_fruit1

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_fruit1)

pvr_fruit1 <- PVR::PVR(pvr_class,
                       trait = scores_fruit1$A1,
                       envVar = env_fruit1$D_lag3,
                       method = "moran")
pvr_fruit1@PVR
pvr_fruit1@PVR$R2
pvr_fruit1@PVR$p

test_p_value <- cbind(scores_fruit1$A1, pvr_fruit1@Selection$Vectors, env_fruit1$D_lag3)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc1y1_daylength_lag3 <- c(list(r.squared = summary(mod)$r.squared,
                               p.value = unname(glance(mod)$p.value)),
                               pvr_fruit1@VarPart)
fruit_pc1y1_daylength_lag3 


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_fruit1)

pvr_fruit1.2 <- PVR::PVR(pvr_class,
                         trait = scores_fruit1$A2,
                         envVar = env_fruit1$D_lag3,
                         method = "moran")
pvr_fruit1.2@PVR
pvr_fruit1.2@PVR$R2
pvr_fruit1.2@PVR$p

test_p_value <- cbind(scores_fruit1$A2, pvr_fruit1.2@Selection$Vectors, env_fruit1$D_lag3)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc2y1_daylength_lag3 <- c(list(r.squared = summary(mod)$r.squared,
                                     p.value = unname(glance(mod)$p.value)),
                                pvr_fruit1.2@VarPart)
fruit_pc2y1_daylength_lag3 


################# FRUTING - YEAR 2 - daylength lag3
env_fruit2
phy_fruit2
scores_fruit2

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_fruit2)

pvr_fruit2 <- PVR::PVR(pvr_class,
                       trait = scores_fruit2$A1,
                       envVar = env_fruit2$D_lag3,
                       method = "moran")
pvr_fruit2@PVR
pvr_fruit2@PVR$R2
pvr_fruit2@PVR$p

test_p_value <- cbind(scores_fruit2$A1, pvr_fruit2@Selection$Vectors, env_fruit2$D_lag3)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc1y2_daylength_lag3 <- c(list(r.squared = summary(mod)$r.squared,
                                     p.value = unname(glance(mod)$p.value)),
                                pvr_fruit2@VarPart)
fruit_pc1y2_daylength_lag3 


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_fruit2)

pvr_fruit2.2 <- PVR::PVR(pvr_class, 
                         trait = scores_fruit2$A2,
                         envVar = env_fruit2$D_lag3,
                         method = "moran")
pvr_fruit2.2@PVR
pvr_fruit2.2@PVR$R2
pvr_fruit2.2@PVR$p

test_p_value <- cbind(scores_fruit2$A2, pvr_fruit2.2@Selection$Vectors, env_fruit2$D_lag3)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc2y2_daylength_lag3 <- c(list(r.squared = summary(mod)$r.squared,
                                     p.value = unname(glance(mod)$p.value)),
                                pvr_fruit2.2@VarPart)
fruit_pc2y2_daylength_lag3 


########## FRUTING - YEAR 1 - CANOPY OPENNESS

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_fruit1)
pvr_fruit1 <- PVR::PVR(pvr_class, 
                       trait = scores_fruit1$A1,
                       envVar = env_fruit1$CO,
                       method = "moran")
pvr_fruit1@PVR
pvr_fruit1@PVR$R2
pvr_fruit1@PVR$p

test_p_value <- cbind(scores_fruit1$A1, pvr_fruit1@Selection$Vectors, env_fruit1$CO)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc1y1_canopyopeness <- c(list(r.squared = summary(mod)$r.squared,
                                     p.value = unname(glance(mod)$p.value)),
                               pvr_fruit1@VarPart)
fruit_pc1y1_canopyopeness 


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_fruit1)
pvr_fruit1.2 <- PVR::PVR(pvr_class, 
                         trait = scores_fruit1$A2,
                         envVar = env_fruit1$CO,
                         method = "moran")
pvr_fruit1.2@PVR
pvr_fruit1.2@PVR$R2
pvr_fruit1.2@PVR$p

test_p_value <- cbind(scores_fruit1$A2, pvr_fruit1.2@Selection$Vectors, env_fruit1$CO)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc2y1_canopyopeness <- c(list(r.squared = summary(mod)$r.squared,
                                    p.value = unname(glance(mod)$p.value)),
                               pvr_fruit1.2@VarPart)
fruit_pc2y1_canopyopeness 


################# FRUTING - YEAR 2 - CANOPY OPENNESS

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_fruit2)
pvr_fruit2 <- PVR::PVR(pvr_class, 
                       trait = scores_fruit2$A1,
                       envVar = env_fruit2$CO,
                       method = "moran")
pvr_fruit2@PVR
pvr_fruit2@PVR$R2
pvr_fruit2@PVR$p

test_p_value <- cbind(scores_fruit2$A1, pvr_fruit2@Selection$Vectors, env_fruit2$CO)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc1y2_canopyopeness <- c(list(r.squared = summary(mod)$r.squared,
                                    p.value = unname(glance(mod)$p.value)),
                               pvr_fruit2@VarPart)
fruit_pc1y2_canopyopeness 


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_fruit2)
pvr_fruit2.2 <- PVR::PVR(pvr_class, 
                         trait = scores_fruit2$A2,
                         envVar = env_fruit2$CO,
                         method = "moran")
pvr_fruit2.2@PVR
pvr_fruit2.2@PVR$R2
pvr_fruit2.2@PVR$p

test_p_value <- cbind(scores_fruit2$A2, pvr_fruit2.2@Selection$Vectors, env_fruit2$CO)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc2y2_canopyopeness <- c(list(r.squared = summary(mod)$r.squared,
                                    p.value = unname(glance(mod)$p.value)),
                               pvr_fruit2.2@VarPart)
fruit_pc2y2_canopyopeness 


########## FRUTING - YEAR 1 - DAYLENGTH LAG1

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_fruit1)
pvr_fruit1 <- PVR::PVR(pvr_class, 
                       trait = scores_fruit1$A1,
                       envVar = env_fruit1$D_lag1,
                       method = "moran")
pvr_fruit1@PVR
pvr_fruit1@PVR$R2
pvr_fruit1@PVR$p

test_p_value <- cbind(scores_fruit1$A1, pvr_fruit1@Selection$Vectors, env_fruit1$D_lag1)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc1y1_daylength_lag1 <- c(list(r.squared = summary(mod)$r.squared,
                                    p.value = unname(glance(mod)$p.value)),
                                pvr_fruit1@VarPart)
fruit_pc1y1_daylength_lag1 


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_fruit1)
pvr_fruit1.2 <- PVR::PVR(pvr_class, 
                         trait = scores_fruit1$A2,
                         envVar = env_fruit1$D_lag1,
                         method = "moran")
pvr_fruit1.2@PVR
pvr_fruit1.2@PVR$R2
pvr_fruit1.2@PVR$p

test_p_value <- cbind(scores_fruit1$A2, pvr_fruit1.2@Selection$Vectors, env_fruit1$D_lag1)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc2y1_daylength_lag1 <- c(list(r.squared = summary(mod)$r.squared,
                                     p.value = unname(glance(mod)$p.value)),
                                pvr_fruit1.2@VarPart)
fruit_pc2y1_daylength_lag1 


################# FRUTING - YEAR 2 - DAYLENGTH LAG1

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_fruit2)
pvr_fruit2 <- PVR::PVR(pvr_class, 
                       trait = scores_fruit2$A1,
                       envVar = env_fruit2$D_lag1,
                       method = "moran")
pvr_fruit2@PVR
pvr_fruit2@PVR$R2
pvr_fruit2@PVR$p

test_p_value <- cbind(scores_fruit2$A1, pvr_fruit2@Selection$Vectors, env_fruit2$D_lag1)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc1y2_daylength_lag1 <- c(list(r.squared = summary(mod)$r.squared,
                                     p.value = unname(glance(mod)$p.value)),
                                pvr_fruit2@VarPart)
fruit_pc1y2_daylength_lag1


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_fruit2)
pvr_fruit2.2 <- PVR::PVR(pvr_class, 
                         trait = scores_fruit2$A2,
                         envVar = env_fruit2$D_lag1,
                         method = "moran")
pvr_fruit2.2@PVR
pvr_fruit2.2@PVR$R2
pvr_fruit2.2@PVR$p

test_p_value <- cbind(scores_fruit2$A2, pvr_fruit2.2@Selection$Vectors, env_fruit2$D_lag1)
colnames(test_p_value) <- c("T", "Phy1", "E")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+E, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc2y2_daylength_lag1 <- c(list(r.squared = summary(mod)$r.squared,
                                     p.value = unname(glance(mod)$p.value)),
                                pvr_fruit2.2@VarPart)
fruit_pc2y2_daylength_lag1


#
# JOINING RESULTS ####

list_pvr <- list(flor_pc1y1_daylength = flor_pc1y1_daylength,
               flor_pc2y1_daylength = flor_pc2y1_daylength,
               flor_pc1y2_daylength = flor_pc1y2_daylength,
               flor_pc2y2_daylength = flor_pc2y2_daylength,
               fruit_pc1y1_canopyopeness = fruit_pc1y1_canopyopeness,
               fruit_pc2y1_canopyopeness = fruit_pc2y1_canopyopeness,
               fruit_pc1y2_canopyopeness = fruit_pc1y2_canopyopeness,
               fruit_pc2y2_canopyopeness = fruit_pc2y2_canopyopeness,
               fruit_pc1y1_daylength_lag1 = fruit_pc1y1_daylength_lag1,
               fruit_pc2y1_daylength_lag1 = fruit_pc2y1_daylength_lag1,
               fruit_pc1y2_daylength_lag1 = fruit_pc1y2_daylength_lag1,
               fruit_pc2y2_daylength_lag1 = fruit_pc2y2_daylength_lag1,
               fruit_pc1y1_daylength_lag3 = fruit_pc1y1_daylength_lag3,
               fruit_pc2y1_daylength_lag3 = fruit_pc2y1_daylength_lag3,
               fruit_pc1y2_daylength_lag3 = fruit_pc1y2_daylength_lag3,
               fruit_pc2y2_daylength_lag3 = fruit_pc2y2_daylength_lag3)
length(list_pvr)
names(list_pvr)

tab_varpart_pvr <- matrix(NA,
                          length(list_pvr),
                          6,
                          dimnames=list(names(list_pvr),
                                        c("r-squared",
                                          "p-value",
                                          "a_environment",
                                          "b_phylogeny_environment",
                                          "c_phylogeny",
                                          "d_unexplained")))

for (i in 1:length(list_pvr)){
  tab_varpart_pvr[i,1]<-list_pvr[i][[1]]$r.squared
  tab_varpart_pvr[i,2]<-list_pvr[i][[1]]$p.value
  tab_varpart_pvr[i,3]<-list_pvr[i][[1]]$a
  tab_varpart_pvr[i,4]<-list_pvr[i][[1]]$b
  tab_varpart_pvr[i,5]<-list_pvr[i][[1]]$c
  tab_varpart_pvr[i,6]<-list_pvr[i][[1]]$d
}

tab_varpart_pvr
round(tab_varpart_pvr, digits = 3)

#