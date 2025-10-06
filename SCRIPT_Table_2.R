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
expl.ax1.flower.y1 <- pcoa$eig[1]/sum(pcoa$eig)
expl.ax2.flower.y1 <- pcoa$eig[2]/sum(pcoa$eig)
expl.ax1.flower.y1
expl.ax2.flower.y1
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
expl.ax1.flower.y2 <- pcoa$eig[1]/sum(pcoa$eig)
expl.ax2.flower.y2 <- pcoa$eig[2]/sum(pcoa$eig)
expl.ax1.flower.y2
expl.ax2.flower.y2
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
expl.ax1.fruit.y1 <- pcoa$eig[1]/sum(pcoa$eig)
expl.ax2.fruit.y1 <- pcoa$eig[2]/sum(pcoa$eig)
expl.ax1.fruit.y1
expl.ax2.fruit.y1
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
expl.ax1.fruit.y2 <- pcoa$eig[1]/sum(pcoa$eig)
expl.ax2.fruit.y2 <- pcoa$eig[2]/sum(pcoa$eig)
expl.ax1.fruit.y2
expl.ax2.fruit.y2
scores_PCoA <- round(pcoa$li, 3)
rownames(scores_PCoA) <- rownames(angles)

# Matching Phylo and Pheno(scores)
phylo_traits <- match.phylo.data(phy, scores_PCoA)
phy_fruit2 <- phylo_traits$phy
scores_fruit2 <- phylo_traits$data
colnames(scores_fruit2) <- colnames(scores_PCoA)
head(scores_fruit2)


#
# PVR ####

############# FLOWERING - YEAR 1
phy_flor1
scores_flower1

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_flor1)
pvr_flor1 <- PVR::PVR(pvr_class,
                      trait = scores_flower1$A1,
                      method = "moran")
pvr_flor1@PVR

test_p_value <- cbind(scores_flower1$A1, pvr_flor1@Selection$Vectors)
colnames(test_p_value) <- c("T", "Phy1")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

flower_pc1y1 <- c(list(axis_expl = expl.ax1.flower.y1,
                       r.squared = summary(mod)$r.squared,
                       p.value = unname(glance(mod)$p.value)))
flower_pc1y1


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_flor1)
pvr_flor1.2 <- PVR::PVR(pvr_class, 
                        trait = scores_flower1$A2,
                        method = "moran")
pvr_flor1.2@PVR

test_p_value <- cbind(scores_flower1$A2, pvr_flor1.2@Selection$Vectors)
colnames(test_p_value) <- c("T", "Phy1", "Phy2")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1+Phy2, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

flower_pc2y1 <- c(list(axis_expl = expl.ax2.flower.y1,
                       r.squared = summary(mod)$r.squared,
                       p.value = unname(glance(mod)$p.value)))
flower_pc2y1


############# FLOWERING YEAR 2
phy_flor2
scores_flower2

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_flor2)
pvr_flor2 <- PVR::PVR(pvr_class, 
                      trait = scores_flower2$A1,
                      method = "moran")
pvr_flor2@PVR

test_p_value <- cbind(scores_flower2$A1, pvr_flor2@Selection$Vectors)
colnames(test_p_value) <- c("T", "Phy1")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

flower_pc1y2 <- c(list(axis_expl = expl.ax1.flower.y2,
                       r.squared = summary(mod)$r.squared,
                       p.value = unname(glance(mod)$p.value)))
flower_pc1y2

# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_flor2)
pvr_flor2.2 <- PVR::PVR(pvr_class, 
                        trait = scores_flower2$A2,
                        method = "moran")
pvr_flor2.2@PVR

test_p_value <- cbind(scores_flower2$A2, pvr_flor2.2@Selection$Vectors)
colnames(test_p_value) <- c("T", "Phy1")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

flower_pc2y2 <- c(list(axis_expl = expl.ax2.flower.y2,
                       r.squared = summary(mod)$r.squared,
                       p.value = unname(glance(mod)$p.value)))
flower_pc2y2

########## FRUTING - YEAR 1
phy_fruit1
scores_fruit1

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_fruit1)
pvr_fruit1 <- PVR::PVR(pvr_class, 
                       trait = scores_fruit1$A1,
                       method = "moran")
pvr_fruit1@PVR

test_p_value <- cbind(scores_fruit1$A1, pvr_fruit1@Selection$Vectors)
colnames(test_p_value) <- c("T", "Phy1")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc1y1 <- c(list(axis_expl = expl.ax1.fruit.y1,
                      r.squared = summary(mod)$r.squared,
                      p.value = unname(glance(mod)$p.value)))
fruit_pc1y1


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_fruit1)
pvr_fruit1.2 <- PVR::PVR(pvr_class, 
                         trait = scores_fruit1$A2,
                         method = "moran")
pvr_fruit1.2@PVR

test_p_value <- cbind(scores_fruit1$A2, pvr_fruit1.2@Selection$Vectors)
colnames(test_p_value) <- c("T", "Phy1")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc2y1 <- c(list(axis_expl = expl.ax2.fruit.y1,
                      r.squared = summary(mod)$r.squared,
                      p.value = unname(glance(mod)$p.value)))
fruit_pc2y1

################# FRUTING - YEAR 2
phy_fruit2
scores_fruit2

# AXIS 1
pvr_class <- PVR::PVRdecomp(phy_fruit2)
pvr_fruit2 <- PVR::PVR(pvr_class, 
                       trait = scores_fruit2$A1,
                       method = "moran")
pvr_fruit2@PVR

test_p_value <- cbind(scores_fruit2$A1, pvr_fruit2@Selection$Vectors)
colnames(test_p_value) <- c("T", "Phy1")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc1y2 <- c(list(axis_expl = expl.ax1.fruit.y2,
                      r.squared = summary(mod)$r.squared,
                      p.value = unname(glance(mod)$p.value)))
fruit_pc1y2


# AXIS 2
pvr_class <- PVR::PVRdecomp(phy_fruit2)
pvr_fruit2.2 <- PVR::PVR(pvr_class, 
                         trait = scores_fruit2$A2,
                         method = "moran")
pvr_fruit2.2@PVR

test_p_value <- cbind(scores_fruit2$A2, pvr_fruit2.2@Selection$Vectors)
colnames(test_p_value) <- c("T", "Phy1")
test_p_value <- as.data.frame(test_p_value)
mod <- lm(T~Phy1, data = test_p_value)
summary(mod)

summary(mod)$r.squared
glance(mod)$p.value

fruit_pc2y2 <- c(list(axis_expl = expl.ax2.fruit.y2,
                      r.squared = summary(mod)$r.squared,
                      p.value = unname(glance(mod)$p.value)))
fruit_pc2y2


#
# JOINING RESULTS ####

list_pvr <- list(flower_pc1y1 = flower_pc1y1,
                 flower_pc2y1 = flower_pc2y1,
                 flower_pc1y2 = flower_pc1y2,
                 flower_pc2y2 = flower_pc2y2,
                 fruit_pc1y1 = fruit_pc1y1,
                 fruit_pc2y1 = fruit_pc2y1,
                 fruit_pc1y2 = fruit_pc1y2,
                 fruit_pc2y2 = fruit_pc2y2)
length(list_pvr)
names(list_pvr)


tab_varpart_pvr <- matrix(NA,
                          length(list_pvr),
                          3,
                          dimnames=list(names(list_pvr),
                                        c("axis.expl",
                                          "r-squared",
                                          "p-value")))

for (i in 1:length(list_pvr)){
  tab_varpart_pvr[i,1]<-list_pvr[i][[1]]$axis_expl
  tab_varpart_pvr[i,2]<-list_pvr[i][[1]]$r.squared
  tab_varpart_pvr[i,3]<-list_pvr[i][[1]]$p.value
}

tab_varpart_pvr
round(tab_varpart_pvr, digits = 3)
