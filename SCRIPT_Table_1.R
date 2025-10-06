rm(list = ls(envir = globalenv()), envir = globalenv())
gc()

list_of_packages <- c("ape", "SYNCSA", "PCPS", "here", "phytools", "dplyr",
                      "car", "MuMIn", "nlme", "readxl")
lapply(list_of_packages, library, character.only = TRUE)

# LOAD DATA
path_datafile <- here('data/Seger_etal_JVS_data.RData')
load(path_datafile)

# ORGANIZING DATA ####

# FLOWERING
flor
w_flor <- tibble::column_to_rownames(flor, var = "species")
head(w_flor)
dim(w_flor)
w_flor <- t(sqrt(w_flor))

# FRUITING
fruit
w_fruto <- tibble::column_to_rownames(fruit, var = "species")
head(w_fruto)
dim(w_fruto)
w_fruto <- t(sqrt(w_fruto))


# PHYLOGENY
tree
is.ultrametric(tree)
phydist <- cophenetic(tree)


# FLOWERING - PLACING SPECIES IN THE SAMER ORDER OF PHYLOGENY
colnames(w_flor) == colnames(phydist)
dados_flor <- organize.syncsa(comm = w_flor, phylodist = phydist)
dados_flor
ls(dados_flor)
dados_flor$community[1:10,1:3]
dados_flor$phylodist[1:3,1:3]
dim(dados_flor$phylodist)
colnames(dados_flor$community) == colnames(dados_flor$phylodist)


# FRUTING - PLACING SPECIES IN THE SAMER ORDER OF PHYLOGENY
colnames(w_fruto) == colnames(phydist)
dados_fruto <- organize.syncsa(comm = w_fruto, phylodist = phydist)
dados_fruto
ls(dados_fruto)
dados_fruto$community[1:10,1:3]
dados_fruto$phylodist[1:3,1:3]
dim(dados_fruto$phylodist)
colnames(dados_fruto$community) == colnames(dados_fruto$phylodist)

# ENVIRONMENTAL VARIABLES
# CO - CANOPY OPENNESS
# D - DAYLENGTH
# T - TEMPERATURE
# R - RAINFALL

env2
env2 <- sqrt(env2)

rownames(dados_flor$community) == rownames(env2)
rownames(dados_fruto$community) == rownames(env2)

# TO RUN THE GLS MODELS, ATTACH THE ENVIRONMENTAL VARIABLES
attach(env2)
T
D
T_lag2
CO
# detach(env)


# AUTOCORRELATION STRUCTURE

time <- seq(1:24) # linear time

# sine and cosine of each month for circular autocorrelation structure
sincos
sincos <- tibble::column_to_rownames(sincos, var = "month")
head(sincos)
sine <- sincos$sine
sine
class(sine)
cosine <- sincos$cosine
cosine


#
# PRINCIPAL COORDINATE OF PHYLOGENETIC STRUCTURE ANALYSIS - FLOWERING ####

pcps_flower <- pcps(dados_flor$community, dados_flor$phylodist, method = "bray")
ls(pcps_flower)
pcps_flower$P[1:3,1:3]
dim(pcps_flower$P)
class(pcps_flower)
pcps_flower$vectors
pcps_flower$values

pcps.v <- as.data.frame(pcps_flower$vectors)
class(pcps.v)
names(pcps.v)

pcps1 <- pcps.v$pcps.1
pcps1
pcps2 <- pcps.v$pcps.2
pcps2


# MODEL SELECTION

# PCPS 1

PCPS1_models <- list(
  PCPS1.1 <- gls(pcps1 ~ D),
  PCPS1.2 <- gls(pcps1 ~ D_lag1),
  PCPS1.3 <- gls(pcps1 ~ D_lag2),
  PCPS1.4 <- gls(pcps1 ~ D_lag3),
  PCPS1.5 <- gls(pcps1 ~ T),
  PCPS1.6 <- gls(pcps1 ~ T_lag1),
  PCPS1.7 <- gls(pcps1 ~ T_lag2),
  PCPS1.8 <- gls(pcps1 ~ T_lag3),
  PCPS1.9 <- gls(pcps1 ~ R),
  PCPS1.10 <- gls(pcps1 ~ R_lag1),
  PCPS1.11 <- gls(pcps1 ~ R_lag2),
  PCPS1.12 <- gls(pcps1 ~ R_lag3),
  PCPS1.13 <- gls(pcps1 ~ R + R_lag1),
  PCPS1.14 <- gls(pcps1 ~ R_lag1 + R_lag2),
  PCPS1.15 <- gls(pcps1 ~ R_lag2 + R_lag3),
  PCPS1.16 <- gls(pcps1 ~ R + T),
  PCPS1.17 <- gls(pcps1 ~ R + T_lag1),
  PCPS1.18 <- gls(pcps1 ~ R + T_lag2),
  PCPS1.19 <- gls(pcps1 ~ R + T_lag3),
  PCPS1.20 <- gls(pcps1 ~ R + D),
  PCPS1.21 <- gls(pcps1 ~ R + D_lag1),
  PCPS1.22 <- gls(pcps1 ~ R + D_lag2),
  PCPS1.23 <- gls(pcps1 ~ R + D_lag3),
  PCPS1.24 <- gls(pcps1 ~ R_lag1 + T),
  PCPS1.25 <- gls(pcps1 ~ R_lag1 + T_lag1),
  PCPS1.26 <- gls(pcps1 ~ R_lag1 + T_lag2),
  PCPS1.27 <- gls(pcps1 ~ R_lag1 + T_lag3),
  PCPS1.28 <- gls(pcps1 ~ R_lag1 + D),
  PCPS1.29 <- gls(pcps1 ~ R_lag1 + D_lag1),
  PCPS1.30 <- gls(pcps1 ~ R_lag1 + D_lag2),
  PCPS1.31 <- gls(pcps1 ~ R_lag1 + D_lag3),
  PCPS1.32 <- gls(pcps1 ~ R_lag2 + D),
  PCPS1.33 <- gls(pcps1 ~ R_lag1 + R_lag2 + D),
  PCPS1.34 <- gls(pcps1 ~ R_lag2 + T),
  PCPS1.35 <- gls(pcps1 ~ R_lag1 + R_lag2 + T),
  PCPS1.36 <- gls(pcps1 ~ CO),
  PCPS1.37 <- gls(pcps1 ~ CO_lag1),
  PCPS1.38 <- gls(pcps1 ~ CO_lag2),
  PCPS1.39 <- gls(pcps1 ~ CO + R),
  PCPS1.40 <- gls(pcps1 ~ CO + R_lag1),
  PCPS1.41 <- gls(pcps1 ~ CO + R_lag1 + R_lag2),
  PCPS1.NULL  <- gls(pcps1 ~ 1))

aicc1<-model.sel(PCPS1_models)
aicc1
dim(aicc1)

sel.1_flor <- aicc1[ ,19:21]
sel.1_flor
sel.1_flor_delta <- subset(sel.1_flor, delta<=2)
sel.1_flor_delta # Daylength is the best model

# TEST OF RESIDUALS NORMALITY
shapiro.test(residuals(PCPS1.1, type="normalized")) # normal residuals

# AUTOCORRELATION FUNCTION ESTIMATION
acf1_n <- acf(residuals(PCPS1.1, type="normalized"), lag = 24) # autocorrelation at lag 2, 4 and 5

# VISUALIZATION BY VARIOGRAM
plot(Variogram(PCPS1.1, 
               form = ~time,
               robust = TRUE,
               resType = "normalized")) # autocorrel.
plot(Variogram(PCPS1.1, 
               form = ~sine+cosine,
               robust = TRUE,
               resType = "normalized")) # autocorrel.

# TEST WHICH TEMPORAL AUTOCORRELATION STRUCTURE CONTROLS AUTOCORREL. IN RESIDUALS
# FIND THE BEST MODEL (LOWEST AIC)
cp <- list(corCAR1, corExp, corRatio, corLin, corGaus, corSpher)
z_circ <- vector('list', length(cp))
for(k in 1:length(cp)) {
  z_circ[[k]] <- gls(pcps1 ~ D,
                     correlation = cp[[k]](form = ~sine + cosine))
}

cp2 <- list(corCAR1, corExp, corRatio, corLin, corGaus, corSpher)
z2 <- vector('list', length(cp2))
for(k in 1:length(cp2)) {
  z2[[k]] <- gls(pcps1 ~ D,
                 correlation = cp2[[k]](form = ~time))
}

res_flor1 <- anova(z_circ[[1]], z_circ[[2]], z_circ[[3]],
                   z_circ[[4]], z_circ[[5]], z_circ[[6]],
                   z2[[1]], z2[[2]], z2[[3]],
                   z2[[4]], z2[[5]], z2[[6]])
res_flor1 %>% 
  arrange(AIC)

shapiro.test(residuals(z2[[6]], type = "normalized")) # normal residuals
acf1_n <- acf(residuals(z2[[6]], type = "normalized"), lag = 24)  # autocorrelation at lag 4 and 5

shapiro.test(residuals(z2[[4]], type = "normalized")) # normal residuals
acf1_n <- acf(residuals(z2[[4]], type = "normalized"), lag = 24) # autocorrelation at lag 5

shapiro.test(residuals(z2[[3]], type = "normalized")) # normal residuals
acf1_n <- acf(residuals(z2[[3]], type = "normalized"), lag = 24) # autocorrelation at lag 5

shapiro.test(residuals(z2[[5]], type = "normalized")) # normal residuals
acf1_n <- acf(residuals(z2[[5]], type = "normalized"), lag = 24) # controlled autocorrelation


summary(z2[[5]])
r.squaredLR(z2[[5]])
anova(z2[[5]])



# PCPS 2

PCPS2_models <- list(
  PCPS2.1 <- gls(pcps2 ~ D),
  PCPS2.2 <- gls(pcps2 ~ D_lag1),
  PCPS2.3 <- gls(pcps2 ~ D_lag2),
  PCPS2.4 <- gls(pcps2 ~ D_lag3),
  PCPS2.5 <- gls(pcps2 ~ T),
  PCPS2.6 <- gls(pcps2 ~ T_lag1),
  PCPS2.7 <- gls(pcps2 ~ T_lag2),
  PCPS2.8 <- gls(pcps2 ~ T_lag3),
  PCPS2.9 <- gls(pcps2 ~ R),
  PCPS2.10 <- gls(pcps2 ~ R_lag1),
  PCPS2.11 <- gls(pcps2 ~ R_lag2),
  PCPS2.12 <- gls(pcps2 ~ R_lag3),
  PCPS2.13 <- gls(pcps2 ~ R + R_lag1),
  PCPS2.14 <- gls(pcps2 ~ R_lag1 + R_lag2),
  PCPS2.15 <- gls(pcps2 ~ R_lag2 + R_lag3),
  PCPS2.16 <- gls(pcps2 ~ R + T),
  PCPS2.17 <- gls(pcps2 ~ R + T_lag1),
  PCPS2.18 <- gls(pcps2 ~ R + T_lag2),
  PCPS2.19 <- gls(pcps2 ~ R + T_lag3),
  PCPS2.20 <- gls(pcps2 ~ R + D),
  PCPS2.21 <- gls(pcps2 ~ R + D_lag1),
  PCPS2.22 <- gls(pcps2 ~ R + D_lag2),
  PCPS2.23 <- gls(pcps2 ~ R + D_lag3),
  PCPS2.24 <- gls(pcps2 ~ R_lag1 + T),
  PCPS2.25 <- gls(pcps2 ~ R_lag1 + T_lag1),
  PCPS2.26 <- gls(pcps2 ~ R_lag1 + T_lag2),
  PCPS2.27 <- gls(pcps2 ~ R_lag1 + T_lag3),
  PCPS2.28 <- gls(pcps2 ~ R_lag1 + D),
  PCPS2.29 <- gls(pcps2 ~ R_lag1 + D_lag1),
  PCPS2.30 <- gls(pcps2 ~ R_lag1 + D_lag2),
  PCPS2.31 <- gls(pcps2 ~ R_lag1 + D_lag3),
  PCPS2.32 <- gls(pcps2 ~ R_lag2 + D),
  PCPS2.33 <- gls(pcps2 ~ R_lag1 + R_lag2 + D),
  PCPS2.34 <- gls(pcps2 ~ R_lag2 + T),
  PCPS2.35 <- gls(pcps2 ~ R_lag1 + R_lag2 + T),
  PCPS2.36 <- gls(pcps2 ~ CO),
  PCPS2.37 <- gls(pcps2 ~ CO_lag1),
  PCPS2.38 <- gls(pcps2 ~ CO_lag2),
  PCPS2.39 <- gls(pcps2 ~ CO + R),
  PCPS2.40 <- gls(pcps2 ~ CO + R_lag1),
  PCPS2.41 <- gls(pcps2 ~ CO + R_lag1 + R_lag2),
  PCPS2.NULL  <- gls(pcps2 ~ 1))

aicc2 <- model.sel(PCPS2_models)
aicc2
dim(aicc2)

sel.2_flor <- aicc2[ ,19:21]
sel.2_flor
sel.2_flor_delta <- subset(sel.2_flor, delta<=2)
sel.2_flor_delta # the null models is the best model


#
# PCPS.sig ANALYSIS - FLOWERING ####

sig_flower_1 <- pcps.sig(dados_flor$community,
                         dados_flor$phylodist,
                         envir = env2,
                         method = "bray",
                         squareroot = TRUE,
                         checkdata = FALSE,
                         formula = pcps.1~D,
                         FUN = FUN.GLS.marginal,
                         correlation = corGaus(form = ~time),
                         choices = 1,
                         runs = 999)
sig_flower_1

r.squaredLR(gls(pcps1~D, correlation = corGaus(form = ~time)))
anova(sig_flower_1$model, type = "marginal")

#
# PRINCIPAL COORDINATE OF PHYLOGENETIC STRUCTURE ANALYSIS - FRUITING ####

pcps_fruit <- pcps(dados_fruto$community, dados_fruto$phylodist, method = "bray")
ls(pcps_fruit)
pcps_fruit$P[1:3,1:3]
dim(pcps_fruit$P)
class(pcps_fruit)
pcps_fruit$vectors
pcps_fruit$values

pcps.v <- as.data.frame(pcps_fruit$vectors)
class(pcps.v)
names(pcps.v)

pcps1 <- pcps.v$pcps.1
pcps1
pcps2 <- pcps.v$pcps.2
pcps2


# MODEL SELECTION

# PCPS 1

PCPS1_models <- list(
  PCPS1.1 <- gls(pcps1 ~ D),
  PCPS1.2 <- gls(pcps1 ~ D_lag1),
  PCPS1.3 <- gls(pcps1 ~ D_lag2),
  PCPS1.4 <- gls(pcps1 ~ D_lag3),
  PCPS1.5 <- gls(pcps1 ~ T),
  PCPS1.6 <- gls(pcps1 ~ T_lag1),
  PCPS1.7 <- gls(pcps1 ~ T_lag2),
  PCPS1.8 <- gls(pcps1 ~ T_lag3),
  PCPS1.9 <- gls(pcps1 ~ R),
  PCPS1.10 <- gls(pcps1 ~ R_lag1),
  PCPS1.11 <- gls(pcps1 ~ R_lag2),
  PCPS1.12 <- gls(pcps1 ~ R_lag3),
  PCPS1.13 <- gls(pcps1 ~ R + R_lag1),
  PCPS1.14 <- gls(pcps1 ~ R_lag1 + R_lag2),
  PCPS1.15 <- gls(pcps1 ~ R_lag2 + R_lag3),
  PCPS1.16 <- gls(pcps1 ~ R + T),
  PCPS1.17 <- gls(pcps1 ~ R + T_lag1),
  PCPS1.18 <- gls(pcps1 ~ R + T_lag2),
  PCPS1.19 <- gls(pcps1 ~ R + T_lag3),
  PCPS1.20 <- gls(pcps1 ~ R + D),
  PCPS1.21 <- gls(pcps1 ~ R + D_lag1),
  PCPS1.22 <- gls(pcps1 ~ R + D_lag2),
  PCPS1.23 <- gls(pcps1 ~ R + D_lag3),
  PCPS1.24 <- gls(pcps1 ~ R_lag1 + T),
  PCPS1.25 <- gls(pcps1 ~ R_lag1 + T_lag1),
  PCPS1.26 <- gls(pcps1 ~ R_lag1 + T_lag2),
  PCPS1.27 <- gls(pcps1 ~ R_lag1 + T_lag3),
  PCPS1.28 <- gls(pcps1 ~ R_lag1 + D),
  PCPS1.29 <- gls(pcps1 ~ R_lag1 + D_lag1),
  PCPS1.30 <- gls(pcps1 ~ R_lag1 + D_lag2),
  PCPS1.31 <- gls(pcps1 ~ R_lag1 + D_lag3),
  PCPS1.32 <- gls(pcps1 ~ R_lag2 + D),
  PCPS1.33 <- gls(pcps1 ~ R_lag1 + R_lag2 + D),
  PCPS1.34 <- gls(pcps1 ~ R_lag2 + T),
  PCPS1.35 <- gls(pcps1 ~ R_lag1 + R_lag2 + T),
  PCPS1.36 <- gls(pcps1 ~ CO),
  PCPS1.37 <- gls(pcps1 ~ CO_lag1),
  PCPS1.38 <- gls(pcps1 ~ CO_lag2),
  PCPS1.39 <- gls(pcps1 ~ CO + R),
  PCPS1.40 <- gls(pcps1 ~ CO + R_lag1),
  PCPS1.41 <- gls(pcps1 ~ CO + R_lag1 + R_lag2),
  PCPS1.NULL  <- gls(pcps1 ~ 1))

aicc1<-model.sel(PCPS1_models)
aicc1
dim(aicc1)

sel.1_fruto <- aicc1[ ,18:21]
sel.1_fruto
sel.1_fruto_delta<-subset(sel.1_fruto, delta<=2)
sel.1_fruto_delta

summary(PCPS1.4) # Daylength lag 3 is the best model


# TEST OF RESIDUALS NORMALITY
shapiro.test(residuals(PCPS1.4, type="normalized")) # normal residuals

# AUTOCORRELATION FUNCTION ESTIMATION
acf1_n <- acf(residuals(PCPS1.4, type="normalized"), lag = 24) # autocorrelation at lag 2

# VISUALIZATION BY VARIOGRAM
plot(Variogram(PCPS1.2, 
               form = ~time,
               robust = TRUE,
               resType = "normalized")) # autocorrel.
plot(Variogram(PCPS1.2, 
               form = ~sine+cosine,
               robust = TRUE,
               resType = "normalized")) # autocorrel.

# TEST WHICH TEMPORAL AUTOCORRELATION STRUCTURE CONTROLS AUTOCORREL. IN RESIDUALS
# FIND THE BEST MODEL (LOWEST AIC)

cp <- list(corCAR1, corExp, corRatio, corLin, corGaus, corSpher)
z_circ <- vector('list', length(cp))
for(k in 1:length(cp)) {
  z_circ[[k]] <- gls(pcps1 ~ D_lag3,
                     correlation =cp[[k]](form = ~sine + cosine))
}

cp2 <- list(corCAR1,  corRatio, corLin, corGaus, corSpher, corExp)
z2 <- vector('list', length(cp2))
for(k in 1:length(cp2)) {
  z2[[k]] <- gls(pcps1 ~ D_lag3,
                 correlation = cp2[[k]](form = ~time))
}

res_fruto1 <- anova(z_circ[[1]], z_circ[[2]], z_circ[[3]],
                    z_circ[[4]], z_circ[[5]], z_circ[[6]],
                    z2[[1]], z2[[2]], z2[[3]],
                    z2[[4]], z2[[5]], z2[[6]])

res_fruto1 %>% 
  arrange(AIC) # corCAR1 and corExp are equal

shapiro.test(residuals(z2[[1]], type = "normalized")) # non-normal residuals
acf1_n <- acf(residuals(z2[[1]], type = "normalized"), lag = 24)  # no autocorrelation


summary(z2[[1]])
r.squaredLR(z2[[1]])
anova(z2[[1]])


# PCPS 2

PCPS2_models <- list(
  PCPS2.1 <- gls(pcps2 ~ D),
  PCPS2.2 <- gls(pcps2 ~ D_lag1),
  PCPS2.3 <- gls(pcps2 ~ D_lag2),
  PCPS2.4 <- gls(pcps2 ~ D_lag3),
  PCPS2.5 <- gls(pcps2 ~ T),
  PCPS2.6 <- gls(pcps2 ~ T_lag1),
  PCPS2.7 <- gls(pcps2 ~ T_lag2),
  PCPS2.8 <- gls(pcps2 ~ T_lag3),
  PCPS2.9 <- gls(pcps2 ~ R),
  PCPS2.10 <- gls(pcps2 ~ R_lag1),
  PCPS2.11 <- gls(pcps2 ~ R_lag2),
  PCPS2.12 <- gls(pcps2 ~ R_lag3),
  PCPS2.13 <- gls(pcps2 ~ R + R_lag1),
  PCPS2.14 <- gls(pcps2 ~ R_lag1 + R_lag2),
  PCPS2.15 <- gls(pcps2 ~ R_lag2 + R_lag3),
  PCPS2.16 <- gls(pcps2 ~ R + T),
  PCPS2.17 <- gls(pcps2 ~ R + T_lag1),
  PCPS2.18 <- gls(pcps2 ~ R + T_lag2),
  PCPS2.19 <- gls(pcps2 ~ R + T_lag3),
  PCPS2.20 <- gls(pcps2 ~ R + D),
  PCPS2.21 <- gls(pcps2 ~ R + D_lag1),
  PCPS2.22 <- gls(pcps2 ~ R + D_lag2),
  PCPS2.23 <- gls(pcps2 ~ R + D_lag3),
  PCPS2.24 <- gls(pcps2 ~ R_lag1 + T),
  PCPS2.25 <- gls(pcps2 ~ R_lag1 + T_lag1),
  PCPS2.26 <- gls(pcps2 ~ R_lag1 + T_lag2),
  PCPS2.27 <- gls(pcps2 ~ R_lag1 + T_lag3),
  PCPS2.28 <- gls(pcps2 ~ R_lag1 + D),
  PCPS2.29 <- gls(pcps2 ~ R_lag1 + D_lag1),
  PCPS2.30 <- gls(pcps2 ~ R_lag1 + D_lag2),
  PCPS2.31 <- gls(pcps2 ~ R_lag1 + D_lag3),
  PCPS2.32 <- gls(pcps2 ~ R_lag2 + D),
  PCPS2.33 <- gls(pcps2 ~ R_lag1 + R_lag2 + D),
  PCPS2.34 <- gls(pcps2 ~ R_lag2 + T),
  PCPS2.35 <- gls(pcps2 ~ R_lag1 + R_lag2 + T),
  PCPS2.36 <- gls(pcps2 ~ CO),
  PCPS2.37 <- gls(pcps2 ~ CO_lag1),
  PCPS2.38 <- gls(pcps2 ~ CO_lag2),
  PCPS2.39 <- gls(pcps2 ~ CO + R),
  PCPS2.40 <- gls(pcps2 ~ CO + R_lag1),
  PCPS2.41 <- gls(pcps2 ~ CO + R_lag1 + R_lag2),
  PCPS2.NULL  <- gls(pcps2 ~ 1))

aicc2 <- model.sel(PCPS2_models)
aicc2
dim(aicc2)

sel.2_fruto <- aicc2[ ,18:21]
sel.2_fruto
sel.2_fruto_delta<-subset(sel.2_fruto, delta<=2)
sel.2_fruto_delta # Canopy Openness and Daylength lag1 are the best models (with deltaAIC <= 2)

summary(PCPS2.36) # CO
summary(PCPS2.2) # D_lag1


# AUTOCORRELATION TESTS FOR MODEL 1 - CANOPY OPENNESS

# TEST OF RESIDUALS NORMALITY
shapiro.test(residuals(PCPS2.36, type="normalized")) # normal residuals

# AUTOCORRELATION FUNCTION ESTIMATION
acf1_n <- acf(residuals(PCPS2.36, type="normalized"), lag = 24) # autocorrelation at lag 2 and 10

# VISUALIZATION BY VARIOGRAM
plot(Variogram(PCPS2.36, 
               form = ~time,
               robust = TRUE,
               resType = "normalized")) # autocorrel.
plot(Variogram(PCPS2.36, 
               form = ~sine+cosine,
               robust = TRUE,
               resType = "normalized")) # autocorrel.

# TEST WHICH TEMPORAL AUTOCORRELATION STRUCTURE CONTROLS AUTOCORREL. IN RESIDUALS
# FIND THE BEST MODEL (LOWEST AIC)

cp <- list(corCAR1, corExp, corRatio,  corGaus, corSpher, corLin) 
z_circ <- vector('list', length(cp))
for(k in 1:length(cp)) {
  z_circ[[k]] <- gls(pcps2 ~ CO,
                     correlation =cp[[k]](form = ~sine + cosine))
}

cp2 <- list(corCAR1, corExp, corRatio, corGaus, corSpher) #corLin gives false convergence and was removed
z2 <- vector('list', length(cp2))
for(k in 1:length(cp2)) {
  z2[[k]] <- gls(pcps2 ~ CO,
                 correlation = cp2[[k]](form = ~time))
}


res_fruto2 <- anova(z_circ[[1]], z_circ[[2]], z_circ[[3]],
                    z_circ[[4]], z_circ[[5]], z_circ[[6]],
                    z2[[1]], z2[[2]], z2[[3]], z2[[4]], z2[[5]]) # z2[[6]],

res_fruto2 %>% 
  arrange(AIC) # corGaus with linear time is the best model


shapiro.test(residuals(z2[[4]], type = "normalized")) # normal residuals
acf1_n <- acf(residuals(z2[[4]], type = "normalized"), lag = 24)  # no autocorrelation


summary(z2[[4]])
r.squaredLR(z2[[4]])
anova(z2[[4]])


# AUTOCORRELATION TESTS FOR MODEL 2

# TEST OF RESIDUALS NORMALITY
shapiro.test(residuals(PCPS2.2, type="normalized")) # normal residuals

# AUTOCORRELATION FUNCTION ESTIMATION
acf1_n <- acf(residuals(PCPS2.2, type="normalized"), lag = 24) # autocorrelation at lag 10

# VISUALIZATION BY VARIOGRAM
plot(Variogram(PCPS2.2, 
               form = ~time,
               robust = TRUE,
               resType = "normalized")) # autocorrel.
plot(Variogram(PCPS2.2, 
               form = ~sine+cosine,
               robust = TRUE,
               resType = "normalized")) # autocorrel.

# TEST WHICH TEMPORAL AUTOCORRELATION STRUCTURE CONTROLS AUTOCORREL. IN RESIDUALS
# FIND THE BEST MODEL (LOWEST AIC)

cp <- list(corCAR1, corExp, corRatio,  corGaus, corSpher, corLin) 
z_circ <- vector('list', length(cp))
for(k in 1:length(cp)) {
  z_circ[[k]] <- gls(pcps2 ~ D_lag1,
                     correlation =cp[[k]](form = ~sine + cosine))
}

cp2 <- list(corCAR1, corExp, corRatio, corGaus, corSpher, corLin)
z2 <- vector('list', length(cp2))
for(k in 1:length(cp2)) {
  z2[[k]] <- gls(pcps2 ~ D_lag1,
                 correlation = cp2[[k]](form = ~time))
}


res_fruto2 <- anova(z_circ[[1]], z_circ[[2]], z_circ[[3]],
                    z_circ[[4]], z_circ[[5]], z_circ[[6]],
                    z2[[1]], z2[[2]], z2[[3]], z2[[4]], z2[[5]], z2[[6]]) # z2[[6]],

res_fruto2 %>% 
  arrange(AIC) # corRatio with circular time was the best model


shapiro.test(residuals(z_circ[[3]], type = "normalized")) # normal residuals

plot(Variogram(z_circ[[3]], 
               form = ~sine+cosine,
               robust = TRUE,
               resType = "normalized")) # no autocorrelation

summary(z_circ[[3]])
r.squaredLR(z_circ[[3]])
anova(z_circ[[3]])

#
# PCPS.sig ANALYSIS - FRUITING ####

sig_fruit_1 <- pcps.sig(dados_fruto$community,
                        dados_fruto$phylodist,
                        envir = env2,
                        method = "bray",
                        squareroot = TRUE,
                        checkdata = FALSE,
                        formula = pcps.1~D_lag3,
                        FUN = FUN.GLS.marginal,
                        correlation = corCAR1(form = ~time),
                        choices = 1,
                        runs = 999)
sig_fruit_1

r.squaredLR(gls(pcps1~D_lag3, correlation = corCAR1(form = ~time)))
anova(sig_fruit_1$model, type = "marginal")


sig_fruit_2_CO <- pcps.sig(dados_fruto$community,
                           dados_fruto$phylodist,
                           envir = env2,
                           method = "bray",
                           squareroot = TRUE,
                           checkdata = FALSE,
                           formula = pcps.2~CO,
                           FUN = FUN.GLS.marginal,
                           correlation = corGaus(form = ~time),
                           choices = 2,
                           runs = 999)
sig_fruit_2_CO

r.squaredLR(gls(pcps2~CO, correlation = corGaus(form = ~time)))
anova(sig_fruit_2_CO$model, type = "marginal")


sig_fruit_2_Dlag1 <- pcps.sig(dados_fruto$community,
                           dados_fruto$phylodist,
                           envir = env2,
                           method = "bray",
                           squareroot = TRUE,
                           checkdata = FALSE,
                           formula = pcps.2~D_lag1,
                           FUN = FUN.GLS.marginal,
                           correlation = corRatio(form = ~sine+cosine),
                           choices = 2,
                           runs = 999)
sig_fruit_2_Dlag1

r.squaredLR(gls(pcps2~D_lag1, correlation = corRatio(form = ~sine+cosine)))
anova(sig_fruit_2_Dlag1$model, type = "marginal")

#