rm(list = ls(envir = globalenv()), envir = globalenv())
gc()

list_of_packages <- c("ggplot2", "dplyr", "stringr", "here", "readxl", "SYNCSA", "PCPS")
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
env

# PRINCIPAL COORDINATE OF PHYLOGENETIC STRUCTURE ANALYSIS - FLOWERING
pcps_flower <- pcps(dados_flor$community, dados_flor$phylodist, method = "bray")
pcps_flower$values

a <- summary(pcps_flower, choices = c(1,2))$scores
ls(a)
pcps_flor <- as.data.frame(a$scores.sites)
pcps_flor
pcps_flor12_sp <- as.data.frame(a$scores.species)
pcps_flor12_sp


# PRINCIPAL COORDINATE OF PHYLOGENETIC STRUCTURE ANALYSIS - FRUITING

pcps_fruit <- pcps(dados_fruto$community, dados_fruto$phylodist, method = "bray")
pcps_fruit$values

a <- summary(pcps_fruit, choices = c(1,2))$scores
ls(a)
pcps_fruto <- as.data.frame(a$scores.sites)
pcps_fruto
pcps_fruto12_sp <- as.data.frame(a$scores.species)
pcps_fruto12_sp

# JOINING DATA
data <- data.frame(pcps1flor = pcps_flor$pcps.1,
                   pcps2flor = pcps_flor$pcps.2,
                   pcps1fruto = pcps_fruto$pcps.1,
                   pcps2fruto = pcps_fruto$pcps.2,
                   day = env$D,
                   day3 = env$D_lag3,
                   CO = env$CO)
head(data) 
data

#
## PCPS FLOWER 1 AND 2 ####

plot1 <- ggplot (data = data, aes(x = pcps1flor, y = pcps2flor)) +
  xlim(NA, 0.225)

plot2 <- 
  plot1 +
  geom_point(aes(size = day, 
                 fill = day), 
             shape = 21, 
             alpha = 0.8, 
             colour="#d95f0e") +
  scale_fill_gradient2(high = "#ffffd4", 
                       mid="#fe9929", 
                       low = "#993404",
                       midpoint = 12, 
                       limits = c(10, 14.1), # limits tem q estar presente e = em ambos scale_
                       guide = "legend") +
  scale_radius(range = c(2.5, 10.5), 
               limits = c(10, 14.1)) +
  theme_bw() +
  geom_vline(xintercept = 0,
             color = "lightgrey",
             size = 0.2) +
  geom_hline(yintercept = 0,
             color = "lightgrey",
             size = 0.2) +
  labs(x ="PCPS 1 (35.9%)",
       y = "PCPS 2 (20.6%)",
       size = "Daylength (h)",
       fill = "Daylength (h)") +
  theme(legend.title = element_text(size = 10.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 1, vjust = 2),
        axis.title.y = element_text(hjust = 1, vjust = -0.6),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())

plot2

plot3 <- plot2 +
  geom_point(data = pcps_flor12_sp, 
             mapping = aes(x = pcps.1, y = pcps.2),
             shape = 4, size = 2)

plot3

# THIS FIGURE WAS EXPORTED AND CLADES' NAMES AND LIMITS WERE INSERTED USING MACBOOK KEYNOTE APP

ggsave(here("figs", "Figure_3.jpg"),
       plot3,  
       dpi = 1200, 
       units = "in", 
       height = 5.95, 
       width = 7.95)


#
## PCPS FRUTO 1 E 2 CANOPY OPENNESS ####

plot1 <- 
  ggplot() +
  theme_bw() +
  geom_vline(xintercept = 0, 
             color = "lightgrey", 
             size = 0.2) +
  geom_hline(yintercept = 0, 
             color = "lightgrey", 
             size = 0.2) +
  geom_point(data, 
             mapping = aes(x = pcps1fruto, 
                           y = pcps2fruto,
                           size = day3,
                           fill = CO),
             shape=21, 
             colour="#d95f0e", 
             alpha = 0.8) +
  labs(x="PCPS 1 (33.7%)",
       y= "PCPS 2 (11.8%)",
       size = "Daylength\n3 month\nlag (h)",
       fill = "Canopy\nopenness\n(%)")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 1, vjust = 2),
        axis.title.y = element_text(hjust = 1, vjust = -0.7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())
plot1


plot2 <- 
  plot1 +
  scale_radius(range = c(2.5, 7.7), limits = c(10, 14.1)) +
  scale_fill_gradient2(high = "#ffffd4", 
                       mid = "#fe9929", 
                       low = "#993404", 
                       limits = c(0, 62.2),
                       midpoint = 30,
                       guide = guide_colorbar(barwidth = 0.9,
                                              barheight = 6,
                                              nbin = 5,
                                              reverse = TRUE)) +
  theme(legend.title = element_text(size = 10.5))

plot2


plot3 <- 
  plot2 +
  geom_point(data = pcps_fruto12_sp, 
             mapping = aes(x = pcps.1, y = pcps.2),
             shape = 4,
             size = 2)

plot3

# THIS FIGURE WAS EXPORTED AND CLADES' NAMES AND LIMITS WERE INSERTED USING MACBOOK KEYNOTE APP

ggsave(here("figs", "Figure_4.jpg"),
       plot3,
       dpi = 1200,
       units = "in",
       height = 5.95,
       width = 7.95)


#