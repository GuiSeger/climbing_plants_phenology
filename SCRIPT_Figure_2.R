rm(list = ls(envir = globalenv()), envir = globalenv())
gc()

BiocManager::install("ggtree")
packageVersion("ggtree")

list_of_packages <- c("tidyverse", "readr", "ggridges", "here", "lubridate", "ggtree",
                      "gridExtra", "grid", "RColorBrewer", "readxl", "phytools")
lapply(list_of_packages, library, character.only = TRUE)

# LOAD DATA
path_datafile <- here('data/Seger_etal_JVS_data.RData')
load(path_datafile)

#
# PHYLOGENETIC TREE ####

tree
position <- which(tree$tip.label == "Orthosia_scoparia")
tree$tip.label[position] <- "Orthosia_subulata"
tree$tip.label
tree$tip.label <- tolower(tree$tip.label)

# REORDERING TIPLABELS TO MODIFY THE TREE
spp <- c("Mikania_campanulata",
         "Mikania_ternata",
         "Mikania_parodii",
         "Mikania_micrantha",
         "Mikania_cordifolia",
         "Mikania_paranensis",
         "Mikania_oreophila",
         "Mikania_laevigata",
         "Mikania_orleansensis",
         "Mikania_involucrata",
         "Mikania_hirsutissima",
         "Mikania_burchellii",
         "Calea_pinnatifida",
         "Baccharis_anomala",
         "Baccharis_trinervis",
         "Lepidaploa_balansae",
         "Piptocarpha_ramboi",
         "Mutisia_campanulata",
         "Mutisia_speciosa",
         "Valeriana_scandens",
         "Manettia_paraguariensis",
         "Manettia_pubescens",
         "Manettia_verticillata",
         "Manettia_tweedieana",
         "Galium_hypocarpium",
         "Oxypetalum_pedicellatum",
         "Oxypetalum_mosenii",
         "Oxypetalum_wightianum",
         "Oxypetalum_pannosum",
         "Orthosia_urceolata",
         "Orthosia_subulata",
         "Orthosia_virgata",
         "Strychnos_brasiliensis",
         "Solanum_flaccidum",
         "Solanum_laxum",
         "Solanum_inodorum",
         "Convolvulus_crenatifolius",
         "Dolichandra_unguiscati",
         "Dolichandra_uncata",
         "Amphilophium_crucigerum",
         "Tanaecium_selloi",
         "Myriopus_paniculatus",
         "Myriopus_breviflorus",
         "Seguieria_aculeata",
         "Tragia_volubilis",
         "Passiflora_caerulea",
         "Passiflora_actinia",
         "Anchietea_pyrifolia",
         "Heteropterys_aenea",
         "Heteropterys_intermedia",
         "Tetrapterys_phlomoides",
         "Cayaponia_pilosa",
         "Cayaponia_palmata",
         "Cayaponia_diversifolia",
         "Apodanthera_laciniosa",
         "Begonia_fruticosa",
         "Rubus_sellowii",
         "Rubus_erythroclados",
         "Celtis_lancifolia",
         "Senegalia_velutina",
         "Senegalia_nitidifolia",
         "Piptadenia_affinis",
         "Lathyrus_nervosus",
         "Canavalia_bonariensis",
         "Dalbergia_frutescens",
         "Urvillea_ulmacea",
         "Serjania_meridionalis",
         "Serjania_glabrata",
         "Tropaeolum_pentaphyllum",
         "Cissus_verticillata",
         "Clematicissus_striata",
         "Clematis_denticulata",
         "Clematis_bonariensis",
         "Cissampelos_pareira",
         "Dioscorea_multiflora",
         "Dioscorea_subhastata",
         "Smilax_cognata",
         "Bomarea_edulis"
)

spp <- tolower(spp)

tree2 <- rotateConstr(tree, spp)
tree3 <- read.tree(text = write.tree(tree2))
tree3$tip.label


#
# RIDGE PLOT OF FLOWERING PHENOLOGY ####

# FLOWERING DATA
flor # matriz bruta # CRIAR rdata

flor <- flor %>% 
  column_to_rownames(var = "species") 

# INDIVIDUALS % WITHIN EACH SPECIES
flor2 <- (t(apply(flor, 1, function(x) x/max(x))))*100
head(flor2)
class(flor2)
flor2 <- data.frame(t(flor2))
str(flor2)
colnames(flor2) == tree3$tip.label
flor2 <- flor2[, tree3$tip.label]
colnames(flor2) == tree3$tip.label

dates<-seq(as.Date("2007-09-01"), as.Date("2009-08-01"), by="months")

flor3 <- cbind(month = dates, flor2)
head(flor3)
class(flor3)
str(flor3)
flor3 <- as_tibble(flor3)


test <- gather(flor3, key = "species", value = "freq", -month)
test
test$species
test
test$freq
glimpse(test)
test <- test %>% map_df(rev)


families2
families2 <- tibble::column_to_rownames(families2, var = "species")
families2 <- families2[unique(test$species), ]
families2 <- data.frame(families2)
colnames(families2) <- "x"
head(families2)

families3 <- rep(families2$x, each = 24)

test <- add_column(test, family = families3)
test

Sys.setlocale(category = "LC_ALL", locale = "english")

# GRAPHIC

# GRAY BACKGROUND AUTUMN AND WINTER

# Autumn at southern hemisphere march 21st to june 20th
# Winter at southern hemisphere june 21st to september 20th

rect_df <- data.frame(x1 = c("2007-08-15", "2008-03-21", "2009-03-21"),
                      x2 = c("2007-09-20", "2008-09-20", "2009-08-15"),
                      y1 = c(-Inf, -Inf, -Inf),
                      y2 = c(Inf, Inf, Inf))
rect_df$x1 <- as.Date(rect_df$x1, format = ("%Y-%m-%d"))
rect_df$x2 <- as.Date(rect_df$x2, format = ("%Y-%m-%d"))
str(rect_df)


colourCount <- length (unique(test$family))
getPalette<-colorRampPalette(brewer.pal(12,"Paired"))


flower_ridge <- ggplot(data = test) +
  geom_rect(data = rect_df,
            aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            alpha = 0.15)+
  geom_density_ridges(aes(month, y = fct_rev(as_factor(species)),
                          height = freq,
                          fill = family),
                      stat = "identity",
                      scale = .8,
                      color = NA)+
  scale_fill_manual(values = getPalette(colourCount)) +
  guides(fill = guide_legend(ncol = 1, title = "Families")) +
  labs(x = "Month", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.key.size = unit(.3, 'cm'),
        legend.title = element_text(size = 8)) +
  scale_x_date(date_minor_breaks = 'month',
               breaks = seq(as.Date("2007-09-01"), as.Date("2009-09-01"), by = "2 months"),
               expand = c(-0.02,0.0),
               labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x),
                                            paste(month(x, label = TRUE), "\n", year(x)),
                                            paste(month(x, label = TRUE))))
flower_ridge


# PHYLOTREE

tree2 <- tree3
tree2$tip.label <- Hmisc::capitalize(tree2$tip.label)

phylotree <- ggtree(tree2, ladderize = F) +
  geom_tiplab(align = T, 
              offset = 50, 
              linesize = 0,
              size = 2.5,
              hjust = 1,
              nudge_y = .2, 
              geom = "label",
              fontface = 3,
              fill = "white",
              label.padding = unit(0.02, "lines"),
              label.size = 0)+
  geom_treescale(width = 10, fontsize = 2.5) +
  xlim(0, 188)

phylotree

# JOINING PHYLOTREE AND FLOWERING RIDGEPLOT

panel <- grid.arrange(phylotree + theme(plot.margin = margin(3, -25, 30, 0)),
                      arrangeGrob(flower_ridge + theme(axis.text.y = element_blank()),
                                  nullGrob(),
                                  heights=c(100, 0.02)),
                      ncol = 2,
                      widths = 2:3)


ggsave(here("figs", "Figure_2A.jpg"), 
       panel,  
       dpi = 600, 
       units="cm", 
       height = 21, 
       width = 29.7)


#
# RIDGE PLOT OF FRUITING PHENOLOGY ####

# FRUITING DATA

fruit

fruit <- fruit %>% 
  column_to_rownames(var = "species") 

fruit2 <- (t(apply(fruit, 1, function(x) x/max(x))))*100
head(fruit2)
class(fruit2)
fruit2 <- data.frame(t(fruit2))
str(fruit2)
colnames(fruit2) == tree3$tip.label
fruit2 <- fruit2[, tree3$tip.label]
colnames(fruit2) == tree3$tip.label

dates<-seq(as.Date("2007-09-01"), as.Date("2009-08-01"), by="months")

fruit3 <- cbind(month = dates, fruit2)
head(fruit3)
class(fruit3)
str(fruit3)
fruit3 <- as_tibble(fruit3)


test <- gather(fruit3, key = "species", value = "freq", -month)
test
test$species
test
test$freq
glimpse(test)
test <- test %>% map_df(rev)

test <- add_column(test, family = families3)
test

Sys.setlocale(category = "LC_ALL", locale = "english")


# GRAPHIC

colourCount <- length (unique(test$family))
getPalette<-colorRampPalette(brewer.pal(12,"Paired"))


fruit_ridge <- ggplot(data = test) +
  geom_rect(data = rect_df,
            aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            alpha = 0.15)+
  geom_density_ridges(aes(month, y = fct_rev(as_factor(species)),
                          height = freq,
                          fill = family),
                      stat = "identity",
                      scale = .8,
                      color = NA)+
  scale_fill_manual(values = getPalette(colourCount)) +
  guides(fill = guide_legend(ncol = 1, title = "Families")) +
  labs(x = "Month", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.key.size = unit(.3, 'cm'),
        legend.title = element_text(size = 8)) +
  scale_x_date(date_minor_breaks = 'month',
               breaks = seq(as.Date("2007-09-01"), as.Date("2009-09-01"), by = "2 months"),
               expand = c(-0.02,0.0),
               labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x),
                                            paste(month(x, label = TRUE), "\n", year(x)),
                                            paste(month(x, label = TRUE))))
fruit_ridge


# JOINING PHYLOTREE AND RIDGEPLOT

panel <- grid.arrange(phylotree + theme(plot.margin = margin(3, -25, 30, 0)),
                      arrangeGrob(fruit_ridge + theme(axis.text.y = element_blank()),
                                  nullGrob(),
                                  heights=c(100, 0.02)),
                      ncol = 2,
                      widths = 2:3)

# THIS FIGURE WAS EXPORTED AND CLADES' NAMES IN THE PHYLOGENY WERE INSERTED USING MACBOOK KEYNOTE APP

ggsave(here("figs", "Figure_2B.jpg"), 
       panel,  
       dpi = 600, 
       units="cm", 
       height = 21, 
       width = 29.7)
