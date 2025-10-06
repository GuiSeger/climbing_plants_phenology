rm(list = ls(envir = globalenv()), envir = globalenv())
gc()

list_of_packages <- c("tidyverse", "reshape", "lubridate", "scales", "gridExtra",
                      "grid", "readxl", "cowplot", "here")

lapply(list_of_packages, library, character.only = TRUE)


# LOAD DATA
path_datafile <- here('data/Seger_etal_JVS_data.RData')
load(path_datafile)

#
# DATA ENTRY ####

dates <- seq(as.Date("2007-09-01"), as.Date("2009-08-01"), by="months")

# READING AND SUMMARIZING FLOWERING

flor

flor2 <- flor %>% 
  column_to_rownames(var = "species") %>%
  setNames(dates) %>% 
  mutate(species = flor$species)
flor2

flor3 <- gather(flor2, "month", "ind", -species)
flor3

flor4 <- flor3 %>% 
  mutate(presence = ind>0) %>% 
  group_by(month) %>% 
  summarise(spp_flor = (length(species[presence])*100/78),
            ind_flor = (sum(ind))*100/3197)
flor4
print(flor4, n = Inf)


# READING AND SUMMARIZING FRUITING

fruit

fruit2 <- fruit %>% 
  column_to_rownames(var = "species") %>% 
  setNames(dates) %>% 
  mutate(species = fruit$species)
fruit2

fruit3 <- gather(fruit2, "month", "ind", -species)
fruit3

fruit4 <- fruit3 %>% 
  mutate(presence = ind>0) %>% 
  group_by(month) %>% 
  summarise(spp_mad = (length(species[presence])*100/78),
            ind_mad = (sum(ind))*100/3197)
fruit4
print(fruit4, n = Inf)

# READING AND SUMMARIZING SPECIES PEAKS PER MONTH

peaks

# JOINING DATA
data <- flor4 %>% 
  left_join(fruit4, by = "month") %>% 
  left_join(peaks, by = "month")
data


# VIEWING ENVIRONMENTAL DATA
env
# CO - CANOPY OPENNESS
# D - DAYLENGTH
# T - TEMPERATURE
# R - RAINFALL

Sys.setlocale(category = "LC_ALL", locale = "english")


#
# GRAPHICS ####

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


# SPECIES % GRAPHIC

graph_spp <- subset(data, select = c(month, spp_flor, spp_mad))
graph_spp <- as.data.frame(graph_spp)
graph_spp <- melt(graph_spp, measure.vars = c("spp_mad", "spp_flor"), id.vars = "month")
levels(graph_spp$variable)[levels(graph_spp$variable) == "spp_flor"] <- "Flowering"
levels(graph_spp$variable)[levels(graph_spp$variable) == "spp_mad"] <- "Fruiting"
graph_spp$month <- as.Date(graph_spp$month)
str(graph_spp)

gg_spp <- ggplot(data = graph_spp)+
  theme_bw() +
  geom_line(aes(x = month, y = value, linetype = variable), size = 1) +
  geom_rect(data = rect_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.15) +
  geom_point(aes(x = month, y = value, shape = variable, color = variable, size = variable), fill = "white") +
  scale_color_manual(values = c("black", "black"), breaks = c("Flowering", "Fruiting")) + 
  scale_shape_manual(values = c(21, 16), breaks = c("Flowering", "Fruiting")) +
  scale_size_manual(values = c(2.6, 2.5), breaks = c("Flowering", "Fruiting")) +
  scale_linetype_discrete(breaks = c("Flowering", "Fruiting")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", margin = margin(r = 5), size = 9.5),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.06, "in"))+
  scale_x_date(breaks = seq(as.Date("2007-09-01"), as.Date("2009-09-01"), by = "2 months"),
               date_minor_breaks = 'month',
               labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE))),
               expand = c(0.0, 0.0))+
  scale_y_continuous(expand = c(0.03, 0.03), limits = c(0, 61), breaks = seq(0, 61, 10))+
  labs(y = "% of species")

gg_spp


# INDIVIDUALS % GRAPHIC

graph_ind <- subset(data, select = c(month, ind_flor, ind_mad))
graph_ind <- as.data.frame(graph_ind)
graph_ind <- melt(graph_ind, measure.vars = c("ind_mad", "ind_flor"), id.vars = "month")
levels(graph_ind$variable)[levels(graph_ind$variable) == "ind_flor"] <- "Flowering"
levels(graph_ind$variable)[levels(graph_ind$variable) == "ind_mad"] <- "Fruiting"
graph_ind$month <- as.Date(graph_ind$month)
str(graph_ind)

gg_ind <- ggplot(data = graph_ind) +
  theme_bw() +
  geom_line(aes(x = month, y = value, linetype = variable), size = 1) +
  geom_rect(data = rect_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.15) +
  geom_point(aes(x = month, y = value, shape = variable, color = variable, size = variable), fill = "white") +
  scale_color_manual(values = c("black", "black")) + 
  scale_shape_manual(values = c(16, 21))+
  scale_size_manual(values = c(2.6, 2.5)) +
  theme(legend.title = element_blank(), legend.key.size =  unit(0.18, "in"),
        legend.position = c(0.193, 0.85), legend.text = element_text(size = 8.5),
        legend.margin = margin(t = -0.23, r = 0.05, b = 0.02, l = 0, unit = 'cm'),
        legend.key.width = unit(0.32, "in"),
        legend.key.height = unit(0.7, 'lines'),
        legend.spacing.x = unit(.1, "cm"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", margin = margin(r = 5), size = 9.5),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.06, "in"))+
  scale_x_date(breaks = seq(as.Date("2007-09-01"), as.Date("2009-09-01"), by = "2 months"),
               date_minor_breaks = 'month',
               labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE))),
               expand = c(0.0, 0.0)) +
  scale_y_continuous(expand = c(0.03, 0.03), limits = c(0, 35), breaks = seq(0, 35, 10)) +
  labs(y = "% of individuals")

gg_ind


# SPECIES PEAKS % GRAPHIC

graph_peaks <- subset(data, select = c(month, peak_flor, peak_fruto))
graph_peaks <- as.data.frame(graph_peaks)
graph_peaks <- melt(graph_peaks, measure.vars = c("peak_fruto", "peak_flor"), id.vars = "month")
levels(graph_peaks$variable)[levels(graph_peaks$variable) == "peak_flor"] <- "Flower"
levels(graph_peaks$variable)[levels(graph_peaks$variable) == "peak_fruto"] <- "Fruit"
graph_peaks$month <- as.Date(graph_peaks$month)
str(graph_peaks)


gg_peak <- ggplot(data = graph_peaks) +
  theme_bw() +
  geom_line(aes(x = month, y = value, linetype = variable), size = 1) +
  geom_rect(data = rect_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.15)+
  geom_point(aes(x = month, y = value, shape = variable, color = variable, size = variable), fill = "white") +
  scale_color_manual(values = c("black", "black")) + 
  scale_shape_manual(values = c(16, 21))+
  scale_size_manual(values = c(2.6, 2.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9.3, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.y = element_text(face = "bold", margin = margin(r = 5), size = 9.5),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.06, "in")) +
  scale_x_date(date_minor_breaks = 'month',
               breaks = seq(as.Date("2007-09-01"), as.Date("2009-09-01"), by = "2 months"),
               expand = c(0.0, 0.0),
               labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE))))+
  scale_y_continuous(expand = c(0.03, 0.03), limits = c(0, 24), breaks = seq(0, 24, 4))+
  labs(y = "% of species peaks")

gg_peak


# CANOPY OPENNESS AND DAYLENGTH GRAPHIC

ylim.prim2 <- c(0, 64)
ylim.sec2 <- c(8, 15)
d <- diff(ylim.prim2)/diff(ylim.sec2)
c <- ylim.prim2[1] - d*ylim.sec2[1] 


gg_clima2 <- ggplot(data = env) +
  theme_bw() +
  geom_rect(data = rect_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.15) +
  geom_bar(aes(x = month, y = CO), position = "stack", stat = "identity",
           fill = "darkorange2", alpha = .8) +
  geom_line(aes(x = month, y = c+D*d), size = .8) +
  geom_point(aes(x = month, y = c+D*d), fill = "white") +
  scale_color_manual(values = c("black", "black")) + 
  scale_shape_manual(values = c(16)) +
  scale_size_manual(values = c(2.6)) +
  theme(legend.position="none", axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8.5, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold",margin = margin(r = 5), size = 8.6),
        axis.title.y.right = element_text(margin = margin(l = 6)), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.06, "in"))+
  scale_x_date(date_minor_breaks = 'month',
               breaks = seq(as.Date("2007-09-01"), as.Date("2009-09-01"), by = "2 months"),
               expand = c(0.0, 0.0),
               labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE)))) +
  scale_y_continuous(name = "Canopy\nopenness (%)", breaks = seq(0, 60, 10),
                     expand = c(0, 0.03), limits = c(0, 64),
                     sec.axis = sec_axis(trans=~(.-c)/d, name = 'Daylength (h)', 
                                         breaks = seq(8, 15, 2)))


gg_clima2


# TEMPERATURE DAY VARIATION GRAPHIC

data_temp
str(data_temp)

data_temp <- data_temp %>% 
  select(c(data, temp_media)) %>% 
  filter(data > "2007-08-31" & data < "2009-08-31")

str(data_temp)

data_temp$data <- as.Date(data_temp$data)
day(data_temp$data) = 1

str(data_temp)

temp<-
  ggplot(data = data_temp) +
  theme_bw() +
  geom_rect(data = rect_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.15) +
  geom_point(aes(x = data, y = temp_media), size = 1.1, alpha = 0.15) + 
  stat_summary(aes(x = data, y = temp_media), fun = mean, geom = "line", size = 0.75) +
  stat_summary(aes(x = data, y = temp_media), fun = mean, geom = "point", size = 1.8, 
               shape = 21, colour = "white", fill = "black", stroke = .3) +
  scale_x_date(breaks = seq(as.Date("2007-09-01"), as.Date("2009-09-01"), by = "2 months"),
               date_minor_breaks = 'month',
               labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE))),
               expand = c(0.0, 0.0)) +
  scale_y_continuous(expand = c(0.03, 0.03), limits = c(0, 24), breaks = seq(0, 24, 4)) +
  labs(y = "Mean\ntemperature (°C)", x = NULL) + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8.4, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold",margin = margin(r = 5), size = 8.2),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.06, "in"))

temp


# RAINFALL VARIATION GRAPHIC

data_rain
str(data_rain)

data_rain <- data_rain %>% 
  filter(data > as.POSIXct("2007-08-31") & data < as.POSIXct("2009-08-31")) %>% 
  replace_na(list(pluv = 0))

str(data_rain)

data_rain$data<- as.Date(data_rain$data)
day(data_rain$data) = 1

str(data_rain)

rain <- ggplot(data = env) +
  theme_bw() +
  geom_rect(data = rect_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.15) +
  geom_bar(aes(x = month, y = R), position = "stack", 
           stat = "identity", fill = "#84C6D8", alpha = .9)+
  geom_point(data = data_rain, aes(x = data, y = pluv), shape = 95, size = 6.2, alpha = .33)+
  scale_x_date(breaks = seq(as.Date("2007-09-01"), as.Date("2009-09-01"), by = "2 months"),
               date_minor_breaks = 'month',
               labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                            paste(month(x, label = TRUE), "\n", year(x)), 
                                            paste(month(x, label = TRUE))),
               expand = c(0.0, 0.0)) +
  scale_y_continuous(expand = c(0, 0.03), limits = c(0, 250), breaks=seq(0, 250, 40)) +
  labs(y = "Rainfall (mm)", x = NULL) + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8.4, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", margin = margin(r = 1.5), size = 8.6),
        axis.title.y.right = element_text(margin = margin(l = 6)), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.06, "in"))

rain


# JOINING GRAPHICS IN A SINGLE FIGURE ####

pheno.plot <- plot_grid(rain, temp, gg_clima2, gg_spp, gg_ind, gg_peak,
                        align = 'v',
                        ncol = 1,
                        rel_heights = c(1.15/10, 1.25/10, 1.2/10, 2./10, 2./10, 2.4/10))

pheno.plot2 <- grid.arrange(pheno.plot,
                            bottom = textGrob("Month-Year",
                                              gp = gpar(fontface = 2, fontsize = 9.5)))

ggsave(here('figs/Figure_1.jpg'),
       pheno.plot2,
       dpi = 900,
       units="cm",
       height = 21,
       width = 14)

