rm(list = ls(envir = globalenv()), envir = globalenv())
gc()

list_of_packages <- c("tidyverse", "reshape", "lubridate", "scales", "gridExtra",
                      "grid", "readxl", "cowplot", "plyr", "here")

lapply(list_of_packages, library, character.only = TRUE)

# LOAD DATA
path_datafile <- here('data/Seger_etal_JVS_data.RData')
load(path_datafile)

dates <- seq(as.Date("2007-09-01"), as.Date("2009-08-01"), by="months")
dates <- c("Espécie", as.character(dates))
dates

fruit3 <- fruit %>% 
  setNames(dates)
fruit3

data <- disp %>% 
  left_join(fruit3, by = "Espécie")  
data


data2 <- data%>%
  pivot_longer(
    cols = '2007-09-01':'2009-08-01',
    names_to = "month",
    values_to = "abundance"
  ) %>%
  mutate(month = as.Date(month, format = "%Y-%m-%d"))
data2

# NUMBER OF SPECIES BY SYNDROME
data3 <- data2%>%
  filter(!abundance == 0)

data4 <- ddply(data3,
               .(Dispersão, month),
               summarize,
               spp = (length(Espécie)*100)/78,
               ind = (sum(abundance)*100)/3197,
               .drop = FALSE)
data4

data5 <- data4 %>%
  pivot_wider(names_from = Dispersão, values_from = c(spp, ind))
str(data5)

data5 <- as.data.frame(data5)

# path_datafile <- here('R/Mean-Modal angles Calculation/Data Modal angle.xlsx')
# 
# peaks <- read_xlsx(path_datafile, sheet = "Peaks modal", range = "AD31:AG55")

peaks_disp
data6 <- cbind(data5, peaks_disp[,2:4])
data6
str(data6)

data <- data6

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


#
# SPECIES % GRAPHIC

graph_spp <- subset(data,select = c(month, spp_A, spp_Z, spp_T))
graph_spp <- as.data.frame(graph_spp)
graph_spp <- melt(graph_spp,measure.vars = c("spp_Z", "spp_A", "spp_T"), id.vars = "month")
levels(graph_spp$variable)[levels(graph_spp$variable) == "spp_A"] <- "Anemo"
levels(graph_spp$variable)[levels(graph_spp$variable) == "spp_Z"] <- "Zoo"
levels(graph_spp$variable)[levels(graph_spp$variable) == "spp_T"] <- "Auto"
str(graph_spp)


gg_spp <- ggplot(data = graph_spp) +
  theme_bw() +
  geom_line(aes(x = month, y = value, linetype = variable, colour = variable), size = 1) +
  geom_rect(data = rect_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = .15) +
  geom_point(aes(x = month, y = value, shape = variable, color = variable, size = variable), fill = "white") +
  scale_color_manual(values = c('black', 'black', 'gray50'), breaks = c("Zoo", "Anemo", "Auto")) +
  scale_shape_manual(values = c(16, 21, 16), breaks = c("Zoo", "Anemo", "Auto")) +
  scale_size_manual(values = c(2.5, 2.6, 2.6), breaks = c("Zoo", "Anemo", "Auto")) +
  scale_linetype_manual(values = c("solid", "dashed", "solid"), breaks = c("Zoo", "Anemo", "Auto")) +
  theme(legend.title = element_blank(), 
        legend.key.size =  unit(0.2, "in"),
        legend.position = c(0.15, 0.81),
        legend.text = element_text(size = 8.2),
        legend.margin = margin(t = -0.23, r = 0.05, b = -0.03, l = -0.01, unit = 'cm'),
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
        axis.ticks.length = unit(0.06, "in")) +
  scale_x_date(breaks = seq(as.Date("2007-09-01"), as.Date("2009-09-01"), by = "2 months"),
               date_minor_breaks = 'month',
               labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE))),
               expand = c(0.0, 0.0)) +
  scale_y_continuous(expand = c(0.03, 0.03), limits = c(0, 32.5), breaks = seq(0, 32.5, 5)) +
  labs(y = " % of species")

gg_spp


#
# INDIVIDUALS % GRAPHIC

graph_ind <- subset(data, select = c(month, ind_A, ind_Z, ind_T))
graph_ind <- as.data.frame(graph_ind)
graph_ind <- melt(graph_ind, measure.vars = c("ind_Z", "ind_A", "ind_T"), id.vars = "month")
levels(graph_ind$variable)[levels(graph_ind$variable)=="ind_A"] <- "Anemo"
levels(graph_ind$variable)[levels(graph_ind$variable)=="ind_Z"] <- "Zoo"
levels(graph_ind$variable)[levels(graph_ind$variable)=="ind_T"] <- "Auto"
str(graph_ind)


gg_ind <- ggplot(data = graph_ind) +
  theme_bw() +
  geom_line(aes(x = month, y = value, linetype = variable, colour = variable), size = 1) +
  geom_rect(data = rect_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.15) +
  geom_point(aes(x = month, y = value, shape = variable, color = variable, size = variable), fill = "white") +
  scale_color_manual(values = c('black', 'black', 'gray50'), breaks = c("Zoo", "Anemo", "Auto")) +
  scale_shape_manual(values = c(16, 21, 16), breaks = c("Zoo", "Anemo", "Auto")) +
  scale_size_manual(values = c(2.5, 2.6, 2.6), breaks = c("Zoo", "Anemo", "Auto")) +
  scale_linetype_manual(values = c("solid", "dashed", "solid"), breaks = c("Zoo", "Anemo", "Auto")) +
  theme(legend.position = "none", 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", margin = margin(r = 5), size = 9.5),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.06, "in")) +
  scale_x_date(breaks = seq(as.Date("2007-09-01"), as.Date("2009-09-01"), by = "2 months"),
               date_minor_breaks = 'month',
               labels = function(x) if_else(is.na(lag(x)) | !year(lag(x)) == year(x), 
                                          paste(month(x, label = TRUE), "\n", year(x)), 
                                          paste(month(x, label = TRUE))),
               expand = c(0.0, 0.0)) +
  scale_y_continuous(expand = c(0.03, 0.03), limits = c(0, 16), breaks = seq(0, 16, 2)) +
  labs(y = " % of individuals")

gg_ind


#
# SPECIES PEAKS % GRAPHIC

graph_peaks <- subset(data, select = c(month, peak_A, peak_Z, peak_T))
graph_peaks <- as.data.frame(graph_peaks)
graph_peaks <- melt(graph_peaks, measure.vars = c("peak_Z","peak_A", "peak_T"), id.vars = "month")
levels(graph_peaks$variable)[levels(graph_peaks$variable) == "peak_A"] <- "Anemo"
levels(graph_peaks$variable)[levels(graph_peaks$variable) == "peak_Z"] <- "Zoo"
levels(graph_peaks$variable)[levels(graph_peaks$variable) == "peak_T"] <- "Auto"
str(graph_peaks)


gg_peak <- ggplot(data = graph_peaks) +
  theme_bw() +
  geom_line(aes(x = month, y = value, linetype = variable, colour = variable), size = 1) +
  scale_linetype_manual(values = c("solid", "dashed", "solid"), breaks = c("Zoo", "Anemo", "Auto")) +
  geom_rect(data = rect_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.15) +
  geom_point(aes(x = month, y = value, shape = variable, color = variable, size = variable), fill = "white") +
  scale_color_manual(values = c('black', 'black', 'gray50'), breaks = c("Zoo", "Anemo", "Auto")) +
  scale_shape_manual(values = c(16, 21, 16), breaks = c("Zoo", "Anemo", "Auto")) +
  scale_size_manual(values = c(2.5, 2.6, 2.6), breaks = c("Zoo", "Anemo", "Auto")) +
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
  scale_y_continuous(expand = c(0.03, 0.03), limits = c(0, 13), breaks = seq(0, 13, 2)) +
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
  geom_line(aes(x = month, y = c+D*d),size = .8) +
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

temp <-
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


#
# JOINING GRAPHICS IN A SINGLE FIGURE ####

disp.plot <- plot_grid(rain, temp, gg_clima2, gg_spp, gg_ind, gg_peak,
                       align = 'v',
                       ncol = 1,
                       rel_heights = c(1.15/10, 1.25/10, 1.2/10, 2./10, 2./10, 2.4/10))

disp.plot2 <- grid.arrange(disp.plot,
                           bottom = textGrob("Month-Year", 
                                             gp = gpar(fontface = 2, fontsize = 9.5)))

ggsave(here('figs/Figure_S7.jpg'),
       disp.plot2,
       dpi = 900,
       units = "cm",
       height = 21,
       width = 14)


#
