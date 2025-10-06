rm(list = ls(envir = globalenv()), envir = globalenv())
gc()

library(here)
library(corrplot)
library(Hmisc)

# LOAD DATA
path_datafile <- here('data/Seger_etal_JVS_data.RData')
load(path_datafile)

# CORRELATION BETWEEN ENVIRONMENTAL VARIABLES
correl.b <- rcorr(as.matrix(env2), type = "pearson")
correl.b

# PLOT CORRELATIONS
# Insignificant correlations are leaved blank
myfun <- function(ff){
  corrplot(ff$r, 
           type = "upper", 
           order = "alphabet",
           addCoef.col = "black", 
           tl.col = "black", 
           tl.srt = 45,
           number.cex = .7,
           p.mat = correl.b$P, 
           sig.level = 0.05, 
           insig = "blank", 
           diag = F)
  recordPlot()
}

ff <- rcorr(as.matrix(env2))
myfun(ff) # shows the plot
corplot <- myfun(ff)

#SAVE THE GRAPHIC
pdf(here("figs", "Figure_S5.pdf"),
    width = 7,
    height = 7)
replayPlot(corplot)
graphics.off()

