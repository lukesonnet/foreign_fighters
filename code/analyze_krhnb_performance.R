#----------
# Authors: George Derpanopoulos and Luke Sonnet
#
# Project: Foreign Fighter Supply using KRHNB
#
# Purpose: Testing KRHNB OOS performance
#----------

## Preamble
rm(list=ls())
setwd('C:/Dropbox/research/foreign_fighters/')
library(foreign)
library(KRLS)
library(randomForest)
library(e1071)
library(MASS)
library(ggplot2)
library(gridExtra)
source("krhnb/krhnb.R")
source("code/testCV.R")
library(latex2exp)
library(extrafont)
#font_install("fontcm")
loadfonts(device="win")
# So windows knows where ghostscript is
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.18/bin/gswin64c.exe")

#load("savedata/krhnb_performance.RData")
#save.image("savedata/krhnb_performance.RData")

## Read data
df <- read.csv("data/foreignFightersImputed.csv", as.is=T)

## Setting up matrices
names(df)
excludeVars <- c("ctryName", "ctryCode", "icsrMean", "lnarea")
X <- df[, !(colnames(df) %in% excludeVars)]
y <- round(df$icsrMean)

#----------
# Foreign Fighters
#----------

## Seeds only necessary because results for random forests fluctuates ever
## so slightly
set.seed(20160415)
foldslcv <- chunk(sample(nrow(X)), 163)
set.seed(20160415)
ff.lcv <- compareKrhnb(y = y, X = X, foldslcv, lambda1 = 0.002, lambda2 = 0.85, hnb = FALSE)
apply(ff.lcv$totalMseDF, 2, function(x) sqrt(mean(x)))
apply(ff.lcv$totalMseDF, 2, function(x) mean(x))
#save(ff.lcv, file = "savedata/krhnb.loocv.RData")
#load("savedata/krhnb.loocv.RData")

## Plot fit
setdiff(names(ff.lcv2$totalPredDF), c("yhat.hnb", "fold"))
plot.df <- reshape2::melt(ff.lcv2$totalPredDF[, setdiff(names(ff.lcv2$totalPredDF), c("yhat.hnb", "fold", "yhat.trunc0.k", "yhat.lm"))],
                          id.vars = c("testy"),
                          variable.name = "method",
                          value.name = "ypred")

table(plot.df$method)

p <- ggplot(plot.df, aes(y = ypred, x = testy, color = method, shape = method)) + 
  geom_jitter(size = 2, alpha = 0.6) + 
  geom_abline(slope = 1, intercept = 0, lty = 2) + 
  ylab(TeX("$\\hat{Y}$")) + xlab(TeX("$Y$")) + 
  scale_color_manual(name = "Method",
                     values = c("yhat.k" = "red",
                                "yhat.knb" = "green",
                                "yhat.rf" = "blue"),
                     labels = c("KRLS", "KRHNB", "RForest")) +
  scale_shape_manual(name = "Method",
                     values = c("yhat.k" = 15,
                                "yhat.knb" = 16,
                                "yhat.rf" = 17),
                     labels = c("KRLS", "KRHNB", "RForest")) +
  theme_bw() + 
  theme(plot.title = element_text( hjust = 0),
        text=element_text(family="CM Roman"),
        legend.position="bottom")
p
ggsave(p, file = "tex/tabs_figs/fullpredict.pdf", width = 6.5, height = 4)
embed_fonts("tex/tabs_figs/fullpredict.pdf")

p2 <- ggplot(plot.df[plot.df$testy < 150, ], aes(y = ypred, x = testy, color = method, shape = method)) + 
  geom_jitter(size = 2, alpha = 0.6) + 
  geom_abline(slope = 1, intercept = 0, lty = 2) + 
  ylim(ggplot_build(p)$panel$ranges[[1]]$y.range) + 
  ylab(TeX("$\\hat{Y}$")) + xlab(TeX("$Y$")) + 
  scale_color_manual(name = "Method",
                     values = c("yhat.k" = "red",
                                "yhat.knb" = "green",
                                "yhat.rf" = "blue"),
                     labels = c("KRLS", "KRHNB", "RForest")) +
  scale_shape_manual(name = "Method",
                     values = c("yhat.k" = 15,
                                "yhat.knb" = 16,
                                "yhat.rf" = 17),
                     labels = c("KRLS", "KRHNB", "RForest")) +
  theme_bw() + 
  theme(text=element_text(family="CM Roman"))
p2
ggsave(p2, file = "tex/tabs_figs/predict150.pdf", width = 6.5, height = 4)
embed_fonts("tex/tabs_figs/predict150.pdf")


## From http://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

## Figure 2
pJoin <- grid.arrange(arrangeGrob(p + theme(legend.position="none"),
                      p2 + theme(legend.position="none"),
                      nrow = 1),
                      g_legend(p), nrow = 2, heights=c(10, 1))
ggsave(pJoin, file = "tex/tabs_figs/predictJoin.pdf", width = 7.5, height = 3.5)
embed_fonts("tex/tabs_figs/predictJoin.pdf")


#----------
# Burgoon 2006
#----------

bdf <- read.dta("data/burgoon/burgoon_2006.dta")
burgoon.nb <- glm.nb(terrorinclead ~ welfarelog + govleft + democ + poplog + 
                       govcap + conflict + tradelog + europe + africa + asia + america + 
                       terrorinc,
                     data = bdf)
summary(burgoon.nb)

X <- burgoon.nb$model[, c("terrorinc", "welfarelog", "govleft", "democ", "poplog", "govcap", "conflict",
                          "tradelog", "europe", "africa", "asia", "america")]
y <- burgoon.nb$model$terrorinclead

set.seed(20160214)
folds <- chunk(sample(nrow(X)), 5)
set.seed(20160214)
b.out <- compareKrhnb(y = y, X = X, folds, lambda1s = c(0.075, 0.1, 0.125),
                      lambda2s = c(3:5), lambdafolds = 10)
apply(b.out$totalMseDF, 2, function(x) sqrt(mean(x)))
apply(b.out$totalMseDF, 2, function(x) mean(x))
#save(b.out, file = "savedata/burgoon.perf.RData")