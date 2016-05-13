#----------
# Authors: George Derpanopoulos and Luke Sonnet
#
# Project: Foreign Fighter Supply using KRHNB
#
# Purpose: Analyze the foreign fighter data
#----------

## Preamble
rm(list=ls())
setwd('C:/Dropbox/research/foreign_fighters/')
## Data cleaning
library(reshape2)
## Model fitting
source("krhnb/krhnb.R")
## Presentation
library(ggplot2)
library(rworldmap)
library(gridExtra)
library(xtable)
library(latex2exp)
library(extrafont)
#font_install("fontcm")
loadfonts(device="win")
# So windows knows where ghostscript is
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.18/bin/gswin64c.exe")

## Function for effect quantiles
bootQuantiles <- function(x) {
  tab <- t(sapply(colnames(x),
                  function(y) quantile(x[, y],
                                       probs = c(0.025, 0.05, 0.5, 0.95, 0.975))))
  tab[order(tab[, "50%"]), ]
}

#----------
# Ingest Data
#----------

## Read data
df <- read.csv("data/foreignFightersImputed.csv", as.is=T)

## Setting up matrices
## Drop area because popdensity = pop / area
excludeVars <- c("ctryName", "ctryCode", "icsrMean", "lnarea")
X <- df[, !(colnames(df) %in% excludeVars)]
y <- round(df$icsrMean)

#----------
# Descriptive plots
#----------

labels <- c("contig" = "Contiguous",
            "sunnipct" = "Sunni Pct",
            "internetUserPer100" = "Internet User per 100",
            "lifeExp" = "Life Expectancy",
            "soc.reg.relig" = "Soc. Reg. Religion",
            "lnrefugeesUnhcr" = "Log Refugees in Country",
            "unempYthTotILO" = "Youth Unemployment",
            "lngdppcPPP" = "Log GDP pc, PPP",
            "migrantStockPctPop" = "Migrants as Pct of Pop",
            "fh_pr" = "Freedom House Pol. Rights",
            "europe" = "Europe",
            "fh_cl" = "Freedom House Civ. Lib.",
            "shiapct" = "Shia Pct",
            "lnrefugeeOrigin" = "Log Refugees from Country",
            "lnpop" = "Log Population",
            "religious.frac" = "Religious Frac.",
            "govt.reg.relig" = "Govt. Reg. Religion",
            "lfpMale" = "Male LFP",
            "govt.fav.relig" = "Govt. Fav. Religion",
            "lnpopDensity" = "Log Pop. Density",
            "lnHomicidesPer100k" = "Log Homicides Per 100k",
            "pctpop15.24" = "Youth Bulge (15-24 pct)",
            "dist" = "Distance (km)")

## Table 2
## Summary statistics
sumStat <- t(apply(cbind(y, X[, c("contig", "europe", "dist", "lnpop",
                                  "lnpopDensity", "pctpop15.24", "sunnipct",
                                  "shiapct", "religious.frac",
                                  "govt.reg.relig", "govt.fav.relig",
                                  "soc.reg.relig", "fh_cl", "fh_pr",
                                  "migrantStockPctPop", "lnrefugeesUnhcr",
                                  "lnrefugeeOrigin", "lngdppcPPP", "lifeExp",
                                  "lfpMale", "unempYthTotILO",
                                  "internetUserPer100", "lnHomicidesPer100k")]),
                   2,
                   function(x) {
  c(mean(x), sd(x), summary(x)[c(1:3, 5, 6)])
}))

labelsY <- c("y" = "Foreign Fighter Supply", labels)

rownames(sumStat) <- labelsY[rownames(sumStat)]
colnames(sumStat)[1:2] <- c("Mean", "Std. Dev.")
print(xtable(sumStat), file = "tex/tabs_figs/sumStat.tex")


## Mapping the outcome

## Add proper iso3 codes
df$isoCode <- df$ctryCode
df$isoCode[df$ctryCode == "ZAR"] <- "COD"
df$isoCode[df$ctryCode == "ROM"] <- "ROU"
df$isoCode[df$ctryCode == "TMP"] <- "TLS"
df$isoCode[df$ctryCode == "WBG"] <- "PSE"
df$isoCode[df$ctryCode == "ZAR"] <- "COD"
df$isoCode[df$ctryCode == "KSV"] <- "XKX"

mapDF <- joinCountryData2Map(df,
                             joinCode = "ISO3",
                             nameJoinColumn = "isoCode",
                             verbose = T)
## Kosovo fails

## Figure 1
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
## Second map is to creat custom ticks on the scale
worldMap2 <- mapCountryData(mapDF,
                            nameColumnToPlot="icsrMean" ,
                            catMethod=seq(0, max(df$icsrMean), by = 500),
                            mapTitle="",
                            addLegend = F)
pdf(file = "tex/tabs_figs/worldMap.pdf", width = 7.5, height = 4.5)
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
## Hack to get the legend looking better
worldMap <- mapCountryData(mapDF,
                           nameColumnToPlot="icsrMean" ,
                           catMethod=seq(0, max(df$icsrMean), by = 1),
                           mapTitle="",
                           addLegend = F,
                           ylim=c(0,50))
do.call(addMapLegend,
        c(worldMap2,
          legendLabels = "all",
          legendShrink = 0.5,
          legendWidth = 0.5,
          legendMar = 3))
do.call(addMapLegend,
        c(worldMap,
          legendLabels = "none",
          legendShrink = 0.5,
          legendWidth = 0.5,
          legendMar = 3))
dev.off()

#----------
# Fit KRHNB
#----------

## Set to false to save time, change to TRUE to re-fit models
fitModels <- FALSE

if (fitModels) {
  ## Seeds not needed here because it uses LOOCV
  # Lambda search manually
  lsearch <- krhnb(X=X, y=y, bootstrap = F,
                   lambda1s = c(0.001, 0.005, .01),
                   lambda2s = c(0.5, 1, 1.5),
                   lambdafolds = 163)
  # lambda1 = 0.001
  # lambda2 = 1
  lsearch2 <- krhnb(X=X, y=y, bootstrap = F,
                    lambda1s = c(0.00075, 0.001, 0.002),
                    lambda2s = c(0.6, 0.75, 0.85, 0.95, 1),
                    lambdafolds = 163)
  #lambda1 = 0.002
  #lambda2 = 0.85
  
  ## Seed from random.org from 1 to .Machine$integer.max
  set.seed(109134803)
  ff.krhnb <- krhnb(X=X, y=y,
                    lambda1 = 0.002, lambda2 = 0.85,
                    bootstrapits = 1000)
  
  save(ff.krhnb, file = "savedata/ff.krhnb.RData")
  
} else{
  load("savedata/ff.krhnb.RData")
}

#----------
# Model output
#----------

## Extract percentile intervals
quantilesAll <- apply(ff.krhnb$bootmfx, 2, quantile, c(0.025, 0.975))
quantilesOnes <- apply(ff.krhnb$bootonesmfx, 2, quantile, c(0.025, 0.975))
quantilesCount <- apply(ff.krhnb$bootcountmfx, 2, quantile, c(0.025, 0.975))

## Create matrix to print table of effects
avgEffectMat <- cbind(c(rbind(plyr::revalue(colnames(ff.krhnb$effectmat), labels),
                              rep("", times = length(labels)))),
                      c(rbind(round(ff.krhnb$meaneffect, 2),
                              paste0("[", round(quantilesAll[1,], 2), ", ", round(quantilesAll[2,], 2), "]"))),
                      c(rbind(round(apply(ff.krhnb$onesmat, 2, mean), 2),
                              paste0("[", round(quantilesOnes[1,], 2), ", ", round(quantilesOnes[2,], 2), "]"))),
                      c(rbind(round(apply(ff.krhnb$countmat, 2, mean), 2),
                              paste0("[", round(quantilesCount[1,], 2), ", ", round(quantilesCount[2,], 2), "]"))))

## Table A1
avgEffectTab <- xtable(avgEffectMat)
print(avgEffectTab, file = "tex/tabs_figs/tabAvgEffect.tex", include.rownames=FALSE)

#----------
# Plotting effects
#----------

## Figure 2
## Build data
plot.df <- data.frame(variable = rep(colnames(ff.krhnb$bootmfx), times = 3),
                      mean = c(ff.krhnb$meaneffect,
                               colMeans(ff.krhnb$onesmat),
                               colMeans(ff.krhnb$countmat)),
                      lwr = c(quantilesAll[1, ],
                              quantilesOnes[1, ],
                              quantilesCount[1, ]),
                      upr = c(quantilesAll[2, ],
                              quantilesOnes[2, ],
                              quantilesCount[2, ]),
                      component = rep(c("Both Components",
                                        "Hurdle Component",
                                        "Count Component"),
                                      each = ncol(ff.krhnb$effectmat)))

## Bootstrap means not far from ML point estimates on full data
## Skewed sampling distributions!
plot(ff.krhnb$meaneffect, apply(ff.krhnb$bootmfx, 2, mean))
abline(0, 1)
plot(colMeans(ff.krhnb$countmat), apply(ff.krhnb$bootcountmfx, 2, mean))
abline(0, 1)

plot.df$variable <- reorder(plot.df$variable,
                            rep(plot.df$mean[1:ncol(ff.krhnb$effectmat)], 3))

## Plot three panels
## Joint effects
pAll <- ggplot(plot.df[plot.df$component == "Both Components",],
               aes(y=mean, x=as.factor(variable), ymin = lwr, ymax = upr)) + 
  geom_hline(yintercept = 0, color = "red", size = 0.6) + 
  geom_point(size = 1.2) + geom_errorbar(width = 0, size = 0.75) +
  scale_x_discrete("", labels = labels) + 
  ylab("Mean Marginal Effect") + 
  ggtitle("Both Components") + 
  coord_flip() +  
  theme_bw() + 
  theme(text=element_text(family="CM Roman"))
pAll

## Hurdle effects
pHurdle <- ggplot(plot.df[plot.df$component == "Hurdle Component",],
                  aes(y=mean, x=as.factor(variable), ymin = lwr, ymax = upr)) + 
  geom_hline(yintercept = 0, color = "red", size = 0.6) + 
  geom_point(size = 1.2) + geom_errorbar(width = 0, size = 0.75) +
  ylab("Mean Marginal Effect") + 
  ggtitle("Hurdle Component") + 
  coord_flip() +
  theme_bw() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(family="CM Roman"))
pHurdle

## Count effects
pCount <- ggplot(plot.df[plot.df$component == "Count Component",],
                 aes(y=mean, x=as.factor(variable), ymin = lwr, ymax = upr)) + 
  geom_hline(yintercept = 0, color = "red", size = 0.6) + 
  geom_point(size = 1.2) + geom_errorbar(width = 0, size = 0.75) +
  ylab("Mean Marginal Effect") + 
  ggtitle("Count Component") + 
  coord_flip() +
  theme_bw() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(family="CM Roman"))
pCount

pJoin <- grid.arrange(pAll, pHurdle, pCount, ncol = 3, widths = c(4.5,2.5,2.5))
ggsave(pJoin, file = "tex/tabs_figs/bootEffects.pdf", width = 10, height = 6)
embed_fonts("tex/tabs_figs/bootEffects.pdf", outfile="tex/tabs_figs/bootEffects.pdf")


#----------
# Heterogeneous Effects
#----------

## Figure 4
## Show pointwise effect histograms
plotHistDF <- melt(ff.krhnb$effectmat[, c("sunnipct", "internetUserPer100",
                                          "govt.reg.relig", "lnHomicidesPer100k")])
plotHistDF$variable <- ifelse(plotHistDF$Var2 == "sunnipct",
                              "Sunni Pct",
                              ifelse(plotHistDF$Var2 == "internetUserPer100",
                                     "Internet User per 100",
                                     ifelse(plotHistDF$Var2 == "govt.reg.relig",
                                            "Govt. Reg. Religion",
                                            "Log Homicides per 100k")))

histPlot <- ggplot(plotHistDF, aes(x = value)) + 
  geom_histogram() +
  facet_wrap( ~ variable, ncol = 2) + 
  ylab("Frequency") +
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        text=element_text(family="CM Roman"))
ggsave(histPlot, file = "tex/tabs_figs/histeffects.pdf", width = 6, height = 4)
embed_fonts("tex/tabs_figs/histeffects.pdf", outfile="tex/tabs_figs/histeffects.pdf")
histPlot

## Figure 5
## Het effects within and outside europe
hetdf1 <- melt(data.frame(europe = ifelse(X[, "europe"] == 1, "Europe", "Elsewhere"),
                          ref = ff.krhnb$effectmat[, "lnrefugeesUnhcr"],
                          int = ff.krhnb$effectmat[, "internetUserPer100"]), id.vars = c("europe"))
hetdf1$variable <- ifelse(hetdf1$variable == "ref", "Log Refugees in Country", "Internet User per 100")
hetplot1 <- ggplot(hetdf1, aes(x=europe, y = value)) +
  geom_boxplot(aes(fill = europe), alpha = 0.5) +
  scale_fill_manual(values = c("Europe" = "darkblue", "Elsewhere" = "lightblue")) +
  ylab("Pointwise Marginal Effects") +
  facet_grid(.~variable) +
  theme_bw() + 
  theme(text=element_text(family="CM Roman"),
        axis.title.x=element_blank(),
        legend.position="None")
hetplot1
ggsave(hetplot1, file = "tex/tabs_figs/hetplot1.pdf", width = 6, height = 3)
embed_fonts("tex/tabs_figs/hetplot1.pdf")

hetplotGovtRegSunni <- qplot(X[, "sunnipct"], ff.krhnb$effectmat[, "govt.reg.relig"]) + 
  xlab("Sunni Pct") + ylab("PWMFX of Govt. Reg. of Religion") +
  stat_smooth(alpha = 0.4) +
  theme_bw() + 
  theme(text=element_text(family="CM Roman"))
hetplotGovtRegSunni
ggsave(hetplotGovtRegSunni, file = "tex/tabs_figs/hetplotGovtRegSunni.pdf", width = 6, height = 3)
embed_fonts("tex/tabs_figs/hetplotGovtRegSunni.pdf")

## Figure A1
## Matrix and plot of all interaction effects
interactMat <- apply(ff.krhnb$effectmat, 2, function(x) {
  lmout <- lm(x ~ scale(as.matrix(X)))
  coef(lmout)[-1]
})

diag(interactMat) <- 0
rownames(interactMat) <- gsub("scale\\(as\\.matrix\\(X\\)\\)", "", rownames(interactMat))
rownames(interactMat) <- plyr::revalue(rownames(interactMat), labels)
colnames(interactMat) <- plyr::revalue(colnames(interactMat), labels)

pdf(file = "tex/tabs_figs/interactPlot.pdf", width = 6, height = 8)
# Drop contig, effects mask all others
corrplot::corrplot(interactMat[rownames(interactMat) != "Contiguous",
                               colnames(interactMat) != "Contiguous"], 
                   is.corr = F, type = "upper", diag = F)
dev.off()

#----------
# Analysis with listwise deletion, and no west bank + gaza
#----------

## Read data
dfRaw <- read.csv("data/foreignFightersPrep.csv", as.is=T)

## Adding Europe and MeNA dummies
dfRaw$europe <- ifelse(dfRaw$ctryName %in% c("Albania", "Armenia", "Austria", "Begium",
                                       "Bosnia and Herzegovina", "Bulgaria",
                                       "Croatia", "Cyprus", "Czech Republic",
                                       "Denmark", "Estonia", "Finland", "France",
                                       "Germany", "Greece", "Hungary", "Iceland",
                                       "Ireland", "Italy", "Latvia", "Lithuania",
                                       "Luxembourg", "Macedonia, FYR,", "Malta",
                                       "Moldova", "Netherlands", "Norway",
                                       "Poland", "Portugal", "Romania", "Serbia",
                                       "Slovak Republic", "Slovenia", "Spain",
                                       "Sweden", "Switzerland", "Ukraine",
                                       "United Kingdom"), 1, 0)

## Return WBG ff to israel
dfRaw$icsrMean[dfList$ctryName == "Israel"] <- sum(dfRaw$icsrMean[dfRaw$ctryName %in% c("Israel", "West Bank and Gaza")])
dfRaw$icsrMean[dfList$ctryName == "West Bank and Gaza"] <- NA 

## Listwise delete
sum(apply(dfRaw, 1, function(x) sum(is.na(x))) ==0) ## 147 complete obs

dfList <- na.omit(dfRaw)

## Setting up matrices
## Drop area because popdensity = pop / area
excludeVars <- c("ctryName", "ctryCode", "icsrMean", "lnarea")
X <- dfList[, !(colnames(dfList) %in% excludeVars)]
y <- round(dfList$icsrMean)

## Switch to keep from refitting each time
fitListwise <- FALSE

if(fitListwise) {
  lsearch <- krhnb(X=X, y=y, bootstrap = F,
                   lambda1s = c(0.001, 0.005, .01, 0.05, 0.1, 0.5),
                   lambda2s = c(0.1, 0.5, 0.75, 1, 1.5, 3),
                   lambdafolds = 147)
  
  lsearch2 <- krhnb(X=X, y=y, bootstrap = F,
                   lambda1s = c(0.00075, 0.001, 0.0015),
                   lambda2s = c(0.3, 0.4, 0.5),
                   lambdafolds = 147)
  
  
  lsearch2 <- krhnb(X=X, y=y, bootstrap = F,
                    lambda1s = 0.001,
                    lambda2s = c(0.39, 0.41),
                    lambdafolds = 147)
  
  ## Seed from random.org from 1 to .Machine$integer.max
  set.seed(205271022)
  ff.krhnb.list <- krhnb(X=X, y=y,
                    lambda1 = 0.01, lambda2 = 0.4,
                    bootstrapits = 1000)

  save(ff.krhnb.list, file = "savedata/ff.krhnb.list.RData")
} else {
  load(file = "savedata/ff.krhnb.list.RData")
}

## Extract percentile intervals
quantilesAll <- apply(ff.krhnb.list$bootmfx, 2, quantile, c(0.025, 0.975))
quantilesOnes <- apply(ff.krhnb.list$bootonesmfx, 2, quantile, c(0.025, 0.975))
quantilesCount <- apply(ff.krhnb.list$bootcountmfx, 2, quantile, c(0.025, 0.975))

## Create matrix to print table of effects
avgEffectMatList <- cbind(c(rbind(plyr::revalue(colnames(ff.krhnb.list$effectmat), labels),
                              rep("", times = length(labels)))),
                      c(rbind(round(ff.krhnb.list$meaneffect, 2),
                              paste0("[", round(quantilesAll[1,], 2), ", ", round(quantilesAll[2,], 2), "]"))),
                      c(rbind(round(apply(ff.krhnb.list$onesmat, 2, mean), 2),
                              paste0("[", round(quantilesOnes[1,], 2), ", ", round(quantilesOnes[2,], 2), "]"))),
                      c(rbind(round(apply(ff.krhnb.list$countmat, 2, mean), 2),
                              paste0("[", round(quantilesCount[1,], 2), ", ", round(quantilesCount[2,], 2), "]"))))

## Table A2
avgEffectTabList <- xtable(avgEffectMatList)
print(avgEffectTabList, file = "tex/tabs_figs/tabAvgEffectList.tex", include.rownames=FALSE)
