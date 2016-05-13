#----------
# Authors: George Derpanopoulos and Luke Sonnet
#
# Project: Foreign Fighter Supply using KRHNB
#
# Purpose: Cleaning and merging data for foreign fighters
#----------

## Preamble
rm(list=ls())
library(foreign)
library(dplyr)
library(countrycode)
library(Amelia)

setwd('C:/Dropbox/Research/foreign_fighters/')

#----------
# World Bank Data
#----------

wb <- read.csv("data/wb.csv", as.is=T)

names(wb) <- c("year", "yearCode", "ctryName", "ctryCode", "internetSubPer100",
               "internetUserPer100", "internetServerPerMil", "area",
               "popDensity", "telephonePer100", "cellPer100",
               "lifeExp", "imr", "gdppcPPP", "pop", "popUrbPct",
               "unempYthMaleILO", "unempYthTotILO",
               "unempTotILO", "lfpMale",
               "lfpTot", "secondaryEnroll", "healthExpGDP", "milExpGDP",
               "secondExpPerStudent", "homicidesPer100k", "battleRelDeaths",
               "milPersonnelOfLabor", "migrantStockPctPop", "refugeeAsylum",
               "refugeeOrigin")

wb <- select(wb, -yearCode)

# Take latest non-missing values from 2010-2012
wbCross <- wb %>%
             group_by(ctryName, ctryCode) %>%
             arrange(year) %>%
             summarise_each(funs( .[max(which(!is.na(.)))] ))
as.data.frame(head(wbCross))

out <- wbCross

#----------
# Foreign Fighters
#----------

ff <- read.csv("data/icsr.csv")
ff <- ff %>%
        filter(icsr.date == 20150126) %>% # just use recent data
        select(countrycode, icsr.low, icsr.high, icsr.mean) %>%
        rename(ctryCode = countrycode,
               icsrLow = icsr.low,
               icsrHigh = icsr.high,
               icsrMean = icsr.mean) %>%
        mutate(icsrMean = ifelse(ctryCode %in% c("WBG", "ISR"), 60, icsrMean))
        # for now we split the Israel/Palestinian territories
        # see appendix where we do not do this

out <- merge(out, ff, all.x=T)

out[is.na(out$icsrMean), c("icsrLow", "icsrHigh", "icsrMean")] <- 0

## Drop Iraq and Syria
out <- filter(out, ctryCode != "IRQ", ctryCode != "SYR")

hist(out$icsrMean, breaks=100)

#----------
# Distance
#----------

cepii <- read.dta("data/dist_cepii.dta")
names(cepii)
cepii <- cepii %>%
           filter(iso_d %in% c("SYR", "IRQ")) %>%
           select(iso_o, contig, contains("dist")) %>%
           rename(ctryCode = iso_o) %>%
           mutate(
             ctryCode = ifelse(ctryCode=="PAL", "WBG", ctryCode),
             ctryCode = ifelse(ctryCode=="YUG", "SRB", ctryCode)
           ) %>%
           group_by(ctryCode) %>%
           summarise(
             contig = max(contig), # we do this to get any contiguity with iraq
             dist = min(dist)      # or syria and the smallest distance
           ) %>%
           rbind(., filter(., ctryCode=="ALB"))

# giving kosovo albania's distance values
cepii[cepii$ctryCode=="ALB",]$ctryCode <- c("KSV", "ALB")
filter(cepii, ctryCode == "ALB" | ctryCode == "KSV")

sort(out$ctryCode)
out[out$ctryCode %in% setdiff(out$ctryCode, cepii$ctryCode), c(1:3, 33)]

# This merge is going to drop a lot of micro states and newly founded states
# Then drop small nations

out <- merge(out, cepii, all=F)
sum(out$icsrMean)
out <- filter(out, pop > 500000)
sum(out$icsrMean)

#----------
# Muslim Population
#----------

mus <- read.csv('data/muslim.csv', as.is=T)

mus[mus$muslimpct=="<0.1",]$muslimpct <- 0.1
mus[mus$muslimpct=="<0.1 ",]$muslimpct <- 0.1
head(mus)
mus$muslimpct <- as.numeric(mus$muslimpct)
mus <- mutate(mus, iso2 = ifelse(is.na(iso2), "NA", iso2))

codes <- countrycode_data %>%
           select(wb, iso2c) %>%
           rename(ctryCode = wb, iso2 = iso2c)


out <- merge(out, codes)
out <- merge(out, select(mus, iso2, contains("muslim")), all.x = T)

out <- mutate(out, muslimpct = ifelse(ctryCode == "KSV", 95.6, muslimpct))
out <- filter(out, !(iso2 == "GB" & muslimpct==0.1)) # fix an error in the merge

#----------
# Freedom House
#----------

fh <- read.csv2('data/qog.txt', as.is = T) %>%
  filter(year == 2012,
         ccodealp != "SDN",
         ccodecow != 260,
         ccodecow != 678) %>%
  select(ccodealp, ccodecow, cname, starts_with('fh_')) %>%
  rename(ctryCode = ccodealp) %>%
  mutate(ctryCode = ifelse(ctryCode == "COD", "ZAR", ctryCode),
         ctryCode = ifelse(ctryCode == "ROU", "ROM", ctryCode))

# manually add kosovo, sudan, serbia, east timor, and west bank gaza
# data from website
newfh <- data.frame(ctryCode = c("KSV", "SDN", "SRB", "TMP", "WBG"),
                    fh_status = c(2, 3, 1, 2, 3),
                    fh_cl = c(4, 7, 2, 4, 6),
                    fh_pr = c(5, 7, 2, 3, 6))
totfh <- fh %>% 
  merge(newfh, all=T) %>%
  select(-starts_with('fh_fot'))

out <- merge(out, totfh, all.x=T)

#----------
# Shia Pct of Muslim
#----------

sh <- read.csv(file = 'data/shia.csv', as.is = T)

out <- merge(out, sh)
out <- mutate(out,
              shiapct = (muslimpct/100) * pctmusShia,
              sunnipct = muslimpct - shiapct)

#----------
# ARDA (Religious Freedom)
#----------

arda <- read.dta("data/arda.DTA")
ar <- arda %>%
  select(ardacode, COUNTR08, GRI_08, GFI_08, MSRI_08) %>%
  rename(
   ctryName = COUNTR08,
   govt.reg.relig = GRI_08,
   govt.fav.relig = GFI_08,
   soc.reg.relig = MSRI_08
  )

ardacode <- read.csv("data/arda_code.csv")
ardacode <- ardacode %>%
  filter(YEAR == "2010") %>%
  select(COWCODE, ARDACODE) %>%
  rename(
    cown = COWCODE,
    ardacode = ARDACODE
  ) %>%
  distinct()

codes2 <- countrycode_data %>%
  select(wb, cown) %>%
  rename(ctryCode = wb) %>%
  filter(!is.na(cown))

ar.m <- merge(ardacode, codes2)
ar2 <- merge(ar, ar.m, all.x = T)
ar2[is.na(ar2$ctryCode),]$ctryName
ar2 <- ar2 %>%
  mutate(
    ctryCode = ifelse(ctryName == "Vietnam", "VNM", ctryCode),
    ctryCode = ifelse(ctryName == "Israeli Occupied Territories (Palestine)", "WBG", ctryCode),
    ctryCode = ifelse(ctryName == "Hong Kong", "HKG", ctryCode),
    ctryCode = ifelse(ctryName == "Macau", "MAC", ctryCode),
    ctryCode = ifelse(ctryName == "Vietnam", "VNM", ctryCode)
  ) %>%
  select(-cown, -ctryName, -ardacode) %>%
  filter(!is.na(ctryCode))

out <- merge(out, ar2, all.x = T)

#----------
# Age Structure
#----------

age <- read.csv("data/age_clean.csv", stringsAsFactors = F)

age[age$country == "Brunei",]$country <- "Brunei Darussalam"
age[age$country == "Congo, Republic of the",]$country <- "Congo, Rep."
age[age$country == "Egypt",]$country <- "Egypt, Arab Rep."
age[age$country == "Hong Kong",]$country <- "Hong Kong SAR, China"
age[age$country == "Iran",]$country <- "Iran, Islamic Rep."
age[age$country == "Kyrgyzstan",]$country <- "Kyrgyz Republic"
age[age$country == "Korea, South",]$country <- "Korea, Rep."
age[age$country == "Laos",]$country <- "Lao PDR"
age[age$country == "Macau",]$country <- "Macao SAR, China"
age[age$country == "Macedonia",]$country <- "Macedonia, FYR"
age[age$country == "Burma",]$country <- "Myanmar"
age[age$country == "Korea, North",]$country <- "Korea, Dem. Rep."
age[age$country == "Russia",]$country <- "Russian Federation"
age[age$country == "Slovakia",]$country <- "Slovak Republic"
age[age$country == "Venezuela",]$country <- "Venezuela, RB"
age[age$country == "Yemen",]$country <- "Yemen, Rep."
age[age$country == "Congo, Democratic Republic of the",]$country <- "Congo, Dem. Rep."
age[age$country == "West Bank",]$country <- "West Bank and Gaza"
age[age$country == "Gaza Strip",]$country <- "West Bank and Gaza"
setdiff(out$ctryName, age$country)

# Average gaza and west bank numbers
age.m <- age %>%
  group_by(country) %>%
  summarise_each(funs(mean(., na.rm = T))) %>%
  rename(ctryName = country)

out <- merge(out, age.m)

#----------
# GTD, not used because of missingness
#----------

ter <- read.csv("data/terror_agg.csv", as.is = T)

setdiff(out$ctryName, ter$country_txt)
setdiff(ter$country_txt, out$ctryName)
ter[ter$country_txt == "Bahamas",]$country_txt <- "Bahamas, The"
ter[ter$country_txt == "Bosnia-Herzegovina",]$country_txt <- "Bosnia and Herzegovina"
ter[ter$country_txt == "Congo (Kinshasa)",]$country_txt <- "Congo, Dem. Rep."
ter[ter$country_txt == "Congo (Brazzaville)",]$country_txt <- "Congo, Rep."
ter[ter$country_txt == "Ivory Coast",]$country_txt <- "Cote d'Ivoire"
ter[ter$country_txt == "Egypt",]$country_txt <- "Egypt, Arab Rep."
ter[ter$country_txt == "Hong Kong",]$country_txt <- "Hong Kong SAR, China"
ter[ter$country_txt == "Iran",]$country_txt <- "Iran, Islamic Rep."
ter[ter$country_txt == "Kyrgyzstan",]$country_txt <- "Kyrgyz Republic"
ter[ter$country_txt == "South Korea",]$country_txt <- "Korea, Rep."
ter[ter$country_txt == "Laos",]$country_txt <- "Lao PDR"
ter[ter$country_txt == "Macedonia",]$country_txt <- "Macedonia, FYR"
ter[ter$country_txt == "Russia",]$country_txt <- "Russian Federation"
ter[ter$country_txt == "Great Britain",]$country_txt <- "United Kingdom"
ter[ter$country_txt == "Northern Ireland",]$country_txt <- "United Kingdom"
ter[ter$country_txt == "Western Sahara",]$country_txt <- "Morocco"
ter[ter$country_txt == "Corsica",]$country_txt <- "France"
ter[ter$country_txt == "Venezuela",]$country_txt <- "Venezuela, RB"
ter[ter$country_txt == "Yemen",]$country_txt <- "Yemen, Rep."
ter[ter$country_txt == "West Bank and Gaza Strip",]$country_txt <- "West Bank and Gaza"

ter.m <- ter %>%
  group_by(country_txt) %>%
  summarise_each(funs(sum(., na.rm = T))) %>%
  rename(ctryName = country_txt)

out <- merge(out, ter.m, all.x = T)
## Fill in NA with 0, bc most of these countries seem to be empty not because 
## data are missing but simply because there have been no events
out <- out %>%
  mutate(n.attacks = ifelse(is.na(n.attacks), 0, n.attacks),
         n.kill = ifelse(is.na(n.kill), 0, n.kill),
         n.attacks.post11 = ifelse(is.na(n.attacks.post11), 0, n.attacks.post11),
         n.kill.post11 = ifelse(is.na(n.kill.post11), 0, n.kill.post11))

#----------
# Fractionalization
#----------

frac <- read.csv("data/fractionalization.csv", stringsAsFactors=F)
head(frac)
frac <- frac[-1,c(1,4,5,6)]

frac[frac == 0] <- NA

frac <- frac %>%
  rename(ctryName = Country)

setdiff(out$ctryName, frac$ctryName)
setdiff(frac$ctryName, out$ctryName)
frac[frac$ctryName == "Bahamas",]$ctryName <- "Bahamas, The"
frac[frac$ctryName == "Brunei",]$ctryName <- "Brunei Darussalam"
frac[frac$ctryName == "Cape Verde",]$ctryName <- "Cabo Verde"
frac[frac$ctryName == "Congo, Dem. Rep. (Zaire)",]$ctryName <- "Congo, Dem. Rep."
frac[frac$ctryName == "Congo",]$ctryName <- "Congo, Rep."
frac[frac$ctryName == "Egypt",]$ctryName <- "Egypt, Arab Rep."
frac[frac$ctryName == "Equatorial Guinea ",]$ctryName <- "Equatorial Guinea"
frac[frac$ctryName == "Gambia, The ",]$ctryName <- "Gambia, The"
frac[frac$ctryName == "Hong Kong",]$ctryName <- "Hong Kong SAR, China"
frac[frac$ctryName == "Iran",]$ctryName <- "Iran, Islamic Rep."
frac[frac$ctryName == "Korea, North",]$ctryName <- "Korea, Dem. Rep."
frac[frac$ctryName == "Korea, South",]$ctryName <- "Korea, Rep."
frac[frac$ctryName == "Kyrgyzstan",]$ctryName <- "Kyrgyz Republic"
frac[frac$ctryName == "Lao People's Dem Rep",]$ctryName <- "Lao PDR"
frac[frac$ctryName == "Macau",]$ctryName <- "Macao SAR, China"
frac[frac$ctryName == "Macedonia (Former Yug. Rep)",]$ctryName <- "Macedonia, FYR"
frac[frac$ctryName == "Myanmar (Burma)",]$ctryName <- "Myanmar"
frac[frac$ctryName == "Serbia/Montenegro (Yugoslavia)",]$ctryName <- "Serbia"
frac[frac$ctryName == "East Timor",]$ctryName <- "Timor-Leste"
frac[frac$ctryName == "Venezuela",]$ctryName <- "Venezuela, RB"
frac[frac$ctryName == "Yemen",]$ctryName <- "Yemen, Rep."
frac[frac$ctryName == "West Bank",]$ctryName <- "West Bank and Gaza"
frac[frac$ctryName == "Gaza Strip",]$ctryName <- "West Bank and Gaza"

frac$religious.frac <- as.numeric(frac$Religion)
frac$ethnic.frac <- as.numeric(frac$Ethnic)
frac$language.frac <- as.numeric(frac$Language)

frac.m <- frac %>%
  group_by(ctryName) %>%
  select(ends_with("frac")) %>%
  summarise_each(funs(mean))

out <- merge(out, frac.m, all.x=T)

#----------
# Add refugees without IDPs
#----------

unhcr <- read.csv("data/unhcr_popstats_export_time_series_2016_05_07_003135.csv",
                  stringsAsFactors = F)
unhcr <- rename(unhcr, country = Country...territory.of.asylum.residence)

unhcr[unhcr$country == "Bolivia (Plurinational State of)",]$country <- "Bolivia"
unhcr[unhcr$country == "Côte d'Ivoire",]$country <- "Cote d'Ivoire"
unhcr[unhcr$country == "Central African Rep.",]$country <- "Central African Republic"
unhcr[unhcr$country == "China, Hong Kong SAR",]$country <- "Hong Kong SAR, China"
unhcr[unhcr$country == "Congo",]$country <- "Congo, Rep."
unhcr[unhcr$country == "Dem. Rep. of the Congo",]$country <- "Congo, Dem. Rep."
unhcr[unhcr$country == "Czech Rep.",]$country <- "Czech Republic"
unhcr[unhcr$country == "Dominican Rep.",]$country <- "Dominican Republic"
unhcr[unhcr$country == "Egypt",]$country <- "Egypt, Arab Rep."
unhcr[unhcr$country == "Gambia",]$country <- "Gambia, The"
unhcr[unhcr$country == "Iran (Islamic Rep. of)",]$country <- "Iran, Islamic Rep."
unhcr[unhcr$country == "Kyrgyzstan",]$country <- "Kyrgyz Republic"
unhcr[unhcr$country == "China, Macao SAR",]$country <- "Macao SAR, China"
unhcr[unhcr$country == "Rep. of Korea",]$country <- "Korea, Rep."
unhcr[unhcr$country == "Rep. of Moldova",]$country <- "Moldova"
unhcr[unhcr$country == "Serbia and Kosovo (S/RES/1244 (1999))",]$country <- "Serbia"
unhcr[unhcr$country == "Slovakia",]$country <- "Slovak Republic"
unhcr[unhcr$country == "South Sudan",]$country <- "Sudan" ## So that it sums to sudan's total
unhcr[unhcr$country == "The former Yugoslav Republic of Macedonia",]$country <- "Macedonia, FYR" 
unhcr[unhcr$country == "United Rep. of Tanzania",]$country <- "Tanzania" 
unhcr[unhcr$country == "United States of America",]$country <- "United States" 
unhcr[unhcr$country == "Venezuela (Bolivarian Republic of)",]$country <- "Venezuela, RB" 
unhcr[unhcr$country == "Yemen",]$country <- "Yemen, Rep." 

unhcrCtry <- unhcr %>%
  group_by(country) %>%
  summarise(refugeesUnhcr = sum(Value)) %>%
  select(country, refugeesUnhcr) %>%
  rename(ctryName = country)

out <- merge(out, unhcrCtry, all.x=T)

## Don't have inequality added bc there is so much missinginess

#----------
# Write Merged File
#----------

write.csv(out, file="data/foreignFighters.csv", row.names=F)

#----------
# Impute Some Data
#----------

sort(apply(out, 2, function(x) sum(is.na(x))))

df <- select(out,
             ctryName, ctryCode, area, popDensity, lifeExp, pop, icsrMean, contig, dist,
             internetUserPer100, migrantStockPctPop,
             homicidesPer100k, refugeeOrigin, refugeesUnhcr, gdppcPPP, lfpMale,
             unempYthTotILO, fh_cl, fh_pr, sunnipct, shiapct, govt.reg.relig,
             govt.fav.relig, soc.reg.relig, pctpop15.24, n.kill,
             religious.frac)

sort(apply(df, 2, function(x) sum(is.na(x))))

cbind(df$ctryName, df$pop, df$icsrMean, apply(df, 1, function(x) sum(is.na(x))))

## Drop Puerto Rico, Macao 7 missing values, non-states
df <- df[!(df$ctryName %in% c("Puerto Rico", "Macao SAR, China")), ]

df2 <- df

## Log transformation will improve the efficiency of estimators as well as
## drastically improve imputation
df2$lnHomicidesPer100k <- log(df$homicidesPer100k)
df2$lnpopDensity <- log(df$popDensity)
df2$lnpop <- log(df$pop)
df2$lnarea <- log(df$area)
df2$lngdppcPPP <- log(df$gdppcPPP)
df2$lnrefugeesUnhcr <- log(df$refugeesUnhcr+1)
df2$lnrefugeeOrigin <- log(df$refugeeOrigin)

df2 <- df2[, !(names(df2) %in% c("homicidesPer100k", "popDensity", "pop",
                                 "area", "gdppcPPP", "refugeesUnhcr",
                                 "refugeeOrigin", "n.kill"))]

write.csv(df2, file="data/foreignFightersPrep.csv", row.names=F)

## BOunds for some variables
bds.names <- c("internetUserPer100", "migrantStockPctPop", 
               "lfpMale", "lnrefugeeOrigin",
               "unempYthTotILO", "fh_cl", "fh_pr", "govt.reg.relig",
               "govt.fav.relig", "soc.reg.relig", "lngdppcPPP", "lnrefugeesUnhcr",
               "religious.frac")
bds.nums <- which( colnames(df2) %in% bds.names )
summary(df2[, bds.nums])
bds.mins <- c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 6, 0, 0)
bds.max <- c(100, 100, 100, 100, 7, 7, 10, 10, 10, 1, 12, 16, 15)

bds <- matrix(ncol = 3,
              c(bds.nums, bds.mins, bds.max))
## Impute 1000 times and take the average
set.seed(20160216)
a <- amelia(df2, idvars = c("ctryName", "ctryCode"),
            m = 1000, bounds = bds, emburn = c(10, 1000), empri = 1, p2s = 0)

## Sum together results across imputed dataset, long way to get mean values
meanA <- Reduce("+", lapply(a$imputations, function(x) x[, !(names(x) %in% c("ctryName", "ctryCode"))])) / length(a$imputations) 
meanDF <- cbind.data.frame(ctryCode = df2$ctryCode, ctryName = df2$ctryName, meanA)

## Adding Europe dummy
meanDF$europe <- ifelse(meanDF$ctryName %in% c("Albania", "Armenia", "Austria", "Begium",
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

write.csv(meanDF, file="data/foreignFightersImputed.csv", row.names=F)
