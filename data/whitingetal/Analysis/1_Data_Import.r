## Packages
packages <- c("arm", "plotrix", "plyr", "car", "AICcmodavg", "grid", "png", "pavo", "boot")

install.packages(packages)

sapply(packages, function(x) library(x, character = TRUE))


### GENERATE ALL DATA FOR ANALYSIS

##--------------------------------------- Upload the data file ------------------------------------##

data_mor <- read.csv("./Data/Ceratophora tennenttii Jan2012.csv", header = TRUE, stringsAsFactors = FALSE)

## Restructure data - ANALYSIS 1
datA1 <- data_mor[,c("ID", "Sex", "HL", "HW", "HH", "Horn.L", "Horn.W", "Crest.H", "Dewlap.H", "SVL", "TotL", "Tail", "Mass")]

datA1$Sex       <- as.factor(datA1$Sex)
datA1$condition <- residuals(lm(log(Mass)~log(SVL), data = datA1))
datA1$LSVL      <- log(datA1$SVL)
datA1$LHorn.L   <- log(datA1$Horn.L)
datA1$LCrest.H  <- log(datA1$Crest.H)
datA1$LDewlap.H <- log(datA1$Dewlap.H)
datA1$LHW       <- log(datA1$HW)
datA1$LHH       <- log(datA1$HH)
datA1$LHL       <- log(datA1$HL)

## Restructure data - ANALYSIS 2
datA3 <- na.omit(data_mor[,c("ID", "Sex", "Max.BF", "TempBF", "HL", "HW", "HH", "Horn.L", "Horn.W", "Crest.H", "Dewlap.H", "SVL", "TotL", "Tail", "Mass")])

## Restructure data - Spec and morphology analysis
#See 2_Spectral_analysis.R file for analysis. To do this use ind_vs_bkg
## So note that we are missing 3 lizards without bite force and then another 3 lizards where we don't have spec data, hence why we have 57 lizards total 
source("./Analysis/2_Spectral_analysis.R")

