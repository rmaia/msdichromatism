#--------------------------------------------------------------------------------------------------
## Title: Behavioural ecology of a Sri Lankan lizard, Ceratophora tennenttii
## Authors: Martin Whiting, Ruchira Somaweera, Daniel Noble	
## Date started: 24/09/2013
#-------------------------------------------------------------------------------------------------
source("./Analysis/Functions/modcheck.R")

##--------------------------------------ANALYSIS 1 ----------------------------------------##
## Purpose: To understand whether males and females are sexually dimorphic

## Relationship between horn length and SVL between the sexes

scatterplotMatrix(~log(HH) + log(SVL) + log(HW)+ log(HL)+ log(Horn.L) +log(Crest.H)+ log(Dewlap.H)|Sex, data = datA1)
scatterplotMatrix(~log(SVL) + log(Horn.L) +log(Crest.H)+ log(Dewlap.H)|Sex, diagonal = "histogram", data = datA1)

#-------------------------------------------SVL-------------------------------------------#
summary(modelSVL1.1 <- lm(SVL~ Sex, data = datA1))

drop1(modelSVL1.1, test = "F")

#----------------------------------------- Mass ------------------------------------------#
summary(modelMass1.1 <- lm(Mass~ Sex, data = datA1))

drop1(modelMass1.1, test = "F")

#--------------------------------------Rostrum lenghth----------------------------------------#

summary(modelA1.1 <- lm(log(Horn.L) ~ log(SVL) + Sex + log(SVL):Sex, data = datA1))

# Standardized coefs, x-meanx/2*sdx
summary(modelA1.1_std <- standardize(modelA1.1_1 <- lm(LHorn.L ~ LSVL + Sex + LSVL:Sex, data = datA1)))  
## Model Checking

modcheck(modelA1.1)
leveneTest(residuals(modelA1.1), datA1$Sex)

## Test the standardization procedure to understand what is happening
#Rescale coeffiecits prior to running model
y <- log(datA1$Horn.L)
datA1$svl_scale <- (log(datA1$SVL)-mean(log(datA1$SVL)))/(2*sd(log(datA1$SVL)))
datA1$sex_scale <- ifelse(datA1$Sex == "f", -0.5, 0.5)

# Run modela nd should match with modelA1.1_11, which is ALMOST does. Maybe related to rounding differences between the two methods

mod <- lm(y ~ svl_scale+sex_scale+svl_scale:sex_scale, data = datA1)


## Compare results to an rlm regression. If large effects on coefficients then we have problems and we can use robust regression methods

summary(modelA1.1rlm <- rlm(LHorn.L ~ LSVL + Sex + LSVL:Sex, data = datA1)) 

boot.huber_HorL <- function(data, indices, maxit = 100){
	data <- data[indices,]
	mod  <- rlm(LHorn.L ~ svl_scale + sex_scale + svl_scale:sex_scale, data = data, maxit = maxit)
	coefficients(mod)
}

LHorL_boot <- boot(datA1, boot.huber_HorL, 1000)
summary(LHorL_boot)

#confidence intervals
for(i in 1:4){
print(boot.ci(LHorL_boot, index = i, type = "norm"))
}

# No really strong changes, suggests that point 52 is not majorly influential on parameter estimates. This is verified with model diagnostics above. 52 is within cooks lines.

## Predict values and plot lines and error bars on plot above

newdata_m      <- data.frame(SVL = datA1[datA1$Sex == "m",]$SVL, Sex = "m")

predict_male   <- predict(modelA1.1, newdata = newdata_m, type = "response", se.fit = TRUE)

newdata_f      <- data.frame(SVL = datA1[datA1$Sex == "f",]$SVL, Sex = "f")

predict_female <- predict(modelA1.1, newdata = newdata_f, type = "response", se.fit = TRUE)


#---------------------------------------Crest.H-----------------------------------#

summary(modelA1.3_1 <- lm(log(Crest.H) ~ log(SVL) + Sex + log(SVL):Sex, data = datA1))

summary(modelA1.3_std <- standardize(modelA1.3_1 <- lm(LCrest.H ~ LSVL + Sex + LSVL:Sex, data = datA1)))

## Model Checking
leveneTest(residuals(modelA1.3_1), datA1$Sex)          ## No evidence for het of variance
check <- modcheck(modelA1.3_1)

## Robust regression because of outliers. Set psi to huber weighting because these are not severe outliers. Huber weights residuals less serverly than bisquare weighting. If results don't differ much from lm method than these outliers are not likely to be very influential on the coefficients. 

summary(modelA1.3_1rlm <- rlm(log(Crest.H) ~ scale(log(SVL)) + Sex + scale(log(SVL)):Sex, data = datA1, psi = psi.huber, scale.est = "proposal 2"))

summary(modelA1.3_1rlm_std <- rlm(LCrest.H ~ svl_scale + sex_scale + svl_scale:sex_scale, data = datA1, psi = psi.huber, scale.est = "proposal 2"))

# Generate confidence intervals

boot.huber_CH <- function(data, indices, maxit = 100){
	data <- data[indices,]
	mod  <- rlm(LCrest.H ~ svl_scale + sex_scale + svl_scale:sex_scale, data = data, maxit = maxit)
	coefficients(mod)
}

LCrest_boot <- boot(datA1, boot.huber_CH, 1000)
summary(LCrest_boot)

par(mfrow = c(4,2))
for(i in 1:4){
plot(LCrest_boot, index = i)
}

#confidence intervals
for(i in 1:4){
print(boot.ci(LCrest_boot, index = i, type = "norm"))
}

## Predict values and plot lines and error bars on plot above

newdata_m_1      <- data.frame(SVL = datA1[datA1$Sex == "m",]$SVL, Sex = "m")

predict_male_1   <- predict(modelA1.3_1rlm, newdata = newdata_m, type = "response", se.fit = TRUE)

newdata_f_1      <- data.frame(SVL = datA1[datA1$Sex == "f",]$SVL, Sex = "f")

predict_female_1 <- predict(modelA1.3_1rlm, newdata = newdata_f, type = "response", se.fit = TRUE)

#-----------------------------------------Dewlap.H---------------------------------------#

summary(modelA1.4_1 <- lm(log(Dewlap.H) ~ log(SVL) + Sex + log(SVL):Sex, data = datA1))

drop1(modelA1.4_1, test = "F")

summary(modelA1.4_2 <- lm(LDewlap.H ~ LSVL + Sex, data = datA1))

summary(modelA1.4_2_std <- standardize(lm(LDewlap.H ~ LSVL + Sex, data = datA1)))

## Model Checking
modcheck(modelA1.4_2_std)                          
leveneTest(residuals(modelA1.4_2_std), datA1$Sex)          

## Robust regression because of outliers. Set psi to huber weighting because these are not severe outliers. Huber weights residuals less serverly than bisquare weighting. If results don't differ much from lm method than these outliers are not likely to be very influential on the coefficients. 

summary(modelA1.4_2rlm <- rlm(log(Dewlap.H) ~ log(SVL) + Sex, data = datA1, psi = psi.huber, scale.est = "proposal 2"))

summary(modelA1.4_2rlm_std <- rlm(LDewlap.H ~ svl_scale + sex_scale, data = datA1, psi = psi.huber, scale.est = "proposal 2"))

# Generate confidence intervals

boot.huber_DH <- function(data, indices, maxit = 100){
	data <- data[indices,]
	mod  <- rlm(LDewlap.H ~ svl_scale + sex_scale, data = data, maxit = maxit)
	coefficients(mod)
}

LDewlap_boot <- boot(datA1, boot.huber_DH, 1000)
summary(LDewlap_boot)

par(mfrow= c(4,2))
for(i in 1:3){
plot(LDewlap_boot, index = i)
}

#confidence intervals
for(i in 1:3){
print(paste("INDEX", i, sep = "_"))
print(boot.ci(LDewlap_boot, index = i, type = "norm"))
}

## Predict values and plot lines and error bars on plot above

newdata_m_2      <- data.frame(SVL = datA1[datA1$Sex == "m",]$SVL, Sex = "m")

predict_male_2   <- predict(modelA1.4_2rlm, newdata = newdata_m, type = "response", se.fit = TRUE)

newdata_f_2      <- data.frame(SVL = datA1[datA1$Sex == "f",]$SVL, Sex = "f")

predict_female_2 <- predict(modelA1.4_2rlm, newdata = newdata_f, type = "response", se.fit = TRUE)


#-----------------------------------------Head width---------------------------------------##

summary(modelA1.5 <- lm(log(HW)~ log(SVL) + Sex + log(SVL):Sex, data = datA1))

summary(modelA1.5_std <- standardize(lm(LHW~ LSVL + Sex + LSVL:Sex, data = datA1)))

## Model Checking
modcheck(modelA1.5)                                  ## Weak evidence for non-linearity
leveneTest(residuals(modelA1.5), datA1$Sex)          ## No evidence for het of variance

## Model predictions

newdata_m.3      <- data.frame(SVL = datA1[datA1$Sex == "m",]$SVL, Sex = "m")

predict_male.3   <- predict(modelA1.5, newdata = newdata_m.3, type = "response", se.fit = TRUE)

newdata_f.3      <- data.frame(SVL = datA1[datA1$Sex == "f",]$SVL, Sex = "f")

predict_female.3 <- predict(modelA1.5, newdata = newdata_f.3, type = "response", se.fit = TRUE)


#-------------------------------------------Head length------------------------------------##

summary(modelA1.6 <- lm(log(HL)~ log(SVL) + Sex + log(SVL):Sex, data = datA1))

drop1(modelA1.6, test = "F")

summary(modelA1.6 <- lm(log(HL)~ log(SVL) + Sex, data = datA1))

summary(modelA1.6_std <- standardize(lm(LHL~ LSVL + Sex, data = datA1)))

## Model Checking
modcheck(modelA1.6)                                  ## Weak evidence for non-linearity
leveneTest(residuals(modelA1.6), datA1$Sex)          ## No evidence for het of variance

## Model predictions

newdata_m.4      <- data.frame(SVL = datA1[datA1$Sex == "m",]$SVL, Sex = "m")

predict_male.4   <- predict(modelA1.6, newdata = newdata_m.4, type = "response", se.fit = TRUE)

newdata_f.4      <- data.frame(SVL = datA1[datA1$Sex == "f",]$SVL, Sex = "f")

predict_female.4 <- predict(modelA1.6, newdata = newdata_f.4, type = "response", se.fit = TRUE)


##---------------------------------------Head Height----------------------------------------#

summary(modelA1.7 <- lm(log(HH)~ log(SVL) + Sex + log(SVL):Sex, data = datA1))

drop1(modelA1.7, test = "F")

summary(modelA1.7 <- lm(log(HH)~ log(SVL) + Sex, data = datA1))

summary(modelA1.7_std <- standardize(lm(LHH~ LSVL + Sex, data = datA1)))

## Model Checking

modcheck(modelA1.7)                           ## Weak evidence for non-linearity
leveneTest(residuals(modelA1.7), datA1$Sex)          ## No evidence for het of variance

## Model predictions

newdata.male7   <- data.frame(SVL = datA1$SVL, Sex = "m")

predict.male7   <- predict(modelA1.7, newdata = newdata.male7, type = "response", se.fit = TRUE)

newdata.female7   <- data.frame(SVL = datA1$SVL, Sex = "f")

predict.female7   <- predict(modelA1.7, newdata = newdata.female7, type = "response", se.fit = TRUE)


#--------------Multivariate analysis of covariance on combined head dimensions--------------##

# Head dimensions. Should not change much. 
summary(mv1 <- lm(cbind(log(HH), log(HW), log(HL))~log(SVL) + Sex + log(SVL):Sex, data = datA1))
Manova(mv1, type = "III")

# All morphological response variables. Just to make sure no type I errors. 
summary(mv2 <- lm(cbind(log(HH), log(HW), log(HL), log(Horn.L), log(Crest.H), log(Dewlap.H))~log(SVL) + log(SVL):Sex+ Sex, data = datA1))

Manova(mv2, type = "III")

##-------------------------------------------ANALYSIS 2-------------------------------------#
## Models

## Bite force, head width, body size and sex.

cor(datA3[c("HW", "HH", "HL", "SVL")]) ## High correlations and may give some trouble modelling. Simplify. HL explains most variance, but barely.

summary(modelA3.1 <- lm(log(Max.BF)~ as.factor(Sex) + SVL + HW + as.factor(Sex):HW + as.factor(Sex):SVL, data = datA3))
summary(modelA3.2 <- lm(log(Max.BF)~ as.factor(Sex) + SVL + HL + as.factor(Sex):HL + as.factor(Sex):SVL, data = datA3))
summary(modelA3.3 <- lm(log(Max.BF)~ as.factor(Sex) + SVL + HH + as.factor(Sex):HH + as.factor(Sex):SVL, data = datA3))

## Use head length - different between the sexes and explains most variance. No interaction with SVL so we can exlcude this.

summary(modelA3.2_std <- standardize(modelA3.2))

drop1(modelA3.2_std, test = "F")

## Model checking
modcheck(modelA3.2)
vif(modelA3.2)
leveneTest(residuals(modelA3.2), datA3$Sex)

## Identifies curvature in SVL and HL - log trans helps a bit but still there. Given the sample size limitations I am hesitatnt to include quadratics in the model, which would improve fit. Comparing the fite for both models to the observed values, shows that the model is still doing a reasonably good job predicting the observed data. So I will still whith what I have and avoid overparameterization.

#Welch's t-test because unequal variance between sexes in Mx.BF - see leveneTest
t.test(log(datA3[datA3$Sex == "m",]$Max.BF),log(datA3[datA3$Sex == "f",]$Max.BF))
#Welch's t-test because unequal variance between sexes condition
datA3$condition <- residuals(lm(log(Mass)~log(SVL), data = datA3))
t.test(datA3[datA3$Sex == "m",]$condition, datA3[datA3$Sex == "f",]$condition)

#-------------------- Spectral and morphology analysis ----------------------#

specdat  <- lapply(ind_vs_bkg_liz, function(x) merge(x, datA3, "ID"))

for(i in 1:4){
names(specdat[[i]])[6] <-  "Sex"
}

condition <- lapply(specdat, function(x) residuals(lm(log(Mass)~log(SVL), data = x)))

for(i in 1:4){
specdat[[i]]$condition <-  condition[[i]]
}

# Mouth and Roof
par(mfrow = c(2,2))
lapply(specdat, function(x) hist(log(x$dS)))

# Are sexes different in dS or dL? Note results may be slightly different here because the sample size of dS changed because of missing data in bite force. Plus this is using lizard visual system and NOT bird, which is in results. 

lapply(specdat, function(x) summary(glm(log(dS) ~ Sex, data = x)))
lapply(specdat, function(x) hist(residuals(glm(log(dS) ~ Sex, data = x))))

lapply(specdat, function(x) wilcox.test(log(x$dS)~x$Sex))
lapply(specdat, function(x) t.test(log(x$dS)~x$Sex))

# Is Horn length correlated with chromatic contrast?
par(mfrow = c(2,2))
lapply(specdat, function(x) plot(log(x$dS) ~ Horn.L, data = x))
lapply(specdat, function(x) cor.test(log(x$dS), x$Horn.L, method = "spearman"))
lapply(specdat, function(x) cor.test(log(x[x$Sex == "M",]$dS), x[x$Sex == "M",]$Horn.L, method = "spearman"))
lapply(specdat, function(x) cor.test(log(x[x$Sex == "F",]$dS), x[x$Sex == "F",]$Horn.L, method = "spearman"))

throatCor <- specdat[["throat"]]

plot(Horn.L ~ log(dS), data = throatCor[throatCor$Sex == "M",])
plot(Horn.L ~ log(dS), data = throatCor[throatCor$Sex == "F",], add = TRUE)

throatdat <- specdat$throat
labialdat <- specdat$labial

## Does body size correlate with labial chromatic contrast?

summary(mod <- glm(SVL ~ dS + Sex, dat = labialdat, family = gaussian))
hist(residuals(mod))

#Sexes Pooled
cor.test(labialdat[labialdat$Sex == "M",]$SVL, labialdat[labialdat$Sex == "M",]$dS, method = "spearman")
cor.test(labialdat[labialdat$Sex == "F",]$SVL, labialdat[labialdat$Sex == "F",]$dS, method = "spearman")
cor.test(labialdat$SVL, labialdat$dS, method = "spearman")
cor.test(labialdat$SVL, labialdat$dS)

#---------------------------------------- Throat - Bite Force -------------------------------##
summary(modthroat_std1 <- standardize(lm(log(Max.BF)~ Sex + TempBF + dS + Horn.L  + SVL + HL  + Sex:HL + Sex:SVL + Sex:dS + Sex:Horn.L, data = throatdat)))

drop1(modthroat_std1, test = "F")

summary(modthroat_std2 <- standardize(lm(log(Max.BF)~ Sex + TempBF + dS + Horn.L  + SVL + HL  + Sex:HL + Sex:SVL + Sex:dS, data = throatdat)))

drop1(modthroat_std2, test = "F")


summary(modthroat_std3 <- standardize(lm(log(Max.BF)~ Sex + TempBF + dS + Horn.L  + SVL  + HL  + Sex:HL + Sex:SVL, data = throatdat)))

drop1(modthroat_std3, test = "F")

## Model validation
modcheck(modthroat_std3)


## Achromatic contrast
summary(modthroat_std1 <- standardize(lm(log(Max.BF)~ Sex + TempBF + dL + Horn.L  + SVL + HL  + Sex:HL + Sex:SVL + Sex:dL + Sex:Horn.L, data = throatdat)))

drop1(modthroat_std1, test = "F")

summary(modthroat_std2 <- standardize(lm(log(Max.BF)~ Sex + TempBF + dL + Horn.L  + SVL + HL  + Sex:HL + Sex:SVL + Sex:dL, data = throatdat)))

drop1(modthroat_std2, test = "F")

summary(modthroat_std3 <- standardize(lm(log(Max.BF)~ Sex + TempBF + dL + Horn.L  + SVL + HL  + Sex:HL + Sex:SVL, data = throatdat)))

drop1(modthroat_std3, test = "F")
#---------------------------------------- Labials - Bite Force -------------------------------##

summary(modlabial_std1 <- standardize(lm(log(Max.BF)~ Sex + TempBF+ dS +  Horn.L  + SVL + HL + Sex:HL + Sex:SVL + Sex:dS + Sex:Horn.L, data = labialdat)))

drop1(modlabial_std1, test = "F")

summary(modlabial_std2 <- standardize(lm(log(Max.BF)~ Sex + TempBF + dS +  Horn.L  + SVL + HL + Sex:HL + Sex:SVL + Sex:dS, data = labialdat)))

drop1(modlabial_std2, test = "F")

summary(modlabial_std3 <- standardize(lm(log(Max.BF)~ Sex + TempBF + dS +  Horn.L  + SVL + HL + Sex:HL + Sex:SVL, data = labialdat)))

drop1(modlabial_std3, test = "F")

## Model validation
modcheck(modlabial_std3)

#Model predictions
summary(modlabial <- lm(log(Max.BF)~ Sex + dS + TempBF + Horn.L  + SVL + HL + Sex:HL + Sex:SVL, data = labialdat))

male_dat          <- data.frame(Sex = "M", dS = seq(min(labialdat$dS), max(labialdat$dS), length.out = 30), TempBF = mean(labialdat$TempBF), SVL = mean(labialdat[labialdat$Sex == "M", "SVL"]), HL = mean(labialdat[labialdat$Sex == "M", "HL"]), Horn.L = mean(labialdat[labialdat$Sex == "M", "Horn.L"]))

female_dat        <- data.frame(Sex = "F", dS = seq(min(labialdat$dS), max(labialdat$dS), length.out = 30), TempBF = mean(labialdat$TempBF), SVL = mean(labialdat[labialdat$Sex == "F", "SVL"]), HL = mean(labialdat[labialdat$Sex == "F", "HL"]), Horn.L = mean(labialdat[labialdat$Sex == "F", "Horn.L"]))

pred_male_labBF   <- predict(modlabial, newdata = male_dat, type = "response", se.fit = TRUE)

pred_female_labBF <- predict(modlabial, newdata = female_dat, type = "response", se.fit = TRUE)

labialdat$fitted <- fitted(modlabial)

# Achromatic contrast
summary(modlabial_stdL1 <- standardize(lm(log(Max.BF)~ Sex + TempBF+ dL +  Horn.L  + SVL + HL + Sex:HL + Sex:SVL + Sex:dL + Sex:Horn.L, data = labialdat)))

drop1(modlabial_stdL1, test = "F")

summary(modlabial_stdL2 <- standardize(lm(log(Max.BF)~ Sex + TempBF + dL +  Horn.L  + SVL + HL + Sex:HL + Sex:SVL + Sex:dL, data = labialdat)))

drop1(modlabial_stdL2, test = "F")

summary(modlabial_stdL3 <- standardize(lm(log(Max.BF)~ Sex + TempBF + dL +  Horn.L  + SVL + HL + Sex:HL + Sex:SVL, data = labialdat)))

drop1(modlabial_stdL3, test = "F")

## Model validation
modcheck(modlabial_stdL3)
#-----------------------------------Condition - Throat---------------------------------------#

#Add condition
throatdat$condition <- resid(lm(log(Mass)~log(SVL), data = throatdat))

summary(modelA2.cond <- standardize(lm(condition ~ as.factor(Sex) + dS + Horn.L + as.factor(Sex):Horn.L + as.factor(Sex):dS, data = throatdat)))

drop1(modelA2.cond, test = "F")

# Final Model and present all effects
summary(modelA2.cond <- lm(condition ~ as.factor(Sex) + dS + Horn.L + as.factor(Sex):dS, data = throatdat))

drop1(modelA2.cond, test = "F")

summary(modelA2.cond <- lm(condition ~ as.factor(Sex) + dS + Horn.L, data = throatdat))

drop1(modelA2.cond, test = "F")

#Standardize effect sizes.
summary(modelA2.cond_std <- standardize(modelA2.cond))

# Model checking 
modcheck(modelA2.cond)
leveneTest(residuals(modelA2.cond), throatdat$Sex)

plot(condition~dS, data = throatdat[throatdat$Sex == "F",], cex = 1.5)
points(condition~dS, data = throatdat[throatdat$Sex == "M",], pch = 16, cex = 1.5)

## Robust regression as there is a major outlier

throatdat$dS_z     <- throatdat$dS-mean(throatdat$dS)/(2*sd(throatdat$dS))
throatdat$Horn.L_z <- throatdat$Horn.L-mean(throatdat$Horn.L)/(2*sd(throatdat$Horn.L))
throatdat$Sex.z    <- as.factor(ifelse(throatdat$Sex == "F", -0.5, 0.5))

summary(modelA2.cond_Rob <- rlm(condition ~ Sex.z + dS_z + Horn.L_z + Sex.z:dS_z, data = throatdat))

summary(modelA2.cond_TEST <- lm(condition ~ Sex.z + dS_z + Horn.L_z + Sex.z:dS_z, data = throatdat))

# Generate confidence intervals

boot.huber_cond <- function(data, indices, maxit = 100){
	data <- data[indices,]
	mod  <- rlm(condition ~ Sex.z + dS_z + Horn.L_z + Sex.z:dS_z, data = data, maxit = maxit)
	coefficients(mod)
}

cond_boot <- boot(throatdat, boot.huber_cond, 1000)
summary(cond_boot)

par(mfrow= c(4,2))
for(i in 1:3){
plot(cond_boot, index = i)
}

#confidence intervals
for(i in 1:5){
print(paste("INDEX", i, sep = "_"))
print(boot.ci(cond_boot, index = i, type = "norm"))
}

summary(boot.ci(cond_boot))

#Achromatic contrast
throatdat$dL_z     <- throatdat$dL-mean(throatdat$dL)/(2*sd(throatdat$dL))

summary(modelA2.condT1 <- standardize(lm(condition ~ as.factor(Sex) + dL + Horn.L + as.factor(Sex):Horn.L + as.factor(Sex):dL, data = throatdat)))

	drop1(modelA2.condT1, test = "F")

# Final Model and present all effects
summary(modelA2.condT2 <- lm(condition ~ as.factor(Sex) + dL + Horn.L + as.factor(Sex):dL, data = throatdat))

	drop1(modelA2.condT2, test = "F")

summary(modelA2.condT3 <- lm(condition ~ as.factor(Sex) + dL + Horn.L, data = throatdat))

	drop1(modelA2.condT3, test = "F")

#Standardize final model
summary(modelA2.condT3 <- standardize(lm(condition ~ as.factor(Sex) + dL + Horn.L, data = throatdat)))
#------------------------------------ Labials - Condition ---------------------------------------#


summary(modelA2.cond_L <- lm(condition ~ as.factor(Sex) + dS + Horn.L + as.factor(Sex):Horn.L + as.factor(Sex):dS, data = labialdat))

drop1(modelA2.cond_L, test = "F")

summary(modelA2.cond_L <- lm(condition ~ as.factor(Sex) + dS + Horn.L + as.factor(Sex):dS, data = labialdat))

drop1(modelA2.cond_L, test = "F")

summary(modelA2.cond_L <- lm(condition ~ as.factor(Sex) + dS + Horn.L, data = labialdat))

t.test(condition~Sex, data = labialdat)

# Standardize model coefficients for final model
summary(modelA1.cond_std_L <- standardize(modelA2.cond_L))


## Model Checking
modcheck(modelA2.cond_L)				 
leveneTest(residuals(modelA2.cond_L), labialdat$Sex)     

## Achromatic contrast

summary(modelA2.cond_L1 <- lm(condition ~ as.factor(Sex) + dL + Horn.L + as.factor(Sex):Horn.L + as.factor(Sex):dL, data = labialdat))

drop1(modelA2.cond_L1, test = "F")

summary(modelA2.cond_L2 <- lm(condition ~ as.factor(Sex) + dL + Horn.L + as.factor(Sex):dL, data = labialdat))

drop1(modelA2.cond_L2, test = "F")

summary(modelA2.cond_L3 <- lm(condition ~ as.factor(Sex) + dL + Horn.L, data = labialdat))

t.test(condition~Sex, data = labialdat)

# Standardize model coefficients for final model
summary(modelA1.cond_std_L3 <- standardize(modelA2.cond_L3))

##----------------------------------------- Generate Tables ----------------------------------------##

source("./Analysis/4_tables.R")

#-------------------------------------------Generate Figures ---------------------------------------##

source("./Analysis/3_Figures.R")

