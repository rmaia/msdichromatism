#----------------------------------------Generate Tables for Manuscript---------------------------------------#

##--------------------------- Table 1 ------------------------------##

## Means and standard deviations for males and females.



data_sum <- ddply(data_mor, .(Sex), summarise, 
	Head.Length = mean(HL), sdHL = sd(HL), nHL = length(HL), 
	Head.width = mean(HW), sdHW = sd(HW), nHW = length(HW), 
	Head.Height = mean(HH), sdHH = sd(HH), nHH = length(HH), 
	Horn.Length = mean(Horn.L), sdHorn.L = sd(Horn.L), nHorn.L = length(Horn.L), 
	mDewlapH = mean(Dewlap.H), sdDewlapH = sd(Dewlap.H), nDewlapH = length(Dewlap.H), 
	mSVL = mean(SVL), sdSVL = sd(SVL), nSVL = length(SVL), 
	mTotL = mean(TotL), sd.TotL = sd(TotL), n.TotL = length(TotL), 
	mMass = mean(Mass), sd.Mass = sd(Mass), n.Mass = length(Mass), 
	mCrest.H = mean(Crest.H), sdCrestH = sd(Crest.H), nCrestH = length(Crest.H) )

setwd("./Analysis/Tables")
write.csv(t(data_sum), file = "Table1.csv")
setwd("../..")

# ------------------------ Table 2 - Crest height, dewlap height and horn length---------------------#

morph_results <- matrix(nrow = 4, ncol = 6)
rownames(morph_results) <- c("Intercept", "log(SVL)", "Sex", "log(SVL)*Sex")
colnames(morph_results) <- rep(c("Est.", "Std.Er"), 3)

# Horn length
morph_results[,1] <- round(as.numeric(modelA1.1_std$coefficients), digits = 2)
morph_results[,2] <- round(as.numeric(coef(summary(modelA1.1_std))[,2]), digits = 3)

# Crest height
morph_results[,3] <- round(as.numeric(modelA1.3_1rlm_std$coefficients), digits = 2)
morph_results[,4] <- round(as.numeric(coef(summary(modelA1.3_1rlm_std))[,2]), digits = 3)

# Dewlap height
morph_results[1:3,5] <- round(as.numeric(modelA1.4_2rlm_std$coefficients), digits = 2)
morph_results[1:3,6] <- round(as.numeric(coef(summary(modelA1.4_2rlm_std))[,2]), digits = 3)

setwd("./Analysis/Tables/")
write.csv(morph_results, file = "Table2.csv")
setwd("../..")

##---------------------------- Create Table 3 - head dimensions models-----------------------------------#

morph_head <- matrix(nrow = 4, ncol = 6)
rownames(morph_head) <- c("Intercept", "log(SVL)", "Sex", "log(SVL)*Sex")
colnames(morph_head) <- rep(c("Est.", "Std.Er"), 3)

## Log head width
morph_head[,1] <- round(as.numeric(modelA1.5_std$coefficients), digits = 2)
morph_head[,2] <- round(as.numeric(coef(summary(modelA1.5_std))[,2]), digits = 3)

## Log head 
morph_head[1:3,3] <- round(as.numeric(modelA1.6_std$coefficients), digits = 2)
morph_head[1:3,4] <- round(as.numeric(coef(summary(modelA1.6_std))[,2]), digits = 3)

## Log head 
morph_head[1:3,5] <- round(as.numeric(modelA1.7_std$coefficients), digits = 2)
morph_head[1:3,6] <- round(as.numeric(coef(summary(modelA1.7_std))[,2]), digits = 3)

setwd("./Analysis/Tables/")
write.csv(morph_head, file = "Table3.csv")
setwd("../..")

#------------------------------ Table 4 - Performance and condition ------------------------------#

table4 <- matrix(nrow = 18, ncol = 4)
rownames(table4)[c(2:9, 11:18)] <- c("Intercept", "Sex (m)", "delta S", "Horn length", "SVL","Head length", "Sex*Head length", "Sex*SVL")
colnames(table4) <- c("Max BF", "Max BF", "Condition", "Condition")
table4[1,] <- rep(c("Est", "Std.Er"), 2)

## Condition
#Throat
table4[1:5, 3] <- round(as.numeric(coefficients(modelA2.cond_std)), digits = 2)
table4[1:5, 4]  <- round(as.numeric(coefficients(summary(modelA2.cond_std))[,2]), digits = 2)

#Labials
table4[11:14, 3] <- round(as.numeric(coefficients(modelA1.cond_std_L)), digits = 2)
table4[11:14, 4]  <- round(as.numeric(coefficients(summary(modelA1.cond_std_L))[,2]), digits = 2)

# Bite force
#Throat
table4[2:9, 1]  <- round(as.numeric(coefficients(modthroat_std)), digits = 2)
table4[2:9, 2]  <- round(as.numeric(coefficients(summary(modthroat_std))[,2]), digits = 2)

#Labials
table4[11:18, 1]  <- round(as.numeric(coefficients(modlabial_std)), digits = 2)
table4[11:18, 2]  <- round(as.numeric(coefficients(summary(modlabial_std))[,2]), digits = 2)

setwd("./Analysis/Tables/")
write.csv(table4, file = "Table4.csv")
setwd("../..")
