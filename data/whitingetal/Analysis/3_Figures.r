
#--------------------------- Figure 2 -------------------------------#

setwd("./Analysis/Figures")

pch <- c(16, 21)
text<- c("Male", "Female")
cex <- 1.2

pdf(file = "Fig2.pdf", width = 6.49, height =  9.45)
par(mfrow=c(3,2), mar = c(4,4,1.5,1), cex.lab = 1.2)

x   <- c(1.6, 1.6)
y   <- c(2.46, 2.4)
plot(log(Horn.L) ~ log(SVL), data = datA1[datA1$Sex == "m",], ylab = "log Rostrum length", pch = 16, cex = cex, xlab = "", xlim = c(1.59, 2.10))
points(log(Horn.L) ~ log(SVL), data = datA1[datA1$Sex == "f",],  pch = 21, cex = cex)
for(i in 1:2){
	points(x[i], y[i], pch =pch[i], cex = cex)
	text(x[i]+(x[i]*0.01), y[i]-(y[i]*0.0055), text[i], adj = c(0,0))
}
lines(predict_male$fit~log(newdata_m$SVL),  lty = 1, col = "black")
lines(spline(predict_male$fit + (1.96*predict_male$se.fit)~log(newdata_m$SVL)),  lty = 2, col = "black")
lines(spline(predict_male$fit - (1.96*predict_male$se.fit)~log(newdata_m$SVL)),  lty = 2, col = "black")
lines(predict_female$fit~log(newdata_f$SVL),  lty = 1, col = "gray75")
lines(spline(predict_female$fit + (1.96*predict_female$se.fit)~log(newdata_f$SVL)),  lty = 2, col = "gray75")
lines(spline(predict_female$fit - (1.96*predict_female$se.fit)~log(newdata_f$SVL)),  lty = 2, col = "gray75")
mtext("A)", adj = -0.20, cex = 1.2)

y <- c(2.05, 1.96)
x   <- c(1.6, 1.6)

plot(log(Crest.H)~log(SVL), data = datA1[datA1$Sex == "m",], pch = 16, ylim = c(0.39, 2.2), cex = 1.2, xlab = "", ylab = "log Crest height", xlim = c(1.59, 2.1))
points(log(Crest.H)~log(SVL), data = datA1[datA1$Sex == "f",], pch = 21, cex = 1.2)
for(i in 1:2){
	points(x[i], y[i], pch =pch[i], cex = cex)
	text(x[i]+(x[i]*0.01), y[i]-(y[i]*0.004), text[i], adj = c(0,0))
}
mtext("B)", adj = -0.20, cex = 1.2)
pred_plot(predict_male_1, newdata_m_1, "SVL")
pred_plot(predict_female_1, newdata_f_1, "SVL", col = "gray75")


y   <- c(2.94, 2.88)
x   <- c(1.6, 1.6)

plot(log(Dewlap.H)~log(SVL), data = datA1[datA1$Sex == "m",], pch = 16, cex = 1.2, ylab = "log Dewlap height",  xlab = "", ylim = c(1.79, 3), xlim = c(1.59, 2.15))
points(log(Dewlap.H)~log(SVL), data = datA1[datA1$Sex == "f",], pch = 21, cex = 1.2)
for(i in 1:2){
	points(x[i], y[i], pch =pch[i], cex = cex)
	text(x[i]+(x[i]*0.01), y[i]-(y[i]*0.0055), text[i], adj = c(0,0))
}
mtext("C)", adj = -0.20, cex = 1.2)
pred_plot(predict_male_2, newdata_m_2, "SVL")
pred_plot(predict_female_2, newdata_f_2, "SVL", col = "gray75")

x   <- c(1.61, 1.61)
y   <- c(2.58 ,2.55)
plot(log(HW) ~ log(SVL), data = datA1[datA1$Sex == "m",], xlim = c(1.599, 2.1), ylab = "log Head width", xlab = "", pch = 16, cex = cex)
points(log(HW) ~ log(SVL), data = datA1[datA1$Sex == "f",], pch = 21, cex = cex)
for(i in 1:2){
	points(x[i], y[i], pch =pch[i], cex = cex)
	text(x[i]+(x[i]*0.01), y[i]-(y[i]*0.0055), text[i], adj = c(0,0))
}
mtext("D)", adj = -0.20, cex = 1.2)
pred_plot(predict_male.3, newdata_m.3, "SVL")
pred_plot(predict_female.3, newdata_f.3, "SVL", col = "gray75")


y   <- c(3.115 ,3.09)

plot(log(HL) ~ log(SVL), data = datA1[datA1$Sex == "m",], xlim = c(1.59, 2.10), ylab = "log Head length", xlab = "", pch = 16, cex = cex)
points(log(HL) ~ log(SVL), data = datA1[datA1$Sex == "f",], pch = 21, cex = cex)
for(i in 1:2){
	points(x[i], y[i], pch =pch[i], cex = cex)
	text(x[i]+(x[i]*0.01), y[i]-(y[i]*0.004), text[i], adj = c(0,0))
}
mtext("E)", adj = -0.20, cex = 1.2)
pred_plot(predict_male.4, newdata_m.4, "SVL")
pred_plot(predict_female.4, newdata_f.4, "SVL", col = "gray75")

y   <- c(2.55 ,2.51)

plot(log(HH) ~ log(SVL), data = datA1[datA1$Sex == "m",], xlim = c(1.59, 2.10), ylab = "log Head depth", xlab = "", pch = 16, cex = cex)
points(log(HH) ~ log(SVL), data = datA1[datA1$Sex == "f",], pch = 21, cex = cex)
for(i in 1:2){
	points(x[i], y[i], pch =pch[i], cex = cex)
	text(x[i]+(x[i]*0.01), y[i]-(y[i]*0.005), text[i], adj = c(0,0))
}
mtext("F)", adj = -0.20, cex = 1.2)
pred_plot(predict.male7, newdata.male7, "SVL")

mtext("log Snout-vent length (SVL)", side = 1, adj =-3, padj = 3)

dev.off()

setwd("../..")
#--------------------------END--------------------------------#


#----------- Figure 3 - Spectral Figure --------------------#

## Need to add background manually using points etc becuase it overrides plot if we use aggplot a second time as is the case with plot function

## Plot mean spectral curves based on grouping variables
y <- c(avg_barkPLot[,2]+2*se_bark, sort(avg_bark[,2]-2*se_bark, decreasing = TRUE))
x <- c(mbkgd_specs_sm[,1], sort(mbkgd_specs_sm[,1], decreasing = TRUE))

y2 <- c((avg_leavesPlot[,2]+(2*se_leaves)), sort((avg_leaves[,2] - (2*se_leaves)), decreasing = TRUE))
x2 <- c(Lbkgd_specs_sm[,1], sort(Lbkgd_specs_sm[,1], decreasing = TRUE))

sex  <- lapply(regions_specs_sm, function(x) substr(colnames(x), 1, 1)[-1])

# Body regions
setwd("pics")

## problems rescaling with quartz. Height and width don't seem to match. Window size needs to be 23 cm wide and 21.5cm height on screen. Converting to inches doesn't work.
quartz()
par(mfrow = c(2,2), mar = c(3,5,1,0.5), mgp = c(3, 0.5, 0))

#Legend
text <- c("Male", "Female", "Bark", "Leaves")
col  <- c("grey50", "black", "black", "gray50")
linetype <- c(1 ,1, 2, 2)

# Position of pictures
yline   <- c(50, 53, 56, 59)
xline   <- rep(300, 4)
x_m     <- c(0.28, 0.78)
x_f     <- c(0.40, 0.90)
image_f   <- c("female_cerato_labial.png", "female_cerato_throat.png")
image_m <- c("male_ceratophora_labial.png", "male_ceratophora_throat.png")

text.mar <- c("A)", "B)")
adjs <- c(-0.18, -0.20)
padjs <- c(0.10, 0.10)
ylabs <- c("", "")

# Plot A and B
for(i in 1:2){
aggplot(regions_specs_sm[[i]], FUN.error = function(x){2*sd(x)/sqrt(length(x))}, by = sex[[i]], lwd = 2, ylim = c(0, 60), lcol = c("black", "gray50"), shadecol = c("gray40", "gray40"), main = "", xlab = "", ylab = ylabs[i], cex.lab = 1.2)

mtext(text.mar[i], side = 3, adj = adjs[i], padj = padjs[i])

#Add female
w = 0.10
fem <- readPNG(image_f[i])
 fem <- as.raster(fem)
 grid.raster(fem, x = x_f[i], y = 0.92, width = w, height = w*0.6)

## Add male
w = 0.11
male <- readPNG(image_m[i])
 male <- as.raster(male)
 grid.raster(male, x = x_m[i], y = 0.92, width = w, height = w*0.6)

#female symbol
text("\\VE", vfont=c("sans serif","bold"), xpd=TRUE, x = 640, y = 60, cex = 1.2)
#male symbol
text("\\MA", vfont=c("sans serif","bold"), xpd=TRUE, x = 485, y = 60, cex = 1.2)

#add background
points(avg_barkPLot[,2]~mbkgd_specs_sm[,1], col = "black", type = "l", lwd = 2, lty = 2)
points(avg_barkPLot[,2]+2*se_bark~mbkgd_specs_sm[,1], col = "black", type = "l", lwd = 1, lty = 3)
points(avg_barkPLot[,2]-2*se_bark~mbkgd_specs_sm[,1], col = "black", type = "l", lwd = 1, lty = 3)

points(avg_leavesPlot[,2]~Lbkgd_specs_sm[,1], col = "gray50", type = "l", lwd = 2, lty = 2)
points(avg_leavesPlot[,2]+2*se_leaves~Lbkgd_specs_sm[,1], col = "gray50", type = "l", lwd = 1, lty = 3)
points(avg_leavesPlot[,2]-2*se_leaves~Lbkgd_specs_sm[,1], col = "gray50", type = "l", lwd = 1, lty = 3)

#legend
for(i in 1:4){
text(xline[i]+75, yline[i], text[i])
arrows(xline[i], yline[i], xline[i]+30, yline[i], length = 0, col = col[i], lty = linetype[i], lwd = 1.5)
}
}

# Mouth
w = 0.07
x_mou      <- c(0.35, 0.85)
text.mar   <- c("C)", "D)")
mouth      <- regions_specs_sm[3:4]
sex_mouth  <- sex[3:4]
mouth_pics <- c("mouth_roof.png", "mouth_tongue.png")

# Plot C and D
for(i in 1:2){
aggplot(mouth[[i]], FUN.error = function(x){2*sd(x)/sqrt(length(x))}, by = sex_mouth[[i]], lwd = 2, xlab = "", ylim = c(0, 60), lcol = c("black", "gray50"), shadecol = c("gray40", "gray40"), ylab = ylabs[i], main = "", cex.lab = 1.2)

mtext(text.mar[i], side = 3, adj = adjs[i], padj = padjs[i])

#Add gaping mouth
mou <- readPNG(mouth_pics[i])
mou <- as.raster(mou)
grid.raster(mou, x = x_mou[i], y = 0.40, width = w, height = w*1.8)

points(avg_barkPLot[,2]~mbkgd_specs_sm[,1], type = "l", lwd = 2, lty = 2)
points(avg_barkPLot[,2]+2*se_bark~mbkgd_specs_sm[,1], type = "l", lwd = 1, lty = 3)
points(avg_barkPLot[,2]-2*se_bark~mbkgd_specs_sm[,1], type = "l", lwd = 1, lty = 3)

## Col = 139, which is a forest green: see colour palette at http://www.statmethods.net/advgraphs/parameters.html
points(avg_leavesPlot[,2]~Lbkgd_specs_sm[,1], type = "l", col = "gray50", lwd = 2, lty = 2)
points(avg_leavesPlot[,2]+2*se_leaves~Lbkgd_specs_sm[,1], col = "gray50", type = "l", lwd = 1, lty = 3)
points(avg_leavesPlot[,2]-2*se_leaves~Lbkgd_specs_sm[,1], col = "gray50", type = "l", lwd = 1, lty = 3)

#legend
for(i in 1:4){
text(xline[i]+75, yline[i], text[i])
arrows(xline[i], yline[i], xline[i]+30, yline[i], length = 0, col = col[i], lty = linetype[i], lwd = 1.5)
}
}

mtext("Wavelength (nm)", side = 1, adj = -1, padj = 2, cex = 1.2)
mtext("Reflectance (%)", side = 2, adj = 1.5, padj = -29.8, cex = 1.2)

setwd("../Analysis/Figures")
quartz.save(file = "Fig3.pdf", type = "pdf")

setwd("../..")
#---------------------------------------------------------------------------------------------------#

#-------------------------------- Figure 5 ----------------------------#
setwd("./Analysis/Figures")

# bird_dat taken from Spectral_analysis.R file.
#JNDs
dat_birdFigJND <- matrix(nrow = 2, ncol =4)
dat_birdFigJND[1:2, 1] <- bird_dat$roof[,2]
dat_birdFigJND[1:2, 2] <- bird_dat$tong[,2]
dat_birdFigJND[1:2, 3] <- bird_dat$roof[,4]
dat_birdFigJND[1:2, 4] <- bird_dat$tong[,4]
colnames(dat_birdFigJND) <- c("Roof", "Tongue", "Roof", "Tongue")
rownames(dat_birdFigJND) <- c("F", "M")

#SEs
dat_birdFigSE <- matrix(nrow = 2, ncol =4)
dat_birdFigSE[1:2, 1] <- bird_dat$roof[,3]
dat_birdFigSE[1:2, 2] <- bird_dat$tong[,3]
dat_birdFigSE[1:2, 3] <- bird_dat$roof[,5]
dat_birdFigSE[1:2, 4] <- bird_dat$tong[,5]
colnames(dat_birdFigSE) <- c("Roof",  "Tongue", "Roof","Tongue")
rownames(dat_birdFigSE) <- c("F", "M")

pdf(file = "Fig5.pdf", height = 6.88, width = 6.88)

JNDBarplot(dat_birdFigJND, error = dat_birdFigSE, pos = 45, ylim = c(0, 50), ylab = "Just Noticable Differences (JNDs)", name = "", xlab = "Mouth patch", col = c("gray", "black"), cex.lab = 1.2) -> bp.out
arrows(x0 = bp.out[1,1], x1 = bp.out[2,2], y0 = 40, y1 = 40, angle = 0, length = 0)
arrows(x0 = bp.out[1,3], x1 = bp.out[2,4], y0 = 40, y1 = 40, angle = 0, length = 0)
text(expression(Delta ~ "S"), x = bp.out[1,1]+(bp.out[2,2]-bp.out[1,1])/2, y = 41)
text(expression(Delta ~ "L"), x = bp.out[2,4] - ((bp.out[2,4] - bp.out[1,3])/2), y = 41)
box()
#barpairs <- bp.out[2,] - ((bp.out[2,]-bp.out[1,])/2)
#text("NS", x = barpairs, y = 35)

dev.off()
setwd("../..")
#---------------------------------------- END ----------------------------#

#-------------------------- Figure 4 ----------------------------#

setwd("./Analysis/Figures")
pdf(file = "Fig4.pdf", height = 6.88, width = 6.88)

plot(fitted ~ dS, data = labialdat, ylab = "Log Maximum Bite Force", xlab = expression(paste("Labial patch ", Delta, S, " (JNDs)")), ylim = c(1.1,2.6), xlim = c(1, 12.5), cex.lab = 1.5, cex = 1.5, omi = c(1,1.5,0.2,0.5))
points(fitted ~ dS, data = labialdat[labialdat$Sex == "M",], pch = 16, cex = 1.5)
text("Males", x = 1+0.5, y = 1+0.16*10, adj = 0)
text("Females", x = 1+0.5, y = 1+0.16*9.5, adj = 0)
points(x = 1, y = 1+0.16*10, pch = 16, cex = 1.2)
points(x = 1, y = 1+0.16*9.5,  cex = 1.2)

#See lines 523-525 for prediction vectors

lines(spline(pred_male_labBF$fit ~ male_dat$dS), lty = 1, lwd = 1.5) 
lines(spline(pred_male_labBF$fit+(1.96*pred_male_labBF$se.fit) ~ male_dat$dS), lty = 2, lwd = 1.5) 
lines(spline(pred_male_labBF$fit-(1.96*pred_male_labBF$se.fit) ~ male_dat$dS), lty = 2, lwd = 1.5) 

lines(spline(pred_female_labBF$fit ~ female_dat$dS), lty = 1, lwd = 1.5, col = "gray") 
lines(spline(pred_female_labBF$fit+(1.96*pred_female_labBF$se.fit) ~ male_dat$dS), lty = 2, lwd = 1.5, col = "gray") 
lines(spline(pred_female_labBF$fit-(1.96*pred_female_labBF$se.fit) ~ male_dat$dS), lty = 2, lwd = 1.5, col = "gray") 

dev.off()

setwd("../..")
