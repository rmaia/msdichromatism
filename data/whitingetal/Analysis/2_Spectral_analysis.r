
## Load functions
source("./Analysis/Functions/func.R")

## Import spec data 
specs <- getspec(where = "Data/Spec Files/Processed", lim = c(300,700), ext = 'txt')

## BODY REGIONS - extract the body regions of interest for the analysis.
## Body region vector

regions    <- c("labial", "throat", "roof", "tong")

## Subset data and store in a list
subsetdat <- list()

for(i in 1:length(regions)){
	subsetdat [[i]] <- specs[,c("wl", grep(regions[i], colnames(specs), value = TRUE))]	
}
names(subsetdat) <- regions

# Body regions - Smooth spectral curves and prepare for plotting and visual modelling. Seems to be a bug in the procspec function. It first needs to smooth then zero the negative numbers, but it seems to be zeroing then smoothing and this results in negative numbers re-introduced. When both arguments are used this should automatically do what I do below to avoid negative values being re-introduced into the data. 

regions_specs_sm <- list()

for(i in 1:length(regions)){
regions_specs_sm[[i]]   <- procspec(subsetdat[[regions[i]]], opt = "smooth" , span = 0.2)
}
names(regions_specs_sm) <- regions

for(i in 1:length(regions)){
regions_specs_sm[[i]]   <- procspec(regions_specs_sm[[regions[i]]], fixneg = 'zero')
}
names(regions_specs_sm) <- regions


### BACKGROUND
bark_backgrd   <- getspec(where = "Data/Spec Files/Background/Bark_Processed")
leaves_backgrd <- getspec(where = "Data/Spec Files/Background/Leaves_Processed")

# Background - Smooth spectral curves and prepare for plotting and visual modelling
mbkgd_specs_sm <- procspec(bark_backgrd, opt = "smooth" , span = 0.2)
Lbkgd_specs_sm <- procspec(leaves_backgrd, opt = "smooth" , span = 0.2)

## Used for Spectral Figure to ensure smoothed plots. 
avg_barkPLot     <- aggspec(mbkgd_specs_sm)
avg_leavesPlot   <- aggspec(Lbkgd_specs_sm)

mbkgd_specs_sm <- procspec(bark_backgrd, fixneg = 'zero')
Lbkgd_specs_sm <- procspec(leaves_backgrd, fixneg = 'zero')

# Bark
avg_bark <- aggspec(mbkgd_specs_sm)
se_bark  <- apply(mbkgd_specs_sm[,2:ncol(mbkgd_specs_sm)], 1, function(x) sd(x)/sqrt(length(x)))

# Leaves
avg_leaves <- aggspec(Lbkgd_specs_sm)
se_leaves  <- apply(Lbkgd_specs_sm[,2:ncol(Lbkgd_specs_sm)], 1, function(x) sd(x)/sqrt(length(x)))

#--------------------------- INDIVIDUAL SPECTRAL ANALYSIS -------------------------------#

## Condense spectral curves across individuals - i.e. average spectra within individuals for body regions

ind_agg <- list()

ind <- lapply(regions_specs_sm, function(x) {substr(names(x), 1, 4)})

for(i in 1:length(regions_specs_sm)){
ind_agg[[i]] <- aggspec(regions_specs_sm[[i]], by = ind[[i]])
}
names(ind_agg) <- regions

# Add individual background - avg_leaves

name <- names(ind_agg)

addbkg <- function(dat){
	dat$Leaves        <- avg_leaves[,2] 
}
 ind_agg <- lapply(ind_agg, function(x) cbind(x, Leaves = avg_leaves[,2]))

## VISUAL MODELING 

# Create a dataframe containing the visual sensitivities of cones from lizards - taken from pg 1336 of Chan et al. 2009. Behavioural Ecology, 20: 1334-1342

# Lizard visual system
liz_vis <- sensmodel(c(360, 440, 493, 571)) 

visual_ind_liz    <- lapply(ind_agg, function(x) vismodel(x, visual = liz_vis, illum= "forestshade", bkg = avg_leaves[,2], relative = FALSE, vonkries = TRUE, qcatch = "fi"))

coldistance_liz   <- lapply(visual_ind_liz, function(x) coldist(x, vis = "tetra", noise = "neural", n1 = 1, n2 = 1, n3 = 3.5, n4 = 6, v = 0.10))

ind_vs_bkg_liz    <- lapply(coldistance_liz, function(x) x[x$patch2 == "Leaves",])

lapply(ind_vs_bkg_liz, tail)

## Now that we have the ds and dl contrasts against the background using a lizard visual system we can test whether "greater contrast against a background based on a lizard visual system" is associated with morphology and performance controlling for sex differences. We are making an important assumption that ds between leaves and lizard patches across lizards can be detected by lizards. In otherwords, the difference between ds for each lizard are discriminable. Pretty big assumption, but seems to be used in the colour literature (See Chan et al. 2009. BE, 20:1334-1342 as example. They implicitly assume that higher ds values of females are more discrimninable as they increase based on measurements across the season.)

# Bird visual system
visual_ind_bird  <- lapply(ind_agg, function(x) vismodel(x, visual = "avg.uv", illum= "forestshade", bkg = avg_leaves[,2], relative = FALSE, vonkries = TRUE, qcatch = "fi"))

coldistance_bird <- lapply(visual_ind_bird, function(x) coldist(x, vis = "tetra", noise = "neural", n1 = 1, n2 = 2, n3 = 3,n4 = 3, v = 0.1))

ind_vs_bkg_bird  <- lapply(coldistance_bird, function(x) x[x$patch2 == "Leaves",])

# Bite force - require datA3

# Add ind ID and Sex to the data frames so that they can be subsetted and grouped as well as matched with IDs

ind_ID_col <- lapply(ind_vs_bkg_liz, function(x) regval(x))
ind_Sex_col <- lapply(ind_vs_bkg_liz, function(x) gsub("([MF])[0-9]+","\\1", x$patch1))

names <- names(ind_vs_bkg_liz) 

for(i in 1:4){
ind_vs_bkg_liz[[i]]$ID  <- ind_ID_col[[i]]
ind_vs_bkg_liz[[i]]$Sex <- ind_Sex_col[[i]]
}

# Same for birds

for(i in 1:4){
ind_vs_bkg_bird[[i]]$ID  <- ind_ID_col[[i]]
ind_vs_bkg_bird[[i]]$Sex <- ind_Sex_col[[i]]
}

## Do males and females differ is dS, LIZARD

liz_dat <- lapply(ind_vs_bkg_liz, function(x) ddply(x, .(Sex), summarise, meandS = mean(dS), se = sd(dS)/sqrt(length(dS))))

# Non- normal so be safe and use non-parametric test
par(mfrow = c(2,2))
lapply(ind_vs_bkg_liz, function(x) hist(log(x$dS)))

#dS
lapply(ind_vs_bkg_liz, function(x) wilcox.test(x$dS ~ x$Sex))
#dL
lapply(ind_vs_bkg_liz, function(x) wilcox.test(x$dL ~ x$Sex)) 

## Do males and females differ is dS, BIRD

bird_dat <- lapply(ind_vs_bkg_bird, function(x) ddply(x, .(Sex), summarise, meandS = mean(dS), seS = sd(dS)/sqrt(length(dS)), meandL = mean(dL), seL = sd(dL)/sqrt(length(dL))))
#dS
lapply(ind_vs_bkg_bird, function(x) wilcox.test(x$dS ~ x$Sex))
#dL
lapply(ind_vs_bkg_bird, function(x) wilcox.test(x$dL ~ x$Sex))

# Does chromatic contrast using a bird visual system differ between the sexes? Is this robust to outliers?

males_out <- c("M034", "M033")

bird_dat <- lapply(ind_vs_bkg_bird, function(x) ddply(x[!x$patch1 %in% males_out,], .(Sex), summarise, meandS = mean(dS), se = sd(dS)/sqrt(length(dS))))

lapply(ind_vs_bkg_bird, function(x) t.test(log(x[!x$patch1 %in% males_out,]$dS) ~ x[!x$patch1 %in% males_out,]$Sex, var.equal = FALSE))
## From this point we can now match the data up using the ind_vs_bkg list of dataframes GO TO 1_Analysis_Script_Noble.R ##