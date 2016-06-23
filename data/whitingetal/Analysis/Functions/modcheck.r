
modcheck <- function(mod, ...){
#Install packages
	packages <- c("car", "plyr")
	sapply(packages, function(x) library(x, character = TRUE))

	checks <- list()

	# Run model checks
		# Residuals
		checks$resid       <- resid(mod)
		checks$ShapiroTest <- (shapiro.test(checks$resid))
	
	# Plots
		quartz()
		par(mfrow = c(2,2))
		hist(checks$resid, main = "Histogram of model residuals", xlab = "Residuals")
		plot(fitted(mod) ~ checks$resid, ylab = "Predicted Y from model", xlab = "Residuals")
		acf(checks$resid, main = "Autocorrelation of Residuals", ...)
	
	# Influence plots
		checks$InfluPlot <- influencePlot(mod)

	quartz()
		checks$resPlot   <- residualPlots(mod)

checks
}
