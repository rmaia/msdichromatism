
## Functions
## problems with text files and need to remove the last line from each file. Function below will do this and re-write the files so that they have the same file names and format.

spec_proc <- function(directory, outputdirect){
files <- list.files(directory)

for(i in 1:length(files)){
		file <- readLines(paste(directory, files[i], sep = "/"), n = 2048)
		setwd(outputdirect)
		write(file, file = files[i])
	}
}

# Use spec_proc function to do this processing
#spec_proc(directory = "~/Dropbox/Behavioural ecology of Ceratophora tennenttii/Data/Spec Files/Individuals_Pooled", "~/Dropbox/Behavioural ecology of Ceratophora tennenttii/Data/Spec Files/Processed")


## Plotting

pred_plot <- function(predictions, newdata, x_var, col = "black"){
lines(predictions$fit~log(newdata[, x_var]),  lty = 1, col = col)
lines(spline(predictions$fit + (1.96* predictions$se.fit)~log(newdata[, x_var])),  lty = 2, col = col)
lines(spline(predictions$fit - (1.96* predictions$se.fit)~log(newdata[, x_var])),  lty = 2, col = col)
}


## PLotting JNDs
JNDBarplot <- function(data, name = "region", error = 1, pos, fontsize = 2, ...){
	barplot(data, beside = TRUE, legend.text = rownames(data), args.legend = list(bty = "n"), ...) -> bp.out
	abline(h = 0)
	arrows(x0 = bp.out, y0 = data, x1 = bp.out, y1 = data+error, length = 0.05, angle = 90)
	arrows(x0 = bp.out, y0 = data, x1 = bp.out, y1 = data-error, length = 0.05, angle = 90)
	text(name, x = max(bp.out)/2, y = pos, cex = fontsize)
	bp.out
}
## Obtaining residuals between two morphological traits to obtain relative values. Y is your response, X your covariate and because morphological traits log = TRUE. But can change to FALSE for raw variables

rel_var <- function(data, Y, X, log = TRUE){
if(log == TRUE){	
	res <- residuals(lm(log(data[,Y])~log(data[,X])))
	data[,paste("rel_", Y, sep ="")] <- res
}else{
	res <- residuals(lm(data[,Y]~data[,X]))
	data[,paste("rel_", Y, sep ="")] <- res
}
return(data)
}

# Need to create an numeric ID vector so we can match IDs up. # Basically below says find M or F (i.e. [MF]) followd by a 0, then extract () subsequent digits from 0-9 until the end. Then substidute only the () part of the regular expression (i.e. use \\1)

regval <- function(data){
ind_num <- gsub("[MF]0([0-9]+)", "\\1", data$patch1)
ind_num_final <- as.numeric(gsub("0([0-9]+)", "\\1", ind_num))

return(ind_num_final)
}
