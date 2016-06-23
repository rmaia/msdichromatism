
getwd()

rm(list= ls())
setwd("../..")

inputFold = "~/Desktop/Lizards/"
outputFold = "~/Desktop/ProcFiles/"

renameFile <- function(inputFold, outputFold){
	folder_name <- list.files(inputFold)
	dir         <- paste0(inputFold,folder_name, "/")

		for(i in 1:length(folder_name)){
			liz       <- folder_name[i]
			directory <- dir[i]
			fileList  <- list.files(directory)
	## write.table can take a file connection in teh "file =" argument. This allows you to modify the path so that files can be written in a different directory. 
			
				setwd(directory)
				for(j in 1:length(fileList)){
					filedat <- read.table(fileList[j])
					write.table(filedat, file = paste0(outputFold, liz, ".", fileList[j]), col.names = FALSE, row.names = FALSE, sep = "\t")
				}
		}
}

renameFile(inputFold, outputFold)

