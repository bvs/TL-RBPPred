### THis function combines features extracted from executable file and seqCTD function ###
combineFeatures <- function(seq,csv) {
seq <- readFASTA(seq)
ctd <- seqCTD(seq)
seqCSV <- read.csv(csv,header=FALSE)
seqCSV <- t(seqCSV[,-1])
mat <- rbind(ctd,seqCSV)
colnames(mat) <- names(seq)
rownames(mat) <- NULL
return(mat)
}
### End of combineFeatures ###
