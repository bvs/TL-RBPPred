### THis function combines features extracted from feature116.exe  and seqCTD.R ###
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
