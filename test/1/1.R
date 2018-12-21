



a <- "E://cibersort//test//sum_513.txt"
Y <- read.table(a, header=T, sep="\t",check.names=F)
str(Y)
data5 <- data2$Genesymbol
data6 <- as.numeric(data5)
str(data7)
data7 <- read.table("./sum_3.txt")
str(y)

B <- as.factor(Y[,'normalized_count'])
A <- as.factor(Y[,'Genesymbol'])
C <- merge(A,B)


colApply <- function(dat, cols = colnames(dat), func = as.vector) { 
  dat[cols] <- lapply(dat[cols], func) 
  return(dat)
} 

 dat <- data.frame(data4, stringsAsFactors = F)
 
 data5 <- colApply(dat, cols = colnames(dat), func = as.factor)

 