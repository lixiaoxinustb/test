setwd("E:/cibersort/test/") 
getwd()
data1<-read.csv("./LM22.csv")
data2<- read.table(file="sum.txt", header=T,blank.lines.skip=F, sep='\t')
name1<-data1$Genesymbol
write.csv(name1,"./name1.csv")
name1<-read.csv("./name1.csv")

data3 <- subset(data2,Genesymbol %in% name1$x) 
write.table(data5,"data5.txt",quote = FALSE,row.names = FALSE,sep="\t")

data9 <- read.table("./o1.txt",header=T)
data10 <- as.vector(data9)
data4 <- as.data.frame(data2)


data8 <- as.vector(data2)
 
write.table(Y,"out.txt",quote = FALSE,row.names = FALSE)

a <- 'E://cibersort//test//sum_3.txt'
Y <- read.table(a, header=T, sep="\t",check.names=F)



