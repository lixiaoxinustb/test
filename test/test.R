data1 <-iris
rownames(data1) <- make.names(data1[,5], unique=TRUE) 

dups <- dim(data1)[1] - length(unique(data1[,5]))
if(dups > 0) {
  warning(paste(dups," duplicated gene symbol(s) found in mixture file!",sep=""))
  rownames(data1) <- make.names(data1[,5], unique=TRUE) 
  ##如果基因标签出现重复，提示警告信息：
}else {rownames(Y) <- Y[,1]}

X <- data.matrix(X) #将数据框转换为数字矩阵，返回值将数据框中的所有变量转换为数字，因素和有序因素被其内部代码所取代
data1 <- data.matrix(data1)
data1<- t(data1)
#order
X <- X[order(rownames(X)),] #按照列名排序
data1 <- data1[order(rownames(data1)),]

P <- perm #number of permutations ##permutations指的是什么

#anti-log if max < 50 in mixture file
if(max(data1) < 50) {data1 <- 2^data1}  #如果样本数据中最大值小于50，所以数值变为2^Y，为什么

#quantile normalization of mixture file分位数归一化混合文件
if(QN == TRUE){
  tmpc <- colnames(data1) #取列名
  tmpr <- rownames(data1)  #取行名
  data1 <- normalize.quantiles(data1) #使用基于分位数的归一化，该函数归一化探针水平强度矩阵，强度矩阵，其中每列对应于芯片，每行是探针
  colnames(data1) <- tmpc
  rownames(data1) <- tmpr
}

#store original mixtures #保存原始混合文件
dataorig <- data1
datamedian <- max(median(dataorig),1)


data1gns <- row.names(data1)

rownames(Y) <- Y[,1]

data1 <- data1[,-5]

data2<-data1[,-5]

data3<- scale(data2)
data4<- scale(data3)

setwd("E:/cibersort/test")
getwd()
X <- read.table("E:/cibersort/test/LM22.txt",header=T,sep="\t",row.names=1,check.names=F) 
Y <- read.table("E:/cibersort/test/sum_3.txt", header=T, sep="\t",check.names=F)
a<- dim(Y)
a
b <- dim(Y)[1]
x <-  c(3:5,11:8,8 + 0:5)
ux <- unique(x)
m0 <- matrix(NA, 4, 0)
rownames(m0)

dups <- dim(Y)[1] - length(unique(Y[,1]))
rownames(Y) <- Y[,1]

Y <- Y[,-1]

DF <- data.frame(a = 1:3, b = letters[10:12],
                 c = seq(as.Date("2004-01-01"), by = "week", len = 3),
                 stringsAsFactors = TRUE)
data.matrix(DF[1:2])
data.matrix(DF)
Y

write.table(predict1, file="predict1.csv", sep="", row.names=T, col.names=T, quote=F)

YX<-Y[YintX,]
data1 <- (data1 - mean(data1)) / sd(as.vector(data1))
as.vector(data1)
pts <- list(x = cars[,1], y = cars[,2])


A <- X
B <- yr
model<-svm(A,B,type="nu-regression",kernel="linear",nu=0.25,scale=F) 

res <- function(i){
  if(i==1){nus <- 0.25}
  if(i==2){nus <- 0.5}
  if(i==3){nus <- 0.75}
  model<-svm(A,B,type="nu-regression",kernel="linear",nu=nus,scale=F)     ###支持向量机，其中参数nu的值用来改变分类算法的准确性，在本算法中，这三个数值都要进行计算
  model
}

out[[1]]$coefs
out[[1]]$SV
weights = t(out[[3]]$coefs) %*% out[[3]]$SV

nusvm[1] 
nusvm[3] <- sqrt((mean((k - yr)^2)))
corrv[3] <- cor(k, yr)

itor <- 2

L= x[,-1]
setwd('E:/cibersort/test')