#dependencies
library(e1071)
library(parallel)
library(preprocessCore)

#Core algorithm 核心算法
CoreAlg <- function(X, y, absolute, abs_method){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)     ###支持向量机，其中参数nu的值用来改变分类算法的准确性，在本算法中，这三个数值都要进行计算
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') 
    out <- mclapply(1:svn_itor, res, mc.cores=1) else ###mclapply是并行计算函数，在windows系统下无法开启多核并行计算，如果不是windows就可以
      out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)  ##### 
  ##### 这里输出的是三个参数下的支持向量模型
  nusvm <- rep(0,svn_itor) ###重复函数,第一个参数重复的内容，第二个参数重复的次数
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV #提取模型的coefs和SV相乘作为权重
    weights[which(weights<0)]<-0  ###把小于0的权重赋值为0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*') ###sweep函数对矩阵进行计算，MARGIN=2代表列，MARGIN=1代表行，此处为在x矩阵按列乘以W
    k <- apply(u, 1, sum) ##对u按行求和
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }  #最后得到 nusvm和corrv，里面分别对应三个模型下的数值
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]] ###以最小的nusvm值挑选对应的模型
  
  #get and normalize coefficients #获得并且归一化系数
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
  if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}



#do permutations 指的是什么？
doPerm <- function(perm, X, Y, absolute, abs_method){
  itor <- 1
  Ylist <- as.list(data.matrix(Y)) ###生成一个列表，奇怪
  dist <- matrix()##一个空的矩阵
  while(itor <= perm){
    #print(itor)
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    ##sample函数是随机抽样函数，第一个参数表示范围，第二个参数表示抽取样本次数，此处是在length(Ylist)这么多样本中随机抽取dim(X)[1]次
    ##从列表中取出这些随机抽取的位置上的值，生成一个数值型向量
    #standardize mixture ### 标准化矩阵
    yr <- (yr - mean(yr)) / sd(yr)
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, absolute, abs_method)
    mix_r <- result$mix_r
    #store correlation  
    if(itor == 1) {dist <- mix_r}else{dist <- rbind(dist, mix_r)}
    itor <- itor + 1
  }#将CoreAlg中得到的mix_r保存在dist中，一共有perm个数值
  newList <- list("dist" = dist)
}



#main function 主要函数
CIBERSORT <- function(sig_matrix, mixture_file, perm=100, QN=TRUE, absolute=FALSE, abs_method='sig.score'){
  
  if(absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score') stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
  
  #read in data
  X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F) 
  #22个免疫细胞参考
  Y <- read.table(mixture_file, header=T, sep="\t",check.names=F)
  #样本
  #to prevent crashing on duplicated gene symbols, add unique numbers to identical names
  #为防止重复基因符号导致崩溃，请为相同的名称添加唯一的数字
  dups <- dim(Y)[1] - length(unique(Y[,1])) #得到基因标签中重复的数量
  ## dim(Y) 结果为行数和列数两个数，dim(Y)[1]取出第一个数为行数
  ##unique()返回向量，数据框或x之类的数组，但删除了重复的元素/行
  ##获取或设置向量（包括列表）和因子的长度
  if(dups > 0) {
    warning(paste(dups," duplicated gene symbol(s) found in mixture file!",sep=""))
    rownames(Y) <- make.names(Y[,1], unique=TRUE) 
    ##如果基因标签出现重复，提示警告信息：
  }else {rownames(Y) <- Y[,1]}
  Y <- Y[,-1] ##此时，带有重复的基因标签列已经换为唯一无重复的基因标签列
  
  X <- data.matrix(X) #将数据框转换为数字矩阵，返回值将数据框中的所有变量转换为数字，因素和有序因素被其内部代码所取代
  Y <- data.matrix(Y)
  
  #order
  X <- X[order(rownames(X)),] #按照列名排序
  Y <- Y[order(rownames(Y)),]
  
  P <- perm #number of permutations ##permutations指的是什么
  
  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}  #如果样本数据中最大值小于50，所以数值变为2^Y，为什么
  
  #quantile normalization of mixture file分位数归一化混合文件
  if(QN == TRUE){
    tmpc <- colnames(Y) #取列名
    tmpr <- rownames(Y)  #取行名
    Y <- normalize.quantiles(Y) #使用基于分位数的归一化，该函数归一化探针水平强度矩阵，强度矩阵，其中每列对应于芯片，每行是探针
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #store original mixtures #保存原始混合文件
  Yorig <- Y
  Ymedian <- max(median(Yorig),1) #求混合文件的中位数，与1相比取最大值赋给Ymedian
  
  #intersect genes  相交基因
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns  #混合文件的基因标签是否包含于参考文件的基因标签中，为真输出TURE，否者输出FALSE
  ###这里损失了三个基因标签 (3-Mar 5-Sep  8-Sep)为什么
  Y <- Y[YintX,]  ###在混合文件中取出包含于参考文件基因标签的行,
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix 标准化签名矩阵
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients相关系数的经验零点分布
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)}
  ##将doPerm中得到的dist排序，默认升序
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures 迭代混合物
  while(itor <= mixtures){
    y <- Y[,itor]
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    #run SVR core algorithm
    result <- CoreAlg(X, y, absolute, abs_method)
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    if(absolute && abs_method == 'sig.score') {
      w <- w * median(Y[,itor]) / Ymedian
    }
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(absolute) out <- c(out, sum(w))
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    itor <- itor + 1
  }
  
  #save results
  write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  if(!absolute){colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")}
  else{colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",paste('Absolute score (',abs_method,')',sep=""))}
  obj
}
