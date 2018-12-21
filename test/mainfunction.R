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