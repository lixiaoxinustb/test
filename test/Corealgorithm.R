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