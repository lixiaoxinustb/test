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