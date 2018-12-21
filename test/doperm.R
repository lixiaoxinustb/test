#do permutations ָ����ʲô��
doPerm <- function(perm, X, Y, absolute, abs_method){
  itor <- 1
  Ylist <- as.list(data.matrix(Y)) ###����һ���б������
  dist <- matrix()##һ���յľ���
  
  while(itor <= perm){
    #print(itor)
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    ##sample���������������������һ��������ʾ��Χ���ڶ���������ʾ��ȡ�����������˴�����length(Ylist)��ô�������������ȡdim(X)[1]��
    ##���б���ȡ����Щ�����ȡ��λ���ϵ�ֵ������һ����ֵ������
    
    #standardize mixture ### ��׼������
    yr <- (yr - mean(yr)) / sd(yr) 
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, absolute, abs_method)
    
    mix_r <- result$mix_r
    
    #store correlation  
    if(itor == 1) {dist <- mix_r}else{dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }#��CoreAlg�еõ���mix_r������dist�У�һ����perm����ֵ
  newList <- list("dist" = dist)
}