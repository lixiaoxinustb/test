CoreAlg <- function(X, y, absolute, abs_method){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)     ###֧�������������в���nu��ֵ�����ı�����㷨��׼ȷ�ԣ��ڱ��㷨�У���������ֵ��Ҫ���м���
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') 
    out <- mclapply(1:svn_itor, res, mc.cores=1) else ###mclapply�ǲ��м��㺯������windowsϵͳ���޷�������˲��м��㣬�������windows�Ϳ���
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)  ##### 
  ##### ��������������������µ�֧������ģ��
  nusvm <- rep(0,svn_itor) ###�ظ�����,��һ�������ظ������ݣ��ڶ��������ظ��Ĵ���
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV #��ȡģ�͵�coefs��SV�����ΪȨ��
    weights[which(weights<0)]<-0  ###��С��0��Ȩ�ظ�ֵΪ0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*') ###sweep�����Ծ�����м��㣬MARGIN=2�����У�MARGIN=1�����У��˴�Ϊ��x�����г���W
    k <- apply(u, 1, sum) ##��u�������
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }  #���õ� nusvm��corrv������ֱ��Ӧ����ģ���µ���ֵ
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]] ###����С��nusvmֵ��ѡ��Ӧ��ģ��
  
  #get and normalize coefficients #��ò��ҹ�һ��ϵ��
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
  if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}