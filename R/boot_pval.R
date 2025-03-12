boot_pval<-function(formula,family, data, B=100, ...)
{
  n <- nrow(data)
  idx <- base::sample(x=1:n, size=n*B,replace=TRUE)
  boot.dat<-array(t(data[idx,]),dim=c(ncol(data),n,B))
  rownames(boot.dat)<-colnames(data)

  vrs<-all.vars(formula)
  n<-nrow(data)
  p<-length(vrs)-1
  ix<-as.numeric(lapply(data,is.numeric))
  pr<-NULL
  for(i in 1:B)
  {
    dat<-data.frame(t(boot.dat[,,i]))
    dat[,ix==1]<-apply(dat[,ix==1],2,as.numeric)
    if (vrs[2]=="."){
      h<-glm(formula, family, dat,...)
    } else {
      h <-gam(formula, family, dat,...)
    }
    pr<-c(pr,sum(residuals(h,type="pearson")^2)/n)
  }
  prob<-mean(1<=pr)
  pv<-2*min(prob,1-prob)
  return(pv)
}
