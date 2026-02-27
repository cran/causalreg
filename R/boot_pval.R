boot_pval<-function(formula, family, data, B=100, ...)
{
  n <- nrow(data)
  vrs<-all.vars(formula)
  pr<-NULL
  for(i in 1:B)
  {
    idx <- base::sample(x=1:n,replace=TRUE)
    dat<-data[idx,]
    if (vrs[2]=="."){
      h<-glm(formula=formula, family=family, data=dat,...)
    } else {
      h <-gam(formula=formula, family=family, data=dat,...)
    }
    pr<-c(pr,sum(residuals(h,type="pearson")^2)/n)
  }
  prob<-mean(1<=pr)
  pv<-2*min(prob,1-prob)
  return(pv)
}
