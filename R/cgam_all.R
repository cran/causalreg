cgam_all<-function(formula, family, data, alpha=0.05,pval.approx=TRUE,B=100,...)
{
  n<-nrow(data)
  vrs<-all.vars(formula)

  dip.name<-formula[[2]]
  if(vrs[2]=="."){
    var.names<-colnames(data)[colnames(data)!=dip.name]
  } else {
    var.names<-attr(terms.formula(formula,data=data),"term.labels")
  }

  p<-length(var.names)


  if (vrs[2]=="."){
    a<-glm(formula=formula, family=family, data=data, ...)
  } else {
    a<-gam(formula=formula, family=family, data=data, ...)
  }
  if(names(a$coefficients)[1]=="(Intercept)"){
    intercept<-TRUE
  } else{
    intercept<-FALSE
  }


  mod.all<-list()
  for (i in 1:p){
    vc <- combn(var.names,i)
    for (j in 1:ncol(vc)){
      if(intercept)
        mod <- as.formula(paste0(dip.name,"~", paste0(vc[,j], collapse = "+")))
      else
        mod <- as.formula(paste0(dip.name,"~", paste0(vc[,j], collapse = "+"),"-1"))
      mod.all <- c(mod.all, mod)
    }
  }

  pearson.all<-NULL
  bic.all<-NULL
  pv.all<-NULL
  for(j in 1:length(mod.all))  {
    if (vrs[2]=="."){
      a.gam<-glm(mod.all[[j]], family, data, ...)
      edf<-length(a.gam$coefficients)
    } else{
      a.gam<-gam(mod.all[[j]], family, data, ...)
      sa.gam<-summary(a.gam)
      edf<-length(sa.gam$p.coeff)+sum(sa.gam$edf)
    }
    pearson.all[j]<- sum(residuals(a.gam,type="pearson")^2)
    if (pval.approx){
      pv.all[j]<-2*min(pchisq(pearson.all[j],n-edf),pchisq(pearson.all[j],n-edf,lower.tail = FALSE))
    } else {
      pv.all[j]<-boot_pval(mod.all[[j]],family,data,B,...)
    }
    bic.all<-c(bic.all,BIC(a.gam))

  }

  if (sum(pv.all>alpha)>0){
    model.potential<-mod.all[pv.all>alpha]
    mod.opt<-model.potential[[which.min(bic.all[pv.all>alpha])]]
    mod.opt<-format(as.formula(deparse(mod.opt),env=NULL))
  } else {
    mod.opt <- "no potential causal model found"
  }


  for(j in 1:length(mod.all))
    mod.all[[j]]<-format(as.formula(deparse(mod.all[[j]]),env=NULL))

  return(list(models=mod.all, pearsonrisk=pearson.all, bic=bic.all, pv=pv.all,model.opt=mod.opt))
}
