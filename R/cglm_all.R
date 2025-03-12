cglm_all<-function(formula, family, data, alpha=0.05,pval.approx=TRUE,B=100,...)
{
  n<-nrow(data)
  a<-glm(formula, family, data)
  if(names(a$coefficients)[1]=="(Intercept)")
  {
    intercept<-TRUE
  } else{
    intercept<-FALSE
  }
  dip.name<-formula[[2]]
  a<-all.vars(formula)
  if(a[2]==".")
    var.names<-colnames(data)[colnames(data)!=dip.name]
  else
    var.names<-attr(terms.formula(formula),"term.labels")

  p<-length(var.names)

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
  pv.all<-NULL
  bic.all<-NULL
  for(j in 1:length(mod.all))  {
    a.glm<-glm(mod.all[[j]], family, data)
    bic.all<-c(bic.all,BIC(a.glm))
    prs<- sum(residuals(a.glm,type="pearson")^2)
    pearson.all<-c(pearson.all,prs/n)
    if (pval.approx){
      p<-length(coef(a.glm))
      pv.all<-c(pv.all,2*min(pchisq(prs,n-p),pchisq(prs,n-p,lower.tail = FALSE)))
    } else {
      pv.all<-c(pv.all,boot_pval(mod.all[[j]],family,data,B,...))
    }
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

  return(list(models=mod.all, pearsonrisk=pearson.all, pv=pv.all, bic=bic.all, model.opt=mod.opt))
}
