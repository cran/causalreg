cglm_all<-function(formula, family, data, alpha=0.05,pval=c("bootstrap","chi-square"),B=100,...)
{
  n<-nrow(data)
  a<-glm(formula, family, data,...)
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
    a.glm<-glm(mod.all[[j]], family=family, data=data)
    bic.all<-c(bic.all,BIC(a.glm))
    prs<- sum(residuals(a.glm,type="pearson")^2)
    pearson.all<-c(pearson.all,prs/n)
    if (pval == "chi-square"){
      p<-length(coef(a.glm))
      pv.all<-c(pv.all,2*min(pchisq(prs,n-p),pchisq(prs,n-p,lower.tail = FALSE)))
    }
    if (pval == "bootstrap") {
      pv.all<-c(pv.all,boot_pval(mod.all[[j]],family=family,data=data,B=B,...))
    }
  }

  if (sum(pv.all>alpha)>0){
    model.potential<-mod.all[pv.all>alpha]
    mod.opt<-model.potential[[which.min(bic.all[pv.all>alpha])]]
    mod.opt<-deparse1(mod.opt)
  } else {
    mod.opt <- "no potential causal model found"
  }


  var.cat<-colnames(data)[sapply(data, function(x) !is.numeric(x))]
  var.bin<-colnames(data)[sapply(data, function(x) length(unique(na.omit(x))) == 2)]
  var.cat<-union(var.cat,var.bin)
  if(mod.opt!="no potential causal model found")
  {
  var.mod<-attr(terms.formula(as.formula(mod.opt),data=data),"term.labels")
  var.noncat <- setdiff(var.mod, var.cat)
  }
  else{
    var.noncat<-var.mod<-var.cat
  }

  if (family == "binomial" & length(var.noncat) < length(var.mod))
  {
    mod.test <- reformulate(var.noncat, response = all.vars(formula(mod.opt))[1])
    fmli<-update.formula(as.formula(mod.opt),mod.test)
    pv<-pv.all[[which(mod.all==fmli)]]
    if(pv>alpha)
    {
      mod.opt<-deparse1(fmli)
    }
    else
    {
      varc<-var.cat[var.cat %in% var.mod]
      for(i in 1:length(varc))
      {
        fmli<-update.formula(mod.opt, paste0("~. -", varc[i]))
        pv<-pv.all[[which(mod.all==fmli)]]
        if(pv>alpha)
        {
          mod.opt<-deparse1(fmli)
        }
      }
    }
  }

  for(j in 1:length(mod.all))
    mod.all[[j]]<-deparse1(mod.all[[j]])

  return(list(models=mod.all, pearsonrisk=pearson.all, pv=pv.all, bic=bic.all, model.opt=mod.opt))
}
