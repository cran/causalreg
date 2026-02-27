cgam_all<-function(formula, family, data, alpha=0.05,pval=c("bootstrap","chi-square"),B=100,...)
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
      a.gam<-glm(mod.all[[j]], family=family, data=data, ...)
      edf<-length(a.gam$coefficients)
    } else{
      a.gam<-gam(mod.all[[j]], family=family, data=data, ...)
      sa.gam<-summary(a.gam)
      edf<-length(sa.gam$p.coeff)+sum(sa.gam$edf)
    }
    if (pval=="chi-square"){
      pearson.all[j]<- sum(residuals(a.gam,type="pearson")^2)
      pv.all[j]<-2*min(pchisq(pearson.all[j],n-edf),pchisq(pearson.all[j],n-edf,lower.tail = FALSE))
    }
    if(pval=="bootstrap") {
      pv.all[j]<-boot_pval(mod.all[[j]],family=family,data=data,B=B,...)
    }
    bic.all<-c(bic.all,BIC(a.gam))

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
    #we assume that no splines are used on categorical variables
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

  return(list(models=mod.all, pearsonrisk=pearson.all, bic=bic.all, pv=pv.all,model.opt=mod.opt))
}
