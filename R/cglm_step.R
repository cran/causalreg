cglm_step<-function(formula, family, data, alpha=0.05, pval=c("bootstrap","chi-square"), B=100,...)
{
  n<-nrow(data)
  a<-glm(formula, family=family, data=data,...)
  if(names(a$coefficients)[1]=="(Intercept)")
  {
    intercept<-TRUE
  } else
  {
    intercept<-FALSE
  }

  dip.name<-formula[[2]]
  a<-all.vars(formula)
  if(a[2]==".")
    var.names<-colnames(data)[colnames(data)!=dip.name]
  else
    var.names<-attr(terms.formula(formula,data=data),"term.labels")



  p<-length(var.names)

  mod.step<-vector(mode="list")
  pv.step<-NULL

  if(intercept)
  {
    mod.step[[1]]<-as.formula(paste0(dip.name,"~","1"))
    pv.step[1]<-0
  } else
  {
    mod.step[[1]]<-as.formula(paste0(dip.name,"~","-1"))
    pv.step[1]<-0
  }


  continue<- TRUE
  var.namesyes<-var.names #for potential inclusion
  var.namesmod<-NULL #already in the model

  j<-0
  while(continue & length(var.namesyes)>0)
  {
    j<-j+1
    #check which variable to add
    pvals<-NULL
    for(i in 1: length(var.namesyes))
    {
      fmli<-update.formula(mod.step[[j]], paste0("~. +", var.namesyes[i]))
      a.glm<-glm(fmli, family=family, data=data,...)
      p.step<-length(coef(a.glm))
      ps<-sum(residuals(a.glm,type="pearson")^2)
      if (pval == "chi-square"){
        pv<-2*min(pchisq(ps,n-p.step),pchisq(ps,n-p.step,lower.tail = FALSE))
      }
      if (pval == "bootstrap") {
        pv<-boot_pval(fmli,family=family, data=data, B=B,...)
      }
      pvals<-c(pvals,pv)
    }
    i.min<-which.max(pvals)
    mod.min<-update.formula(mod.step[[j]], paste0("~. +", var.namesyes[i.min]))
    #rearrange all variables names
    a<-attr(terms.formula(mod.min),"term.labels")
    mod.min<-reformulate(sort(a),response=dip.name,intercept=intercept)
    a.glm<-glm(mod.min, family=family, data=data,...)
    p.step<-length(coef(a.glm))
    ps<- sum(residuals(a.glm,type="pearson")^2)
    if(pval == "chi-square"){
      pv.min<-  2*min(pchisq(ps,n-p.step),pchisq(ps,n-p.step,lower.tail = FALSE))}
    if (pval == "bootstrap"){
      pv.min<-boot_pval(mod.min, family=family, data=data, B=B,...)
    }
    pval.diff<-pv.min-pv.step[j]
    continue <- (pval.diff>0)|(pv.min>alpha)
    if(continue)
    {
      mod.step[[j+1]]<-mod.min
      pv.step[j+1]<-pv.min
      var.namesmod<-c(var.namesmod,var.namesyes[i.min])
      var.namesyes<-var.namesyes[-i.min]
    }
  }

  model.current<-mod.step[[length(mod.step)]]
  current.glm<-glm(model.current, family=family, data=data,...)
  bic.current<-BIC(current.glm)
  j<-j-1
  improve = TRUE
  while (improve){
    j<-j+1
    if(length(var.namesmod) >= 1)
    {
      bics<-NULL
      for(i in 1:length(var.namesmod))
      {
        fmli<-update.formula(model.current, paste0("~. -", var.namesmod[i]))
        bics[i]<-BIC(glm(fmli, family=family, data=data,...))
      }
      i.min<-which.min(bics)
      mod.min<-update.formula(model.current, paste0("~. -", var.namesmod[i.min]))
      bic.min<-BIC(glm(mod.min, family=family, data=data,...))
      if(bic.min < bic.current)
      {
        model.current<-mod.min
        #rearrange all variables names
        a<-attr(terms.formula(model.current),"term.labels")
        if(length(a)>0)
        {
          model.current<-reformulate(sort(a),response=dip.name,intercept=intercept)
        }
        bic.current<-bic.min
        var.namesmod<-var.namesmod[-i.min]
        mod.step[[j+1]]<-model.current
      }
      else
      {
        improve <- FALSE
      }
    }
    else
      improve<-FALSE
  }

  mod.opt<-mod.step[[length(mod.step)]]

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
    a.glm<-glm(fmli, family=family, data=data,...)
    p.step<-length(coef(a.glm))
    ps<-sum(residuals(a.glm,type="pearson")^2)
    if (pval == "chi-square"){
      pv<-2*min(pchisq(ps,n-p.step),pchisq(ps,n-p.step,lower.tail = FALSE))
    }
    if (pval == "bootstrap") {
      pv<-boot_pval(fmli,family=family, data=data, B=B,...)
    }
    if(pv>alpha)
    {
      mod.opt<-fmli
    }
    else
    {
      varc<-var.cat[var.cat %in% var.mod]
      for(i in 1:length(varc))
      {
        fmli<-update.formula(mod.opt, paste0("~. -", varc[i]))
        a.glm<-glm(fmli, family=family, data=data,...)
        p.step<-length(coef(a.glm))
        ps<-sum(residuals(a.glm,type="pearson")^2)
        if (pval == "chi-square"){
          pv<-2*min(pchisq(ps,n-p.step),pchisq(ps,n-p.step,lower.tail = FALSE))
        }
        if (pval == "bootstrap") {
          pv<-boot_pval(fmli,family=family, data=data, B=B,...)
        }
        if(pv>alpha)
        {
          mod.opt<-fmli
        }
      }
    }
  if(mod.opt!=mod.step[[length(mod.step)]])
  {
    mod.step[[length(mod.step)+1]]<-mod.opt
  }
  }

  for(j in 1:length(mod.step))
    mod.step[[j]]<-deparse1(mod.step[[j]])
  return(list(models=mod.step, model.opt=mod.step[[length(mod.step)]]))
}
