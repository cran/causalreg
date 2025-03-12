cgam_step<-function(formula, family, data, alpha=0.05,pval.approx=TRUE, B=100,...)
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
    a<-glm(formula, family, data, ...)
  } else {
    a<-gam(formula, family, data, ...)
  }

  if(names(a$coefficients)[1]=="(Intercept)"){
    intercept<-TRUE
  } else{
    intercept<-FALSE
  }

  mod.step<-vector(mode="list")
  pv.step<-NULL

  if(intercept)
  {
    mod.step[[1]]<-as.formula(paste0(dip.name,"~","1"))
    pv.step[1]<-0

  }
  else
  {
    mod.step[[1]]<-as.formula(paste0(dip.name,"~","-1"))
    pv.step[1]<-0
  }
    var.namesyes<-var.names #for potential inclusion
    var.namesmod<-NULL #already in the model

  continue<- TRUE
  j<-0
  while(continue & length(var.namesyes)>0)
  {
    j<-j+1
    #check which variable to add
    pvals<-NULL
    for(i in 1: length(var.namesyes))
    {
      fmli<-update.formula(mod.step[[j]], paste0("~. +", var.namesyes[i]))
      if (vrs[2]=="."){
        a.gam<-glm(fmli, family, data,...)
        edf<-length(a.gam$coefficients)
      } else{
        a.gam<-gam(fmli, family, data, ...)
        sa.gam<-summary(a.gam)
        edf<-length(sa.gam$p.coeff)+sum(sa.gam$edf)
      }
      ps<- sum(residuals(a.gam,type="pearson")^2)
      if (pval.approx){
        pv<-2*min(pchisq(ps,n-edf),pchisq(ps,n-edf,lower.tail = FALSE))
      } else {
        pv<-boot_pval(fmli,family,data,B,...)
      }
      pvals<-c(pvals,pv)
    }
    i.min<-which.max(pvals)
    mod.min<-update.formula(mod.step[[j]], paste0("~. +", var.namesyes[i.min]))
    #rearrange all variables names
    a<-attr(terms.formula(mod.min),"term.labels")
    mod.min<-reformulate(sort(a),response=dip.name,intercept=intercept)

    if (vrs[2]=="."){
      a.gam<-glm(mod.min, family, data, ...)
      edf<-length(a.gam$coefficients)
    } else{
      a.gam<-gam(mod.min, family, data, ...)
      sa.gam<-summary(a.gam)
      edf<-length(sa.gam$p.coeff)+sum(sa.gam$edf)
    }
    ps<- sum(residuals(a.gam,type="pearson")^2)
    if (pval.approx){
      pv.min<-2*min(pchisq(ps,n-edf),pchisq(ps,n-edf,lower.tail = FALSE))
    } else {
      pv.min<-boot_pval(fmli,family,data,B,...)
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
  if (vrs[2]=="."){
    current.gam<-glm(model.current, family, data, ...)
  } else{
    current.gam<-gam(model.current, family, data, ...)
  }
  bic.current<-BIC(current.gam)
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
        if (vrs[2]=="."){
          bics[i]<-BIC(glm(fmli, family, data, ...))
        } else{
          bics[i]<-BIC(gam(fmli, family, data, ...))
        }
      }
      i.min<-which.min(bics)
      mod.min<-update.formula(model.current, paste0("~. -", var.namesmod[i.min]))
      if (vrs[2]=="."){
        bic.min<-BIC(glm(mod.min, family, data, ...))
      } else{
        bic.min<-BIC(gam(mod.min, family, data, ...))
      }
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

  for(j in 1:length(mod.step))
    mod.step[[j]]<-format(as.formula(deparse(mod.step[[j]]),env=NULL))

  mod.opt<-mod.step[[length(mod.step)]]
  return(list(models=mod.step, model.opt=mod.opt))
}
