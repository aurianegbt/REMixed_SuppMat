## This function computes the mean estimate, relative bias, empirical and estimated sd, 95% coverage on the parameter from a scenario for the non-significant biomarkers when wrongly selected.
## For biomarkers are not selected, standard error isn't computed. So coverage are computed on the set on replicates where it is selected. (usefull for noise table)

compute_noise_table <-
  function(parameterValue,genes,trueValueDF,Nbr_genes,method,replicates=1:200){

    if(Nbr_genes==20){
      signif=5
    }else{
      signif=10
    }

    parameterValue = parameterValue[parameterValue$Nbr_genes==Nbr_genes &
                                      parameterValue$method==method ,]

    genes_aux <- genes[genes$Nbr_genes==Nbr_genes &
                         genes$method==method &
                         genes$genes >signif,]

    trueValue = setNames(trueValueDF[trueValueDF$Nbr_genes==Nbr_genes,"trueValue"],trueValueDF[trueValueDF$Nbr_genes==Nbr_genes,"parameter"])

    genes_selected0 <- genes_aux[genes_aux$selected,c("model","parameter")]
    genes_selected <- rbind(genes_selected0,genes_selected0 %>% mutate(parameter=str_replace_all(genes_selected0$parameter,"alpha_1","alpha_0")))
    genes_selected <- rbind(genes_selected,genes_selected0 %>% mutate(parameter=str_replace_all(genes_selected0$parameter,"alpha_1","sigma_G")))

    parameterValue_aux <- merge(genes_selected,parameterValue,by=c("model","parameter"))

    mean_Estimation = tapply(parameterValue_aux$values,parameterValue_aux$parameter,mean)
    relative_bias = as.numeric((mean_Estimation - trueValue[names(mean_Estimation)])/trueValue[names(mean_Estimation)])*100
    empirical_sd = tapply(parameterValue_aux$values,parameterValue_aux$parameter,sd)
    estimated_sd = tapply(parameterValue_aux$se,parameterValue_aux$parameter,FUN=function(x){mean(x,na.rm=TRUE)})
    count = sapply(names(mean_Estimation),FUN=function(p){nrow(parameterValue_aux[parameterValue_aux$parameter==p,])})

    coverHess = setNames(rep(0,length(mean_Estimation)),names(mean_Estimation))
    nCount = setNames(rep(0,length(mean_Estimation)),names(mean_Estimation))
    for(p in names(mean_Estimation)){
      for (j in unique(parameterValue_aux[parameterValue_aux$parameter==p,"model"])){
        sdhat = parameterValue_aux[parameterValue_aux$model==j & parameterValue_aux$parameter==p,"se"]
        if(!is.na(sdhat)){
          nCount[p] <- nCount[p] + 1
          lw = parameterValue_aux[parameterValue_aux$parameter==p & parameterValue_aux$model==j,"values"] - 1.96*sdhat
          up =parameterValue_aux[parameterValue_aux$parameter==p & parameterValue_aux$model==j,"values"] + 1.96*sdhat
          bool = as.numeric(trueValue[p] > lw & trueValue[p] < up)
          coverHess[p] = coverHess[p] + bool
        }
      }
    }
    coverHess <- coverHess/nCount

    res =data.frame(parameter = names(mean_Estimation),
                    trueValue = unname(as.numeric(trueValue[names(mean_Estimation)])),
                    mean_Estimation= unname(mean_Estimation),
                    relative_bias = unname(relative_bias),
                    empirical_sd = unname(empirical_sd[names(mean_Estimation)]),
                    estimated_sd = unname(estimated_sd[names(mean_Estimation)]),
                    count = unname(count[names(mean_Estimation)]),
                    coverHess = unname(coverHess[names(mean_Estimation)])*100,
                    tab="noise")

    return(res)
  }
