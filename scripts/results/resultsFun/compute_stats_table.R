## This function computes the mean estimate, relative bias, empirical and estimated sd, 95% coverage on the parameter from a scenario
## For biomarkers are not selected, standard error isn't computed. So coverage are computed on the set on replicates where it is selected. (usefull for noise table)

compute_stats_table <-
  function(parameterValue,trueValueDF,Nbr_genes,method,replicates=1:200){


    parameterValue = parameterValue[parameterValue$Nbr_genes==Nbr_genes &
                                      parameterValue$method==method,]

    # Estimation accuracy for all parameters throughout all replicates
    trueValue = setNames(trueValueDF[trueValueDF$Nbr_genes==Nbr_genes,"trueValue"],trueValueDF[trueValueDF$Nbr_genes==Nbr_genes,"parameter"])

    mean_Estimation = tapply(parameterValue$values,parameterValue$parameter,mean)
    relative_bias = as.numeric((mean_Estimation - trueValue[names(mean_Estimation)])/trueValue[names(mean_Estimation)])*100
    empirical_sd = tapply(parameterValue$values,parameterValue$parameter,sd)
    estimated_sd = tapply(parameterValue$se,parameterValue$parameter,FUN=function(x){mean(x,na.rm=TRUE)})
    countNA = tapply(parameterValue$se,parameterValue$parameter,FUN=function(x){sum(is.na(x))})

    coverHess = setNames(rep(0,length(mean_Estimation)),names(mean_Estimation))
    nCount = setNames(rep(length(replicates),length(mean_Estimation)),names(mean_Estimation))
    for (j in replicates){
      sdhat = setNames(parameterValue[parameterValue$model==j,"se"],parameterValue[parameterValue$model==j,"parameter"])[names(mean_Estimation)]
      for(par in names(sdhat)){
        if(is.na(sdhat[par])){
          nCount[par] <- nCount[par] -1
        }else{
          lw = parameterValue[parameterValue$parameter==par & parameterValue$model==j,"values"] - 1.96*sdhat[par]
          up =parameterValue[parameterValue$parameter==par & parameterValue$model==j,"values"] + 1.96*sdhat[par]
          bool = as.numeric(trueValue[par] > lw & trueValue[par] < up)
          coverHess[par] = coverHess[par] + bool
        }
      }
    }
    coverHess <- coverHess/nCount

    res = data.frame(parameter = names(mean_Estimation),
                     trueValue = unname(as.numeric(trueValue[names(mean_Estimation)])),
                     mean_Estimation= unname(mean_Estimation),
                     relative_bias = unname(relative_bias),
                     empirical_sd = unname(empirical_sd[names(mean_Estimation)]),
                     estimated_sd = unname(estimated_sd[names(mean_Estimation)]),
                     count = nb_replicates-unname(countNA[names(mean_Estimation)]),
                     coverHess = unname(coverHess[names(mean_Estimation)])*100,
                     tab="all")

  return(res)
}
