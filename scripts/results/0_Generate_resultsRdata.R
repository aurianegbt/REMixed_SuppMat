## SimulationResults not available on this repository github
## Code cannot be run
## Only for visibility/reproductibility

library(dplyr)
library(stringr)

nb_replicates <- 200
Nbr_genes = c(20,50)
p.max = 0.05

resultsStats <- data.frame()
parameterValue <- data.frame()
genes <- data.frame()

for(pjt in project){
  for(nb in Nbr_genes){
    signif = ifelse(nb==20,5,10)
    for(arr in 1:nb_replicates){
      tryCatch({
        if(nb==20){
          load(paste0("outputs/simulationResults/simulationVaccinology_20/resREMIX_",arr,".RData"))
        }else{
          load(paste0("outputs/simulationResults/simulationVaccinology_50/resREMIX_final_",arr,".RData"))
        }


        resultsStats <- rbind(resultsStats,data.frame(model=arr,
                                                      method ="REMix",
                                                      Nbr_genes = nb,
                                                      lambda.min=res$info$lambda,
                                                      time = res$finalRes$time[["elapsed"]],
                                                      FP = sum(res$finalRes$saemBeforeTest$param[paste0("alpha_1",(signif+1):nb,"_pop")]!=0),
                                                      FN = sum(res$finalRes$saemBeforeTest$param[paste0("alpha_1",1:signif,"_pop")]==0),
                                                      TN = sum(res$finalRes$saemBeforeTest$param[paste0("alpha_1",(signif+1):nb,"_pop")]==0),
                                                      TP = sum(res$finalRes$saemBeforeTest$param[paste0("alpha_1",1:signif,"_pop")]!=0)))

        resultsStats <- rbind(resultsStats,data.frame(model=arr,
                                                      method = "REMixWald",
                                                      Nbr_genes = nb,
                                                      lambda.min=res$info$lambda,
                                                      time = res$finalRes$time[["elapsed"]],
                                                      FP = sum(res$finalRes$alpha[(signif+1):nb]!=0),
                                                      FN = sum(res$finalRes$alpha[1:signif]==0),
                                                      TN = sum(res$finalRes$alpha[(signif+1):nb]==0),
                                                      TP = sum(res$finalRes$alpha[1:signif]!=0)))

        parameterValue <- rbind(parameterValue,merge(data.frame(model=arr,
                                                                method="REMix",
                                                                Nbr_genes = nb,
                                                                parameter =names(res$finalRes$saemBeforeTest$param),
                                                                values = res$finalRes$saemBeforeTest$param),res$finalRes$saemBeforeTest$standardError$stochasticApproximation[,-3],by="parameter",all.x = TRUE))

        parameterValue <- rbind(parameterValue,merge(data.frame(model=arr,
                                                                method="REMixWald",
                                                                Nbr_genes = nb,
                                                                parameter =names(res$finalRes$param),
                                                                values = res$finalRes$param),res$finalRes$standardError$stochasticApproximation[,-3],by="parameter",all.x = TRUE))

        genes <- rbind(genes,
                       data.frame(model=arr,
                                  method="REMix",
                                  Nbr_genes = nb,
                                  parameter = paste0("alpha_1",1:nb),
                                  genes = 1:nb,
                                  selected = unname(res$finalRes$saemBeforeTest$param[paste0("alpha_1",1:nb,"_pop")]!=0)))

        genes <- rbind(genes,
                       data.frame(model=arr,
                                  method="REMixWald",
                                  Nbr_genes = nb,
                                  parameter = paste0("alpha_1",1:nb),
                                  genes = 1:nb,
                                  selected = unname(res$finalRes$alpha!=0)))

        rm("res")
        if(nb==20){
          rm("resCV")
        }
      },error=function(e){print(arr)})
    }
  }
}


replicates = unique(genes$model)
precision = abs(floor(log10(1/length(replicates))))


resultsStats = resultsStats %>% mutate(FPR=FP/(FP+TN),
                                       FNR = FN/(FN+TP),
                                       FDR = FP/(TP+FP),
                                       F1_score = 2*TP/(2*TP+FP+FN))

parameterValue[str_detect(parameterValue$parameter,"alpha_"),"parameter"] <- str_remove_all(parameterValue[str_detect(parameterValue$parameter,"alpha_"),"parameter"],"_pop")

parameterValue[str_detect(parameterValue$parameter,"ay"),"parameter"] <- str_replace_all(parameterValue[str_detect(parameterValue$parameter,"ay"),"parameter"],"ay","sigma_")

# TO EDIT FOR NEW PROJECT

trueValueDF = data.frame()

trueValue = read.csv(paste0("data/simulationSetup/simulX/simulationVaccinology_20/simulation/Simulation/populationParameters.txt"))
trueValue <- trueValue[trueValue$rep==1,-c(1,69)]
trueValueDF<- rbind(trueValueDF,cbind(data.frame(parameter=names(trueValue),trueValue=as.numeric(unname(trueValue))),Nbr_genes=20))

if(50 %in% Nbr_genes){
  trueValue = read.csv(paste0("data/simulationSetup/simulX/simulationVaccinology_50/simulation/Simulation/populationParameters.txt"))
  trueValue <- trueValue[trueValue$rep==1,-c(1,159)]
  trueValueDF<- rbind(trueValueDF,cbind(data.frame(parameter=names(trueValue),trueValue=as.numeric(unname(trueValue))),Nbr_genes=50))
}



parameterValue <- merge(parameterValue,trueValueDF,by=c("parameter","Nbr_genes","Proj"),all.x=TRUE)

save(resultsStats,parameterValue,trueValueDF,genes,file=paste0("outputs/simulationResults/simulationResults.RData"))



# For Initialisation  -----------------------------------------------------
results = data.frame()
for(nb in Nbr_genes){
  signif = ifelse(nb==20,5,10)
  for(rep in 1:nb_replicates){
    load(paste0("outputs/simulationResults/initialization/",rep,nb,"_gene.RData"))

    for(i in 1:length(res)){
      results <- rbind(results,
                       data.frame(
                         Nbr_genes=nb,
                         rep=rep,
                         LL=-2*res[[i]]$LL["OFV"],
                         BIC = res[[i]]$L["BIC"],
                         BICc = res[[i]]$L["BICc"],
                         genes_1 = res[[i]]$genes[1],
                         genes_2=res[[i]]$genes[2],
                         alpha_gene_1 = res[[i]]$parameters[[paste0("alpha_1",res[[i]]$genes[1],"_pop")]],
                         alpga_gene_2 = res[[i]]$parameters[[paste0("alpha_1",res[[i]]$genes[2],"_pop")]],
                         delta_AB_pop = res[[i]]$parameters[["delta_AB_pop"]],
                         phi_S_pop = res[[i]]$parameters[["phi_S_pop"]],
                         delta_S_pop = res[[i]]$parameters[["delta_S_pop"]],
                         row.names = 1))
    }


    # View(results %>% arrange(desc(LL)))
  }
}


results = results %>% mutate(replicates = paste0("replicate ",rep))

results$replicates <- factor(results$replicates,levels=paste0("replicate ",1:200))

max_genes = results %>% group_by(Nbr_genes,replicates,Proj) %>%
  summarize(LL=max(LL)) %>% merge(results[,c("replicates","LL","genes_1","genes_2","delta_S_pop","phi_S_pop","delta_AB_pop")],by=c("replicates","LL"))

save(results,max_genes,file=paste0("outputs/simulationResults/InitialisationIllustration.RData"))
