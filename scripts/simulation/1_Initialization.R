# empty the memory
rm(list=ls())
dir <- function(x){if(!dir.exists(x)){dir.create(x)}}


# Load Functions and softwares
library(lixoftConnectors)
library(foreach)
library(dplyr)
initializeLixoftConnectors()
library(REMixed)

# Arguments of batch
Nb_rep <- 1
Nbr_genes <- 20

temporaryDirectory <- paste0("tmp",Nbr_genes,"_",Nb_rep,"/")
dir(temporaryDirectory)

# Create project
sim <- read.csv(paste0("data/simulationFiles/simulationVaccinology_",Nbr_genes,"/data.txt"))
sim <- sim %>% filter(rep==Nb_rep) %>% select(ID,TIME,obs,obsid)

write.csv(sim,file=paste0(temporaryDirectory,"data.txt"),quote = F,row.names = F)

file.copy(from=paste0("data/simulationFiles/simulationVaccinology_",Nbr_genes,"/model.txt"),
          to = paste0(temporaryDirectory,"model.txt"))

newProject(modelFile = paste0(temporaryDirectory,"model.txt"),
           data = list(dataFile=paste0(temporaryDirectory,"data.txt"),headerTypes=c("id","time","observation","obsid")))

pset <- list(nbsmoothingiterations=100,
             nbexploratoryiterations=250,
             simulatedannealing=T, smoothingautostop=T,exploratoryautostop=T)
pop.set <- lixoftConnectors::getPopulationParameterEstimationSettings()
pop.set <- modifyList(pop.set, pset[intersect(names(pset), names(pop.set))])

setPopulationParameterEstimationSettings(pop.set)

setMapping(append(list(list(data="yAB",prediction="tAB",model="yAB")),
                  lapply(1:Nbr_genes,FUN=function(n){list(data=paste0("yG",n),prediction=paste0("G",n),model=paste0("yG",n))})))

setErrorModel(yAB="constant")

setPopulationParameterInformation(delta_S_pop=list(initialValue=0.01),
                                  delta_AB_pop = list(initialValue=0.01))

for(i in 1:Nbr_genes){
  cmd=paste0("setIndividualParameterDistribution(alpha_0",i,"='normal')")
  eval(parse(text = cmd))
  cmd = paste0("setIndividualParameterVariability(alpha_0",i,"=FALSE)")
  eval(parse(text = cmd))
  cmd=paste0("setIndividualParameterDistribution(alpha_1",i,"='normal')")
  eval(parse(text = cmd))
  cmd = paste0("setErrorModel(yG",i,"='constant')")
  eval(parse(text = cmd))
  cmd = paste0("setIndividualParameterVariability(alpha_1",i,"=FALSE)")
  eval(parse(text = cmd))
}

saveProject(paste0(temporaryDirectory,"project.mlxtran"))

# MLE estimation for non kept genes
MLE = list()
for(k in 1:Nbr_genes){
  Y = sim %>% filter(obsid==paste0("yG",k))
  Y <- Y$obs
  MLE[[k]] <- MASS::fitdistr(Y, "normal")
}
rm("Y")

# group creation
Nbr_grp = Nbr_genes/2
genes = 1:Nbr_genes
sizes = rep(Nbr_genes%/% Nbr_grp,Nbr_grp)
if((Nbr_genes %% Nbr_grp)!=0)
  sizes[1:(Nbr_genes %% Nbr_grp)] <- sizes[1:(Nbr_genes %% Nbr_grp)] + 1
shuffled = sample(genes)

split_list <- split(shuffled, rep(1:Nbr_grp, sizes))

res = lapply(split_list,FUN=function(genes){

  cat("--------------------------------------------------\n")
  cat("Starting with genes",paste0(genes,collapse=", "),"\n")

  NEWtemporaryDirectory <- paste0(temporaryDirectory,"/tmp_",paste0(genes,collapse="_"),"/")

  if(!dir.exists(NEWtemporaryDirectory)){
    dir(NEWtemporaryDirectory)
  }

  loadProject(paste0(temporaryDirectory,"project.mlxtran"))
  saveProject(paste0(NEWtemporaryDirectory,"project.mlxtran"))

  for(i in 1:Nbr_genes){
    if(!(i %in% genes)){
      cmd = paste0("setPopulationParameterInformation(alpha_1",i,"_pop=list(initialValue=0,method='FIXED'))")
      eval(parse(text = cmd))
      cmd = paste0("setPopulationParameterInformation(alpha_0",i,"_pop=list(initialValue=",MLE[[i]]$estimate[["mean"]],",method='FIXED'))")
      eval(parse(text = cmd))
      cmd = paste0("setPopulationParameterInformation(ayG",i,"=list(initialValue=",MLE[[i]]$estimate[["sd"]],",method='FIXED'))")
      eval(parse(text = cmd))
    }else{
      cmd = paste0("setPopulationParameterInformation(alpha_1",i,"_pop=list(initialValue=1))")
      eval(parse(text = cmd))
    }
  }

  saveProject(paste0(NEWtemporaryDirectory,"project.mlxtran"))

  runPopulationParameterEstimation()
  runLogLikelihoodEstimation()

  saveProject(paste0(NEWtemporaryDirectory,"project.mlxtran"))

  cat("Estimated Parameters :\n")
  print(getEstimatedPopulationParameters()[c("delta_S_pop","phi_S_pop","delta_AB_pop",paste0("alpha_1",genes,"_pop"),"omega_delta_S","omega_phi_S","omega_delta_AB")])
  cat("Estimated logLikelihood :\n")
  print(getEstimatedLogLikelihood()$importanceSampling)


  return(list(genes=genes,
              project=paste0(NEWtemporaryDirectory,"project.mlxtran"),
              parameters = getEstimatedPopulationParameters(),
              LL=getEstimatedLogLikelihood()$importanceSampling))
})

LL = sapply(res,FUN=function(r){r$LL[["OFV"]]})

keep = which.min(LL)

loadProject(res[[keep]]$project)
saveProject(paste0(temporaryDirectory,"project.mlxtran"))

runConditionalDistributionSampling()

saveProject(paste0(temporaryDirectory,"project.mlxtran"))

save(MLE,res,file=paste0("outputs/initialization/",Nb_rep,Nbr_genes,"_gene.RData"))

for(r in res){
  unlink(paste0(temporaryDirectory,"/tmp_",paste0(r$genes,collapse="_")),force=TRUE,recursive = TRUE)
}

