# empty the memory
rm(list=ls())
dir <- function(x){if(!dir.exists(x)){dir.create(x)}}

# Load Functions and softwares
suppressWarnings({
  suppressMessages(library(lixoftConnectors))
  library(foreach)
  suppressMessages(initializeLixoftConnectors(software = "monolix",
                                              path = "/cm/shared/dev/modules/generic/apps/tools/monolix/2023R1bis/",
                                              force = TRUE))})

library(REMixed)
library(dplyr)

# Arguments of batch
arr <- as.numeric(slurmR::Slurm_env(x='SLURM_ARRAY_TASK_ID')) # nbr Init
args = commandArgs(trailingOnly=TRUE)
s <- as.numeric(args[1])
seed <- c(17102501,25041012,16021506,30110607,12082607,07022503,20012305,180645,19820556,51472189)[s]


Nbr_genes = 34

cat("===========================================================================\n
- - - - - - - - - - - - - - Work Information - - - - - - - - - - - - - - -\n")
cat("Working on Rinchai data - Initialization\n")
cat("\n===========================================================================\n")

pathToResults ="outputs/applicationResults/Rinchai/initializationDelay"
dir(pathToResults)
pathToResults = paste0(pathToResults,"/seed_",s)
dir(pathToResults)
pathToResults <- paste0(pathToResults,"/init_",arr,".RData")
temporaryDirectory <- paste0("/beegfs/agabaut/tmpRinchai_delay_s",s,"_",arr,"/")
dir(temporaryDirectory)

file.copy(from = "data/applicationFiles/Rinchai/model/model_for_REMix_delay.txt",
          to=paste0(temporaryDirectory,"model.txt"),
          overwrite = TRUE)

file.copy(from = "data/applicationFiles/Rinchai/data/data_NORM.txt",
          to=paste0(temporaryDirectory,"data.txt"),
          overwrite = TRUE)

newProject(data=list(dataFile=paste0(temporaryDirectory,"data.txt"),
                     headerTypes=c("id","time","observation","obsid")),
           modelFile=paste0(temporaryDirectory,"model.txt"))

saveProject(paste0(temporaryDirectory,"project.mlxtran"))

setMapping(append(list(list(data="yAB",prediction="tAb",model="yAB")),
                  lapply(1:Nbr_genes,FUN=function(n){list(data=paste0("yG",n),prediction=paste0("G",n),model=paste0("yG",n))})))

setErrorModel(yAB="constant")

for(i in 1:Nbr_genes){
  cmd=paste0("setIndividualParameterDistribution(alpha_0",i,"='normal')")
  eval(parse(text = cmd))
  cmd = paste0("setIndividualParameterVariability(alpha_0",i,"=FALSE)")
  eval(parse(text = cmd))
  cmd=paste0("setIndividualParameterDistribution(alpha_1",i,"='normal')")
  eval(parse(text = cmd))
  cmd = paste0("setIndividualParameterVariability(alpha_1",i,"=FALSE)")
  eval(parse(text = cmd))
  cmd = paste0("setErrorModel(yG",i,"='constant')")
  eval(parse(text = cmd))
}

setPopulationParameterInformation(delta_V_pop=list(initialValue=0.46,method="FIXED"),
                                  td_pop=list(initialValue=1,method="FIXED"),
                                  fM1_pop = list(initialValue=1,method="FIXED"),
                                  fM2_pop = list(initialValue=4.5),
                                  delta_Ab_pop=list(initialValue=0.03,method="FIXED"),
                                  k_pop=list(initialValue=0.25,method="FIXED"),
                                  delta_S_pop=list(initialValue=0.01,method="FIXED"))

setIndividualParameterVariability(delta_V=FALSE,
                                  td = FALSE,
                                  fM1 = FALSE,
                                  fM2 = TRUE,
                                  theta = TRUE,
                                  delta_S=FALSE,
                                  delta_Ab=FALSE,
                                  k=FALSE)

saveProject(paste0(temporaryDirectory,"project.mlxtran"))


set.seed(seed)
genes = 1:Nbr_genes
Nbr_grp = Nbr_genes %/% 2
sizes = rep(Nbr_genes%/% Nbr_grp,Nbr_grp)
if((Nbr_genes %% Nbr_grp)!=0)
  sizes[1:(Nbr_genes %% Nbr_grp)] <- sizes[1:(Nbr_genes %% Nbr_grp)] + 1
shuffled = sample(genes)

split_list <- split(shuffled, rep(1:Nbr_grp, sizes))

genes <- split_list[[arr]]


dynFUN <- function(t,y,parms){
  y <- c(S=y[["S"]],AB=y[["AB"]],F1=0,F2=0,F3=0,F4=0,F5=0,F6=0)
  out <- deSolve::ode(y,t,
                      function(t,y,parms){
                        with(as.list(c(y, parms)), {

                          t_0 = 0
                          t_1 = 21

                          if(t < t_1){
                            C = fM1
                            ttilde = t_0
                          }else{
                            C = fM2
                            ttilde = t_1
                          }

                          dS = C*exp(-delta_V*(t-ttilde))-k*S
                          dF1 = k*S-k*F1
                          dF2 = k*F1-k*F2
                          dF3 = k*F2-k*F3
                          dF4 = k*F3-k*F4
                          dF5 = k*F4-k*F5
                          dF6 = k*F5-delta_S*F6
                          dAB = theta * F6 - delta_Ab * AB

                          list(c(dS,dAB,dF1,dF2,dF3,dF4,dF5,dF6))
                        })
                      },parms)

  out <- out[,c("time","S","AB")]

  td = parms[["td"]]

  out[,"S"] <- sapply(out[,"time"],function(t){
    if(t<td){
      return(0)
    }else if(t < 21+td){
      C=parms[["fM1"]]
      return(C*exp(-parms[["delta_V"]]*(t-td)))
    }else{
      C=parms[["fM2"]]
      return(C*exp(-parms[["delta_V"]]*(t-21-td)))
    }
  })
  class(out) <- c("deSolve","matrix")
  return(out)
}

ObsModel.transfo = list(S=list(AB=log10),linkS="yAB",R=rep(list(S=function(x){x}),Nbr_genes),linkR = paste0("yG",1:Nbr_genes))

alpha = list(alpha0=setNames(paste0("alpha_0",1:Nbr_genes),paste0("yG",1:Nbr_genes)),alpha1=setNames(paste0("alpha_1",1:Nbr_genes),paste0("yG",1:Nbr_genes)))

sim <- read.csv("data/applicationFiles/Rinchai/data/data_NORM.txt")

MLE = list()
for(k in 1:Nbr_genes){
  Y = sim %>% filter(yobs==paste0("yG",k))
  Y <- Y$obs
  MLE[[k]] <- MASS::fitdistr(Y, "normal")
}
rm("Y")

for(i in 1:Nbr_genes){
  cmd=paste0("setIndividualParameterDistribution(alpha_0",i,"='normal')")
  eval(parse(text = cmd))
  cmd = paste0("setIndividualParameterVariability(alpha_0",i,"=FALSE)")
  eval(parse(text = cmd))
  cmd = paste0("setPopulationParameterInformation(alpha_0",i,"_pop=list(initialValue=MLE[[i]]$estimate['mean']))")
  eval(parse(text = cmd))
  cmd=paste0("setIndividualParameterDistribution(alpha_1",i,"='normal')")
  eval(parse(text = cmd))
  cmd = paste0("setIndividualParameterVariability(alpha_1",i,"=FALSE)")
  eval(parse(text = cmd))
  cmd = paste0("setErrorModel(yG",i,"='constant')")
  eval(parse(text = cmd))
  cmd = paste0("setPopulationParameterInformation(",getContinuousObservationModel()$parameters[[paste0('yG',i)]],"=list(initialValue=MLE[[i]]$estimate['sd']))")
  eval(parse(text = cmd))
}

for(i in 1:Nbr_genes){
  if(!(i %in% genes)){
    cmd = paste0("setPopulationParameterInformation(alpha_1",i,"_pop=list(initialValue=0,method='FIXED'))")
    eval(parse(text = cmd))
    cmd = paste0("setPopulationParameterInformation(alpha_0",i,"_pop=list(method='FIXED'))")
    eval(parse(text = cmd))
    cmd = paste0("setPopulationParameterInformation(ayG",i,"=list(method='FIXED'))")
    eval(parse(text = cmd))
  }else{
    cmd = paste0("setPopulationParameterInformation(alpha_1",i,"_pop=list(initialValue=1))")
    eval(parse(text = cmd))
    cmd = paste0("setPopulationParameterInformation(alpha_0",i,"_pop=list(method='MLE'))")
    eval(parse(text = cmd))
    cmd = paste0("setPopulationParameterInformation(ayG",i,"=list(method='MLE'))")
    eval(parse(text = cmd))

  }
}
saveProject(paste0(temporaryDirectory,"project.mlxtran"))

runPopulationParameterEstimation()
runLogLikelihoodEstimation()

saveProject(paste0(temporaryDirectory,"project.mlxtran"))

cat("Estimated Parameters :\n")
print(getEstimatedPopulationParameters())
cat("Estimated logLikelihood :\n")
print(getEstimatedLogLikelihood()$importanceSampling)


res <- list(genes=genes,
            project=temporaryDirectory,
            parameters = getEstimatedPopulationParameters(),
            LL=getEstimatedLogLikelihood()$importanceSampling)


save(res,file=pathToResults)
