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
arr <- as.numeric(slurmR::Slurm_env(x='SLURM_ARRAY_TASK_ID')) # rank lambda

Nbr_genes = 34
seed = 7
grp = 7

cat("===========================================================================\n
- - - - - - - - - - - - - - Work Information - - - - - - - - - - - - - - -\n")
cat("Working on Rinchai data \n")
cat("Run n°",arr,".\n")
cat("\n===========================================================================\n")

pathToResults ="outputs/applicationResults/Rinchai/runDelay"
dir(pathToResults)
pathToResults <- paste0(pathToResults,"/res_",arr,".RData")
temporaryDirectory <- paste0("/beegfs/agabaut/tmpRinchai_delay_s",seed,"_",grp,"/")
if(!dir.exists(temporaryDirectory)){
  stop("Error, Init hasn't been done.")
}
if(!file.exists(paste0(temporaryDirectory,"project.mlxtran"))){
  stop("Error, can't find mlxtran file.")
}

project = paste0(temporaryDirectory,"project.mlxtran")


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

loadProject(paste0(temporaryDirectory,"project.mlxtran"))
# if(!getLaunchedTasks()$conditionalDistributionSampling){
#   runConditionalDistributionSampling()
#   saveProject()
# }
#
# data <- readMLX(paste0(temporaryDirectory,"project.mlxtran"),
#                 ObsModel.transfo,
#                 alpha)
#
# lambda_max = REMixed:::lambda.max(dynFUN,y =c(S=5,AB=0.5), data=data)
#
if(!file.exists(paste0("/beegfs/agabaut/tmpRinchai_delay_s",seed,"_",grp,"_",arr,"/project/remix/Build_1.mlxtran"))){
  dir(paste0("/beegfs/agabaut/tmpRinchai_delay_s",seed,"_",grp,"_",arr))
  saveProject(paste0("/beegfs/agabaut/tmpRinchai_delay_s",seed,"_",grp,"_",arr,"/project.mlxtran"))
  project <- paste0("/beegfs/agabaut/tmpRinchai_delay_s",seed,"_",grp,"_",arr,"/project.mlxtran")
}else{
  project <-paste0("/beegfs/agabaut/tmpRinchai_delay_s",seed,"_",grp,"_",arr,"/project/remix/Build_1.mlxtran")
}
#
# print(lambda_max)
#
# stop()

lambda_max = 2000

lambda <- lambda_max * 0.001**((30:1)/30)[arr]

resCV = cv.remix(project = project,
                 dynFUN = dynFUN,
                 y = c(S=5,AB=0.5),
                 ObsModel.transfo = ObsModel.transfo,
                 alpha = alpha,
                 selfInit = TRUE,
                 eps1=10**(-2),
                 eps2=0.5,
                 lambda_max = lambda,
                 alambda=1,
                 nlambda=1,
                 unlinkBuildProject = FALSE)

save(resCV,file=pathToResults)
