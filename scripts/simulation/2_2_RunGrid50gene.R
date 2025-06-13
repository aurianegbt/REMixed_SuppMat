## This scrip only to compute for a lambda rank arr in the grid in order to parallelize later on every point of the grid
# empty the memory
rm(list=ls())
dir <- function(x){if(!dir.exists(x)){dir.create(x)}}


# Load Functions and softwares
library(lixoftConnectors)
library(foreach)
library(dplyr)
initializeLixoftConnectors()
library(REMix)

# Arguments of batch
Nb_rep <- 1
arr_lam <- 1
Nbr_genes <- 50

if(!file.exists(paste0("outputs/simulationResults/simulationVaccinology_",Nbr_genes,"/resREMIX_0_",Nb_rep,".RData"))){
  stop("error, first lambda estimation does not exists.")
}else{
  load(paste0("outputs/simulationResults/simulationVaccinology_",Nbr_genes,"/resREMIX_0_",Nb_rep,".RData"))
}

lambda_max <- resCV$lambda/0.001
lambda <- lambda_max * 0.001**((15:1)/15)[arr_lam+1]

temporaryDirectory <- paste0("tmp",Nbr_genes,"_",Nb_rep,"_",arr_lam,"/")
dir(temporaryDirectory)
pathToResults = paste0("outputs/simulationResults/simulationVaccinology_",Nbr_genes,"/resREMIX_",arr_lam,"_",Nb_rep,".RData")

project = paste0(temporaryDirectory,"project.mlxtran")
dynFUN <- function(t,y,parms){
  #phi_S,phi_L,delta_Ab,delta_S=0.23,delta_L=0.000316
  if(!setequal(names(y),c("AB","S"))){
    bonus = setdiff(names(y),c("AB","S"))
    malus = setdiff(c("AB","S"),names(y))
    if(length(bonus)!=0 & length(malus) !=0){
      stop(paste0("Missing initial condition for ",malus," and ",bonus," isn't in the model."))
    }else if(length(bonus)!=0){
      stop(paste0(bonus," isn't dynamic of the model."))
    }else if(length(minus) !=0){
      stop(paste0("Missing initial condition for ",malus))
    }
  }
  if(!setequal(names(parms),c("phi_S","delta_AB","delta_S"))){
    stop(paste0("Missing parmeters ",setdiff(c("phi_S","delta_AB","delta_S",),names(parms))," and ",setdiff(names(parms),c("phi_S","delta_AB","delta_S"))," isn't in the model."))
  }
  parms <- c(parms,S0=y[["S"]])
  out <- deSolve::ode(y["AB"],t,
                      function(t,y,parms){
                        with(as.list(c(y, parms)), {
                          dAB = phi_S * S0 * exp(-delta_S*t) - delta_AB * AB
                          list(c(dAB))
                        })
                      },parms)
  out <- cbind(out,S=y[["S"]]*exp(-parms[["delta_S"]]*out[,"time"]))
  return(out)
}
alpha0 = setNames(paste0("alpha_0",1:Nbr_genes),paste0("yG",1:Nbr_genes))

resCV = cv.remix(project = project,
                 dynFUN = dynFUN,
                 y = c(S=5,AB=1000),
                 ObsModel.transfo = list(S=list(AB=log10),linkS="yAB",R=rep(list(S=function(x){x}),Nbr_genes),linkR = paste0("yG",1:Nbr_genes)),
                 alpha = list(alpha0=alpha0,alpha1=setNames(paste0("alpha_1",1:Nbr_genes),paste0("yG",1:Nbr_genes))),
                 selfInit = TRUE,
                 eps1=10**(-2),
                 eps2=1,
                 lambda_max = lambda,
                 alambda=1,
                 nlambda=1)

save(resCV,file=pathToResults)
