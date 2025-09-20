library(REMixed)

load(paste0("~/Travail/REMix_PLAFRIM/outputs/applicationResults/Rinchai/runDelay/res_1.RData"))

cv.remix <- resCV

for(arr in 2:30){
  if(file.exists(paste0("~/Travail/REMix_PLAFRIM/outputs/applicationResults/Rinchai/runDelay/res_",arr,".RData"))){
    load(paste0("~/Travail/REMix_PLAFRIM/outputs/applicationResults/Rinchai/runDelay/res_",arr,".RData"))


    if(length(resCV$lambda)!=0){
      cv.remix$lambda <- c(cv.remix$lambda,resCV$lambda)
      cv.remix$LL <- c(cv.remix$LL,resCV$LL)
      cv.remix$LL.pen <- c(cv.remix$LL.pen,resCV$LL.pen)

      cv.remix$res <- append(cv.remix$res,resCV[["res"]])
      cv.remix$outputs <- append(cv.remix$outputs,resCV[["outputs"]])

    }
    rm("resCV")
  }
}

class(cv.remix) <- "cvRemix"

Nbr_genes=34


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


final = retrieveBest(cv.remix)
final$info$project <- '/beegfs/agabaut/tmpRinchai_delay_s7_7/project.mlxtran'


finalTEST = computeFinalTest(final,dynFUN,y=c(S=5,AB=0.5),ObsModel.transfo)

save(finalTEST,file="outputs/applicationResults/Rinchai/finalSEL_test.RData")
