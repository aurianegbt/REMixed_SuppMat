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

library(REMixed)
library(dplyr)
library(ggplot2)
load("data/applicationFiles/Rinchai/correspondance.RData")

lab_map <- setNames(
  correspondance_function$yobs,
  paste0("alpha_1", sub("^yG", "", correspondance_function$yG))
)

plotIC(cv.remix,dismin=FALSE)+geom_line(color="black")+xlab("Penalty parameter")+ylab("BICc")

ggsave(filename="outputs/figures/finalFigures/IC_app.eps",device="eps",height=5,width=6)

plotCalibration(cv.remix)+theme(legend.position="right")+
  scale_colour_discrete(
    name   = "Function",
    labels = function(v) unname(lab_map[v])
  )+xlab("Penalty parameter")+ylab("Biomarkers contribution")

ggsave(filename="outputs/figures/finalFigures/Calibration.eps",device="eps",height=6,width=12)

# #
# selection = data.frame()
# for(i in 1:length(cv.remix$res)){
#   selection <- rbind(selection,unname(extract(cv.remix,n = i)$finalRes$alpha!=0))
# }
# colnames(selection) <- names(extract(cv.remix,n = i)$finalRes$alpha)
#
# selection <- (t(sapply(selection,function(x){ifelse(x,"X","")})))
# colnames(selection) <- round(cv.remix$lambda,digits=2)
# View(selection)
#
# parameterValue = data.frame()
# for(i in 1:length(cv.remix$res)){
#   aux = extract(cv.remix,n = i)
#   parameterValue <- rbind(parameterValue,aux$finalRes$param[c("fM2_pop","theta_pop")])
#   colnames(parameterValue) = c("fM2_pop","theta_pop")
# }
# parameterValue <- cbind(lambda=cv.remix$lambda,parameterValue)
