library(stringr)
library(ggplot2)
library(lixoftConnectors)
initializeLixoftConnectors(software = "simulx",force=T)


## Arguments

Nb_genes = 20  # 20 -> 5 ; 50 -> 10
Nb_rep = 200
Nb_genes_signif = ifelse(Nb_genes==20,5,10)
Nbr_ind = 25

## Create Simulation Model

dir <- function(d,unlink=FALSE,recursive=TRUE){if(!dir.exists(d)){dir.create(d)}else if(unlink){unlink(d,recursive=recursive);dir.create(d)}}
set.seed(25041710)

Normal = readLines("data/Normal.txt")

line_input = str_sub(Normal[10],end=-2)
line_output = str_sub(Normal[22],end=-2)

input = character()
output = character()
def_genes = c()
obs_genes = c()

for(i in 1:Nb_genes){
  input = paste0(input,",alpha_0",i,",alpha_1",i,",sigma_G",i)
  output = paste0(output,",G",i)

  def_genes = c(def_genes,paste0("G",i,"= alpha_0",i," + alpha_1",i," * S"))
  obs_genes = c(obs_genes,paste0("yG",i," = {distribution=normal,prediction=G",i,",errorModel=constant(sigma_G",i,")}"))
}
input <- paste0(input,"}")
output <- paste0(output,"}")

NewModel <- c(Normal[1:9],paste0(line_input,input),Normal[11:20],def_genes,Normal[20:21],paste0(line_output,output),Normal[23:26],obs_genes)

dir(paste0("data/simulationSetup/simulX/simulationVaccinology_",Nb_genes,"/"))
writeLines(NewModel,paste0("data/simulationSetup/simulX/simulationVaccinology_",Nb_genes,"/model_for_simulx.txt"))

# Create Simulation

newProject(modelFile = paste0("data/simulationSetup/simulX/simulationVaccinology_",Nb_genes,"/model_for_simulx.txt"))

outputElement = data.frame(delta_S_pop=0.23,
                           phi_S_pop=611,
                           delta_AB_pop=0.025,
                           sigma_AB=0.11,
                           omega_phi_S = 0.60,
                           omega_delta_AB = 0.3,
                           omega_delta_S = 0.1)

alpha_1 = c(round(sample(c(-1,1),size = Nb_genes_signif,replace=TRUE)*rgamma(Nb_genes_signif,2,3),digits=2),rep(0,Nb_genes-Nb_genes_signif))
alpha_0 = round(rnorm(Nb_genes),digits=2)
sigma = round(abs(rnorm(Nb_genes,mean=0.3,sd=0.05)),digits=2)

dir(paste0("outputs/figures/simulationPlot/simulationVaccinology_",Nb_genes,"/"))
summary(sigma)
summary(alpha_0)
boxplot(sigma)
print(alpha_1[1:Nb_genes_signif])


for(i in 1:Nb_genes){
  auxDF = data.frame(sigma[i],alpha_1[i],alpha_0[i])
  colnames(auxDF) <- c(paste0("sigma_G",i),
                       paste0("alpha_1",i),
                       paste0("alpha_0",i))

  outputElement <- cbind(outputElement,auxDF)
}

definePopulationElement(name="PopParameters",
                        element = outputElement)

defineOutputElement(name=paste0("yAB"),element = list(output="yAB",data=data.frame(time=c(0,7,21,123,180,300))))
for(i in 1:Nb_genes){
  defineOutputElement(name=paste0("yG",i),element = list(output=paste0("yG",i),data=data.frame(time=c(1:21))))
}

setGroupSize("simulationGroup1", Nbr_ind)
setGroupElement("simulationGroup1", elements = c("PopParameters","yAB",paste0("yG",1:Nb_genes)))

setNbReplicates(Nb_rep)

runSimulation()
saveProject(paste0("data/simulationSetup/simulX/simulationVaccinology_",Nb_genes,"/simulation.smlx"))

sim <- getSimulationResults(rep = 1)

exportProject(settings = list(targetSoftware = "monolix",filesNextToProject = T,dataFilePath = "data.txt", modelFileName = "model.txt"),force = T)
saveProject(paste0("data/simulationSetup/simulX/simulationVaccinology_",Nb_genes,"/monolix.mlxtran"))
unlink(paste0("data/simulationSetup/simulX/simulationVaccinology_",Nb_genes,"/monolix.mlxtran"),force=TRUE)
unlink(paste0("data/simulationSetup/simulX/simulationVaccinology_",Nb_genes,"/monolix.mlxproperties"),force=TRUE)


ggplot(sim$res$yG1,aes(x=time,group=id,color=as.factor(id),y=yG1))+geom_line(alpha=0.6)+geom_point(alpha=0.6)+theme(legend.position="none") + scale_fill_brewer()
ggsave(filename = paste0("outputs/figures/simulationPlot/simulationVaccinology_",Nb_genes,"/yG1.png"))

ggplot(sim$res$yG2,aes(x=time,group=id,color=as.factor(id),y=yG2))+geom_line(alpha=0.6)+geom_point(alpha=0.6)+theme(legend.position="none") + scale_fill_brewer()
ggsave(filename = paste0("outputs/figures/simulationPlot/simulationVaccinology_",Nb_genes,"/yG2.png"))

ggplot(sim$res$yG3,aes(x=time,group=id,color=as.factor(id),y=yG3))+geom_line(alpha=0.6)+geom_point(alpha=0.6)+theme(legend.position="none") + scale_fill_brewer()
ggsave(filename = paste0("outputs/figures/simulationPlot/simulationVaccinology_",Nb_genes,"/yG3.png"))

ggplot(sim$res$yG4,aes(x=time,group=id,color=as.factor(id),y=yG4))+geom_line(alpha=0.6)+geom_point(alpha=0.6)+theme(legend.position="none") + scale_fill_brewer()
ggsave(filename = paste0("outputs/figures/simulationPlot/simulationVaccinology_",Nb_genes,"/yG4.png"))


ggplot(sim$res$yG5,aes(x=time,group=id,color=as.factor(id),y=yG5))+geom_line(alpha=0.6)+geom_point(alpha=0.6)+theme(legend.position="none") + scale_fill_brewer()
ggsave(filename = paste0("outputs/figures/simulationPlot/simulationVaccinology_",Nb_genes,"/yG5.png"))


ggplot(sim$res$yG11,aes(x=time,group=id,color=as.factor(id),y=yG11))+geom_line(alpha=0.6)+geom_point(alpha=0.6)+theme(legend.position="none") + scale_fill_brewer()
ggsave(filename = paste0("outputs/figures/simulationPlot/simulationVaccinology_",Nb_genes,"/yG11.png"))

ggplot(sim$res$yAB,aes(x=time,group=id,color=as.factor(id),y=yAB))+geom_line(alpha=0.6)+geom_point(alpha=0.6)+theme(legend.position="none") + scale_fill_brewer()
ggsave(filename = paste0("outputs/figures/simulationPlot/simulationVaccinology_",Nb_genes,"/yAB.png"))

dir.create(paste0("data/simulationFiles/simulationVaccinology_",Nb_genes))
file.copy(from=paste0("data/simulationSetup/simulX/simulationVaccinology_",Nb_genes,"/data.txt"),
          to = paste0("data/simulationFiles/simulationVaccinology_",Nb_genes,"/data.txt"),overwrite = TRUE)
file.copy(from=paste0("data/simulationSetup/simulX/simulationVaccinology_",Nb_genes,"/model.txt"),
          to = paste0("data/simulationFiles/simulationVaccinology_",Nb_genes,"/model.txt"),overwrite = TRUE)

