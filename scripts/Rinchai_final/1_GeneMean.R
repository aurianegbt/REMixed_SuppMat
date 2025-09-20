library(dplyr)
library(ggplot2)

patientToKeep = data.frame(PZB = paste0("PZB",c(4:16,19,25)),
                           CVX=paste0("CVX",c(4:16,19,25)))
timeInDays = setNames(c(0:9,14,21:30,35),
                      c(paste0("C1D0",1:9),"C1D10","C1D14",paste0("C2D0",1:9),"C2D10","C2D14"))

sample_info = read.csv("data/applicationFiles/Rinchai/raw/Data_2025-05-07_10-33-35.csv")
sample_info$sample.id <- unname(sapply(sample_info$Specimen.ID,FUN=function(x){stringr::str_split(x,"_")[[1]][3]}))
sample_info$Participant.ID <- unname(sapply(sample_info$Specimen.ID,FUN=function(x){stringr::str_split(x,"_")[[1]][1]}))
sample_info$time.ID <- unname(sapply(sample_info$Specimen.ID,FUN=function(x){stringr::str_split(x,"_")[[1]][2]}))

annotations = read.csv("data/applicationFiles/Rinchai/annotations.txt")
load("data/applicationFiles/Rinchai/data/data_norm.RData")
load("data/applicationFiles/modules_Chaussabel_Gen3_gmt.RData")

# GENE MEAN REPRESENTATIVE ------------------------------------------------
data <- merge(annotations,log_counts,by.x="ensembl_gene_id",by.y="ID")
data <- data[data$hgnc_symbol!="",-1]
colnames(data)[1] <- "ID"


genefunction = unique(modules_Chaussabel_Gen3_gmt$genesets.descriptions)
representativeGenes = data.frame()
# genefunction = "Interferon"

for(gf in genefunction){
  print(paste0("starting function :",gf))

  genes = unique(unlist(modules_Chaussabel_Gen3_gmt$genesets[modules_Chaussabel_Gen3_gmt$genesets.descriptions==gf]))

  genesDF = data[data$ID %in% genes,]

  data_genes = data.frame()
  for(g in 1:nrow(genesDF)){
    aux = data.frame(obs=as.numeric(genesDF[g,-1]),
                     yobs=genesDF$ID[g],
                     sampleID = colnames(genesDF[,-1]),
                     id=sapply(colnames(genesDF[,-1]),function(sid){sample_info[sample_info$sample.id==paste0("0",sid),"Participant.ID"]},USE.NAMES = FALSE),
                     time=sapply(colnames(genesDF[,-1]),function(sid){unname(timeInDays[stringr::str_split(sample_info[sample_info$sample.id==paste0("0",sid),"Specimen.ID"],"_")[[1]][2]])},USE.NAMES = FALSE))

    aux <- merge(aux,patientToKeep,by.x="id",by.y="CVX")
    aux <- aux %>% select(-id) %>% rename(id=PZB)

    data_genes = rbind(data_genes,aux)
  }
  save(data_genes,file=paste0("data/applicationFiles/Rinchai/data_gene_function/genes_",stringr::str_replace_all(stringr::str_replace_all(gf,"/","_")," ","_"),".RData"))

  representativeGenes = rbind(representativeGenes,cbind(data_genes %>%
                                                          group_by(id,time) %>%
                                                          summarise(obs=mean(obs)),
                                                        yobs=gf))

}

save(representativeGenes,file="data/applicationFiles/Rinchai/data/genes.RData")

#### NORMALIZATION

load("data/applicationFiles/Rinchai/data/genes.RData")

correspondance_function = data.frame(yG=paste0("yG",1:length(unique(representativeGenes$yobs))),yobs=unique(representativeGenes$yobs))
correspondance_patientID = data.frame(patient=1:length(unique(data$id)),id=unique(data$id))

save(correspondance_function,correspondance_patientID,file="data/applicationFiles/Rinchai/correspondance.RData")

for(gf in unique(representativeGenes$yobs)){
  aux = representativeGenes %>% filter(yobs==gf) %>% filter(time==0)

  mGi0 = mean(aux$obs)
  sGi0 = sd(aux$obs)

  representativeGenes[representativeGenes$yobs==gf,"obs"] <- (representativeGenes[representativeGenes$yobs==gf,"obs"]-mGi0)/sGi0

  data_genes = representativeGenes %>% filter(yobs==gf)

  ggplot(data_genes,aes(x=time,y=obs,group=interaction(yobs, id),color=id))+
    geom_line(alpha=0.005)+
    stat_summary(fun = mean, geom = "point", aes(group = id), size = 1) +  # Moyenne par temps
    stat_summary(fun = mean, geom = "line", aes(group = id), linewidth = 0.7) +
    # stat_summary(fun.data = function(y) {
    #   data.frame(ymin = mean(y) - sd(y), ymax = mean(y) + sd(y), y = mean(y))
    # }, geom = "errorbar", aes(group = yobs),linewidth=0.7) +
    scale_color_viridis_d(option = "plasma", direction = -1)+
    theme (text = element_text(size=28),
           axis.text.y = element_text(hjust = 1, size=20),
           axis.title= element_text (size=20),
           axis.text.x = element_text(angle = 45, hjust = 1, size=22),
           legend.text = element_text(size=12),
           legend.title =element_text(size=20),
           strip.text = element_text(size=20, face="bold"),
           plot.background = element_rect(fill='transparent', color=NA),
           legend.background = element_rect(fill='transparent'),
           legend.box.background = element_rect(fill='transparent',color=NA))+
    ylab ("Gene expressions")+
    xlab("time (in days)") +
    ggtitle(stringr::str_replace_all(gf,"/"," "))

  ggsave(file=paste0("outputs/figures/applicationPlot/Rinchai/genefunctionTraj/meanGene_by_id_",stringr::str_replace_all(stringr::str_replace_all(gf,"/","_")," ","_"),"_NORM.png"),
         height = 2000,width=2000,unit="px")
}

load("data/applicationFiles/Rinchai/correspondance.RData")
data <- read.csv("data/applicationFiles/Rinchai/data/igGDataRinchai.txt")


data_genes <- representativeGenes %>%
  left_join(correspondance_function, by = c("yobs"="yobs")) %>%
  mutate(yobs = yG) %>%
  select(-yG)

data_genes <- data_genes %>%
  left_join(correspondance_patientID, by = c("id" = "id"))  %>%
  as.data.frame() %>%
  select(patient,time,obs,yobs) %>%
  rename(id=patient)

data <- data %>%
  left_join(correspondance_patientID, by = c("id" = "id"))  %>%
  select(-id) %>%
  rename(id=patient) %>%
  rename(yobs=yObs) %>%
  select(id,time,obs,yobs)

data_all <- rbind(data,data_genes)

write.csv(data_all,file="data/applicationFiles/Rinchai/data/data_NORM.txt",quote = F,row.names = F)


