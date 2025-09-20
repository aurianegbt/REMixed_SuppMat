## COnduct  selection on Interferon
library(dplyr)
library(ggplot2)
library(readxl)

patientToKeep = data.frame(PZB = paste0("PZB",c(4:16,19,25)),
                           CVX=paste0("CVX",c(4:16,19,25)))
timeInDays = setNames(c(0:9,14,21:30,35),
                      c(paste0("C1D0",1:9),"C1D10","C1D14",paste0("C2D0",1:9),"C2D10","C2D14"))

sample_info = read.csv("data/applicationFiles/Rinchai/sample_info.txt")

# Antibodies data ---------------------------------------------------------

sciadv_abp9961_supplemental_file_s1 <- read_excel("data/applicationFiles/Rinchai/raw/sciadv.abp9961_supplemental_file_s1.xls")
colnames(sciadv_abp9961_supplemental_file_s1) <- paste0(colnames(sciadv_abp9961_supplemental_file_s1),c(c("",""),paste0(" ",sciadv_abp9961_supplemental_file_s1[1,][-c(1,2)])))


data <-  sciadv_abp9961_supplemental_file_s1[-1,] %>%
  filter(Patient %in% patientToKeep$PZB) %>%
  select(Patient,`Day post-vaccination`,`SARS-Cov2 S1...13 Total IgG`) %>%
  mutate(yObs="yAB") %>%
  rename(time="Day post-vaccination") %>%
  rename(obs = "SARS-Cov2 S1...13 Total IgG") %>%
  rename(id = "Patient") %>%
  mutate(obs = as.numeric(obs)) %>%
  mutate(time = as.numeric(time)) %>%
  mutate(obs = log10(obs))

rm(sciadv_abp9961_supplemental_file_s1)

ggplot(data,aes(x=time,y=obs,group=id,color=id))+geom_line()+ylab("log10(Ab)")+
  scale_color_viridis_d(name="Patient ID",option = "plasma", direction = -1)

ggsave(file="outputs/figures/applicationPlot/Rinchai/antibodyTrajectories.png",height=1500,width=2000,units = "px")

head(data)
write.csv(data,file="data/applicationFiles/data/igGDataRinchai.txt",quote = F,row.names = F)



# SAMPLE INFO -------------------------------------------------------------


sample_info = read.csv("data/applicationFiles/Rinchai/raw/Data_2025-05-07_10-33-35.csv")
sample_info <- cbind(sample.id = sapply(stringr::str_split(sample_info$Specimen.ID,"_"),FUN=function(x){x[3]}),sample_info)

write.csv(sample_info,file="data/applicationFiles/Rinchai/sample_info.txt")
# which_to_delete = colnames(count_data)[which(count_data$samples$lib.size < 80000)]
# View(sample_info[sample_info$sample.id %in% paste0("0",which_to_delete),])


# RAW COUNT MATRIX --------------------------------------------------------

covax_prime = read.csv("data/applicationFiles/Rinchai/raw/GSE190001_COVAX_raw_count_PRIME.txt",sep="\t")
colnames(covax_prime) <- stringr::str_remove_all(colnames(covax_prime),"X")
covax_boost = read.csv("data/applicationFiles/Rinchai/raw/GSE190001_COVAX_raw_count_BOOST.txt",sep="\t")
colnames(covax_boost) <- stringr::str_remove_all(colnames(covax_boost),"X")

# all(covax_prime[,"ID"]==covax_boost[,"ID"])
covax = cbind(covax_prime,covax_boost[,-1])

count_data <- edgeR::DGEList(counts = covax[, -1])
count_data <- edgeR::calcNormFactors(count_data, method = "TMMwsp")

log_counts <- mapply(function(v, lib_size, norm_factor) {
  log2((v + 0.5) / (lib_size * norm_factor + 1) * 1e6)
},
as.data.frame(count_data$counts),
count_data$samples$lib.size,
count_data$samples$norm.factors,
SIMPLIFY = TRUE
)

sample_sums <- colSums(2^log_counts)

barplot(sample_sums,
        las = 2,
        main = "Somme des valeurs des gènes par sample",
        xlab = "Samples",
        ylab = "Somme des valeurs",
        cex.names = 0.6,
        col=as.factor(sample_info$DOSE))

# summary(sample_sums[sample_info$DOSE=="BOOST"])
library(ggplot2)

df_to_plot <- data.frame(gene=names(sample_sums),sums=unname(sample_sums))

ggplot(df_to_plot,aes(x=gene,y=sums))+geom_col() +
  xlab("Samples") + ylab("Sums") +
  theme(axis.text.x=element_blank())

ggsave(file="outputs/figures/applicationPlot/Rinchai/barplots_norm_genes.png")

log_counts <- cbind(ID=covax[,1],log_counts)

save(log_counts,file="data/applicationFiles/Rinchai/data/data_norm.RData")

# GENE SYMBOL -------------------------------------------------------------
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_ids <- data_norm[,"ID"]

annotations <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = rownames(log_counts),
                     mart = ensembl)



write.csv(annotations,file="data/applicationFiles/Rinchai/annotations.txt",quote=F,row.names=F)

