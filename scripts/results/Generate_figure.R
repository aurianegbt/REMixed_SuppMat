library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggpattern)
library(ggh4x)
dir <- function(d){if(!dir.exists(d)){dir.create(d)}}
# Save format
PNG <- TRUE
JPEG <- TRUE
EPS <- TRUE

saveFigures <- function(pathToSave,plot=last_plot(),height=NA,width=NA,PNG=TRUE,JPEG=TRUE,EPS=TRUE,dpi=600){
  if(PNG){
    ggsave(paste0(pathToSave,".png"),plot=plot,
           height = height,width=width,dpi=dpi, units = "px",device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(pathToSave,".jpeg"),plot=plot,
           height = height,width=width,dpi=dpi, units = "px",device=grDevices::jpeg)
  }
  if(EPS){
    ggsave(paste0(pathToSave,".eps"),plot=plot,
           height = height/dpi,width=width/dpi,dpi=dpi, units = "in",device=cairo_ps)
  }
}

pathToResults = paste0("outputs/figures/finalFigures/")
dir(pathToResults)
bb = theme(plot.background = element_rect(linewidth=0.7,color="black"))
Nbr_genes = c(20,50)

# Graphs Figure 1 -  Variable Selection Perfomance ------------------------
load(paste0("outputs/simulationResults/simulationResults.RData"))

replicates = unique(genes$model)
nb_replicates = length(replicates)
method=c("REMix","REMixWald")

# Selection frequency
df_to_plot <- genes[genes$selected,]
df_to_plot <- df_to_plot %>% mutate(Nbr_genes = paste0(df_to_plot$Nbr_genes," genes"))

valueToPrint <- data.frame()
for(m in method){
  for(n in Nbr_genes){
    signif = ifelse(n==20,5,10)
    aux = df_to_plot[df_to_plot$Nbr_genes==paste0(n," genes") & df_to_plot$method==m,]
    valueToPrint <- rbind(valueToPrint,data.frame(
      Nbr_genes = paste0(n," genes"),
      method=m,
      signif=signif,
      text = paste0(round(mean(table(aux[aux$genes %in% (signif+1):n,"parameter"])/nb_replicates)*100,digits=1),"%"), # fréquence de sélection moyenne des gènes le long des réplicats pour cette méthode et ce nombre de gène
      coory=mean(table(aux[aux$genes %in% (signif+1):n,"parameter"])/nb_replicates)*200,
      coorx=as.character(n)
    ))
  }
}


errorStatsParCov <-  resultsStats[,c("Nbr_genes","method","model","FP","FN","TP","TN","FPR","FNR","F1_score")]

df = data.frame(Nbr_genes = rep(Nbr_genes,length(method)),
                Method=rep(method,each=length(Nbr_genes)),
                FPR_mean = sapply(split(errorStatsParCov$FPR,
                                        errorStatsParCov[,c("Nbr_genes","method")]),FUN= mean),

                FPR_median = sapply(split(errorStatsParCov$FPR,
                                          errorStatsParCov[,c("Nbr_genes","method")]),FUN = median),

                FPR_sd = sapply(split(errorStatsParCov$FPR,
                                      errorStatsParCov[,c("Nbr_genes","method")]),FUN = sd),

                FPR_q975 = sapply(split(errorStatsParCov$FPR,
                                        errorStatsParCov[,c("Nbr_genes","method")]),
                                  FUN = function(x){return(quantile(x,0.975))}),

                FPR_q25 = sapply(split(errorStatsParCov$FPR,
                                       errorStatsParCov[,c("Nbr_genes","method")]),
                                 FUN = function(x){quantile(x,0.025)}),

                FNR_mean = sapply(split(errorStatsParCov$FNR,
                                        errorStatsParCov[,c("Nbr_genes","method")]),FUN = mean),

                FNR_median = sapply(split(errorStatsParCov$FNR,
                                          errorStatsParCov[,c("Nbr_genes","method")]),FUN = median),

                FNR_sd = sapply(split(errorStatsParCov$FNR,
                                      errorStatsParCov[,c("Nbr_genes","method")]),FUN = sd),

                FNR_q975 =  sapply(split(errorStatsParCov$FNR,
                                         errorStatsParCov[,c("Nbr_genes","method")]),
                                   FUN = function(x){quantile(x,0.975)}),

                FNR_q25 = sapply(split(errorStatsParCov$FNR,
                                       errorStatsParCov[,c("Nbr_genes","method")]),
                                 FUN = function(x){quantile(x,0.025)}),

                F1_score_mean = sapply(split(errorStatsParCov$F1_score,
                                             errorStatsParCov[,c("Nbr_genes","method")]),FUN = mean),

                F1_score_median = sapply(split(errorStatsParCov$F1_score,
                                               errorStatsParCov[,c("Nbr_genes","method")]),FUN = median),

                F1_score_sd = sapply(split(errorStatsParCov$F1_score,
                                           errorStatsParCov[,c("Nbr_genes","method")]),FUN = sd),

                F1_score_q975 =  sapply(split(errorStatsParCov$F1_score,
                                              errorStatsParCov[,c("Nbr_genes","method")]),
                                        FUN = function(x){quantile(x,0.975)}),

                F1_score_q25 = sapply(split(errorStatsParCov$F1_score,
                                            errorStatsParCov[,c("Nbr_genes","method")]),
                                      FUN = function(x){quantile(x,0.025)}))



percent <- function(x, digits = 1, format = "f", ...) {      # Create user-defined function
  paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}

CB <- function(df,char){
  lb = df[,paste0(char,"_q25")]
  ub = df[,paste0(char,"_q975")]
  return(paste0("[",sapply(lb,FUN=function(x){percent(max(0,x))}),
                ";",sapply(ub,FUN=function(x){percent(min(1,x))}),"]"))}

df <- cbind(df, FPR_CB = CB(df,"FPR"),FNR_CB = CB(df,"FNR"),F1_score_CB=CB(df,"F1_score"))

df[,c(3,4,8,9,13,14)] <- apply(df[,c(3,4,8,9,13,14)] ,2,FUN=function(x){round(x,3)})

valueScores <- data.frame()
for(m in method){
  for(n in Nbr_genes){
    aux <- df[df$Method==m & df$Nbr_genes==n,]
    valueScores <- rbind(valueScores,data.frame(
      Nbr_genes = paste0(n," genes"),
      method=m,
      x = ifelse(n==20,"6","11"),
      text = paste0("FPR : ",aux$FPR_mean*100,"% ",aux$FPR_CB,"\n",
                    "FNR : ",aux$FNR_mean*100,"% ",aux$FNR_CB,"\n",
                    "F1-score : ",aux$F1_score_mean*100,"% ",aux$F1_score_CB)
    ))
  }
}

fillscale = c()
if(20 %in% Nbr_genes){
  fillscale <- c(fillscale,c(rep("darkred",5),rep("grey",15)))
}
if(50 %in% Nbr_genes){
  fillscale <- c(fillscale,c(rep("darkred",10),rep("grey",40)))
}


selectPropREMix <-
  ggplot(df_to_plot %>% filter(method=="REMix"),aes(x=factor(genes)))+geom_bar(fill=fillscale)+
  ggh4x::facet_grid2(Nbr_genes~.,scales="free_x",independent="x")+
  ylab("Proportion") +
  xlab("Genes") +
  scale_y_continuous(limits = c(0,200),
                     breaks = seq(0,200,200/4),
                     labels = scales::percent(seq(0,1,0.25)))+
  geom_text(data=valueToPrint%>% filter(method=="REMix"),aes(y=coory,x=coorx,label=text),hjust=1,vjust=-1,size=5)+
  geom_text(data=valueScores%>% filter(method=="REMix"),aes(label=text,x=x),y=200,hjust=0,vjust=1)+
  geom_segment(data=valueToPrint%>% filter(method=="REMix"),aes(x = signif+1,xend=coorx,y=coory,yend=coory),color="grey16",linetype="dashed") +
  theme(axis.title=element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=7),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA))

saveFigures(paste0(pathToResults,"selectionProportionREMix"),
            height = 2500,width = 6000,PNG = PNG,JPEG=JPEG,EPS=EPS)


selectPropREMixWald <-
  ggplot(df_to_plot %>% filter(method=="REMixWald"),aes(x=factor(genes)))+geom_bar(fill=fillscale)+
  ggh4x::facet_grid2(Nbr_genes~.,scales="free_x",independent="x")+
  ylab("Proportion") +
  xlab("Genes") +
  scale_y_continuous(limits = c(0,200),
                     breaks = seq(0,200,200/4),
                     labels = scales::percent(seq(0,1,0.25)))+
  geom_text(data=valueToPrint%>% filter(method=="REMixWald"),aes(y=coory,x=coorx,label=text),hjust=1,vjust=-1,size=5)+
  geom_text(data=valueScores%>% filter(method=="REMixWald"),aes(label=text,x=x),y=200,hjust=0,vjust=1)+
  geom_segment(data=valueToPrint%>% filter(method=="REMixWald"),aes(x = signif+1,xend=coorx,y=coory,yend=coory),color="grey16",linetype="dashed") +
  theme(axis.title=element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=7),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA))

saveFigures(paste0(pathToResults,"selectionProportionREMixWald"),
            height = 2500,width = 6000,PNG = PNG,JPEG=JPEG,EPS=EPS)

# Exact Model

value = data.frame()
for(m in method){
  for(nb in Nbr_genes){
    aux = resultsStats[resultsStats$Nbr_genes==nb  & resultsStats$method==m & resultsStats$FN==0,]

    value <- rbind(value,data.frame(Nbr_genes=paste0(nb," genes"),
                                    method=m,
                                    Type=c("Exact",rep("Overselection",length(unique(aux$FP))-1)),
                                    Value=sapply(sort(unique(aux$FP)),
                                                 FUN=function(t){
                                                   nrow(aux[aux$FP<=t,])
                                                 })/length(replicates)*100,
                                    FP = sort(unique(aux$FP)),
                                    ExactNB = sapply(sort(unique(aux$FP)),
                                                     FUN=function(t){
                                                       nrow(aux[aux$FP==t,])
                                                     })/length(replicates)*100))
  }
}


compStatsREMix <-
  ggplot(mapping=aes(x=Nbr_genes, y=Value,pattern=Type)) +
  geom_col(data=value[value$Type=="Overselection" & value$method=="REMix",], aes(group = Nbr_genes), width = 0.7, position = position_dodge(width = 0.7),color="transparent")  +
  geom_col(data=value[value$method=="REMix",],aes(group=Nbr_genes,fill=factor(FP),y=ExactNB),width=0.7)+
  geom_bar_pattern(data=value[value$Type=="Exact" & value$method=="REMix",], position="dodge",width=0.7,color="transparent",pattern_fill="transparent",fill=NA,pattern_density=0.05,pattern_spacing=0.025,stat='identity') +
  labs(x="", y = "Proportion") +
  scale_fill_manual(values=setNames(colorRampPalette(c("navajowhite", "darkred"))(length(unique(value$FP))),sort(unique(value$FP))),name="False Positives")+
  scale_y_continuous(limits = c(0,100),
                     breaks = seq(0,100,100/4),
                     labels = scales::percent(seq(0,1,0.25)))+
    # scale_fill_gradient(low="navajowhite",high="darkred",name="False Positives")+
  theme(axis.title=element_text(size=16),
        legend.position = "bottom",
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=7),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box = "vertical"
        )+
    guides(fill = guide_legend(override.aes = list(colour = NA), keywidth = 0.5, keyheight = 1),
           pattern = guide_legend(override.aes = list(fill = "transparent", color = "transparent")))

saveFigures(paste0(pathToResults,"ComparisonStatsREMix"),
            height = 3000,width = 3000,PNG = PNG,JPEG=JPEG,EPS=EPS)


compStatsREMixWald <-
  ggplot(mapping = aes(x = Nbr_genes, y = Value, pattern = Type)) +
  geom_col(data = value[value$Type == "Overselection" & value$method == "REMixWald", ], aes(group = Nbr_genes),width = 0.7,position = position_dodge(width = 0.7), color = "transparent") +
  geom_col(data = value[value$method == "REMixWald", ], aes(group = Nbr_genes, fill = factor(FP), y = ExactNB),width = 0.7) +
  geom_bar_pattern(data = value[value$Type == "Exact" & value$method == "REMixWald", ],
                   position = "dodge",
                   width = 0.7,
                   color = "transparent",
                   pattern_fill = "transparent",
                   fill = NA,
                   pattern_density = 0.05,
                   pattern_spacing = 0.025,
                   stat = 'identity') +
  labs(x = "", y = "Proportion") +
  scale_fill_manual(values = setNames(colorRampPalette(c("navajowhite", "darkred"))(length(unique(value$FP))),
                                      sort(unique(value$FP))),name = "False Positives") +
  theme(axis.title=element_text(size=16),
        legend.position = "bottom",
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=7),
        strip.text = element_text(size = 11),
        legend.box = "vertical",
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill = "transparent", color = NA))+
  guides(fill = guide_legend(override.aes = list(colour = NA), keywidth = 0.5, keyheight = 1),
         pattern = guide_legend(override.aes = list(fill = "transparent", color = "transparent")))+
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 100 / 4),
                     labels = scales::percent(seq(0, 1, 0.25)))


saveFigures(paste0(pathToResults,"ComparisonStatsREMixWald"),
            height = 3000,width = 3000,PNG = PNG,JPEG=JPEG,EPS=EPS)

## Arrange

ggarrange(compStatsREMix+bb,
          selectPropREMix+bb,
          nrow=1,widths=c(1.2,1.8),labels=c("A1","A2"))

saveFigures(paste0(pathToResults,"Figure1R"),
            height = 3000,width = 7000,PNG = PNG,JPEG=JPEG,EPS=EPS)

ggarrange(compStatsREMixWald + bb,
          selectPropREMixWald + bb,
          nrow = 1, widths = c(1.2, 2.8), labels = c("B1", "B2"))

saveFigures(paste0(pathToResults,"Figure1RW"),
            height = 3000,width = 7000,PNG = PNG,JPEG=JPEG,EPS=EPS)

# Graphs Figure 2 - Initialisation  ---------------------------------------
load(paste0("outputs/simulationResults/InitialisationIllustration.RData"))

results_plot = results %>% group_by(Nbr_genes,genes_1,genes_2) %>% summarize(LL=mean(LL))

# initLL20 <-
#   ggplot(results_plot%>% filter(Nbr_genes==20),aes(x=genes_1,y=genes_2,color=LL))+
#   geom_point()+
#   # facet_grid(Nbr_genes~.)+
#   scale_color_gradient2(high="darkred",mid="tomato",low="navajowhite",midpoint=sum(range((results%>% filter(Nbr_genes==20))$LL))/2) +
#   geom_segment(x=5.5,xend=5.5,y=5.5,yend=20,color="black",linetype = "dashed") +
#   geom_segment(x=5.5,xend=20,y=5.5,yend=5.5,color="black",linetype = "dashed") +
#   geom_segment(x=5.5,xend=5.5,y=0,yend=5.5,color="darkred",linetype = "dashed") +
#   geom_segment(x=0,xend=5.5,y=5.5,yend=5.5,color="darkred",linetype = "dashed") +
#   geom_point(data = max_genes%>% filter(Nbr_genes==20), aes(color=LL,x = genes_1, y = genes_2, shape = "Log-likelihood\nmaximum"),size = 4,inherit.aes = FALSE) +
#   scale_shape_manual(name = "",values = c("Log-likelihood\nmaximum" = 10)) +
#   guides(shape = guide_legend()) +
#   labs(x="2nd gene",y="1st gene", color = "Log-Likelihood")+
#   theme(axis.title=element_text(size=12),
#         # legend.position="bottom",
#         strip.text = element_text(size = 11),
#         plot.background = element_rect(fill='transparent', color=NA),
#         legend.key = element_rect(fill = "transparent", color = NA),
#         legend.background = element_rect(fill = "transparent", color = NA))

initLL <-
  ggplot(results_plot%>% filter(Nbr_genes==50),aes(x=genes_1,y=genes_2,color=LL))+
  geom_point(size=0.7)+
  # facet_grid(Nbr_genes~.)+
  scale_color_gradient2(high="darkred",mid="tomato",low="navajowhite",midpoint=sum(range((results%>% filter(Nbr_genes==50))$LL))/2) +
  geom_segment(x=10.5,xend=10.5,y=10.5,yend=50,color="black",linetype = "dashed") +
  geom_segment(x=10.5,xend=50,y=10.5,yend=10.5,color="black",linetype = "dashed") +
  geom_segment(x=10.5,xend=10.5,y=0,yend=10.5,color="darkred",linetype = "dashed") +
  geom_segment(x=0,xend=10.5,y=10.5,yend=10.5,color="darkred",linetype = "dashed") +
  geom_point(data = max_genes%>% filter(Nbr_genes==50), aes(color=LL,x = genes_1, y = genes_2, shape = "Log-likelihood\nmaximum over\npairs"),size = 2,inherit.aes = FALSE) +
  scale_shape_manual(name = "",values = c("Log-likelihood\nmaximum over\npairs" = 10)) +
  guides(shape = guide_legend()) +
  labs(x="2nd gene",y="1st gene", color = "Log-Likelihood")+
  theme(axis.title=element_text(size=12),
        # legend.position="bottom",
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA))

aux = results %>% filter(rep==1,Nbr_genes==50)
df_to_plot = aux %>%
  mutate(genes=sapply(1:nrow(aux),FUN=function(i){
    paste0("(",paste0(range(aux$genes_1[i],aux$genes_2[i]),collapse=";"),")")
  })) %>%
  mutate(true = sapply(1:nrow(aux),FUN=function(i){
    if(aux$genes_1[i] %in% 1:10 && aux$genes_2[i] %in% 1:10){
      "Both"
    }else if(aux$genes_1[i] %in% 1:10 || aux$genes_2[i] %in% 1:10){
      "One"
    }else{
      "None"
    }
  })) %>%
  arrange(desc(LL))

df_to_plot$genes <- factor(df_to_plot$genes,levels=df_to_plot$genes)
df_to_plot$true <- factor(df_to_plot$true,levels=c("None","One","Both"))

init1 <- ggplot(df_to_plot, aes(x = genes, y = LL, shape = true, color = true)) +
  geom_point(size = 3) +
  ylab("Log-Likelihood") +
  xlab("Genes pair") +
  theme(
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
    strip.text = element_text(size = 11),
    plot.background = element_rect(fill = 'transparent', color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  scale_color_manual(
    name = "Contains\nTrue Positives",
    values = c("Both" = "dodgerblue", "One" = "darkseagreen", "None" = "darkred")
  ) +
  scale_shape_manual(
    name = "Contains\nTrue Positives",
    values = c("Both" = 16, "One" = 17, "None" = 15) # Example shape codes
  )


ggarrange(initLL+bb,
          init1+bb,
          ncol=1,heights = c(1.25,1),
          labels = c("A","B"))

saveFigures(paste0(pathToResults,"Figure2"),
            height = 3000,width = 3000,PNG = PNG,JPEG=JPEG,EPS=EPS)

### FIgure SUppMat

initLL <-
  ggplot(results_plot%>% filter(Nbr_genes==20),aes(x=genes_1,y=genes_2,color=LL))+
  geom_point(size=0.7)+
  # facet_grid(Nbr_genes~.)+
  scale_color_gradient2(high="darkred",mid="tomato",low="navajowhite",midpoint=sum(range((results%>% filter(Nbr_genes==20))$LL))/2) +
  geom_segment(x=5.5,xend=5.5,y=5.5,yend=20,color="black",linetype = "dashed") +
  geom_segment(x=5.5,xend=20,y=5.5,yend=5.5,color="black",linetype = "dashed") +
  geom_segment(x=5.5,xend=5.5,y=0,yend=5.5,color="darkred",linetype = "dashed") +
  geom_segment(x=0,xend=5.5,y=5.5,yend=5.5,color="darkred",linetype = "dashed") +
  geom_point(data = max_genes%>% filter(Nbr_genes==20), aes(color=LL,x = genes_1, y = genes_2, shape = "Log-likelihood\nmaximum over\npairs"),size = 2,inherit.aes = FALSE) +
  scale_shape_manual(name = "",values = c("Log-likelihood\nmaximum over\npairs" = 10)) +
  guides(shape = guide_legend()) +
  labs(x="2nd gene",y="1st gene", color = "Log-Likelihood")+
  theme(axis.title=element_text(size=12),
        # legend.position="bottom",
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +bb

saveFigures(paste0(pathToResults,"Figure2Supp"),
            height = 1800,width = 3000,PNG = PNG,JPEG=JPEG,EPS=EPS)

