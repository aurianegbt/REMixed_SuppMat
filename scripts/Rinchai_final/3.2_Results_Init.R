library(REMixed)
library(ggplot2)
library(dplyr)

load("data/applicationFiles/Rinchai/correspondance.RData")

results = data.frame()
ranked= data.frame()

for(seed in 1:10){
  res_init <- list()
  for(i in 1:17){
    load(paste0("~/Travail/REMix_PLAFRIM/outputs/applicationResults/Rinchai/initializationDelay/seed_",seed,"/init_",i,".RData"))

    res$genes <- correspondance_function[correspondance_function$yG %in% paste0("yG",res$genes),"yobs"]

    res_init <- append(res_init,list(res))
  }

  aux <- data.frame(seed=seed,
                    arr = 1:17,
                    genes1 = unname(sapply(res_init,function(i){i$genes[1]})),
                    genes2 = unname(sapply(res_init,function(i){i$genes[2]})),
                    genes3 = unname(sapply(res_init,function(i){i$genes[3]})),
                    genes = unname(sapply(res_init,function(i){paste0("(",paste0(sort(i$genes),collapse=","),")")})),
                    LL = sapply(res_init,function(i){-1/2*i$LL[["OFV"]]}),
                    alpha_1G1 = sapply(res_init,function(i){i$parameters[paste0("alpha_1",rownames(correspondance_function[correspondance_function$yobs==i$genes[1],]),"_pop")]}),
                    alpha_1_G2 = sapply(res_init,function(i){i$parameters[paste0("alpha_1",rownames(correspondance_function[correspondance_function$yobs==i$genes[2],]),"_pop")]}),
                    alpha_1_G3 = sapply(res_init,function(i){ifelse(!is.na(i$genes[3]),i$parameters[paste0("alpha_1",rownames(correspondance_function[correspondance_function$yobs==i$genes[3],]),"_pop")],0)}),
                    delta_V = sapply(res_init,function(i){i$parameters["delta_V_pop"]}),
                    theta = sapply(res_init,function(i){i$parameters["theta_pop"]}),
                    fm2 = sapply(res_init,function(i){i$parameters["fM2_pop"]}),
                    td = sapply(res_init,function(i){i$parameters["td_pop"]}))  %>%
    arrange(desc(LL))

  aux <- cbind(aux,rank=1:17)

  results <- rbind(results,
                   aux)

  ranked = rbind(ranked,
                 rbind(aux[,c("seed","genes1","rank")] %>% rename(genes=genes1),
                       rbind(aux[,c("seed","genes2","rank")]%>% rename(genes=genes2),
                             aux[,c("seed","genes3","rank")]%>% rename(genes=genes3))))

  ranked <- ranked[!is.na(ranked$genes),]


}

library(ggplot2)

df_to_plot <- results %>% arrange(desc(LL))
#
# ggplot(df_to_plot,aes(x=factor(genes,levels = unique(genes)),y=LL))+geom_point()+theme(axis.text.x = element_text(angle=90))
#
# ggplot(df_to_plot[1:20,],aes(x=factor(genes,levels = unique(genes)),y=LL))+geom_point()+theme(axis.text.x = element_text(angle=90))

df_to_plot_sorted = data.frame()
already_explore = c()
for(g in unique(sort(correspondance_function$yobs))){
  aux1 = df_to_plot %>% filter(genes1 %in% g)
  aux2 = df_to_plot %>% filter(genes2 %in% g)

  aux = rbind(aux1,(aux2 %>% rename(genes_2=genes1) %>% rename(genes1 = genes2) %>% rename(genes2=genes_2))[,colnames(aux1)]) %>% arrange(genes2) %>% select(-arr,-seed,-rank) %>% unique()

  already_explore= c(already_explore,g)

  aux = aux %>% filter(!(genes2 %in% already_explore))

  df_to_plot_sorted = rbind(df_to_plot_sorted,aux)
}


ggplot(df_to_plot_sorted,aes(x=factor(genes1,levels=sort(unique(genes1),decreasing = TRUE)),
                      y=factor(genes2,levels=sort(unique(genes2),decreasing = TRUE)),
                      color=LL))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.4),
        axis.text.y = element_text(vjust = 0.4))+
  geom_point(size=3)+
  scale_color_gradient2(high="darkred",mid="tomato",low="navajowhite",midpoint=sum(range(df_to_plot$LL))/2)+
  geom_point(data = df_to_plot[1,,drop=FALSE], aes(color=LL,x = genes1, y = genes2, shape = "Log-likelihood\nmaximum over\npairs"),size = 5,inherit.aes = FALSE) +
  scale_shape_manual(name = "",values = c("Log-likelihood\nmaximum over\npairs" = 10)) +
  guides(shape = guide_legend()) +
  labs(x="2nd gene",y="1st gene", color = "Log-Likelihood")+
  theme(axis.title=element_text(size=12),
        # legend.position="bottom",
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA))

ggsave(filename="outputs/figures/finalFigures/Initialization_app.eps",device="eps",height=6,width=8)

