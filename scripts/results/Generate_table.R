library(stringr)
library(dplyr)
library(huxtable)
library(xtable)
library(flextable)
library(webshot2)
library(chromote)
Nbr_genes =c(20,50)

options(chromote.chrome_args = c("--headless=new"))
dir <- function(d){if(!dir.exists(d)){dir.create(d)}}

tex2_names = c(delta_AB_pop = "\\delta_{Ab}",
               delta_S_pop = "\\delta_S",
               phi_S_pop = "\\varphi_S\\times10^{-2}",
               omega_delta_AB = "\\omega_{\\delta_{Ab}}",
               omega_delta_S = "\\omega_{\\delta_S}",
               omega_phi_S = "\\omega_{\\varphi_S}",
               setNames(paste0("\\alpha_{1",1:50,"}"),paste0("alpha_1",1:50)),
               setNames(paste0("\\alpha_{0",1:50,"}"),paste0("alpha_0",1:50)),
               sigma_AB = "\\sigma_{Ab}",
               setNames(paste0("\\sigma_{G",1:50,"}"),paste0("sigma_G",1:50)))
options(scipen=999)

source("scripts/results/resultsFun/format_latex_table.R", echo=FALSE)


# Load results ------------------------------------------------------------
load(paste0("outputs/simulationResults/simulationResults.RData"))

replicates = unique(genes$model)
nb_replicates = length(replicates)
precision = abs(floor(log10(1/length(replicates))))
method=c("REMix","REMixWald")

# Table Results -----------------------------------------------------------
dir("outputs/figures/")
pathToResults = "outputs/figures/finalTables"
dir(pathToResults)


for(ng in Nbr_genes){
  for(meth in method){
    stats = compute_stats_table(parameterValue,trueValueDF,ng,meth)

    signif <- ifelse(ng == 20, 5, 10)
    keep <- c("delta_AB_pop", "delta_S_pop", "phi_S_pop", "omega_delta_AB", "omega_delta_S", "omega_phi_S",
              paste0("alpha_1", 1:signif), paste0("alpha_0", 1:signif), "sigma_AB", paste0("sigma_G", 1:signif))
    # all(stats[stats$parameter %in% keep, ] %>% select(count)==200)
    df <- stats[stats$parameter %in% keep, ] %>% select(-count)
    # df_all <- stats %>% select(-count)

    df[df$parameter == "phi_S_pop", c("trueValue", "mean_Estimation", "empirical_sd", "estimated_sd")] <-
      df[df$parameter == "phi_S_pop", c("trueValue", "mean_Estimation", "empirical_sd", "estimated_sd")] * 0.01

    df <-
      format_df(df)
    # df_all <- format_df(df_all)

    lines = to_latex_table(df,keep,tex2_names,
                           caption=paste0("Parameter estimation accuracy for population and biomarker-specific parameters (",ng,"-gene scenario, 200 replicates,",ifelse(meth=="REMix"," before test)"," after test)"),"."))

    export_to_png(lines,paste0(pathToResults,"/",ng,"_",meth,".png"))

    writeLines(lines,con=paste0(pathToResults,"/",ng,"_",meth,".txt"))

    stats = compute_noise_table(parameterValue,genes,trueValueDF,ng,meth) %>%
      select(parameter,trueValue,mean_Estimation,empirical_sd,estimated_sd,count)

    df <- format_noise_df(stats)

    lines = noise_to_latex_table(df,tex2_names,
                                 caption=paste0("Parameter estimation accuracy for noisy biomarker-specific parameters (",ng,"-gene scenario, 200 replicates,",ifelse(meth=="REMix"," before test)"," after test)"),"."))
    writeLines(lines,con=paste0(pathToResults,"/",ng,"_",meth,"_noise.txt"))
    export_to_png(lines,paste0(pathToResults,"/",ng,"_",meth,"_noise.png"))
  }
}

