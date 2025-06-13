# These functions allows to save in txt file the stats table computed in a nice formatted way

# PAS NOISE ---------------------------------------------------------------

compute_stats_table <-
  function(parameterValue,trueValueDF,Proj,Nbr_genes,method,replicates=1:200){

    parameterValue = parameterValue[parameterValue$Nbr_genes==Nbr_genes &
                                      parameterValue$method==method &
                                      parameterValue$Proj==Proj,]

    # Estimation accuracy for all parameters throughout all replicates
    trueValue = setNames(trueValueDF[trueValueDF$Nbr_genes==Nbr_genes,"trueValue"],trueValueDF[trueValueDF$Nbr_genes==Nbr_genes,"parameter"])

    mean_Estimation = tapply(parameterValue$values,parameterValue$parameter,mean)
    relative_bias = as.numeric((mean_Estimation - trueValue[names(mean_Estimation)])/trueValue[names(mean_Estimation)])*100
    empirical_sd = tapply(parameterValue$values,parameterValue$parameter,sd)
    estimated_sd = tapply(parameterValue$se,parameterValue$parameter,FUN=function(x){mean(x,na.rm=TRUE)})
    countNA = tapply(parameterValue$se,parameterValue$parameter,FUN=function(x){sum(is.na(x))})

    coverHess = setNames(rep(0,length(mean_Estimation)),names(mean_Estimation))
    nCount = setNames(rep(length(replicates),length(mean_Estimation)),names(mean_Estimation))
    for (j in replicates){
      sdhat = setNames(parameterValue[parameterValue$model==j,"se"],parameterValue[parameterValue$model==j,"parameter"])[names(mean_Estimation)]
      for(par in names(sdhat)){
        if(is.na(sdhat[par])){
          nCount[par] <- nCount[par] -1
        }else{
          lw = parameterValue[parameterValue$parameter==par & parameterValue$model==j,"values"] - 1.96*sdhat[par]
          up =parameterValue[parameterValue$parameter==par & parameterValue$model==j,"values"] + 1.96*sdhat[par]
          bool = as.numeric(trueValue[par] > lw & trueValue[par] < up)
          coverHess[par] = coverHess[par] + bool
        }
      }
    }
    coverHess <- coverHess/nCount

    res = data.frame(parameter = names(mean_Estimation),
                     trueValue = unname(as.numeric(trueValue[names(mean_Estimation)])),
                     mean_Estimation= unname(mean_Estimation),
                     relative_bias = unname(relative_bias),
                     empirical_sd = unname(empirical_sd[names(mean_Estimation)]),
                     estimated_sd = unname(estimated_sd[names(mean_Estimation)]),
                     count = nb_replicates-unname(countNA[names(mean_Estimation)]),
                     coverHess = unname(coverHess[names(mean_Estimation)])*100,
                     tab="all")

    return(res)
  }


format_num <- function(x, digits,sd=FALSE, cover = FALSE) {
  sapply(x, function(xi) {
    if (is.na(xi)) return(xi)

    xi_str <- format(round(xi, digits), nsmall = digits)

    if (!sd && !cover && xi > 0)
      xi_str=paste0("~", xi_str)

    return(xi_str)
  })
}

format_df <- function(df) {
  df <- df %>%
    mutate(
      trueValue = format_num(trueValue, 3),
      mean_Estimation = format_num(mean_Estimation, 3),
      relative_bias = format_num(relative_bias, 2),
      empirical_sd = format_num(empirical_sd, 4,sd=TRUE),
      estimated_sd = format_num(estimated_sd, 4,sd=TRUE),
      coverHess = format_num(coverHess,1,cover=TRUE)
    )
  return(df)
}

to_latex_table <- function(df,keep,tex2_names,caption=""){
  df <- df %>% select(-tab) %>% filter(parameter %in% keep)
  df$parameter <- factor(df$parameter,levels=keep)
  df <- df %>% arrange(parameter)

  colnames(df) <- c("parameter","Target Value","Mean Estimate","Relative Bias (\\%)","Empirical SD","Estimated SD","Estimated Cover (\\%)")

  latex_names <- paste0("$", unname(tex2_names[keep]), "$")
  hux <- huxtable(cbind(Parameter = latex_names, df %>% select(-parameter)), add_rownames = FALSE)[-1, ]
  hux <- set_escape_contents(hux, FALSE)
  tex <- xtable(hux, caption = caption)
  lines = capture.output(print(tex, sanitize.text.function = identity, include.rownames = FALSE))[-c(1, 2)]

  header1 <- "    & Target & Mean & Relative & Empirical & Estimated & Estimated \\\\"
  header2 <- "Parameter & Value & Estimate & Bias (\\%) & SD & SD & Cover (\\%) \\\\ \\hline"

  # Remplacement de l’en-tête généré automatiquement
  lines <- c(lines[1:4], header1, header2, lines[6:length(lines)])
  return(lines)
}



# NOISE -------------------------------------------------------------------

compute_noise_table <-
  function(parameterValue,genes,trueValueDF,Proj,Nbr_genes,method,replicates=1:200){

    if(Nbr_genes==20){
      signif=5
    }else{
      signif=10
    }

    parameterValue = parameterValue[parameterValue$Nbr_genes==Nbr_genes &
                                      parameterValue$method==method &
                                      parameterValue$Proj==Proj,]

    genes_aux <- genes[genes$Nbr_genes==Nbr_genes &
                         genes$method==method &
                         genes$Proj==Proj &
                         genes$genes >signif,]

    trueValue = setNames(trueValueDF[trueValueDF$Nbr_genes==Nbr_genes,"trueValue"],trueValueDF[trueValueDF$Nbr_genes==Nbr_genes,"parameter"])

    genes_selected0 <- genes_aux[genes_aux$selected,c("model","parameter")]
    genes_selected <- rbind(genes_selected0,genes_selected0 %>% mutate(parameter=str_replace_all(genes_selected0$parameter,"alpha_1","alpha_0")))
    genes_selected <- rbind(genes_selected,genes_selected0 %>% mutate(parameter=str_replace_all(genes_selected0$parameter,"alpha_1","sigma_G")))

    parameterValue_aux <- merge(genes_selected,parameterValue,by=c("model","parameter"))

    mean_Estimation = tapply(parameterValue_aux$values,parameterValue_aux$parameter,mean)
    relative_bias = as.numeric((mean_Estimation - trueValue[names(mean_Estimation)])/trueValue[names(mean_Estimation)])*100
    empirical_sd = tapply(parameterValue_aux$values,parameterValue_aux$parameter,sd)
    estimated_sd = tapply(parameterValue_aux$se,parameterValue_aux$parameter,FUN=function(x){mean(x,na.rm=TRUE)})
    count = sapply(names(mean_Estimation),FUN=function(p){nrow(parameterValue_aux[parameterValue_aux$parameter==p,])})

    coverHess = setNames(rep(0,length(mean_Estimation)),names(mean_Estimation))
    nCount = setNames(rep(0,length(mean_Estimation)),names(mean_Estimation))
    for(p in names(mean_Estimation)){
      for (j in unique(parameterValue_aux[parameterValue_aux$parameter==p,"model"])){
        sdhat = parameterValue_aux[parameterValue_aux$model==j & parameterValue_aux$parameter==p,"se"]
        if(!is.na(sdhat)){
          nCount[p] <- nCount[p] + 1
          lw = parameterValue_aux[parameterValue_aux$parameter==p & parameterValue_aux$model==j,"values"] - 1.96*sdhat
          up =parameterValue_aux[parameterValue_aux$parameter==p & parameterValue_aux$model==j,"values"] + 1.96*sdhat
          bool = as.numeric(trueValue[p] > lw & trueValue[p] < up)
          coverHess[p] = coverHess[p] + bool
        }
      }
    }
    coverHess <- coverHess/nCount

    res =data.frame(parameter = names(mean_Estimation),
                    trueValue = unname(as.numeric(trueValue[names(mean_Estimation)])),
                    mean_Estimation= unname(mean_Estimation),
                    relative_bias = unname(relative_bias),
                    empirical_sd = unname(empirical_sd[names(mean_Estimation)]),
                    estimated_sd = unname(estimated_sd[names(mean_Estimation)]),
                    count = unname(count[names(mean_Estimation)])/length(replicates)*100,
                    coverHess = unname(coverHess[names(mean_Estimation)])*100,
                    tab="noise")

    return(res)
  }

# Fonction de formatage numérique pour ce type de tableau
format_noise_num <- function(x, digits, sd = FALSE) {
  sapply(x, function(xi) {
    if (is.na(xi)) return(xi)

    xi_str <- format(round(xi, digits), nsmall = digits)

    if (!sd && as.numeric(xi) >= 0) xi_str <- paste0("~", xi_str)

    return(xi_str)
  })
}

# Fonction pour homogénéiser visuellement la colonne "Number of Selection"
format_selection_freq <- function(x) {
  x <- format(round(x, 1), nsmall = 1)

  max_len <- max(nchar(x))
  return(sapply(x, function(val) {
    val_str <- as.character(val)
    padding <- paste(rep("~", max_len - nchar(val_str)), collapse = "")
    paste0(padding, val_str)
  },USE.NAMES = FALSE))
}

# Fonction de formatage du tableau noise
format_noise_df <- function(df) {
  df <- df %>%
    mutate(
      trueValue = format_noise_num(trueValue, 3),
      mean_Estimation = format_noise_num(mean_Estimation, 3),
      empirical_sd = format_noise_num(empirical_sd, 4, sd = TRUE),
      estimated_sd = format_noise_num(estimated_sd, 4, sd = TRUE),
      count = format_selection_freq(count)
    )
  return(df)
}

# Fonction latex pour tableau avec gènes groupés
noise_to_latex_table <- function(df,tex2_names,caption="") {
  df <- df %>% arrange(parameter)
  g <- as.integer(str_extract(df$parameter, "(?<=_.).+"))
  latex_names <- paste0("$", tex2_names[df$parameter], "$")

  df_out <- data.frame(Gene = paste0("G", g),
                       Parameter = latex_names,
                       df[, c("trueValue", "mean_Estimation", "empirical_sd", "estimated_sd", "count")])

  df_out$Gene <- factor(df_out$Gene,levels=paste0("G",sort(unique(g))))

  hux <- huxtable(df_out %>% arrange(Gene,Parameter), add_rownames = FALSE)[-1, ]
  hux <- set_escape_contents(hux, FALSE)
  tex <- xtable(hux, caption = caption)
  lines = capture.output(print(tex, sanitize.text.function = identity, include.rownames = FALSE))[-c(1, 2)]

  header1 <- "     &           & Target & Mean     & Empirical & Estimated & Selection \\\\"
  header2 <- "Gene & Parameter & Value  & Estimate & SD        & SD        & Frequency  (\\%) \\\\ \\hline"

  lines <- c(lines[1:4], header1, header2, lines[6:length(lines)])

  multicol = seq(8,8+length(unique(g))*3-3,length.out=length(unique(g)))
  toerase = sort(union(seq(8,8+length(unique(g))*3-3,length.out=length(unique(g)))+1,seq(8,8+length(unique(g))*3-3,length.out=length(unique(g)))+2))
  for(l in multicol){
    newline = str_split(lines[l],"&")[[1]]
    newline[1] <- paste0("\\multirow{3}{*}{",str_remove_all(newline[1]," "),"} ")
    newline[7] <- paste0("\\multirow{3}{*}{",str_remove_all(newline[7]," \\\\\\\\"),"}"," \\\\")
    lines[l] <- paste0(newline,collapse="&")
  }
  for(l in toerase){
    newline = str_split(lines[l],"&")[[1]]
    newline[1] <- " "
    newline[7] <- " \\\\"
    lines[l] <- paste0(newline,collapse="&")
  }

  before=lines[1:7]
  middle=lines[8:max(toerase)]
  after =lines[(max(toerase)+1):length(lines)]


  middle_modified <- character(0)

  for(i in seq(1,length(middle),3)){
    end <- min(i+2,length(middle))
    middle_modified <- c(middle_modified,middle[i:end])
    if(end<length(middle)){
      middle_modified <- c(middle_modified,"\\hline")
    }
  }

  lines <- c(before, middle_modified, after)

  return(lines)
}


# For all -----------------------------------------------------------------

export_to_png <- function(latex_lines,pathToSave){
  temp_dir <- dirname(pathToSave)
  tex_path <- file.path(temp_dir, "temp_table.tex")
  pdf_path <- file.path(temp_dir, "temp_table.pdf")

  latex_doc <- c(
    "\\documentclass{article}",
    "\\usepackage{booktabs}",
    "\\usepackage[margin=1in]{geometry}",
    "\\usepackage{multirow}",
    "\\usepackage{multicol}",
    "\\begin{document}",
    latex_lines,
    "\\end{document}"
  )

  writeLines(latex_doc, tex_path)
  tinytex::latexmk(tex_path)

  img <- magick::image_read_pdf(pdf_path, density = 300)
  img <- magick::image_trim(img)
  magick::image_write(img, pathToSave, format = "png")

  unlink(tex_path)
  unlink(pdf_path)
}
