data <- read.csv("data/applicationFiles/Rinchai/data/data.txt")
nbgenes <- length(unique(data$yobs))-1

linesmodel <- readLines("data/applicationFiles/Rinchai/model/model_blank_delay.txt")

input_lines <- stringr::str_sub(linesmodel[2],end=-2)
newinput_lines <- paste0(input_lines,",",paste0(c("alpha_0","alpha_1"),rep(1:nbgenes,each=2),collapse=","),"}")

output_lines <- stringr::str_sub(linesmodel[49],end=-2)
newoutput_lines <- paste0(output_lines,",",paste0("G",1:nbgenes,collapse=","),"}")

add_gene <- paste0("G",1:nbgenes," = alpha_0",1:nbgenes," + alpha_1",1:nbgenes," * Cp*exp(-delta_V*(t-ttildep-td))")

model <- c(linesmodel[1],newinput_lines, linesmodel[3:47],add_gene,linesmodel[47:48],newoutput_lines)

writeLines(model,con="data/applicationFiles/Rinchai/model/model_for_REMix_delay.txt")



