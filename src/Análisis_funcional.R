library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)

up_go <- read.csv("/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/anotacion/data/DAVIDChartReport_upregulated_BP.csv", header = TRUE)
down_go <- read.csv("/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/anotacion/data/DAVIDChartReport_downregulated_BP.csv", header = TRUE)

