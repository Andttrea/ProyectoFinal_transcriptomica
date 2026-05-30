library(dplyr)

path_archivo <- "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/count_table.tsv"

metadata_table <- read.csv("/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/metadata/SraRunTable.csv",check.names = FALSE)

metadata <- metadata_table %>% select(srr_id = "Run", tratamiento = "treatment", geo_id = "Sample Name") %>% mutate(condicion = tratamiento)

tabla_counts <- read.table(path_archivo, header = TRUE, sep = "\t", check.names = FALSE)

srr_names <- colnames(tabla_counts)[-1]

condition_names <- metadata$condicion[match(srr_names, metadata$srr_id)]

colnames(tabla_counts) <- c("Geneid", make.unique(condition_names, sep = "_rep"))

write.table(tabla_counts,"/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/condition_table.tsv", sep="\t", quote=F, row.names=F)