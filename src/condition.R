library(dplyr)

path_archivo <- "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/counts_table.tsv"

metadata_table <- read.csv("/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/metadata/SraRunTable.csv",check.names = FALSE)

metadata <- metadata_table %>% select(srr_id = "Run", age = "AGE", sex = "sex") %>% mutate(condicion = paste0(age, "_", sex))

tabla_counts <- read.table(path_archivo, header = TRUE, sep = "\t", check.names = FALSE)

srr_names <- colnames(tabla_counts)[-1]

condition_names <- metadata$condicion[match(srr_names, metadata$srr_id)]

condition_names_clean <- gsub(" ", "_", condition_names)

# Reemplazar las barras "/" por guiones bajos "_"
condition_names_clean <- gsub("/", "_", condition_names_clean)

# Eliminar los puntos "." (los reemplazamos por texto vacío "")
condition_names_clean <- gsub("\\.", "", condition_names_clean)

colnames(tabla_counts) <- c("Geneid", make.unique(condition_names_clean, sep = "_rep"))

write.table(tabla_counts,"/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/condition_table.tsv", sep="\t", quote=F, row.names=F)
