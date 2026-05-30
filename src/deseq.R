# Cargamos las librerias 
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(edgeR)
library(grid)

# Definimos rutas de archivos y directorios de salida 
gene_name_map <- read.delim("/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/GENCODE/G_ID_NAME.tsv", header = FALSE, row.names = 1)
counts_file <- "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/condition_table.tsv"
out_dir <- "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/deseq/"
fig_dir <- paste0(out_dir, "figuras/")

# Leemos la tabla de conteos 
counts_table <- read.delim(counts_file,header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Definimos los factores
age <- factor(c(rep("embryonic_day_145", 6), rep("embryonic_day_115", 6)), levels = c("embryonic_day_115", "embryonic_day_145"))
sex <- factor(
  c(rep("female", 3), rep("male", 3),
    rep("female", 3), rep("male", 3)))

sample_names <- colnames(counts_table)

# Generamos la metadata
metadata <- data.frame(sample_names, age, sex)
meta_data <- metadata %>% remove_rownames %>% column_to_rownames(var = "sample_names")
# Checaremos que los nombres de la tabla de metadatos estan en el mismo orden que las columnas de los conteos
k <- all(colnames(counts_table) == rownames(meta_data))
print(paste0("Los nombres de las columnas de counts_table y los nombres de fila de meta_data coinciden: ", k))

# Hacemos nuestro diseño experimental
dds <- DESeqDataSetFromMatrix(countData = round(counts_table), colData = meta_data, design = ~ 0 + age + sex)
design <- model.matrix( ~ 0 + age + sex)

keep <- filterByExpr(dds, design) 
suma_keep <- sum(keep)
suma_keep
dds <- dds[keep,]

# Visualizamos los datos con PCA
vsd <- vst(dds)
plotPCA(vsd, intgroup = c("age")) + aes(shape=sex) + theme_classic(base_size=18, base_line_size = 1)
ggsave(paste0(fig_dir, "PCA_plot.png"), plot = PCA_plot, width = 12, height = 7)

# Realizamos el análisis de DESeq2
dds <- DESeq(dds)
resultsNames(dds)

# Hacemos el contraste 
contrast <- makeContrasts(E145_vs_E115 = ageembryonic_day_145 - ageembryonic_day_115, levels = design)
# Aplicamos el contraste a los resultados de DESeq2
res <- results(dds, contrast=contrast[,"E145_vs_E115"])
res$Gene_name <- gene_name_map[rownames(res),] 

# Definimos FDR y LFC
FDR <- 0.05
LFC <- 1

up <- (res$log2FoldChange > LFC) & (res$padj < FDR) 
up[is.na(up)] <- FALSE
cat ("Upregulated: ", sum(up), "\n")

down <- (res$log2FoldChange < -LFC) & (res$padj < FDR)
down[is.na(down)] <- FALSE
cat ("Downregulated: ", sum(down), "\n")

















annotation <- read.delim("/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/lenght_table.tsv", row.names = 1)
# Calculamos TPM  para el background de genes y guardamos la tabla de TPM log2 transformada
mcols(dds)$basepairs = annotation[rownames(dds),]
log2_fpkm = log2(fpkm(dds)+0.1) 
fpkm2tpm_log2 <- function(fpkm) { fpkm - log2(sum(2^fpkm)) + log2(1e6) } 
log2_tpm = apply(log2_fpkm, 2, fpkm2tpm_log2) 
gene_names = gene_name_map[rownames(log2_tpm),]
write.table(cbind(gene_names, log2_tpm), paste0(out_dir, "TPM_log2-table.txt"), sep="\t", quote=FALSE)

