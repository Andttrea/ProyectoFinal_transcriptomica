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
FDR <- 0.01
LFC <- 1

up <- (res$log2FoldChange > LFC) & (res$padj < FDR) 
up[is.na(up)] <- FALSE
cat ("Upregulated: ", sum(up), "\n")
# Output: Upregulated:  2301

down <- (res$log2FoldChange < -LFC) & (res$padj < FDR)
down[is.na(down)] <- FALSE
cat ("Downregulated: ", sum(down), "\n")
# Output: Downregulated:  972

# Guardamos la tabla de resultados
write.table(res[up,], paste0(out_dir, "deseq-DEG_up_0.01.txt"), sep="\t", quote=FALSE, row.names=TRUE)
write.table(res[down,], paste0(out_dir, "deseq-DEG_down_0.01.txt"), sep="\t", quote=FALSE, row.names=TRUE)

# Convertimos resultados a data frame
res_df <- as.data.frame(res)

# Quitamos genes con padj NA
res_df <- res_df[!is.na(res_df$padj), ]

# Asignamos los colores para las categorías
vpcolors <- c("NO" = "gray", "DOWN" = "#790ebc", "UP" = "#558207")

# Creamos la columna DE
res_df$DE <- "NO"
res_df$DE[res_df$log2FoldChange > LFC & res_df$padj < FDR] <- "UP"
res_df$DE[res_df$log2FoldChange < -LFC & res_df$padj < FDR] <- "DOWN"

# La convertimos en factor para controlar el orden de la leyenda
res_df$DE <- factor(res_df$DE, levels = c("DOWN", "NO", "UP"))

# Creamos la gráfica
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), col = DE)) +
  geom_point(alpha = 0.4, size = 1.5) +
  labs(
    title = "Volcano plot: E14.5 vs E11.5",
    x = "log2 Fold Change",
    y = "-log10 FDR"
  ) +
  scale_color_manual(values = vpcolors) +
  geom_vline(xintercept = c(-LFC, LFC), col = "black", linetype = "longdash") +
  geom_hline(yintercept = -log10(FDR), col = "black", linetype = "longdash") +
  theme_classic(base_size = 15, base_line_size = 1)

# Mostrar la gráfica
volcano_plot

# Guardar la gráfica
ggsave(
  paste0(fig_dir, "volcano_plot_E145_vs_E115.png"),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Ahora haremos los heatmaps

# Filtramos solo genes significativos
significant <- res[up | down, ]

# Ordenamos por padj
significant_order <- significant[order(significant$padj), ]

# Tomamos los 2000 genes más significativos
top_genes <- head(rownames(significant_order), n = 2000)

zscore_t <-  t(scale(t(log2_tpm[top_genes,]))) 

# Ordenamos columnas por edad y luego por sexo
orden_columnas <- order(meta_data$age, meta_data$sex)
zscore_significant <- zscore_t[, orden_columnas]

# Usamos los nombres actuales de tus columnas como etiquetas
sample_labels <- c(
  "E14.5_F", "E14.5_F_1", "E14.5_F_2",
  "E14.5_M", "E14.5_M_1", "E14.5_M_2",
  "E11.5_F", "E11.5_F_1", "E11.5_F_2",
  "E11.5_M", "E11.5_M_1", "E11.5_M_2"
)

sample_labels <- sample_labels[orden_columnas]

# Creamos el heatmap
heatmap_plot <- Heatmap(
  zscore_significant,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  name = "Z-score",
  km = 2,
  column_title = "Heatmap de los genes diferencialmente expresados",
  column_names_gp = gpar(
    col = "black",
    fontsize = 10,
    fontface = "bold"
  ),
  top_annotation = HeatmapAnnotation(
    Edad = meta_data$age[orden_columnas],
    col = list(
      Edad = c(
        "embryonic_day_115" = "#58b4e1",
        "embryonic_day_145" = "#a473d1"
      )
    )
  )
)
heatmap_plot

# Guardamos la imagen
png(
  paste0(fig_dir, "heatmap_top_genes.png"),
  width = 11.25,
  height = 7.5,
  res = 300,
  units = "in"
)

draw(heatmap_plot)
dev.off()


# Para el top 20

# Seleccionamos los 20 genes más significativos
top_20_ids <- head(rownames(significant_order), n = 20)

zscore_top <- t(scale(t(log2_tpm[top_20_ids, ])))

zscore_top_ordenado <- zscore_top[, orden_columnas]


heatmap_top20 <- Heatmap(zscore_top_ordenado,
  cluster_rows = T,
  cluster_columns = F,
  row_labels = gene_name_map[rownames(zscore_top_ordenado), ],
  name = "Z-score",
  km = 2,
  column_title = "Top 20 Genes Significativos",
  column_names_gp = gpar(
  col = "black",
  fontsize = 10,
  fontface = "bold"
  ),
  row_names_gp = gpar(
    fontsize = 10,
    fontface = "italic" # Los nombres de genes siempre van en cursivas
    ),
   top_annotation = HeatmapAnnotation(
    Edad = meta_data$age[orden_columnas],
    col = list(
      Edad = c(
        "embryonic_day_115" = "#58b4e1",
        "embryonic_day_145" = "#a473d1")
  )))


heatmap_top20

# ---------------- V2 Heatmap -----------------------------------------------------------------------------------------#
# Etiquetas de columnas usando los nombres actuales de la tabla
sample_labels_top20 <- colnames(zscore_top20)[orden_columnas]

# Obtenemos nombres de genes
gene_labels <- gene_name_map[rownames(zscore_top20_ordenado), 1]
gene_labels <- as.character(gene_labels)

# Si algún gen no tiene nombre, usamos el Ensembl ID
gene_labels[is.na(gene_labels) | gene_labels == ""] <- rownames(zscore_top20_ordenado)[is.na(gene_labels) | gene_labels == ""]

# Creamos el heatmap top 20
heatmap_top20 <- Heatmap(
  zscore_top20_ordenado,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_labels = gene_labels,
  column_labels = sample_labels_top20,
  name = "Z-score",
  km = 2,
  column_title = "Top 20 genes diferencialmente expresados",
  column_names_gp = gpar(
    col = "black",
    fontsize = 10,
    fontface = "bold"
  ),
  row_names_gp = gpar(
    fontsize = 10,
    fontface = "italic"
  ),
  top_annotation = HeatmapAnnotation(
    Edad = meta_data$age[orden_columnas],
    col = list(
      Edad = c(
        "embryonic_day_115" = "#58b4e1",
        "embryonic_day_145" = "#a473d1"
      )
    )
  )
)

heatmap_top20
# ---------------- V2 Heatmap -----------------------------------------------------------------------------------------#

# Guardamos el heatmap top 20
png(
  paste0(fig_dir, "heatmap_top_20_genes"),
  width = 11.25,
  height = 7.5,
  res = 300,
  units = "in"
)

draw(heatmap_top20)
dev.off()


# --------------------------------------------------- Para calcular TPM -----------------------------------------------------------------#

annotation <- read.delim("/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/lenght_table.tsv", row.names = 1)
# Calculamos TPM  para el background de genes y guardamos la tabla de TPM log2 transformada
mcols(dds)$basepairs = annotation[rownames(dds),]
log2_fpkm = log2(fpkm(dds)+0.1) 
fpkm2tpm_log2 <- function(fpkm) { fpkm - log2(sum(2^fpkm)) + log2(1e6) } 
log2_tpm = apply(log2_fpkm, 2, fpkm2tpm_log2) 
gene_names = gene_name_map[rownames(log2_tpm),]
write.table(cbind(gene_names, log2_tpm), paste0(out_dir, "TPM_log2-table.txt"), sep="\t", quote=FALSE)

# --------------------------------------------------- Para calcular TPM -----------------------------------------------------------------#

saveRDS(res, file = "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/deseq/res.rds")


# Seleccion de genes para STRING

up_string <- as.data.frame(res) %>% filter(padj < 0.01 & log2FoldChange > 1)
up_order <- up_string[order(up_string$padj), ]
string <- head(rownames(up_order), n = 1000)
write.table(
  string,
  file = "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/deseq/genes_UP_top1000_STRING.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
