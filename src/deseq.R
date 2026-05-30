# Cargamos las librerias 
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(edgeR)
library(grid)

# Definimos rutas de archivos y directorios de salida 
counts_file <- "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/condition_table.tsv"
out_dir <- "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/deseq/"
fig_dir <- paste0(out_dir, "figuras/")

# Leemos la tabla de conteos 
counts_table <- read.delim(counts_file,header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Generamos factores
treatment <- factor(c(rep("7mM_glucose",5),rep("5_5p5mM_glucose",5)), levels = c("5_5p5mM_glucose", "7mM_glucose"))

sample_names <- colnames(counts_table)

# Generamos la metadata
metadata <- data.frame(sample_names, treatment)
meta_data <- metadata %>% remove_rownames %>% column_to_rownames(var = "sample_names")
# Checaremos que los nombres de la tabla de metadatos estan en el mismo orden que las columnas de los conteos
k <- all(colnames(counts_table) == rownames(meta_data))
print(paste0("Los nombres de las columnas de counts_table y los nombres de fila de meta_data coinciden: ", k))

# Hacemos el análisis de expresión diferencial 
dds <- DESeqDataSetFromMatrix(countData = round(counts_table), colData = meta_data, design =  ~ 0 + treatment)
design <- model.matrix(~ 0 + treatment)

keep <- filterByExpr(dds, design) 
suma_keep <- sum(keep)
dds <- dds[keep,]

vsd <- vst(dds)
plotPCA(vsd, intgroup = "treatment") + theme_classic(base_size=25, base_line_size = 1)

# Usar la primera columna como rownames
gene_ids <- counts_table[[1]]
gene_counts <- counts_table[, -1]

rownames(gene_counts) <- gene_ids

# Asegurar que los conteos sean enteros
gene_counts <- round(as.matrix(gene_counts))

# Crear metadata a partir de los nombres de las columnas
sample_names_original <- colnames(gene_counts)

# Quitar sufijos _rep1, _rep2, etc. para recuperar condición
condition_raw <- gsub("_rep[0-9]+$", "", sample_names_original)

# Convertir nombres de condición a nombres seguros
condition <- case_when(
  condition_raw == "7 mM glucose" ~ "glucose_7mM",
  condition_raw == "5/5.5 mM glucose" ~ "normo_5_5p5mM",
  TRUE ~ condition_raw
)

# Filtrar solo las condiciones que nos interesan
keep_samples <- condition %in% c("normo_5_5p5mM", "glucose_7mM")

gene_counts <- gene_counts[, keep_samples]
condition <- condition[keep_samples]
sample_names_original <- sample_names_original[keep_samples]

# Crear nombres de muestra seguros y únicos
sample_names <- make.unique(condition, sep = "_rep")

colnames(gene_counts) <- sample_names

meta_data <- data.frame(
  sample_name = sample_names,
  original_name = sample_names_original,
  condition = factor(condition, levels = c("normo_5_5p5mM", "glucose_7mM"))
)

meta_data <- meta_data %>%
  remove_rownames() %>%
  column_to_rownames(var = "sample_name")

# Verificar que el orden coincida
print(meta_data)
print(table(meta_data$condition))

k <- all(colnames(gene_counts) == rownames(meta_data))
print(paste0("Los nombres de gene_counts y meta_data coinciden: ", k))

if (!k) {
  stop("Error: los nombres de columnas de conteos no coinciden con la metadata.")
}

# Guardar metadata usada
write.table(
  meta_data,
  paste0(out_dir, "metadata_deseq2.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

# ============================
# Crear objeto DESeq2
# ============================

dds <- DESeqDataSetFromMatrix(
  countData = gene_counts,
  colData = meta_data,
  design = ~ condition
)

# ============================
# Filtrar genes de baja expresión
# ============================

design_matrix <- model.matrix(~ condition, data = meta_data)

keep <- filterByExpr(counts(dds), design = design_matrix)

suma_keep <- sum(keep)

write.table(
  suma_keep,
  paste0(out_dir, "genes_filtrados.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

dds <- dds[keep, ]

# ============================
# Transformación VST y PCA
# ============================

vsd <- vst(dds, blind = FALSE)

pca_plot <- plotPCA(vsd, intgroup = "condition") +
  theme_classic(base_size = 18) +
  labs(
    title = "PCA: 7 mM glucose vs 5/5.5 mM glucose",
    color = "Condition"
  )

ggsave(
  paste0(fig_dir, "PCA_plot.png"),
  plot = pca_plot,
  width = 8,
  height = 6
)

# Guardar matriz VST
vst_matrix <- assay(vsd)

write.table(
  vst_matrix,
  paste0(out_dir, "vst_matrix.tsv"),
  sep = "\t",
  quote = FALSE
)

# ============================
# Correr DESeq2
# ============================

dds <- DESeq(dds)

saveRDS(dds, paste0(out_dir, "dds_DESeq2.rds"))

# Contraste:
# log2FC positivo = mayor expresión en 7 mM glucose
# log2FC negativo = mayor expresión en 5/5.5 mM glucose

res <- results(
  dds,
  contrast = c("condition", "glucose_7mM", "normo_5_5p5mM")
)

res <- as.data.frame(res)
res$gene_id <- rownames(res)

# Ordenar por padj
res <- res %>%
  arrange(padj)

write.table(
  res,
  paste0(out_dir, "DESeq2_results_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(res, paste0(out_dir, "DESeq2_results_all.rds"))

# ============================
# Genes significativos
# ============================

FDR <- 0.05
LFC <- 0.5

res$DE <- "NO"
res$DE[!is.na(res$padj) & res$padj < FDR & res$log2FoldChange > LFC] <- "UP"
res$DE[!is.na(res$padj) & res$padj < FDR & res$log2FoldChange < -LFC] <- "DOWN"

up <- res %>% filter(DE == "UP")
down <- res %>% filter(DE == "DOWN")
deg <- res %>% filter(DE %in% c("UP", "DOWN"))

write.table(
  nrow(up),
  paste0(out_dir, "genes_up_regulados.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  nrow(down),
  paste0(out_dir, "genes_down_regulados.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  up,
  paste0(out_dir, "DESeq2_DEG_up.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  down,
  paste0(out_dir, "DESeq2_DEG_down.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  deg,
  paste0(out_dir, "DESeq2_DEG_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================
# Volcano plot
# ============================

res$minus_log10_padj <- -log10(res$padj)

vpcolors <- c("NO" = "gray", "DOWN" = "#2C7BB6", "UP" = "#D7191C")

volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = minus_log10_padj, color = DE)) +
  geom_point(alpha = 0.5, size = 1.4) +
  scale_color_manual(values = vpcolors) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "longdash", color = "black") +
  geom_hline(yintercept = -log10(FDR), linetype = "longdash", color = "black") +
  theme_classic(base_size = 15) +
  labs(
    title = "Volcano plot: 7 mM glucose vs 5/5.5 mM glucose",
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = "DE"
  )

ggsave(
  paste0(fig_dir, "volcano_plot.png"),
  plot = volcano_plot,
  width = 8,
  height = 6
)

# ============================
# Heatmap de genes significativos
# ============================

if (nrow(deg) > 1) {
  
  # Tomar hasta 2000 genes significativos ordenados por padj
  top_genes <- head(deg$gene_id, n = min(2000, nrow(deg)))
  
  zscore_t <- t(scale(t(vst_matrix[top_genes, , drop = FALSE])))
  
  orden_columnas <- order(meta_data$condition)
  zscore_significant <- zscore_t[, orden_columnas, drop = FALSE]
  
  heatmap_plot <- Heatmap(
    zscore_significant,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    name = "Z-score",
    km = 2,
    column_title = "Heatmap de genes diferencialmente expresados",
    column_names_gp = gpar(
      col = "black",
      fontsize = 10,
      fontface = "bold"
    ),
    top_annotation = HeatmapAnnotation(
      Condition = meta_data$condition[orden_columnas],
      col = list(
        Condition = c(
          "normo_5_5p5mM" = "#4DAF4A",
          "glucose_7mM" = "#E41A1C"
        )
      )
    )
  )
  
  png(
    paste0(fig_dir, "heatmap_DEG_genes.png"),
    width = 11.25,
    height = 7.5,
    res = 300,
    units = "in"
  )
  draw(heatmap_plot)
  dev.off()
}

# ============================
# Heatmap top 20
# ============================

if (nrow(deg) >= 2) {
  
  top_20_ids <- head(deg$gene_id, n = min(20, nrow(deg)))
  
  zscore_top <- t(scale(t(vst_matrix[top_20_ids, , drop = FALSE])))
  
  orden_columnas <- order(meta_data$condition)
  zscore_top_ordenado <- zscore_top[, orden_columnas, drop = FALSE]
  
  heatmap_top20 <- Heatmap(
    zscore_top_ordenado,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    name = "Z-score",
    column_title = "Top 20 genes significativos",
    column_names_gp = gpar(
      col = "black",
      fontsize = 10,
      fontface = "bold"
    ),
    row_names_gp = gpar(
      fontsize = 8,
      fontface = "italic"
    ),
    top_annotation = HeatmapAnnotation(
      Condition = meta_data$condition[orden_columnas],
      col = list(
        Condition = c(
          "normo_5_5p5mM" = "#4DAF4A",
          "glucose_7mM" = "#E41A1C"
        )
      )
    )
  )
  
  png(
    paste0(fig_dir, "heatmap_top_20_genes.png"),
    width = 11.25,
    height = 7.5,
    res = 300,
    units = "in"
  )
  draw(heatmap_top20)
  dev.off()
}

cat("Análisis DESeq2 terminado.\n")
cat("Genes filtrados retenidos:", suma_keep, "\n")
cat("Genes UP en 7 mM glucose:", nrow(up), "\n")
cat("Genes DOWN en 7 mM glucose:", nrow(down), "\n")
cat("Resultados guardados en:", out_dir, "\n")

