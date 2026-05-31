# Cargar librerías necesarias
library(DESeq2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot) 
library(ggplot2)    

# Definimos rutas de entrada y salida
input_file <- "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/deseq/res.rds"
output_dir <- "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/anotacion/GSEA"

# Cargamos los resultados de DESeq2
res <- readRDS(file = input_file)

# Convertir a data.frame por si viene como objeto DESeqResults
res_df <- as.data.frame(res)

# Preparamos  la lista de genes (gene_list) y filtramos los valores NA en la columna 'stat'. 
# NOTA: DESeq2 suele generar NAs por conteos bajos o valores atípicos, y si pasamos un NA a clusterProfiler, la función fallará
res_df <- res_df[!is.na(res_df$stat), ]

# Ordenamos de mayor a menor basado en el estadístico (stat)
# Esto es esencial para GSEA, ya que evalúa si los genes de una vía se acumulan en los extremos (up o down regulados)
res_df <- res_df[order(-res_df$stat), ]

# Extraer el vector numérico y asignarle los nombres de los genes (Ensembl IDs)
gene_list <- res_df$stat
names(gene_list) <- rownames(res_df)

# 4. Ejecutar el Gene Set Enrichment Analysis (GSEA)
gse <- gseGO(
  geneList     = gene_list,
  ont          = "BP",            # Biological Process
  keyType      = "ENSEMBL",
  OrgDb        = org.Mm.eg.db,    # Base de datos para Mus musculus
  eps          = 1e-300,
  pvalueCutoff = 0.05             # Umbral de significancia
)

write.table(gse[ , c(1:10)], file ="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/anotacion/GSEA/gse_results.txt", sep = "\t")

# Generamos el gseaplot2 múltiple
# NOTA: Seleccionamos los geneSetID c(1, 43, 392) porque, en conjunto, capturan la fisiopatología 
# de la sarcopenia en la transición de 9 a 24 meses. Se eligieron estas 3 vías ortogonales para evitar redundancia:
# - ID 1 (extracellular matrix organization): Indica el colapso del andamiaje e integridad estructural del tejido.
# - ID 43 (muscle cell proliferation): Refleja la pérdida de la capacidad regenerativa de los mioblastos.
# - ID 392 (response to hypoxia): Muestra la incapacidad del músculo viejo para adaptarse al estrés celular y metabólico.
# Juntas, estas vías (con enriquecimiento negativo) cuentan la historia del declive funcional del músculo.

# Extraemos el nombre del pathway para el título
top_pathway_title <- gse$Description[1] 
p1 <- gseaplot2(gse, geneSetID = c("GO:0006260","GO:0099504", "GO:0016055"), title = "PATHWAYS UP-REGULATED Y DOWN-REGULATED", pvalue_table = TRUE)
"GO:0099504"
p1
# Guardar el plot de GSEA en tu carpeta
ggsave(
  filename = file.path(output_dir, "GSEA.png"),
  plot = p_hibrido,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white" 
)

# Generar el gráfico clásico para tu primera vía (ej. Replicación de ADN)
p_clasico_1 <- gseaplot(
  gse, 
  geneSetID = "GO:0060078", # Cambia por tu ID exacto
  by = "all",               # Esto asegura que salgan los 3 paneles (curva, líneas, ranking)
  title = "DNA replication (E11.5)"
)

print(p_clasico_1)


vias_seleccionadas <- c("GO:0006260", "GO:0099504", "GO:0016055", "GO:0060078") 

# Paleta limpia: Verde clásico de GSEA, Rojo suave y Azul claro
colores_limpios <- c("#1ab4a7", "#8c1782", "#4c822e", "#dca235")

# Generar el gráfico de 3 vías pulido
p_hibrido <- gseaplot2(
  gse, 
  geneSetID = vias_seleccionadas,
  color = colores_limpios,
  pvalue_table = FALSE,           # Quitamos la tabla para limpiar el fondo
  base_size = 20,                 # Letra un poco más grande
  rel_heights = c(1.5, 0.4, 0.6)  # Hacemos las líneas de "hits" un poco más delgadas para que respire el gráfico
)

print(p_hibrido)
