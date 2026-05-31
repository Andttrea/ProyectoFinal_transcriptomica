library(dplyr)
library(ggplot2)
library(stringr)

generar_go_dotplot <- function(file_path, output_path, plot_title, top_n = 15) {
  
  # Leer archivo
  datos_go <- read.csv(file_path, header = TRUE)
  
  # Procesar datos
  datos_plot <- datos_go %>%
    mutate(
      log10_fdr = -log10(FDR),
      Term = str_replace(Term, "GOTERM_BP_DIRECT~", ""),
      Term = str_replace(Term, "GOTERM_CC_DIRECT~", ""),
      Term = str_wrap(Term, width = 40)
    ) %>%
    arrange(FDR) %>%
    slice_head(n = top_n)
  
  # Hacer gráfica
  p <- ggplot(datos_plot, aes(x = Fold.Enrichment, y = reorder(Term, Fold.Enrichment))) +
    geom_point(aes(size = Count, color = log10_fdr)) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = plot_title,
      x = "Fold Enrichment",
      y = "Términos de GO",
      size = "Conteo de genes",
      color = "-log10(FDR)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(face = "bold")
    )
  
  # Guardar
  ggsave(output_path, plot = p, width = 9, height = 7, dpi = 300)
  
  return(p)
}

generar_go_dotplot(
  file_path = file.path(base, "data", "DAVIDChartReport_downregulated_BP.csv"),
  output_path = file.path(base, "GO_BP_down.png"),
  plot_title = "PLOT DE PROCESOS CELULARES EN GENES DOWN"
)

generar_go_dotplot(
  file_path = file.path(base, "data", "DAVIDChartReport_upregulated_BP.csv"),
  output_path = file.path(base, "GO_BP_up.png"),
  plot_title = "PLOT DE PROCESOS CELULARES EN GENES UP"
)

generar_go_dotplot(
  file_path = file.path(base, "data", "DAVIDChartReport_downregulated_CC.csv"),
  output_path = file.path(base, "GO_CC_down.png"),
  plot_title = "PLOT DE COMPONENTES CELULARES DOWNREGULADOS"
)

generar_go_dotplot( file_path = file.path(base, "data", "DAVIDChartReport_upregulated_CC.csv"), output_path = file.path(base, "GO_CC_up.png"), plot_title = "PLOT DE COMPONENTES CELULARES UPREGULADOS" )

