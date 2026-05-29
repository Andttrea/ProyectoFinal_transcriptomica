Number of input reads\|Uniquely mapped reads %\|% of reads mapped to multiple loci\|% of reads unmapped: too short\|Mismatch rate per base

| SRR | Number of input reads | Uniquely mapped reads % | % of reads mapped to multiple loci | % of reads unmapped: too short | Mismatch rate per base |
|---|---:|---:|---:|---:|---:|
| SRR32154135 | 18137105 | 84.46% | 7.45% | 7.74% | 0.43% |
| SRR32154136 | 19300857 | 82.67% | 7.09% | 9.91% | 0.49% |
| SRR32154137 | 24275773 | 86.10% | 8.70% | 5.09% | 0.25% |
| SRR32154138 | 20734715 | 85.34% | 8.71% | 5.80% | 0.26% |
| SRR32154139 | 21662597 | 85.51% | 9.56% | 4.72% | 0.22% |
| SRR32154150 | 23500033 | 84.19% | 7.99% | 7.43% | 0.40% |
| SRR32154151 | 21893699 | 83.74% | 7.33% | 8.58% | 0.42% |
| SRR32154152 | 24594715 | 80.73% | 9.19% | 9.81% | 0.31% |
| SRR32154153 | 22409418 | 85.67% | 8.62% | 5.57% | 0.27% |
| SRR32154154 | 21874899 | 86.54% | 9.07% | 4.16% | 0.21% |

No se aplicó un filtrado adicional por MAPQ > 10 debido a que las métricas de alineamiento obtenidas con STAR indicaron una calidad adecuada y consistente entre las muestras. En todas las bibliotecas, el porcentaje de lecturas mapeadas de manera única fue alto, entre 80.73% y 86.54%, mientras que la proporción de lecturas multimapeadas se mantuvo en un rango moderado, entre 7.09% y 9.56%. Además, la tasa de mismatch por base fue baja en todas las muestras, entre 0.21% y 0.49%, lo que sugiere una buena concordancia entre las lecturas y el genoma de referencia. Por lo tanto, filtrar por MAPQ podría eliminar lecturas potencialmente informativas sin aportar una mejora sustancial al conjunto de datos. En consecuencia, se decidió conservar los BAM generados por STAR sin filtrado adicional y continuar con el conteo de lecturas por gen, controlando posteriormente la asignación de lecturas ambiguas durante el análisis de cuantificación.