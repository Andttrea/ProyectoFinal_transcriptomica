| SRR | Number of input reads | Uniquely mapped reads % | % of reads mapped to multiple loci | % of reads unmapped: too short | Mismatch rate per base |
|---|---:|---:|---:|---:|---:|
| SRR27534519 | 29387207 | 92.18% | 5.59% | 1.94% | 0.14% |
| SRR27534520 | 35961680 | 92.11% | 5.69% | 1.93% | 0.14% |
| SRR27534521 | 31749298 | 92.19% | 5.76% | 1.79% | 0.12% |
| SRR27534522 | 27247824 | 92.78% | 5.70% | 1.21% | 0.18% |
| SRR27534523 | 28735943 | 92.73% | 5.77% | 1.19% | 0.18% |
| SRR27534524 | 23008311 | 92.46% | 6.02% | 1.23% | 0.18% |
| SRR27534525 | 28274148 | 91.13% | 6.72% | 1.84% | 0.17% |
| SRR27534526 | 26699471 | 91.46% | 6.59% | 1.66% | 0.17% |
| SRR27534527 | 28399978 | 90.91% | 6.85% | 1.98% | 0.19% |
| SRR27534528 | 29280831 | 92.00% | 6.43% | 1.28% | 0.14% |
| SRR27534529 | 30838207 | 92.13% | 6.36% | 1.22% | 0.13% |
| SRR27534530 | 28290961 | 91.98% | 6.49% | 1.25% | 0.13% |

Los archivos BAM generados por STAR fueron producidos directamente en formato ordenado por coordenada mediante la opción SortedByCoordinate, por lo que no fue necesario aplicar un ordenamiento adicional con samtools sort. Posteriormente, los archivos fueron indexados con samtools index para generar archivos .bai, necesarios para el acceso eficiente a los alineamientos y su posible visualización en herramientas como IGV.c

En el módulo de niveles de duplicación, las muestras no presentaron advertencias en FastQC. Aunque en RNA-seq es esperable observar cierta proporción de lecturas duplicadas debido a la alta abundancia de algunos transcritos, los resultados no sugieren una duplicación excesiva a nivel de control de calidad inicial. Además, debido a que FastQC evalúa duplicación de secuencias antes del alineamiento, este módulo no permite distinguir completamente entre duplicación técnica por PCR y duplicación biológica asociada a genes altamente expresados.

La cuantificación de lecturas por gen se realizó con featureCounts utilizando la anotación en formato GTF. Se especificó -t exon para contar lecturas asignadas a regiones exónicas y -g gene_id para agrupar los conteos a nivel de gen. Debido a que los datos son paired-end, se usaron -p y --countReadPairs para contar fragmentos en lugar de lecturas individuales. Además, se empleó -s 2 porque las bibliotecas fueron preparadas con TruSeq Stranded mRNA, protocolo direccional de orientación reversa. La opción -B permitió contar únicamente fragmentos en los que ambas lecturas del par estuvieran alineadas correctamente.

La cuantificación con featureCounts mostró una alta proporción de fragmentos asignados a genes en todas las muestras, con valores entre aproximadamente 18.9 y 30.1 millones de fragmentos asignados. Las principales categorías de fragmentos no asignados correspondieron a lecturas multimapeadas y lecturas alineadas fuera de regiones anotadas como exónicas. La categoría de ambigüedad fue baja en comparación con los fragmentos asignados, lo que sugiere una asignación mayoritariamente clara a nivel de gen. Debido a que las lecturas multimapeadas pueden introducir incertidumbre en la cuantificación, se mantuvo un enfoque conservador y no se incluyeron en los conteos finales para el análisis de expresión diferencial.

| SRR | Genes asignados | Genes no asignados | % asignado | % no asignado |
|---|---:|---:|---:|---:|
| SRR27534519 | 24251569 | 6867489 | 77.93% | 22.07% |
| SRR27534520 | 30111086 | 8058471 | 78.89% | 21.11% |
| SRR27534521 | 26687449 | 7078947 | 79.04% | 20.96% |
| SRR27534522 | 22193158 | 6862381 | 76.38% | 23.62% |
| SRR27534523 | 23414995 | 7288006 | 76.26% | 23.74% |
| SRR27534524 | 18957642 | 5692174 | 76.91% | 23.09% |
| SRR27534525 | 23581676 | 6909347 | 77.34% | 22.66% |
| SRR27534526 | 22404491 | 6385784 | 77.82% | 22.18% |
| SRR27534527 | 23753909 | 6872943 | 77.56% | 22.44% |
| SRR27534528 | 24921983 | 6706367 | 78.80% | 21.20% |
| SRR27534529 | 26424602 | 6866913 | 79.37% | 20.63% |
| SRR27534530 | 24124433 | 6460700 | 78.88% | 21.12% |

Se evaluó el efecto de utilizar la opción --largestOverlap en featureCounts. Esta opción incrementó ligeramente el número de fragmentos asignados y redujo la categoría Unassigned_Ambiguity; sin embargo, el cambio fue mínimo en relación con el total de fragmentos asignados. Por esta razón, se decidió emplear una cuantificación conservadora sin --largestOverlap, evitando forzar la asignación de lecturas que podían solaparse con más de una característica genómica.

| SRR | Genes asignados | Genes no asignados | % asignado | % no asignado |
|---|---:|---:|---:|---:|
| SRR27534519 | 24306959 | 6812099 | 78.11% | 21.89% |
| SRR27534520 | 30179155 | 7990402 | 79.07% | 20.93% |
| SRR27534521 | 26752451 | 7013945 | 79.23% | 20.77% |
| SRR27534522 | 22240434 | 6815105 | 76.54% | 23.46% |
| SRR27534523 | 23466631 | 7236370 | 76.43% | 23.57% |
| SRR27534524 | 18999553 | 5650263 | 77.08% | 22.92% |
| SRR27534525 | 23633764 | 6857259 | 77.51% | 22.49% |
| SRR27534526 | 22453467 | 6336808 | 77.99% | 22.01% |
| SRR27534527 | 23806569 | 6820283 | 77.73% | 22.27% |
| SRR27534528 | 24981809 | 6646541 | 78.99% | 21.01% |
| SRR27534529 | 26488182 | 6803333 | 79.56% | 20.44% |
| SRR27534530 | 24183959 | 6401174 | 79.07% | 20.93% |


