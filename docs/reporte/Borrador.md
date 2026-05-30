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

Para determinar la orientación de la librería, se evaluó el conteo con featureCounts usando -s 0, -s 1 y -s 2. La opción -s 0 produjo el mayor número de lecturas asignadas a genes en todas las muestras, mientras que las opciones stranded (-s 1 y -s 2) redujeron considerablemente el número de lecturas asignadas. Por lo tanto, las bibliotecas se trataron como no orientadas (unstranded) y se utilizó -s 0 para generar la matriz final de conteos.

v 1 -> sin --countReadPairs 
v test 

v 2 -> s 0
v 3 -> s 1
v 4 -> s 2
v 5 -> -t gene
v 6 -> no -B -C

Para evaluar si los parámetros estrictos de conteo paired-end estaban causando la baja asignación en algunas muestras, se comparó featureCounts con y sin -B y -C. La eliminación de estos parámetros produjo cambios mínimos en el número de fragmentos asignados, por lo que se concluyó que la baja asignación no se debía a filtros de pares alineados o quiméricos. En cambio, la prueba con -t gene aumentó notablemente el número de fragmentos asignados en las muestras problemáticas, lo que sugiere que una proporción importante de sus lecturas cae en regiones génicas no exónicas. Por ello, se mantuvo el conteo estándar por exones agrupados por gen (-t exon -g gene_id) para la matriz final.v

Sí. Para tu proyecto actual, puedes justificar el uso de **DESeq2** sobre **edgeR** de una forma bastante sólida, especialmente porque el propio artículo original usó **DESeq2** para el análisis diferencial de los datos RNA-seq. En el paper de placenta, después de alinear con STAR y contar con featureCounts, los autores indican que usaron **DESeq2 sobre conteos crudos** para realizar el análisis de expresión diferencial, incorporando además la variable de paciente en el diseño para manejar la variabilidad entre muestras pareadas. 

Podrías escribir algo así:

> Para el análisis de expresión diferencial se decidió utilizar **DESeq2**, ya que es una herramienta ampliamente empleada para datos de RNA-seq basados en conteos crudos y modela la expresión génica mediante una distribución binomial negativa, lo cual permite considerar la variabilidad biológica entre réplicas. Además, DESeq2 estima factores de normalización para corregir diferencias en profundidad de secuenciación y aplica una estimación moderada de la dispersión, lo que resulta útil cuando se trabaja con un número limitado de muestras, como ocurre en este experimento. Esta elección también permite mantener coherencia metodológica con el estudio original, en el cual los autores procesaron los datos con fastp, alinearon contra GRCh38/hg38 usando STAR, cuantificaron con featureCounts y realizaron el análisis diferencial con DESeq2 sobre conteos crudos. En comparación con edgeR, que también es una alternativa válida basada en binomial negativa, DESeq2 ofrece un flujo más directo para normalización, transformación de varianza, visualización mediante PCA y extracción de contrastes, por lo que fue seleccionado como método principal para identificar genes diferencialmente expresados.

Y si quieres mencionar tu experiencia previa comparando ambos métodos, puedes agregar:

> En análisis previos se observó que tanto DESeq2 como edgeR conservaron la estructura biológica principal de los datos; sin embargo, DESeq2 tendió a mostrar agrupamientos ligeramente más compactos en PCA y una detección más amplia de genes diferencialmente expresados, mientras que edgeR se comportó de forma más selectiva. Por ello, edgeR puede utilizarse como método complementario o de validación, pero DESeq2 se eligió como enfoque principal por su consistencia con el artículo original, su facilidad de integración con matrices de conteo generadas por featureCounts y su desempeño robusto en diseños con pocas réplicas. En tu reporte anterior, justamente comparaste ambos enfoques y describiste que DESeq2 y edgeR mantuvieron la misma estructura general, aunque DESeq2 mostró clusters ligeramente más compactos y edgeR mayor dispersión. 

Una versión más corta para tu metodología sería:

> Se utilizó **DESeq2** para el análisis de expresión diferencial porque permite trabajar directamente con conteos crudos generados por featureCounts, modela la variabilidad de RNA-seq mediante una distribución binomial negativa y realiza una estimación moderada de la dispersión, lo cual es adecuado para experimentos con un número limitado de réplicas. Además, su uso mantiene la coherencia con el estudio original, donde los autores emplearon STAR, featureCounts y DESeq2 para identificar genes diferencialmente expresados en explantes placentarios humanos expuestos a distintas concentraciones de glucosa. Aunque edgeR también es una herramienta válida, DESeq2 fue seleccionado como método principal por su flujo integrado de normalización, transformación de varianza, PCA y análisis diferencial.

| SRR | Assigned | Unassigned | % Assigned | % Unassigned |
|---|---:|---:|---:|---:|
| SRR32154135 | 3072222 | 16401075 | 15.78% | 84.22% |
| SRR32154136 | 3211620 | 16876600 | 15.99% | 84.01% |
| SRR32154137 | 18653790 | 9445143 | 66.39% | 33.61% |
| SRR32154138 | 13676099 | 9595659 | 58.77% | 41.23% |
| SRR32154139 | 16671629 | 8431921 | 66.41% | 33.59% |
| SRR32154150 | 4313937 | 21480086 | 16.72% | 83.28% |
| SRR32154151 | 4239563 | 19166354 | 18.11% | 81.89% |
| SRR32154152 | 16583791 | 12469062 | 57.08% | 42.92% |
| SRR32154153 | 14557038 | 10506392 | 58.08% | 41.92% |
| SRR32154154 | 16921140 | 8188579 | 67.39% | 32.61% |
| Promedio |  -  |  -  | 44.07% | 55.93% |
