#!/bin/bash 

# Archivo de anotación GTF
GTF_FILE="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/GENCODE/gencode.v49.chr_patch_hapl_scaff.annotation.gtf"

# Carpeta donde están los BAM generados por STAR
BAM_DIR="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star"


# Crear carpeta de salida
mkdir -p "$OUT_DIR"

echo "Iniciando conteo con featureCounts..."

featureCounts -o "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/v6/counts_matrix_v6.txt" \
    -T 10 \
    -a "$GTF_FILE" \
    -t exon \
    -g gene_id \
    -s 0 \
    --largestOverlap \
    --countReadPairs \
    -p \
    "$BAM_DIR"/*Aligned.sortedByCoord.out.bam

echo "Conteo completado"