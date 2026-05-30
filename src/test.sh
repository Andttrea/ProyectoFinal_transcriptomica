#!/bin/bash 

# Asignamos la variable para la ruta del archivo GTF y psara la ruta de los archivos BAM generados por STAR
GTF_FILE="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/GENCODE/gencode.vM10.chr_patch_hapl_scaff.annotation.gtf"
BAM="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/"

# Iniiciamos el featureCounts
#Flags: 
# -T: número de hilos a usar, en este caso 8 para optimizar el proceso
# -a: especifica el archivo de anotación GTF o GFF, que es necesario para asignar las lecturas a las características genómicas correctas
# --largestOverlap: asigna una lectura a la característica con la que tiene la mayor superposición
# -p: indica que los datos son de secuenciación paired-end, lo que permite a featureCounts manejar correctamente las lecturas emparejadas
# -t exon: especifica que solo se contarán las lecturas que se alineen a los exones, lo que es común en análisis de expresión génica
# -g gene_id: indica que las lecturas se contarán a nivel de gen,
# -s 0: especifica que los datos no son strand-specific
# --countReadPairs: cuenta las lecturas emparejadas como una sola unidad
# -B: requiere que ambas lecturas de un par estén correctamente alineadas para ser contadas. mejora la precisión 
# -C: evita contar las lecturas que elimina pares quiméricos o inconsistentes.
# "$BAM": inidiica el archivo BAM de entrada que se va a procesar
featureCounts -o "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/v3/counts_matrix.txt" \
    -T 10 \
    -a "$GTF_FILE" \
    --largestOverlap \
    -t exon \
    -g gene_id \
    -s 2 \
    --countReadPairs \
    -p \
    -C \
    -B \
    "$BAM"/*Aligned.sortedByCoord.out.bam
