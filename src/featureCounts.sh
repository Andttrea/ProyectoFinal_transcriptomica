#!/bin/bash 

# Asignamos la variable para la ruta del archivo GTF y psara la ruta de los archivos BAM generados por STAR
GTF_FILE="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/GENCODE/gencode.v49.chr_patch_hapl_scaff.annotation.gtf"
BAM="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/"

# Iniiciamos el featureCounts
#Flags: 
# -T: número de hilos a usar, en este caso 8 para optimizar el proceso
# -a: especifica el archivo de anotación GTF o GFF, que es necesario para asignar las lecturas a las características genómicas correctas
# --largestOverlap: asigna una lectura a la característica con la que tiene la mayor superposición
# -p: indica que los datos son de secuenciación paired-end, lo que permite a featureCounts manejar correctamente las lecturas emparejadas
# -B: requiere que ambas lecturas de un par estén correctamente alineadas para ser contadas. mejora la precisión 
# -C: evita contar las lecturas que se alinean a múltiples ubicaciones en el genoma, lo que ayuda a reducir el ruido en los datos de conteo
# "$BAM": inidiica el archivo BAM de entrada que se va a procesar
featureCounts -o "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/featurecounts/counts_matrix.txt" \
-T 10 \
-a "$GTF_FILE" \
--largestOverlap \
-p \
-C \
-B \
"$BAM"/*Aligned.sortedByCoord.out.bam
