#!/bin/bash
# Extraemos el nombre de la base de la muestra
for R1 in /export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/raw/*_1.fastq; do
base=$(basename "$R1" _1.fastq)
# Definimos el archivo Read2 correspondiente
R2="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/raw/${base}_2.fastq"
# Definimos los archivos de salida
out_R1="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/trimmed/cutadapt/${base}_1_trimmed.fastq"
out_R2="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/trimmed/cutadapt/${base}_2_trimmed.fastq"

ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" 

# Ejecutamos fastp
echo "Ejecutando cutadapt para $base"
# Flags:
# -j : número de hilos a utilizar para el procesamiento, en este caso 10 para acelerar el proceso de recorte
# -a : secuencia del adaptador a recortar en las lecturas del primer par, en este caso la secuencia del adaptador para R1
# -A : secuencia del adaptador a recortar en las lecturas del segundo par, en este caso la secuencia del adaptador para R2
# -a "A{15}" : recorta secuencias de adenina de longitud 15 o más en las lecturas del primer par
# -A "A{15}" : recorta secuencias de adenina de longitud 15 o más en las lecturas del segundo par
# -a "G{15}" : recorta secuencias de guanina de longitud 15 o más en las lecturas del primer par
# -A "G{15}" : recorta secuencias de guanina de longitud 15 o más en las lecturas del segundo par
# -m : longitud mínima de las lecturas después del recorte, en este caso 50 para eliminar secuencias muy cortas que podrían introducir ruido en el alineamiento
# -u : número de bases a recortar desde el inicio de las lecturas del primer par, en este caso 12 para ser laxos y evitar problemas de calidad al inicio de las lecturas
# -U : número de bases a recortar desde el inicio de las lecturas del segundo par, en este caso 12 para ser laxos y evitar problemas de calidad al inicio de las
# -o : archivo de salida para las lecturas del primer par, en este caso el archivo out
# -p : archivo de salida para las lecturas del segundo par, en este caso el archivo out_R2
# "$R1" "$R2" : archivos de entrada para las lecturas del primer y segundo par, en este caso los archivos R1 y R2 que corresponden a la muestra actual
cutadapt \
-j 12 \
-a "$ADAPTER_R1" \
-A "$ADAPTER_R2" \
-a "A{15}" \
-A "A{15}" \
-a "G{15}" \
-A "G{15}" \
-u 12 -U 12 \
-m 50 \
-o "$out_R1" \
-p "$out_R2" \
"$R1" "$R2"  > "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/trimmed/cutadapt/${base}_cutadapt.log"

done
echo "Limpieza completada para todas las muestras"
