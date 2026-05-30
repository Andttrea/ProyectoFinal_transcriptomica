#!/bin/bash
# Extraemos el nombre de la base de la muestra
for R1 in /export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/fastq/*_1.fastq; do
base=$(basename "$R1" _1.fastq)
# Definimos el archivo Read2 correspondiente
R2="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/fastq/${base}_2.fastq"
# Definimos los archivos de salida
out_R1="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/trimmed/${base}_1_trimmed.fastq"
out_R2="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/trimmed/${base}_2_trimmed.fastq"
echo "Procesando $base"
# Ejecutamos fastp
echo "Ejecutando fastp para $base"
# Flags:
# -i : archivo de entrada para las lecturas del primer par, en este caso el archivo R1 que corresponde a la muestra actual# -I : archivo de entrada para las lecturas del segundo par, en este caso el
#archivo R2 que corresponde a la muestra actual
# -o : archivo de salida para las lecturas del primer par, en este caso el
#archivo out_R1 que corresponde a la muestra actual
# -O : archivo de salida para las lecturas del segundo par, en este caso el
#archivo out_R2 que corresponde a la muestra actual
# -f : número de bases a recortar desde el inicio de las lecturas, en este
#caso 12 para ser laxos y evitar problemas de calidad al inicio de las
#lecturas
# -F : número de bases a recortar desde el inicio de las lecturas del
#segundo par, en este caso 12 para ser laxos y evitar problemas de calidad al
#inicio de las lecturas
# el cúal observamos en la calidad de las lecturas.
# -l : longitud mínima de las lecturas después del recorte, en este caso 50 para eliminar secuencias muy cortas que podrían introducir ruido en el alineamiento
# -w = 10 : número de hilos a utilizar para el procesamiento, en este caso 10 para acelerar el proceso de recorte
# -h : archivo de salida para el reporte HTML generado por fastp.
# --n_base_limit : número máximo de bases ambiguas (N) permitidas en una lectura después del recorte, en este caso 5 para eliminar lecturas con demasiadas bases ambiguas que podrían afectar la calidad del alineamiento
fastp -i "$R1" -I "$R2" \
-o "$out_R1" -O "$out_R2" \
-f 12 -F 12 \
-l 50 \
--n_base_limit 5 \
-w 12 \
-h "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/trimmed/${base}_fastp_report.html"

done
echo "Limpieza completada para todas las muestras"