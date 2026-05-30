#!/bin/bash
# Asignamos las variables para nuestro indice de referencia
index="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/GENCODE/index/"
Star_time_file="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/star_time.txt"

# Iniciamos nuestro ciclo
for R1 in /export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/trimmed/*_1_trimmed.fastq; do
# Asignamos el nombre de la muestra a partir del nombre del archivo
    base=$(basename $R1 _1_trimmed.fastq)
    # Definimos donde se encuentran los archivos R2
    R2="/export/storage/users/andreavg/ProyectoFinal_transcriptomica/data/trimmed/${base}_2_trimmed.fastq"

    echo "Procesando muestra: $base"
    #Flags:
    # --runThreadN : número de hilos a usar, en este caso usaremos 8 para
    #manejar la carga de trabajo de manera eficiente
    # --genomeDir : especifica el directorio del índice de referencia, guardado
    #en la variable "index"
    # --readFilesIn : es el archivo de entrada para las lecturas del primer par
    #y el segundo par, en este caso los archivos R1 y R2 que corresponden a la
    #muestra actual
    # --readFilesCommand : especifica el comando para descomprimir los archivos
    #de entrada, en este caso usamos "zcat" para manejar archivos comprimidos
    # --outSAMunmapped None : evita que se escriban las lecturas no alineadas en
    #el archivo de salida, para reducir el tamaño del archivo resultante
    # --outSAMtype SAM : especifica el formato de salida para el alineamiento,
    #en este caso SAM para mantener la consistencia con los resultados de hisat2
    # --outFileName Prefix : especifica el prefijo para los archivos de salida
    # time se utiliza para medir el tiempo de ejecución del comando STAR y se
    #redirige la salida de error estándar, que es donde se imprime el tiempo, al
    #archivo correspondiente para single-end.
    # usando las {} se asegura que el tiempo se mida solo para el comando STAR y
    #no para otros comandos que puedan estar en el script.
    { time STAR --runThreadN 14 \
    --genomeDir "$index" \
    --readFilesIn "$R1" "$R2" \
    --outSAMunmapped None \
    --outReadsUnmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "/export/storage/users/andreavg/ProyectoFinal_transcriptomica/results/star/${base}" ; } 2>> "$Star_time_file"

done
echo "Alineamiento y indexación completados para todas las muestras"