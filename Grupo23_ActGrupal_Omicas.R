# Script para resolución actividad grupal Omicas y Secuenciación

setwd("~/Documentos/MasterBioinformatica/1Q/Secuenciacion_Omicas/Actividades/Actividad grupal/Resolucion_Act_Grupal/Documentos_R_ActGrupal_Omicas")

# Cargamos el diseño experimental y la traducción de los transcritos a genes
f <- list.files(pattern="Design.csv", ignore.case=TRUE)[1]
sample_info <- read.csv(f, fileEncoding="UTF-8-BOM", check.names=FALSE)

tx2gene <- read.delim("Transcrito_a_Gen.tsv", fileEncoding="UTF-8-BOM", check.names=FALSE)

colnames(tx2gene) = c("TXNAME", "GENEID")

# Definimos rutas a los ficheros de cuantificación de Salmon
files = file.path("Salmon", sample_info$Sample, "quant.sf")
names(files) = sample_info$Sample

files

BiocManager::install("tximport")
library(tximport)

# Leemos los datos de expresión con tximport
txi = tximport(files, type = "salmon", tx2gene = tx2gene)

# Exploración de la matriz
txi$counts

## ESTO LO PUSO EL PROFESOR ASI PORQUE LOS DOCUMENTOS ESTAN DENTRO DE DETERMINADAS CARPETAS, PERO TAL VEZ ES MAS FACIL IMPORTARLO UNO A UNO Y YA ESTA