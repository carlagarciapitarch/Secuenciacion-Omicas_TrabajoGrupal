#instalacion de edgeR

BiocManager::install("edgeR")
BiocManager::install("tximport")
library(edgeR)
library(tximport)
setwd("C:/Users/alexi/Desktop/MASTER/MASTER/omica/ACT2/quantsfs")

#Generación de la matriz de conteos

#CARGA DEL ARCHIVO DE MAPEO 
tx2gene <- read.table("Transcrito_a_Gen.tsv", header = FALSE, sep = "\t") 
colnames(tx2gene) <- c("TXNAME", "GENEID")

# 3. DEFINICIÓN DE MUESTRAS GRUPOS Y ARCHIVOS DE ENTRADA

sampleName <- c("Marge", "Patty", "Selma", "Bart", "Lisa", "Maggie")
group <- factor(c("Obeso2", "Obeso2", "Obeso2", "Normopeso", "Normopeso", "Normopeso"))

# Lista de los 6 archivos .sf
files <- c(
  "Marge.sf", 
  "Patty.sf", 
  "Selma.sf", 
  "Bart.sf", 
  "Lisa.sf", 
  "Maggie.sf"
)

# Creación del data.frame de muestras
samples <- data.frame(sampleName, group, files)
colnames(samples)[2] <- "condition" 

# Asignar los nombres de las muestras a la lista de archivos
names(files) <- samples$sampleName


#IMPORTACIÓNDE DATOS CON tx2gene

txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM") 

#OBTENCIÓN DE LA MATRIZ DE CONTEOS BRUTOS (EdgeR)

# EdgeR usa números enteros.
counts_matrix <- round(txi$counts)

# Asignar los nombres de las muestras a las columnas
colnames(counts_matrix) <- samples$sampleName

###################################
#Ahora usamos EdgeR para analizar y comparar Obeso 2 vs Normopeso
###################################

#CREACIÓN DEL OBJETO DGEList
dge <- DGEList(counts = counts_matrix, group = group)

#FILTRADO Y NORMALIZACIÓN

#Eliminamos genes con baja expresión 

keep <- filterByExpr(dge) 
dge <- dge[keep, , keep.lib.sizes=FALSE]
print(paste("Número de genes restantes después del filtrado:", dim(dge)[1]))

#Normalización TMM 
dge <- calcNormFactors(dge) 

#Diseño de la matriz del modelo 
design <- model.matrix(~group)

#Estimar la dispersión (variabilidad biológica y técnica)
dge <- estimateDisp(dge, design) 
fit <- glmQLFit(dge, design)

#PRUEBA DE EXPRESIÓN DIFERENCIAL 

# El contraste es Obeso2 (coeficiente 2) vs Normopeso (coeficiente 1, la referencia)
# El resultado de la prueba QLF indica qué genes cambiaron significativamente.
qlf <- glmQLFTest(fit, coef = 2) 

#OBTENCIÓN DE LA TABLA DE RESULTADOS

# Obtener la tabla de resultados completa (incluyendo logFC, PValue y FDR)
results_table <- topTags(qlf, n = Inf)$table

#IDENTIFICACIÓN Y EXPORTACIÓN DE GENES DIFERENCIALMENTE EXPRESADOS 

# Definir umbrales de significancia y Fold Change 
FDR_threshold <- 0.05
logFC_threshold <- 1 

# Identificar los Genes Diferencialmente Expresados
degs <- results_table[results_table$FDR < FDR_threshold & abs(results_table$logFC) > logFC_threshold, ]

#EXPORTACIÓN PARA TAREA DE VISUALIZACIÓN Y ENRIQUECIMIENTO
write.csv(results_table, "Obeso2_vs_Normopeso_Results_Completa.csv")
write.csv(degs, "Obeso2_vs_Normopeso_DEGs_Lista.csv")

print(paste("Se identificaron", dim(degs)[1], "genes diferencialmente expresados (FDR < 0.05 y |logFC| > 1)."))

##########
#Visualización
##########
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)
#GENERACIÓN DEL VOLCANO PLOT

titulo_plot <- "Expresión Diferencial: Obeso 2 vs Normopeso"
logFC_threshold <- 1 
FDR_threshold <- 0.05 

EnhancedVolcano(
  results_table,
  lab = rownames(results_table), 
  x = 'logFC',
  y = 'FDR',
  title = titulo_plot,
  pCutoff = FDR_threshold,
  FCcutoff = logFC_threshold,
  pointSize = 1.5,
  labSize = 3.0,
  colAlpha = 0.5,
  legendPosition = 'right',
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  colConnectors = 'grey30',
  # Personalizar colores
  col = c('grey30', 'grey30', 'grey30', 'red3'), 
  caption = paste0('FC cutoff: ', logFC_threshold, ' | FDR cutoff: ', FDR_threshold)
)

#HEATMAP

#Obtener los conteos normalizados (logaritmo de CPM)
logCPM <- cpm(dge, prior.count=2, log=TRUE)

# 2.2. Seleccionar solo Genes Diferencialmente Expresados 
degs_genes <- rownames(degs) 
logCPM_degs <- logCPM[degs_genes, ]
logCPM_scaled <- t(scale(t(logCPM_degs)))

# 3. CREACIÓN DEL HEATMAP

annotation_col <- data.frame(Group = factor(samples$condition))
rownames(annotation_col) <- samples$sampleName

pheatmap(
  logCPM_scaled, 
  cluster_rows = TRUE, 
  show_rownames = TRUE, # Si hay muchos genes, oculte los nombres de fila
  cluster_cols = TRUE, 
  annotation_col = annotation_col,
  main = "Heatmap de Genes Diferencialmente Expresados (Obeso 2 vs Normopeso)",
  fontsize_row = 5,
  fontsize_col = 10,
  filename = "Heatmap_DEGs_Obeso2_vs_Normopeso.png" 
)
