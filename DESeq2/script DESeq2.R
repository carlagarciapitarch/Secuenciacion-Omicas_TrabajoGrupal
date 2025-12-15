BiocManager::install("DESeq2")
BiocManager::install("tximport")
install.packages("ggrepel")
library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(ggrepel)

getwd()
samples <- data.frame(
  sample = c("BART", "LISA", "MAGGIE", 
             "PATTY", "SELMA", "MARGE"),
  condition = factor(c("normopeso", "normopeso", "normopeso", 
                       "obeso", "obeso", "obeso")),
  row.names = c("BART", "LISA", "MAGGIE", 
                "PATTY", "SELMA", "MARGE")
)

#definir las rutas de los archivos
files <- c(
  BART = "BARTquant.sf",
  LISA = "LISAquant.sf",
  MAGGIE = "MAGGIEquant.sf",
  PATTY = "PATTYquant.sf",
  SELMA = "SELMA.sf",
  MARGE = "MARGE.sf"
)
all(file.exists(files))

#hay que obtener el archivo tx2gene
BiocManager::install("biomaRt")
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
tx2gene <- getBM(attributes = c("refseq_mrna", "ensembl_gene_id", "external_gene_name"),
                 mart = ensembl)
tx2gene <- tx2gene[tx2gene$refseq_mrna != "", c("refseq_mrna", "external_gene_name")]
colnames(tx2gene) <- c("TXNAME", "GENEID")
write.csv(tx2gene, "tx2gene.csv", row.names = FALSE)
head(tx2gene)
#eliminar filas con genes vacíos
tx2gene <- tx2gene[tx2gene$GENEID != "", ]

#guardar archivo limpio
write.csv(tx2gene, "tx2gene.csv", row.names = FALSE)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                ignoreTxVersion = TRUE)
dds <- DESeqDataSetFromTximport(txi, colData = samples, 
                                design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#establecer el nivel de referencia (normopeso como control, para ver qué genes se expresan más en el grupo obeso)
dds$condition <- relevel(dds$condition, ref = "normopeso")

#ahora sí, se hace el análisis de la expresión diferencial
#ejecutar DESeq2
dds <- DESeq(dds)
#obtener resultados
res <- results(dds, contrast = c("condition", "obeso", "normopeso"))
#resumen de resultados
summary(res)
#ordenar por p-valor ajustado
resOrdered <- res[order(res$padj),]
#ver los genes más significativos
head(resOrdered, 20)

#voy a hacer un análisis de los resultados previo a la visualización de patrones
#número de genes diferencialmente expresados (p-valor ajustado < 0.05)
sum(res$padj < 0.05, na.rm = TRUE)

#genes sobreexpresados en obesos (log2FC > 1, padj < 0.05)
#si Log2FC > 1, el gen se expresa al menos 2 veces más en obesos
#si Log2FC < -1, el gen se expresa al menos 2 veces menos en obesos (o sea, más en normopeso)
genes_up <- sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm = TRUE)
print(paste("Genes sobreexpresados:", genes_up))

#genes subexpresados en obesos (log2FC < -1, padj < 0.05)
genes_down <- sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm = TRUE)
print(paste("Genes subexpresados:", genes_down))
#no hay ninguno

#se transforman los datos para la visualización
rld <- rlog(dds, blind = FALSE)

#hacemos volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                              ifelse(res_df$log2FoldChange > 1, "Sobreexpresado", "Subexpresado"),
                              "No significativo")

#volcano plot con etiquetas para los genes más significativos
top_genes_labels <- res_df[order(res_df$padj), ][1:20, ]  # Top 20 genes

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Sobreexpresado" = "purple", 
                                "Subexpresado" = "blue", 
                                "No significativo" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes_labels,
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 20,
                  box.padding = 0.5) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Expresión Diferencial",
       subtitle = "Genes sobreexpresados en el grupo Obeso 2 en referencia al grupo Normopeso",
       x = "Log2 Fold Change",
       y = "-Log10 P-valor ajustado",
       color = "Estado")

ggsave("volcano_plot.pdf", width = 12, height = 10)

#heatmap solo de los genes diferencialmente expresados
#primero se seleccionan los genes significativos (padj < 0.05)
sig_genes <- rownames(res[which(res$padj < 0.05),])
#heatmap de los top 50 genes más significativos
top_genes <- head(order(res$padj), 50)
mat <- assay(rld)[top_genes,]
mat <- mat - rowMeans(mat)  # Centrar por filas

#anotación de columnas
annotation_col <- data.frame(
  Condicion = samples$condition,
  row.names = rownames(samples)
)
ann_colors <- list(
  Condicion = c(obeso = "blue", normopeso = "purple")
)

#heatmap
pdf("heatmap_top50.pdf", width = 10, height = 12)
pheatmap(mat, 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "purple"))(100),
         main = "Top 50 Genes Diferencialmente Expresados")
dev.off()

#heatmap de todos los genes significativos
if(length(sig_genes) > 0 & length(sig_genes) <= 200) {
  mat_all <- assay(rld)[sig_genes,]
  mat_all <- mat_all - rowMeans(mat_all)
  
  pdf("heatmap_all_significant.pdf", width = 10, height = 14)
  pheatmap(mat_all, 
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           cluster_rows = TRUE, 
           cluster_cols = TRUE,
           show_rownames = FALSE,
           show_colnames = TRUE,
           color = colorRampPalette(c("blue", "white", "purple"))(100),
           main = "Genes Significativos (padj < 0.05)")
  dev.off()
}


#por buscar más patrones, hago un PCA
#para ver agrupamiento de muestras
pdf("PCA_plot.pdf", width = 8, height = 6)
plotPCA(rld, intgroup = "condition") +
  theme_minimal() +
  ggtitle("PCA - Obeso vs Normopeso")
dev.off()

#PCA personalizado
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = name), size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% varianza")) +
  ylab(paste0("PC2: ", percentVar[2], "% varianza")) +
  theme_minimal() +
  scale_color_manual(values = c("purple", "blue")) +
  ggtitle("Análisis de Componentes Principales")

ggsave("PCA_plot_custom.pdf", width = 10, height = 8)


#guardo resultados
#tabla completa de resultados
write.csv(as.data.frame(resOrdered), 
          file = "resultados_DESeq2_completos.csv")
#genes significativos
sig_results <- subset(res, padj < 0.05)
write.csv(as.data.frame(sig_results), 
          file = "genes_significativos.csv")
#genes sobreexpresados
up_regulated <- subset(res, padj < 0.05 & log2FoldChange > 1)
write.csv(as.data.frame(up_regulated), 
          file = "genes_sobreexpresados.csv")

#y como resumen final
cat("\n=== RESUMEN DEL ANÁLISIS ===\n")
cat("Total de genes analizados:", nrow(res), "\n")
cat("Genes significativos (p-valor ajustado < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("Genes sobreexpresados (FC > 2, p-valor ajustado < 0.05):", genes_up, "\n")
cat("Genes subexpresados (FC < 0.5, p-valor ajustado < 0.05):", genes_down, "\n")

