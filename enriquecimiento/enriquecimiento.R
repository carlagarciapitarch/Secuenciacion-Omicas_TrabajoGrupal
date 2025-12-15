#Establecer el directorio de trabajo y la instalacion de paquetes 
setwd("C:/Users/angel/Desktop/Master/Secuenciación y Ómicas de Próxima Generación/actividad 2/mubio03_act2/TallerGrupal_Ficheros/Salmon")
BiocManager::install("clusterProfiler",update=TRUE,ask=FALSE)
BiocManager::install("org.Hs.eg.db",update=TRUE,ask=FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

#---------PASO 1:  Preparacion y Mapeo de la Lista de Genes----------------
#Definir la lista de genes expresados por Símbolo.Estos genes provienen del analisis de expresion diferencial( con EdgeR)
genes_significativos<- c("SIM1","UBR2","MC4R","FTO","NTRK2")

#Pasar de Simbolo a ENTREz_IDs ya que el analisis de enriquecimiento requiere IDs de Entrez

entrez_ids<- mapIds(x=org.Hs.eg.db,
                    keys=genes_significativos,
                    column="ENTREZID",
                    keytype= "SYMBOL",
                    multiVals="first")

entrez_ids<- na.omit(entrez_ids)
cat("Entrez IDs listos para el análisis:", entrez_ids, "\n")


#----------------PASO 2: Analisis de Enriquecimiento Funcional (GO)-------------
# - ont = "ALL": Analiza los tres dominios de la ontología (BP: Proceso Biológico, CC: Componente Celular, MF: Función Molecular).
# - pAdjustMethod = "BH": Utiliza el método Benjamini-Hochberg (BH) para corregir el p-valor (FDR).
# - pvalueCutoff = 0.05: Solo se consideran significativos los términos con p-valor ajustado (FDR) menor a 0.05.
# - readable = TRUE: Convierte los IDs de Entrez de los genes enriquecidos a sus Símbolos originales en la tabla de resultados.

ego <- enrichGO(gene = entrez_ids, 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENTREZID",
                ont = "ALL", 
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)

#Reporta el numero de terminos GO que cumplieron con los umbrales de significancia.
resultados_df <- as.data.frame(ego)
print(paste("Número de términos GO significativos encontrados:", nrow(resultados_df)))


#---------------PASO 3: Analisis de Enriquecimiento de Rutas(KEGG)------------
# Función enrichKEGG: Realiza el análisis de enriquecimiento para rutas de la base de datos KEGG.
# - organism = 'hsa': Especifica el organismo (hsa = Homo sapiens).
# - pvalueCutoff y qvalueCutoff: Umbrales de significancia para p-valor y q-valor (FDR).
ekegg <- enrichKEGG(gene = entrez_ids, 
                    organism = 'hsa', 
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)

#Reporta el conteo de rutas encontradas 
resultados_kegg_df <- as.data.frame(ekegg)
print(paste("Número de rutas KEGG significativas encontradas:", nrow(resultados_kegg_df)))
#Dado que la lista de rutas significativas es pequeña es comun no encontrar rutas significativas por lo que 
#se repite el analisis KEGG con un umbral de p-valor mas bajo que permite visualizar las rutas mas relevantes
#para estos genes.
ekegg_relajado <- enrichKEGG(gene = entrez_ids, 
                             organism = 'hsa', 
                             pvalueCutoff = 1.0, 
                             qvalueCutoff = 1.0)


#-------------PASO 4: Visualizacion --------------------
#Grafico KEGG

dotplot(ekegg_relajado, showCategory = 5) + 
  theme_bw() + 
  labs(
    title = "Rutas KEGG más enriquecidas (Umbral relajado)",
    x = "Proporción Génica (Gene Ratio)",
    colour = "p.ajustado",
    size = "Conteo de Genes"
  ) +
  theme(
    axis.text.y = element_text(size = 10), 
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right" 
  )

#Grafico GO

dotplot(ego, showCategory = 10) + 
  facet_grid(ONTOLOGY ~ ., scale = "free") + 
  theme_minimal() + 
  scale_colour_gradient(low = "#4daf4a", high = "#e41a1c") +
  labs(
    title = "Términos GO Enriquecidos por Dominio Funcional",
    x = "Proporción Génica (Gene Ratio)",
    colour = "p.ajustado" 
  ) +
  theme(
    axis.title.y = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold", size = 11),
    axis.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



