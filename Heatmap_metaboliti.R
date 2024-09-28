library(stringr)
library(openxlsx)
library(readxl)
library(pheatmap)
library(DAAG)
library(factoextra)
library(lawstat)
library(broom)
library(AICcmodavg)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(rmarkdown)
setwd('C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/Analisi metabolomica/script_analisi_metabolomica/')

rm(list=ls())

#Importo il dataset
matrice_metaboliti <- read.xlsx("matrice_metaboliti.xlsx",rowNames =  TRUE)
matrice_metaboliti <- matrice_metaboliti[,order(colnames(matrice_metaboliti))]

#Trasformazione logaritmica
matrice_metaboliti <- log2(matrice_metaboliti +1)
samples_df <- read.xlsx('Samples_data_phyloseq.xlsx',rowNames = TRUE)
rownames(samples_df) <- samples_df$label

samples_df <- samples_df[order(row.names(samples_df)),]

#Rimuovo i metaboliti che sono 0 in tutte le condizioni
d <- apply(matrice_metaboliti,1,mean)
ind <- which(d == 0)
matrice_metaboliti <- matrice_metaboliti[-ind,] 
rm(ind)
rm(d)


ind <- grep(pattern = "X BL", rownames(samples_df))
samples_df <- samples_df[-ind,]

ind <- grep(pattern = "X EP", rownames(samples_df))
samples_df <- samples_df[-ind,]



  
  
### Uso la funzione grep: cosi in base a cosa do in input posso selezionare le
# colonne che voglio. 
# Mi basta commentare queste 6 righe codice per fare la heatmap di tutti i campioni

#ind1 <- grep("LrH",colnames(matrice_metaboliti))
#matrice_metaboliti <- matrice_metaboliti[,-ind1]
#ind2 <- grep("Lrt",colnames(matrice_metaboliti))
#matrice_metaboliti <- matrice_metaboliti[,-ind2]
#ind3 <- grep("CT",colnames(matrice_metaboliti))
#matrice_metaboliti <- matrice_metaboliti[,-ind3]
#ind4 <- grep("BL",colnames(matrice_metaboliti))
#matrice_metaboliti <- matrice_metaboliti[,-ind4]

# farlo se voglio fare 'heatmap solo per quelli stat. significativi per fermentazione
#matrice_metaboliti <- t(matrix_metabo_analyst)

#Cutoff sugli zeri
#x <- apply(matrice_metaboliti,1,mean)
#ind <- which(x == 0)
#matrice_metaboliti <- matrice_metaboliti[-ind,]


#seleziono i 4 metaboliti espressi diff. tra t0 e t6
ind <- grep(pattern = "p-Cresol", rownames(matrice_metaboliti))
ind2 <- grep(pattern = "Hexanal",rownames(matrice_metaboliti))
ind3 <- grep(pattern = "Decanal",rownames(matrice_metaboliti))
ind4 <- grep(pattern = "2-Hexanone",rownames(matrice_metaboliti))
matrice_metaboliti <- matrice_metaboliti[c(ind,ind2,ind3,ind4),]
#Annotazioni per CTRL, LR, LRH
ctrl <- grep("Ct",colnames(matrice_metaboliti))
ctrl[] <- "CTRL"
lrh <- grep("LrH",colnames(matrice_metaboliti))
lrh[] <- "LrH"
lr <- grep("Lrt",colnames(matrice_metaboliti))
lr[] <- "Lr"
samp <- c(ctrl,lrh,lr)
annotation_campioni <- data.frame(Sampletype = samp)
rownames(annotation_campioni) <- colnames(matrice_metaboliti)

#Annotazioni per T0 e T6
test1 <- grepl("0",colnames(matrice_metaboliti))
tempi <- ifelse(test1,"T0","T6")
annotation_tempi <- data.frame(Tempi = tempi)
rownames(annotation_tempi) <- colnames(matrice_metaboliti)

#Annotazioni per BL ed EP
test2 <- grepl("BL",colnames(matrice_metaboliti))
var <- ifelse(test2,"Inizio ferm.","Fine ferm.")
annotation_ferm <- data.frame(Fasi = var)
rownames(annotation_ferm) <- colnames(matrice_metaboliti)

annotation <- cbind(annotation_campioni,annotation_ferm,annotation_tempi)

### COLORI ANNOTAZIONI
#"mediumseagreen","yellow3" per tempi
#"turquoise","violet" per fasi
# "wheat3","deepskyblue2","firebrick2" per tipo di campione

#Colori delle variabili
colori <- c("turquoise","violet")
names(colori) <- unique(var)
ann_colori_fasi <- list(Fasi = colori)

colori_tempi <- c("mediumseagreen","yellow3" )
names(colori_tempi) <- unique(tempi)
ann_colori_tempi <- list(Tempi = colori_tempi)

colori_campioni <- c("wheat3","deepskyblue2","firebrick2")
names(colori_campioni) <- unique(samp)
ann_colori_campioni <- list(Sampletype = colori_campioni)

annotation_colors <- c(ann_colori_fasi,ann_colori_tempi,ann_colori_campioni)  


#Heatmap
pheatmap(matrice_metaboliti,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         clustering_method = "complete",
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = annotation,
         annotation_colors = annotation_colors,
         annotation_names_col = T,
         show_rownames = F,fontsize_row = 9,
         show_colnames = T,angle_col = 45,
         cutree_rows = 2,
         cutree_cols = 2,
         scale = "row",cellwidth = 6,cellheight = 2,
         border_color = NA,
         color = colorRampPalette(colors = c("blue","blue3","black","yellow3","yellow"))(100))

#pheatmap 4 metaboliti tempo t0-t6
pheatmap(matrice_metaboliti,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         clustering_method = "complete",
         
         annotation_col = annotation_tempi,
         annotation_names_col = T,
         show_rownames = T,fontsize_row = 9,
         show_colnames = T,angle_col = 45,
         cutree_rows = 2,
         cutree_cols = 2,
         scale = "row",
         border_color = NA )

