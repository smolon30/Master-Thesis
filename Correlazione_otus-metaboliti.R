library(psych)
library(stringr)
library(openxlsx)
library(readxl)
library(reshape2)
setwd('C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/Analisi metabolomica/script_analisi_metabolomica/')
rm(list = ls())

#Importo le ASV significative
matrice_ASV <- read.xlsx("C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/ASV_correlazione.xlsx")

#Importo i metaboliti significativi
matrice_metaboliti <- read.xlsx("C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/metaboliti_correlazione.xlsx")

matrice_ASV_phylum <- read.xlsx("C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/Matrice_phylum_cytoscape.xlsx")

matrice_ASV <- matrice_ASV[,-c(16,17,19,21)] #tengo fino a famiglia
rownames(matrice_ASV_phylum) <- matrice_ASV_phylum$ID
rownames(matrice_ASV) <- matrice_ASV$ID
rownames(matrice_metaboliti) <- matrice_metaboliti$ID


#Correlazione di Pearson
pearson_corr <- corr.test(x = matrice_ASV[,-c(1:3)], y = matrice_metaboliti[,-c(1:3)],method = 'pearson',use = "complete",adjust = "BH")

#Creo matrice correlazione
matrice_correlazione <- merge(melt(pearson_corr$r, value.name = 'corr'),melt(pearson_corr$p, value.name = 'p-value'),by = c('Var1','Var2'))

#Seleziono le correlazioni significative
ind <- which(matrice_correlazione$`p-value`<= 0.051)
matrice_correlazione_signif <- matrice_correlazione[ind,]
#matrice_correlazione_signif <- matrice_correlazione_signif[-c(1:21),]

#write.xlsx(x = matrice_correlazione_signif,file = 'C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/matrice_correlazione_cytoscape_finoaL5.xlsx', colNames = T, rowNames = T)
