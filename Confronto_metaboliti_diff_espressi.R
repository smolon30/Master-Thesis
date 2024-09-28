####################

## STUDIO DISTRIBUZIONE DATI SULLA MATRICE IN CUI HO FATTO LA MEDIA TRA CAMPIONE E DUPLICATO

library(stringr)
library(openxlsx)
library(readxl)
library(ggplot2)
setwd('C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Analisi metabolomica/')
rm(list=ls())
par(mar = c(5,10,2,2))

caciotta_campioni <- read.xlsx('matrice_metaboliti.xlsx',rowNames = TRUE)
caciotta_campioni <- log2(caciotta_campioni + 1)
# voglio fare la media tra ogni campione e il suo duplicato biologico, e poi usare questi 
# dati per fare le analisi.


soglia_fc <- 2 #voglio tenermi solo i metaboliti che stanno aumentando o si stanno riducendo piu di 1.5 volte rispetto al normale, 
#cioè in scala logaritmica quelli che aumentano o si riducono di log2(1.5) volte


#Fold Change metaboliti tra inizio e fine fermentazione
ind  <- grep(pattern = "BL", colnames(caciotta_campioni))
medie_BL <- caciotta_campioni[,ind]
ind2 <- grep(pattern = "EP",x = colnames(caciotta_campioni))
medie_EP <- caciotta_campioni[,ind2]

#calcolo fold change
logFC = rowMeans(medie_EP) - rowMeans(medie_BL)


hist(logFC, main = "Frequenza di distribuzione (logaritmica) del FC", breaks = 50, xlab = "LOG(FC)", ylab="Frequenza", col = "red" )
abline(v=c(-log2(soglia_fc), log2(soglia_fc)),lty = 2,lwd = 4, col = 'grey')
# quelli all'esterno delle linee tratteggiate sono variati almeno di 1.5 volte rispetto al valore iniziale.

ind <- which(abs(logFC) < log2(soglia_fc))
caciotta_campioni1 <- caciotta_campioni[-ind,] # n.b da qui leggere i nomi dei metaboliti che stanno variando
logFC <- logFC[-ind]
rm(ind)

#Barplot
X <- sort(logFC,decreasing = FALSE)
barplot(X, horiz = TRUE,cex.axis = 1.5,col = "steelblue" ,border = "black", las = 1,cex.names = 0.8, xlab = "LOG(FC)",
        font.axis = 3)
# QUELLI che nel grafico sono maggiori di 0 stanno aumentando nella condizione LRH rispetto a CTRL, gli altri stanno diminuendo

