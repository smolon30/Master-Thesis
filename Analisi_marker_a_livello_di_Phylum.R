library(readr)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(openxlsx)
library(vegan)
library(ggsignif)
library(dplyr)
library(ggpubr)
library(Matrix)
library(reshape2)
library(vegan)
library(moments)
library(ape)
library(metagenomeSeq)
library(metagMisc)
library(MASS)
library(lefser)
library(SummarizedExperiment)
library(microbiomeMarker)
library(rmarkdown)
setwd("C:/Users/asus 147456/Desktop/UNIVERSITA/1.TESI MAGISTRALE/Analisi batterica/Script_analisi_batterica/")

rm(list=ls())
## Importo i dati
ASV_table <- read_tsv(file = 'ASVtable.tsv')
taxonomy <- read_tsv(file = 'taxonomy.tsv')
taxonomy <- taxonomy[,-3]

samples_df <- read.xlsx('Samples_data_phyloseq.xlsx',rowNames = TRUE)
samples_df$Samples <- NA
samples_df$Samples <- rownames(samples_df)
samples_df <- samples_df[order(row.names(samples_df)),]



## Elimino i campioni 25 e 26 che sono quelli di controllo della fermentazione

ind <- grep(pattern = "25_S25",colnames(ASV_table))
ASV_table <- ASV_table[,-ind]

ind <- grep("26_S26",colnames(ASV_table))
ASV_table <- ASV_table[,-ind]

ind <- grep(pattern = "25_S25", rownames(samples_df))
samples_df <- samples_df[-ind,]

ind <- grep(pattern = "26_S26", rownames(samples_df))
samples_df <- samples_df[-ind,]


## Filtraggio ASV dove ci sono tutti zeri

ind2 <- 0
num_rows <- nrow(ASV_table)
for (i in 1:num_rows){
  ind <- which(ASV_table[i,] == 0)
  if (length(ind) == 24)(
    ind2[i] <- i)
  
}

ind2 <- which(ind2 != "NA")
ind2 <- ind2[-1]

ASV_table <- ASV_table[-ind2,]
taxonomy <- taxonomy[-ind2,]

## Divido la tabella della Tassonomia 

taxonomy <- separate(taxonomy,col = 2,into = c("Kingdom","Phylum","Class","Order",
                                               "Family","Genus","Species"), sep = " ")
n <- nrow(ASV_table)
m <- ncol(taxonomy)
for (i in 1:n){
  taxonomy[i,2] <- sub(';','',taxonomy[i,2])
}

for (i in 1:n){
  taxonomy[i,3] <- sub(';','',taxonomy[i,3])
}


for (i in 1:n){
  taxonomy[i,4] <- sub(';','',taxonomy[i,4])
}

for (i in 1:n){
  taxonomy[i,5] <- sub(';','',taxonomy[i,5])
}

for (i in 1:n){
  taxonomy[i,6] <- sub(';','',taxonomy[i,6])
}



for (i in 1:n){
  taxonomy[i,7] <- sub(';','',taxonomy[i,7])
}



# Rinomino i taxa_names per fare il test di mann whitney (commentare il codice se non serve)


## Disegno la curva di rarefazione
#raremax <- min(colSums(ASV_table[,-1]))
#rarecurve(t(ASV_table[,-1]),step = 20,sample = raremax,
#          col ="navyblue",cex = 0.4,)
##Species_rare <- rarefy(ASV_table[,-1],raremax)
#N_species <- c(1:nrow(ASV_table[,-1]))

# drarefy ci dice la probabilita che le specie siano presenti in una "comunita"
# di dimensione "sample"

### Quante specie elimino: grafico (?)

#plot(N_species, Species_rare, xlab = "Observed Number of Species", ylab = "Rarefied No. of Species")



######## PACCHETTO PHYLOSEQ

## Per prima cosa devo creare l'oggetto phyloseq, a partire dalla ASV_table e dalla 
# taxonomy table.

## Devo impostare la prima colonna come rowNames
ASV_table <- as.matrix(ASV_table)
rownames(ASV_table) <- ASV_table[,1]
ASV_table <- ASV_table[,-1]
x <- apply(ASV_table,1,as.numeric)
x <- t(x)
colnames(x) <- colnames(ASV_table)

taxonomy <- as.matrix(taxonomy)
rownames(taxonomy) <- taxonomy[,1]
taxonomy <- taxonomy[,-1]

# CREO gli "oggetti" da mettere nell'"oggetto phyloseq"
# 1. Taxonomy
taxonomy_phyloseq <- tax_table(taxonomy)

# 2. ASVtable
ASV_table_phyloseq <- otu_table(x,taxa_are_rows = T)

# 3. Samplesdata
samples <- sample_data(samples_df)



#Creo l'oggetto di classe "phyloseq"
Global_data <- phyloseq(ASV_table_phyloseq,taxonomy_phyloseq,samples)

# Ora posso fare tutte le analisi statistiche con questo. Per prima cosa devo
# fare la rarefazione (devo trovare la funzione).

# Tentativo 1. Normalizzazione con la mediana

#total = median(sample_sums(Global_data))
#standf = function(x, t=total) round(t * (x / sum(x)))
#z.rarefied = transform_sample_counts(Global_data, standf)

#Normalizzazione con rarefy_even_depth
Global_data.rarefied = rarefy_even_depth(Global_data, sample.size=min(sample_sums(Global_data)),
                                         replace=F)


## TEST DI MANN-WHITNEY

# Normalizzazione CSS. prima di fare il confronto tra gruppi (non per la beta diversity)
# - Non devo usare la matrice rarefatta.
#"CSS": cumulative sum scaling, calculates scaling factors as the cumulative sum of metabolites abundances up to a data-derived threshold.

# Uso pacchetto metagMisc

#Normalizzazione CSS
physeq_css <- phyloseq_transform_css(Global_data)
physeq_css_relative <- transform_sample_counts(physeq_css,function(x) x / sum(x) )

write.xlsx(x = otu_table(physeq_css_relative),file = "C:/Users/asus 147456/Desktop")
write.xlsx(x = tax_table(physeq_css_relative),file = "C:/Users/asus 147456/Desktop")

#Collassamento a livello di Phylum
physeq_css_collapse <- tax_glom(physeq_css_relative, taxrank = "Phylum")


x <- apply(otu_table(physeq_css_collapse),MARGIN = 1,FUN = mean)



#Assegno i rownames
table_otu <- as.matrix(otu_table(physeq_css_collapse))
table_otu <- as.data.frame(table_otu)
table_tax <- as.matrix(tax_table(physeq_css_collapse))
table_tax <- as.data.frame(table_tax)
rownames(table_tax) <- table_tax[,2]
rownames(table_otu) <- table_tax[,2]

# 3. Devo filtrare le ASV che non sono presenti in almeno il 25%: significa che se ho 0 in almeno 
# il 75% dei campioni elimino quella ASV
#ind2 <- 0
#num_rows <- nrow(otu_table(physeq_rel_css))
#for (i in 1:num_rows){
# ind <- which(otu_table(physeq_rel_css)[i,] == 0)
#if (length(ind) > 21)(
#ind2[i] <- i)
### con la soglia = 21 elimino sono campioni che non sono presenti in almeno il 25.5% 
# dei campioni
#}

#ind2 <- which(ind2 != "NA")
#ind2 <- ind2[-1]
#otu_table(physeq_rel_css) <- otu_table(physeq_rel_css)[-ind2,]

#Ricreo oggetto phyloseq e lo nomino "ps" 
table_tax_ps <- tax_table(as.matrix(table_tax))
table_otu_ps <- otu_table(as.matrix(table_otu),taxa_are_rows = T)
ps <- phyloseq(table_otu_ps,table_tax_ps,samples)

## FILTRO PERCENTUALE OTUS

#filter <- phyloseq::genefilter_sample(ps,filterfun_sample(function(x) x > 0), A = 0.1*nsamples(ps)) #con 0.1 seleziono le asv che sono presenti in almeno il 10% dei campioni. 
# questo perche A e il numero di campioni tale che se nella asv ci sono piu campioni maggiori>0 di A, mi tengo quel campione.
# restituisce un vettore logico con gli indici delle ASV che voglio tenere.

#ps_filtered <- prune_taxa(filter,ps)
#con prune_taxa tengo solo le ASV che ho selezionato con il filtro, ed elimino tutte le altre


# 4. Applico il test di MANN-WHITNEY per ogni riga

# TEST DI MANN-WHITNEY:

# - Ipotesi nulla : la differenza tra le mediane di coppie di osservazioni e zero.
# - Ipotesi alternativa : la differenza tra le mediane di coppie di osservazioni e diversa da zero.

# Si usa quando la distribuzione dei dati non e gaussiana.

# prendendo un campione a caso, si vede che i dati non sono distribuiti normalmente
# (ci sono tanti zeri, lo riverifico dopo il filtraggio)

# Si usa wilcox.test per fare il test di Mann Whitney, che e la stessa cosa del test di wilcoxon
# questa funzione la posso usare solo tra due vettori.



## DEVO ESTRARMI NUOVAMENTE LE TABLE SE VOGLIO FARE IL TEST SULLA MATRICE FILTRATA.


#Separo i gruppi, cosi da poter fare il test di mann whitney
CTRL <- as.matrix(table_otu[,samples_df["SampleType"] == "CTRL"])
LR <- as.matrix(table_otu[,samples_df["SampleType"] == "LR"])
LRH <- as.matrix(table_otu[,samples_df["SampleType"] == "LRH"])

BL <- as.matrix(table_otu[,samples_df["Phase"] == "BL"])
EP <- as.matrix(table_otu[,samples_df["Phase"] == "EP"])

T0 <- as.matrix(table_otu[,samples_df["Time"] == "T0"])
T6 <- as.matrix(table_otu[,samples_df["Time"] == "T6"])

#Calcolo medie
#Medie delle abbondanze relative
media_nei_CTRL <- rowMeans(CTRL)
media_nei_LR <- rowMeans(LR)
media_nei_LRH <- rowMeans(LRH)
media_nei_BL <- rowMeans(BL)
media_nei_EP <- rowMeans(EP)
media_nei_T0 <- rowMeans(T0)
media_nei_T6 <- rowMeans(T6)


#Applico il test di mann whitney per tutte le asv: CTRL vs LR
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- t.test(CTRL[i,],LR[i,])$p.value
}
p_value

#Corrego con metodo BH
p_value_adjust <- p.adjust(p_value, 'BH')

#Creo vettore da importare su excel
v1 <- cbind(media_nei_CTRL,media_nei_LR, p_value,p_value_adjust)
rownames(v1) <- rownames(table_otu)
v1 <- data.frame(v1)
#write.xlsx(x = v1,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = v1,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)


#test di mann whitney : CTRL vs LRH
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- wilcox.test(CTRL[i,],LRH[i,])$p.value
}

p_value
p_value_adjust <- p.adjust(p_value, 'BH')
v2 <- cbind(media_nei_CTRL,media_nei_LRH, p_value,p_value_adjust)
rownames(v2) <- rownames(table_otu)

v2 <- data.frame(v2)
#write.xlsx(x = v2,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = v2,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)

# test di mann whitney: LR vs LRH
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- wilcox.test(LR[i,],LRH[i,])$p.value
}

p_value
p_value_adjust <- p.adjust(p_value, 'BH')

v3 <- cbind(media_nei_LR,media_nei_LRH, p_value,p_value_adjust)
rownames(v3) <- rownames(table_otu)
v3 <- data.frame(v3)
#write.xlsx(x = v3,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = v3,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)

# test di mann whitney: EP vs BL
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- wilcox.test(BL[i,],EP[i,])$p.value
}

p_value
p_value_adjust <- p.adjust(p_value, 'BH')

v4 <- cbind(media_nei_BL,media_nei_EP, p_value,p_value_adjust)
rownames(v4) <- rownames(table_otu)

v4 <- data.frame(v4)
#write.xlsx(x = v4,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = v4,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)

ind <- which(v4[,4] < 0.05)
v4_signif <- v4[ind,]
v4_signif[,1] <- sort(v4_signif[,1], decreasing = FALSE)
v4_signif[,2] <- sort(v4_signif[,2], decreasing = FALSE)
#write.xlsx(x = v4_signif,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)

colnames(table_otu) <- samples_df$label
# matrice per cytoscape 
table_otu <- table_otu[-c(3,5),]
x = data.frame(t(table_otu))
write.xlsx(x = x,file = 'C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/Matrice_phylum_cytoscape.xlsx',rowNames = T, colNames = T)
# test di mann whitney: T0 vs T6
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- wilcox.test(T0[i,],T6[i,])$p.value
}

p_value
p_value_adjust <- p.adjust(p_value, 'BH')

v5 <- cbind(media_nei_T0,media_nei_T6, p_value,p_value_adjust)
#rownames(v5) <- rownames(table_otu)

#v5 <- data.frame(v5)
#write.xlsx(x = v5,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = v5,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)




### USO PACCHETTO microbiomeMarker PER FARE I TEST STATISTICI

#ps_ctrl_lrh <- phyloseq::subset_samples(ps_genus_lefse,
#                                        SampleType %in% c('CTRL','LRH') )


#results_lefse <- run_lefse(ps = ps_ctrl_lrh,group = 'SampleType')

#Calcolo risultato LefSe
results_lefse_Phase <- run_lefse(ps = ps,group = 'Phase')
#results_lefse_Time <- run_lefse(ps = ps,group = 'Time')

marker_phase <- marker_table(results_lefse_Phase)
marker_phase <- as.data.frame(marker_phase)

#marker_time <- marker_table(results_lefse_Time)
#marker_time <- as.data.frame(marker_time)
#write.xlsx(x = marker_phase,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = marker_time, file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
## il p-value che ottengo con run_lefse e il p_value del test di kruskal wallis


#LEfSe (Linear discriminant analysis Effect Size) determina le caratteristiche 
#(organismi, cladi, unita tassonomiche operative, geni o funzioni) che hanno maggiori 
#probabilita di spiegare le differenze tra le classi accoppiando test standard per 
#la significativita statistica con test aggiuntivi che codificano la coerenza biologica 
#e la rilevanza degli effetti

#Uso LefSe solo a livello di Phylum

mm_lefse <- run_lefse(ps = ps, group = "Phase",taxa_rank = "Phylum")

marker_phase <- marker_table(mm_lefse)
