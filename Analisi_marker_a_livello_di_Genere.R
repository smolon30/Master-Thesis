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
library(tidyr)
library(pheatmap)
library(rmarkdown)
setwd("C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/Analisi batterica/Script_analisi_batterica/")

rm(list=ls())
soglia_fc <- 1.5
soglia_pval <- 0.05
## Importo i dati
ASV_table <- read_tsv(file = 'ASVtable.tsv')
taxonomy <- read_tsv(file = 'taxonomy.tsv')
taxonomy <- taxonomy[,-3]

samples_df <- read.xlsx('Samples_data_phyloseq.xlsx',rowNames = TRUE)
samples_df$Samples <- NA
samples_df$Samples <- rownames(samples_df)
samples_df <- samples_df[order(row.names(samples_df)),]

colnames(ASV_table)[2:27] <- samples_df$label
rownames(samples_df) <-samples_df$label
## Elimino i campioni 25 e 26 che sono quelli di controllo della fermentazione

ind <- grep(pattern = "X",colnames(ASV_table))
ASV_table <- ASV_table[,-ind]



ind <- grep(pattern = "X", rownames(samples_df))
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




# Normalizzazione CSS. prima di fare il confronto tra gruppi (non per la beta diversity)
# - Non devo usare la matrice rarefatta.
#"CSS": cumulative sum scaling, calculates scaling factors as the cumulative sum of metabolites abundances up to a data-derived threshold.

# Uso pacchetto metagMisc

#Normalizzazione CSS
physeq_css <- phyloseq_transform_css(Global_data)

#Trasformazione in abbondanze relative
physeq_css_relative <- transform_sample_counts(physeq_css,function(x) x / sum(x) )

#plot_heatmap(physeq_css_relative)



library(microbiomeMarker)
#Collassamento a livello di genere o di phylum (o quello che voglio)
physeq_css_collapse <- tax_glom(physeq_css_relative, taxrank = "Genus")

#Creo una nuova colonna che mi dia la risoluzione a livello inferiore dove non ho il genere.
table_tax <- as.matrix(tax_table(physeq_css_collapse))
table_tax <- as.data.frame(table_tax)

ind <- which(table_tax[,6] == 'g__') 
table_tax[ind,8] <- paste(table_tax[ind,6],table_tax[ind,5])
table_tax[ind,6] <- table_tax[ind,8]

ind <- which(table_tax[,8] == 'g__ f__')
table_tax[ind,9] <- paste(table_tax[ind,8],table_tax[ind,4])
table_tax[ind,6] <- table_tax[ind,9]

ind <- which(table_tax[,9] == 'g__ f__ o__') 
table_tax[ind,10] <- paste(table_tax[ind,9],table_tax[ind,3])
table_tax[ind,6] <- table_tax[ind,10]

ind <- which(table_tax[,10] == 'g__ f__ o__ c__') 
table_tax[ind,11] <- paste(table_tax[ind,10],table_tax[ind,2])
table_tax[ind,6] <- table_tax[ind,11]

table_tax <- table_tax[,-(8:11)]

ind <- which(table_tax[,6] == 'g__Clostridium')
table_tax[ind,6] <- c("f_Peptostreptococcaceae_g_Clostridium", "f_Clostridiaceae_g_Clostridium")


#Assegno i Rownames
table_otu <- as.matrix(otu_table(physeq_css_collapse))
table_otu <- as.data.frame(table_otu)
rownames(table_tax) <- table_tax[,6]
rownames(table_otu) <- table_tax[,6]

ind1 <- which(rownames(table_otu) == 'g__[Ruminococcus]')
ind2 <- which(rownames(table_otu) == 'g__Ruminococcus')

vector <- 0
for (i in 1:ncol(table_otu)){
  vector[i] <- sum(table_otu[ind1,i],table_otu[ind2,i])
}


table_tax_genus <- as.data.frame(table_tax[,6],row.names = table_tax[,6])
colnames(table_tax_genus)[1] <- 'Genus'
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

#Ricreo oggetto phyloseq
table_tax_ps <- tax_table(as.matrix(table_tax))
table_otu_ps <- otu_table(as.matrix(table_otu),taxa_are_rows = T)
ps <- phyloseq(table_otu_ps,table_tax_ps,samples)

table_tax_ps_genus <- tax_table(as.matrix(table_tax_genus))

### ps e IL MIO NUOVO OGGETTO PHYLOSEQ.

## FILTRO PERCENTUALE OTUS

#filter <- phyloseq::genefilter_sample(ps,filterfun_sample(function(x) x > 0), A = 0.1*nsamples(ps)) #con 0.1 seleziono le asv che sono presenti in almeno il 10% dei campioni. 
# questo perche A possiede un numero di campioni tale che se nella asv ci sono piu campioni maggiori>0 di A, mi tengo quel campione.
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

#Medie delle abbondanze relative
media_nei_CTRL <- rowMeans(CTRL)
media_nei_LR <- rowMeans(LR)
media_nei_LRH <- rowMeans(LRH)
media_nei_BL <- rowMeans(BL)
media_nei_EP <- rowMeans(EP)
media_nei_T0 <- rowMeans(T0)
media_nei_T6 <- rowMeans(T6)

#Applico il test di mann whitney per tutte le ASV: CTRL vs LR
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- t.test(CTRL[i,],LR[i,])$p.value
}
p_value
p_value_adjust <- p.adjust(p_value, 'BH')
v1 <- cbind(p_value,p_value_adjust)
rownames(v1) <- rownames(table_otu)
v1 <- data.frame(v1)
#write.xlsx(x = v1,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = v1,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)


#Test di mann whitney : CTRL vs LRH
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- wilcox.test(CTRL[i,],LRH[i,])$p.value
}
p_value
p_value_adjust <- p.adjust(p_value, 'BH')
v2 <- cbind(p_value,p_value_adjust)
rownames(v2) <- rownames(table_otu)
v2 <- data.frame(v2)
#write.xlsx(x = v2,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = v2,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)


#Test di mann whitney: LR vs LRH
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- wilcox.test(LR[i,],LRH[i,])$p.value
}
p_value
p_value_adjust <- p.adjust(p_value, 'BH')
v3 <- cbind(p_value,p_value_adjust)
rownames(v3) <- rownames(table_otu)
v3 <- data.frame(v3)
#write.xlsx(x = v3,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = v3,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)


#Test di mann whitney: EP vs BL
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- wilcox.test(BL[i,],EP[i,])$p.value
}
p_value
p_value_adjust_phase <- p.adjust(p_value, 'BH')
v4 <- cbind(p_value,p_value_adjust_phase)
rownames(v4) <- rownames(table_otu)
v4 <- data.frame(v4)
#write.xlsx(x = v4,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = v4,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)
ind <- which(v4[,1] < 0.05)
v4_signif <- v4[ind,]
v4_signif[,1] <- sort(v4_signif[,1], decreasing = FALSE)
v4_signif[,2] <- sort(v4_signif[,2], decreasing = FALSE)
#write.xlsx(x = v4_signif,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)

#Test di mann whitney: T0 vs T6
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- wilcox.test(T0[i,],T6[i,])$p.value
}
p_value
p_value_adjust <- p.adjust(p_value, 'BH')
v5 <- cbind(p_value,p_value_adjust)
rownames(v5) <- rownames(table_otu)
v5 <- data.frame(v5)
#write.xlsx(x = v5,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = v5,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney_Phylum',colNames = TRUE,rowNames = TRUE)


### USO PACCHETTO microbiomeMarker PER FARE I TEST STATISTICI



results_lefse_Phase <- run_lefse(ps = ps,group = 'Phase')
results_lefse_Time <- run_lefse(ps = ps,group = 'Time')

marker_phase <- marker_table(results_lefse_Phase)
marker_phase <- as.data.frame(marker_phase)

marker_time <- marker_table(results_lefse_Time)
marker_time <- as.data.frame(marker_time)
#write.xlsx(x = marker_phase,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#write.xlsx(x = marker_time, file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
## il p-value che ottengo con run_lefse e il p_value del test di kruskal wallis

#write.xlsx(table_otu,file = 'C:/Users/asus 147456/Desktop/TESI MAGISTRALE/Metaboliti formaggini/File_excel_tabelle_mann_whitney',colNames = TRUE,rowNames = TRUE)
#LEfSe (Linear discriminant analysis Effect Size) determina le caratteristiche 
#(organismi, cladi, unita tassonomiche operative, geni o funzioni) che hanno maggiori 
#probabilita di spiegare le differenze tra le classi accoppiando test standard 
#per la significativita statistica con test aggiuntivi che codificano la coerenza biologica
#e la rilevanza degli effetti

#File per Cytoscape
colnames(table_otu) <- samples_df$label
otu_genere_cytoscape <- as.data.frame(t(table_otu))
#write.xlsx(x = otu_genere_cytoscape,file = "C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/otu_genere_cytoscape.xlsx")


#Selezioni le ASV statisticamente significative (per fase)

table_otu$pvalue <- p_value_adjust_phase
ind <- which(table_otu$pvalue < 0.05)
table_otu_signif <- table_otu[ind,]

#Plot dell'abbondanza dei marker batterici per fase
plot_abundance(results_lefse_Phase, group = 'Phase')

#Creazione file per correlazione e Cytoscape (ASV)
table_otu_signif_traslata <- as.data.frame(t(table_otu_signif))
table_otu_signif_traslata <- table_otu_signif_traslata[-25,] #tolgo il p-value


table_otu_signif_traslata$Group <- samples_df$Phase
table_otu_signif_traslata$Sample_ID <- samples_df$Samples

table_otu_signif_traslata <- table_otu_signif_traslata %>% relocate(Group)
table_otu_signif_traslata <- table_otu_signif_traslata %>% relocate(Sample_ID)

#write.xlsx(x = table_otu_signif_traslata,file = "C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/ASV_correlazione.xlsx",rowNames = T,colNames = T)

# plot abundance 
colnames(table_otu) <- samples_df$Samples
table_otu_signif <- table_otu[ind,]
table_otu_signif_ps <- table_otu_signif[,-c(25:27)]
table_otu_signif_ps <- otu_table(as.matrix(table_otu_signif_ps),taxa_are_rows = T)

table_tax[,1:5] <- NA
table_tax_signif_ps <- tax_table(as.matrix(table_tax))

#ps_phase_signif <- phyloseq(table_otu_signif_ps,table_tax_signif_ps,samples_df)

#lefse_signif <- run_lefse(ps = ps_phase_signif, group = 'Phase')#
plot_abundance(results_lefse_Phase, group = 'SampleType',markers = c('k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae',
                                                         'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia',
                                                         'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__ f__Ruminococcaceae',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Roseburia',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|o__Clostridiales_f__|g__ f__ o__Clostridiales',
                                                         'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella',
                                                         'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__[Paraprevotellaceae]',
                                                         'k__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__RF32|o__RF32_f__|g__ f__ o__RF32',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Anaerostipes',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiaceae|g__ f__Clostridiaceae',
                                                         'k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales',
                                                         'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__[Barnesiellaceae]',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Coprococcus',
                                                         'k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus',
                                                         'k__Bacteria|p__Firmicutes|c__Negativicutes|o__Selemonadales|o__Selemonadales_f__Veillonellaceae|g__Megasphaera',
                                                         'k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia',
                                                         'k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Alcaligenaceae|g__Sutterella',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Dorea',
                                                         'k__Bacteria|p__Firmicutes|c__Negativicutes|o__Selemonadales|o__Selemonadales_f__',
                                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Peptostreptococcaceae|f_Peptostreptococcaceae_g_Clostridium',
                                                         'k__Bacteria|p__Proteobacteria|c__Deltaproteobacteria|o__Desulfovibrionales|f__Desulfovibrionaceae|g__Bilophila',
                                                         'k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Erwinia',
                                                         'k__Bacteria|p__Proteobacteria|c__Deltaproteobacteria|o__Desulfovibrionales|f__Desulfovibrionaceae'))

plot_ef_bar(results_lefse_Phase, markers = c('k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae',
                                            'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia',
                                            'k__Bacteria|p__Firmicutes',
                                            'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__ f__Ruminococcaceae',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Roseburia',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|o__Clostridiales_f__|g__ f__ o__Clostridiales',
                                            'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella',
                                            'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__[Paraprevotellaceae]',
                                            'k__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__RF32|o__RF32_f__|g__ f__ o__RF32',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Anaerostipes',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiaceae|g__ f__Clostridiaceae',
                                            'k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales',
                                            'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__[Barnesiellaceae]',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Coprococcus',
                                            'k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus',
                                            'k__Bacteria|p__Firmicutes|c__Negativicutes|o__Selemonadales|o__Selemonadales_f__Veillonellaceae|g__Megasphaera',
                                            'k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia',
                                            'k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Alcaligenaceae|g__Sutterella',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Dorea',
                                            'k__Bacteria|p__Firmicutes|c__Negativicutes|o__Selemonadales|o__Selemonadales_f__',
                                            'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Peptostreptococcaceae|f_Peptostreptococcaceae_g_Clostridium',
                                            'k__Bacteria|p__Proteobacteria|c__Deltaproteobacteria|o__Desulfovibrionales|f__Desulfovibrionaceae|g__Bilophila',
                                            'k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Erwinia',
                                            'k__Bacteria|p__Proteobacteria|c__Deltaproteobacteria|o__Desulfovibrionales|f__Desulfovibrionaceae'))




#uso lefse solo a livello di genus 
mm_lefse <- run_lefse(ps,group = 'Phase',taxa_rank = 'Genus' )

#lefse a livello di phylum
phylum_lefse <- run_lefse(ps, group = 'Phase',taxa_rank = 'Phylum')
write.xlsx(x = marker_table(mm_lefse),file = "C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/genus_lefse.xlsx",rowNames = T,colNames = T)

#Plot abbondanza marker batterici solo a livello di Genus  
plot_abundance(mm = mm_lefse,group = 'Phase')

#Plot dell'effect size dei marker batterici
plot_ef_bar(mm_lefse)

#Plot heatmap dei marker batterici
plot_heatmap(mm_lefse, group = 'Phase',cluster_sample = TRUE,sample_label = T)





#Heatmap dei soli marker batterici trovati con la fermentazione
colnames(table_otu) <- samples_df$label
table_otu_signif <- table_otu[ind,]
table_otu_signif <- table_otu_signif[,-c(24:27)]
pheatmap(table_otu_signif,scale = 'row' ,clustering_method  = 'complete', 
         angle_col = 45,
         fontsize_row  = 8,
         fontsize_col = 8,cutree_rows = 2,cutree_cols = 2,
         color = colorRampPalette(colors = c("blue","blue3","black","yellow3","yellow"))(100),
         cellwidth = 12,
         cellheight = 12)



#Volcano plot dei batterici statisticamente significativi con Test di Mann-Whitney

#Calcolo del Fold Change
ind1 <- grep( pattern = "BL",x = colnames(table_otu))
ind2 <- grep( pattern = "EP", x = colnames(table_otu))
table_otu <- log2(table_otu +1)
x <- rowMeans(table_otu[,ind1])
y <- rowMeans(table_otu[,ind2])
log2FC <- y-x

#Calcolo dei batteri stat. significativi
p_value <- 0
for (i in 1:nrow(table_otu)){
  
  p_value[i] <- wilcox.test(BL[i,],EP[i,])$p.value
}
p_value
p_value_adjust <- p.adjust(p_value, 'BH')
y <- -log10(p_value_adjust)
table_otu$Log2FC <- log2FC
table_otu$pvalue <- p_value_adjust

#Dico alla matrice quali sono up-regolati
table_otu$diffexpressed <- "NO"
#Chiamo UP quelli con log2FC > 0.6 e pvalue < 0.05 
table_otu$diffexpressed[table_otu$log2FC > log2(soglia_fc) & table_otu$pvalue < 0.05] <- "UP"
#Chiamo DOWN quelli con log2FC < -0.6 e pvalue < 0.05
table_otu$diffexpressed[table_otu$log2FC < -log2(soglia_fc) & table_otu$pvalue < 0.05] <- "DOWN"

#Creo vettore colori
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

#Creo vettore Label
table_otu$delabel <- NA
table_otu$delabel[table_otu$diffexpressed != "NO"] <- rownames(table_otu)[table_otu$diffexpressed != "NO"]

#Volcano Plot
ggplot(data=table_otu, aes(x=Log2FC, y= -log10(pvalue), col = diffexpressed, label = delabel)) + 
  geom_vline(xintercept=c(-log2(soglia_fc), log2(soglia_fc)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  theme_minimal()+
  geom_point()+
  scale_color_manual(values = mycolors)



  
  
