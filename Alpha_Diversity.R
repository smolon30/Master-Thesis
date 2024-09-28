## ALPHA DIVERSITY
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
library(broom)
library(AICcmodavg)
library(ranacapa)
library(rmarkdown)
setwd("path")

rm(list=ls())

#Importo i dati
ASV_table <- read_tsv(file = 'ASVtable.tsv')
taxonomy <- read_tsv(file = 'taxonomy.tsv')
taxonomy <- taxonomy[,-3]
samples_df <- read.xlsx('Samples_data_phyloseq.xlsx',rowNames = TRUE)
samples_df$Samples <- NA
samples_df$Samples <- rownames(samples_df)
samples_df <- samples_df[order(row.names(samples_df)),]
#write.xlsx(x = ASV_table,file = 'path',colNames = TRUE,rowNames = TRUE)



#colnames(ASV_table)[2:27] <- samples_df$label
#ind <- grep(pattern = "BL",x = colnames(ASV_table))
#ASV_table[,ind] <- ASV_table[,ind] + ASV_table[,19]

#ind <- grep(pattern = "EP", x = colnames(ASV_table))
#ASV_table[,ind] <- ASV_table[,ind] + ASV_table[,20]


#matrice_controllo_ferm <- matrice_controllo_ferm[,-c(19,20)]



#Elimino i campioni 25 e 26 che sono quelli di controllo della fermentazione

ind <- grep(pattern = "25_S25",colnames(ASV_table))
ASV_table <- ASV_table[,-ind]

ind <- grep("26_S26",colnames(ASV_table))
ASV_table <- ASV_table[,-ind]

ind <- grep(pattern = "25_S25", rownames(samples_df))
samples_df <- samples_df[-ind,]

ind <- grep(pattern = "26_S26", rownames(samples_df))
samples_df <- samples_df[-ind,]




# FILTRO CAMPIONI: PRENDO SOLO EP oppure BL , oppure solo i T0 o solo i T6
#( non runnare se si vogliono considerare tutti i campioni!)
# elimino quelli che prendo con grep

#ind1 <- grep(pattern = "BL",samples_df$label)
#ind1 <- ind1 +1

##ind <- grep(pattern = "BL",samples_df$label)


#ASV_table <- ASV_table[,-ind1]
#samples_df <- samples_df[-ind,]


## Filtraggio ASV dove ci sono tutti zeri (mettere 24 quando considero tutti i campioni)
#Filtraggio ASV che sono "0" in tutti i campioni
ind2 <- 0
num_rows <- nrow(ASV_table)
for (i in 1:num_rows){
  ind <- which(ASV_table[i,] == 0)
  if (length(ind) == 24)(
    ind2[i] <- i)
  
}
ind2 <- which(ind2 != "NA")
ind2 <- ind2[-1]

#ASV_table <- ASV_table[-ind2,]
#taxonomy <- taxonomy[-ind2,]

#Divido la tabella della Tassonomia 
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



## Voglio scrivere un ciclo for che elimini tutti i ;

for (i in 1:n){
  for (j in 1:(m-1))
    taxonomy[i,j+1] <- sub(';','',taxonomy[i,j+1])
}

# Rinomino i taxa_names per fare il test di mann whitney (commentare il codice se non serve)


## Disegno la curva di rarefazione
#raremax <- min(colSums(ASV_table[,-1]))
#rarecurve(t(ASV_table[,-1]),step = 20,sample = raremax,
#          col ="navyblue",cex = 0.4,)
##Species_rare <- rarefy(ASV_table[,-1],raremax)
#N_species <- c(1:nrow(ASV_table[,-1]))

# di dimensione "sample"

### Quante specie elimino: grafico (?)

#plot(N_species, Species_rare, xlab = "Observed Number of Species", ylab = "Rarefied No. of Species")



######## PACCHETTO PHYLOSEQ

## Per prima cosa devo creare l'oggetto phyloseq, a partire dalla ASV_table e dalla 
# taxonomy table.

#Imposto la prima colonna come RowNames
ASV_table <- as.matrix(ASV_table)
rownames(ASV_table) <- ASV_table[,1]
ASV_table <- ASV_table[,-1]
x <- apply(ASV_table,1,as.numeric)
x <- t(x)
colnames(x) <- colnames(ASV_table)

taxonomy <- as.matrix(taxonomy)
rownames(taxonomy) <- taxonomy[,1]
taxonomy <- taxonomy[,-1]

#Creo gli "oggetti" da mettere nell'"oggetto phyloseq"
# 1. Taxonomy
taxonomy_phyloseq <- tax_table(taxonomy)

# 2. ASVtable
ASV_table_phyloseq <- otu_table(x,taxa_are_rows = T)

# 3. Sampledata
samples <- sample_data(samples_df)



#Creo l'oggetto di classe "phyloseq"
Global_data <- phyloseq(ASV_table_phyloseq,taxonomy_phyloseq,samples)

# Ora posso fare tutte le analisi statistiche con questo. Per prima cosa devo
# fare la rarefazione (devo trovare la funzione).


#Normalizzazione con rarefy_even_depth
Global_data.rarefied = rarefy_even_depth(Global_data, sample.size=min(sample_sums(Global_data)),
                                         replace=F)
write.xlsx(x = otu_table(Global_data.rarefied),file = "C:/Users/asus 147456/Desktop")
write.xlsx(x = tax_table(Global_data.rarefied),file = "C:/Users/asus 147456/Desktop")
#Curva di rarefazione
#p <- ggrare(Global_data, step = 1, color = "SampleType", label = "Sample", se = FALSE)+
  theme_bw()

### ALPHA DIVERSITY

#lista indici per alpha diversity
#c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher").

plot_richness(Global_data.rarefied,x = "SampleType", measures=c("Chao1"),color = 'SampleType')
plot_richness(Global_data.rarefied,x = "Phase", measures=c("Simpson"),color = 'Phase')
plot_richness(Global_data.rarefied,x = "Time", measures=c("Shannon"),color = 'Time')

## Oltre al plot devo anche creare la matrice con i valori di alpha diversity.
# Fare:
#1. Vettore misure alpha diversity.
#2. curve di rarefazione.
#3. Boxplot alpha diversity tra i diversi SampleType.
#4. Bar graphs. devo riportare le abbondanze relative, io invece ho quelle assolute.


#Shannon diversity
Shannon_diversity <- estimate_richness(Global_data.rarefied,split = T,measures = "Shannon")


#Test di Shapiro-Wilk
shapiro.test(Shannon_diversity$Shannon)

# Shapiro-Wilk normality test
#
# data:  Shannon_diversity$Shannon
# W = 0.94103, p-value = 0.172

## La distribuzione e gaussiana
#boxplot(CTRL,LR,LRH,col = c("olivedrab2","deepskyblue2","firebrick2"),main = "Shannon Diversity", xlab = "SampleType",
#ylab = "Shannon Index", names = c("CTRL","LR","LRH"), border = "black")

Shannon_diversity$SampleType <- NA
Shannon_diversity$Phase <- NA
Shannon_diversity$Time <- NA

# voglio fare solo ctrl vs tutti
#ind <- which(samples_df$SampleType == "LRH")
#samples_df$SampleType <- "LR"

#Assegnazione variabili
Shannon_diversity$SampleType[which(rownames(Shannon_diversity) %in% rownames(samples_df))] <- samples_df$SampleType 
Shannon_diversity$Phase[which(rownames(Shannon_diversity) %in% rownames(samples_df))] <- samples_df$Phase
Shannon_diversity$Time[which(rownames(Shannon_diversity) %in% rownames(samples_df))] <- samples_df$Time




#Plot Shannon diversity
ggplot(Shannon_diversity, aes(x = SampleType, y = Shannon)) + 
  geom_boxplot(fill = c("wheat","deepskyblue2","firebrick2"),width = 0.5,) + 
  geom_signif(comparisons = list(c("CTRL", "LR")), map_signif_level=TRUE, test = 't.test') +
  geom_signif(comparisons = list(c("CTRL", "LRH")),map_signif_level=TRUE,y_position = 5.3, test = 't.test')+
  geom_signif(comparisons = list(c("LR", "LRH")),map_signif_level=TRUE,y_position = 5.4, test = 't.test')+
  theme_bw() +
  ggtitle("Shannon Diversity") +
  theme(text = element_text(size = 18),aspect.ratio = 1)

# t-test Shannon Diversity: CTRL VS LR e CTRL vs LRH
CTRL <- Shannon_diversity[samples_df["SampleType"] == "CTRL",]
LR <- Shannon_diversity[samples_df["SampleType"]== "LR",]
LRH <- Shannon_diversity[samples_df["SampleType"]== "LRH",]
t.test(CTRL$Shannon,LR$Shannon)
t.test(CTRL$Shannon,LRH$Shannon)

#anova tra i 3 gruppi
one.way.anova.Shannon <- aov(Shannon ~ SampleType , data = Shannon_diversity )
summary(one.way.anova.Shannon)

#Plot Shannon diversity: BL vs EP
ggplot(Shannon_diversity, aes(x = Phase, y = Shannon)) + 
  geom_boxplot(fill = c("turquoise","violet"),width = 0.5) + 
  geom_signif(comparisons = list(c("BL","EP")), map_signif_level=TRUE, test = 't.test') +
  theme_bw() +
  ggtitle("Shannon Diversity") +
  theme(text = element_text(size = 18),aspect.ratio = 1.5)


# t-test fermentazioni

BL <- Shannon_diversity[samples_df["Phase"] == "BL",]
EP <- Shannon_diversity[samples_df["Phase"] == "EP",]

t.test(BL$Shannon,EP$Shannon)



#Plot Shannon diversity: T0 vs T6
ggplot(Shannon_diversity, aes(x = Time, y = Shannon)) + 
  geom_boxplot(fill = c("mediumseagreen","khaki1"),width = 0.5) + 
  geom_signif(comparisons = list(c("T0","T6")), map_signif_level=TRUE, test = 't.test') +
  theme_bw() +
  ggtitle("Shannon Diversity") +
  theme(text = element_text(size = 22),aspect.ratio = 1.5)

# t-test per tempi

T0 <- Shannon_diversity[samples_df["Time"] == "T0",]
T6 <- Shannon_diversity[samples_df["Time"] == "T6",]

t.test(T0$Shannon,T6$Shannon)




#Simpson diversity
Simpson_diversity <- estimate_richness(Global_data.rarefied,measures = "Simpson")



#Test di shapiro-wilk
shapiro.test(Simpson_diversity$Simpson)
Simpson_diversity$SampleType <- NA
Simpson_diversity$Phase <- NA
Simpson_diversity$Time <- NA

Simpson_diversity$SampleType[which(rownames(Simpson_diversity) %in% rownames(samples_df))] <- samples_df$SampleType 
Simpson_diversity$Phase[which(rownames(Simpson_diversity) %in% rownames(samples_df))] <- samples_df$Phase
Simpson_diversity$Time[which(rownames(Simpson_diversity) %in% rownames(samples_df))] <- samples_df$Time


#La distribuzione e gaussiana:faccio il t-test.

#boxplot(CTRL,LR,LRH,col = c("olivedrab2","deepskyblue2","firebrick2"),main = "Simpson's Diversity", xlab = "SampleType",
#ylab = "Simpson Index", names = c("CTRL","LR","LRH"), border = "black")

#Plot Simpson Diversity: CTRL vs LR vs LRH
ggplot(Simpson_diversity, aes(x = SampleType, y = Simpson)) + 
  geom_boxplot(fill = c("wheat","deepskyblue2","firebrick2"),width = 0.5) + 
  geom_signif(comparisons = list(c("CTRL", "LR")), map_signif_level=TRUE,test = "t.test") +
  geom_signif(comparisons = list(c("CTRL", "LRH")),map_signif_level=TRUE,test = "t.test",y_position = 0.995) +
  geom_signif(comparisons = list(c("LR", "LRH")),map_signif_level=TRUE,y_position = 1, test = 't.test')+
  theme_bw() +
  ggtitle("Simpson Diversity")+
  theme(text = element_text(size = 18),aspect.ratio = 1.5)
  
  

CTRL <- Simpson_diversity[samples_df["SampleType"] == "CTRL",]
LR <- Simpson_diversity[samples_df["SampleType"]== "LR",]
LRH <- Simpson_diversity[samples_df["SampleType"]== "LRH",]

t.test(CTRL$Simpson,LR$Simpson)
t.test(CTRL$Simpson,LRH$Simpson)

one.way.anova.simpson <- aov(Simpson ~ SampleType , data = Simpson_diversity )
summary(one.way.anova.simpson)

#Plot Simpson Diversity per fase di fermentazione
ggplot(Simpson_diversity, aes(x = Phase, y = Simpson)) + 
  geom_boxplot(fill = c("turquoise","violet"),width = 0.5) + 
  geom_signif(comparisons = list(c("BL", "EP")), map_signif_level=TRUE,test = "t.test") +
  theme_bw() +
  ggtitle("Simpson Diversity")+
  theme(text = element_text(size = 18),aspect.ratio = 1.5)


# t-test per fase di fermentazione

BL <- Simpson_diversity[samples_df["Phase"] == "BL",]
EP <- Simpson_diversity[samples_df["Phase"] == "EP",]

t.test(BL$Simpson,EP$Simpson)



#Plot per tempo di conservazione
ggplot(Simpson_diversity, aes(x = Time, y = Simpson)) + 
  geom_boxplot(fill = c("mediumseagreen","khaki1"),width = 0.5) + 
  geom_signif(comparisons = list(c("T0", "T6")), map_signif_level=TRUE,test = "t.test") +
  theme_bw() +
  ggtitle("Simpson Diversity")+
  theme(text = element_text(size = 22),aspect.ratio = 1.5)

# t-test per tempi

T0 <- Simpson_diversity[samples_df["Time"] == "T0",]
T6 <- Simpson_diversity[samples_df["Time"] == "T6",]

t.test(T0$Simpson,T6$Simpson)


#Chao1 Diversity
Chao1_diversity <- estimate_richness(Global_data.rarefied,measures = "Chao1")
#Chao1_diversity <- Chao1_diversity[,-2] # tolgo la seconda colonna che indica la standard 
# deviation. La prima colonna invece mi da l'indice di Chao1.

#boxplot(CTRL,LR,LRH,col = c("olivedrab2","deepskyblue2","firebrick2"),main = "Chao1 Diversity", xlab = "SampleType",
#ylab = "Chao1 Index", names = c("CTRL","LR","LRH"), border = "black")
Chao1_diversity$SampleType <- NA
Chao1_diversity$Phase <- NA
Chao1_diversity$Time <- NA

rownames(Chao1_diversity) <- rownames(Simpson_diversity)

Chao1_diversity$SampleType[which(rownames(Chao1_diversity) %in% rownames(samples_df))] <- samples_df$SampleType
Chao1_diversity$Phase[which(rownames(Chao1_diversity) %in% rownames(samples_df))] <- samples_df$Phase 
Chao1_diversity$Time[which(rownames(Chao1_diversity) %in% rownames(samples_df))] <- samples_df$Time 


shapiro.test(Chao1_diversity$Chao1)

#Plot Chao1 Diversity: CTRL vs LR vs LRH
ggplot(Chao1_diversity, aes(x = SampleType, y = Chao1)) + 
  geom_boxplot(fill = c("wheat","deepskyblue2","firebrick2"),width = 0.5) + 
  geom_signif(comparisons = list(c("CTRL", "LR")), map_signif_level=TRUE,test = "wilcox.test") +
  geom_signif(comparisons = list(c("CTRL", "LRH")),map_signif_level=TRUE,test = "wilcox.test",y_position = 300) +
  geom_signif(comparisons = list(c("LR", "LRH")),map_signif_level=TRUE,y_position = 330, test = 'wilcox.test')+
  theme_bw() +
  ggtitle("Chao1 Diversity")+
  theme(text = element_text(size = 18), aspect.ratio = 1.5)

CTRL <- Chao1_diversity[samples_df["SampleType"] == "CTRL",]
LR <- Chao1_diversity[samples_df["SampleType"]== "LR",]
LRH <- Chao1_diversity[samples_df["SampleType"]== "LRH",]

t.test(CTRL$Chao1,LR$Chao1)
t.test(CTRL$Chao1,LRH$Chao1)

one.way.anova.chao1 <- aov(Chao1 ~ SampleType , data = Chao1_diversity )
summary(one.way.anova.chao1)

#Plot Chao1 Diversity per fase di fermentazione
ggplot(Chao1_diversity, aes(x = Phase, y = Chao1)) + 
  geom_boxplot(fill = c("turquoise","violet"),width = 0.5) + 
  geom_signif(comparisons = list(c("EP", "BL")), map_signif_level=TRUE,test = "t.test") +
  theme_bw() +
  ggtitle("Chao1 Diversity")+
  theme(text = element_text(size = 18), aspect.ratio = 1.5)


# t-test per fase di fermentazione

BL <- Chao1_diversity[samples_df["Phase"] == "BL",]
EP <- Chao1_diversity[samples_df["Phase"] == "EP",]

t.test(BL$Chao1,EP$Chao1)


#Plot Chao1 per tempi
ggplot(Chao1_diversity, aes(x = Time, y = Chao1)) + 
  geom_boxplot(fill = c("mediumseagreen","khaki1"),width = 0.5) + 
  geom_signif(comparisons = list(c("T0", "T6")), map_signif_level=TRUE,test = "t.test") +
  theme_bw() +
  ggtitle("Chao1 Diversity")+
  theme(text = element_text(size = 22),aspect.ratio = 1.5)

# t-test per tempi

T0 <- Chao1_diversity[samples_df["Time"] == "T0",]
T6 <- Chao1_diversity[samples_df["Time"] == "T6",]

t.test(T0$Chao1,T6$Chao1)

