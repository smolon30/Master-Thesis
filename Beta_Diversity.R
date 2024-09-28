library(tidyr)
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
setwd("C:/Users/asus 147456/Desktop/1.TESI MAGISTRALE/Analisi batterica/Script_analisi_batterica/")

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

# FILTRO CAMPIONI: prendo solo EP o BL, T0 o T6

#ind1 <- grep(pattern = "0",samples_df$label)
#ind1 <- ind1 +1

#ind <- grep(pattern = "0",samples_df$label)


#ASV_table <- ASV_table[,-ind1]
#samples_df <- samples_df[-ind,]



## Filtraggio ASV dove ci sono tutti zeri


ind2 <- 0
num_rows <- nrow(ASV_table)
for (i in 1:num_rows){
  ind3 <- which(ASV_table[i,] == 0)
  if (length(ind3) == 24)(
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
  taxonomy[i,7] <- sub(']','',taxonomy[i,7])
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

# drarefy ci dice la probabilit che le specie siano presenti in una "comunita"
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


### BETA DIVERSITY: metodo Bray-Curtis: 
# B = 0 : indica che i due campioni non differiscono, hanno la stessa composizione
# B = 1 : indica che i due campioni hanno il massimo della diversit. 

# 1. Trasformo in abbondanze relative e creo un nuovo oggetto phyloseq.

#relative_phyloseq <- transform_sample_counts(Global_data.rarefied,function(x) x / sum(x))


# 2. CalcolO la distanza Bray Curtis tra i campioni e converto il risultato ad una matrice:
# questa "distanza" e compresa tra 0 e 1 e indica la dissimilarita tra due campioni 
# ( 0 = campioni uguali;)
# ( 1 = massimo della diversita.)

bray_matrix <- phyloseq::distance(Global_data.rarefied, method = "bray") # uso questa funzione

bray_matrix <- as.matrix(bray_matrix)

# Questa matrice rappresenta la dissimilarita (distanza) tra ogni campione.

# Noi vogliamo generare un boxplot rappresentando i gruppi separati, quindi dobbiamo 
# filtrare questa matrice.
# Con il codice seguente otteniamo la matrice df.bray nella quale 
# le distanze sono separate per gruppi.

sample_data(Global_data.rarefied)$SampleType <- factor((sample_data(Global_data.rarefied)$SampleType), levels=c("CTRL","LR","LRH"))

sub_dist <- list()
groups_all <- sample_data(Global_data.rarefied)$SampleType #creo una variabile che mi identifichi i sampletype

# levels(groups_all) va da 1 a 3(CTRL,LR,LRH), allora nel ciclo group va da 1 a 3
for (group in levels(groups_all)) { 
  row_group <- which(groups_all == group) #creo la variabile row_group che mi restituisce gli indici di dove sono gli elementi di quel determinato level.
  sample_group <- sample_names(Global_data.rarefied)[row_group] #creo una variabile alla quale assegno i nomi dei campioni, corrispondenti agli indici dei rowgroup
  sub_dist[[group]] <- bray_matrix[ sample_group, sample_group]# creo la matrice sub_dist, estraendo i valori da bray_matrix
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
} #con il ciclo for creo la variabile sub_dist che e una lista di 3 matrici rappresentanti ognuna le distanze
# tra campioni dello stesso tipo

braygroups<- melt(sub_dist) # melt mette insieme (fonde) sub_dist(che e una lista di 3) e crea un dataframe
df.bray <- braygroups[complete.cases(braygroups), ] #elimina le righe dove non ci sono valori, e quindi lascia soltanto le righe dove non ci sono NA e Nan.
df.bray$L1 <- factor(df.bray$L1, levels=names(sub_dist)) #creo una colonna L1 che mi tiene conto del tipo di campione.

# 84 variabili perche ci sono 28 coppie possibili per CTRL, 28 per LR e 28 per LRH. 
# Quindi il totale e 84 distanze. Ora le vogliamo visualizzare con un boxplot

ggplot(df.bray, aes(x=L1, y=value, colour = L1)) +
  ggtitle("Bray-Curtis Diversity") +
  geom_jitter(width = 0,) +
  geom_boxplot(alpha = 0.6,width = 0.4,borders = c("wheat","deepskyblue2","firebrick2")) +  
  theme(legend.position="none") +
  ylab("Bray-Curtis distance") +
  xlab("SampleType") +
  geom_signif(comparisons = list(c("CTRL", "LR")), map_signif_level=TRUE,test = "wilcox.test",colour = 'black') +
  geom_signif(comparisons = list(c("CTRL", "LRH")),map_signif_level=TRUE,test = "wilcox.test", y_position = 1.04,colour = 'black') +
  theme_classic()


### Altro grafico: PCoA plot

#Creo variabile con metodo Bray-Curtis
ord = ordinate(Global_data.rarefied, method="PCoA", distance = "bray")
# ordinate e una funzione di phyloseq che ordina i dati in base al metodo che scegliamo
# ci dobbiamo creare l'ord che contiene i nostri 24 campioni.

#PCoA con Bray-Curtis: CTRL vs LR vs LRH
plot_ordination(Global_data.rarefied, ord, color = "SampleType") + 
  geom_text(label = samples_df$label[1:24], size = 6, nudge_y = 0.03)+
  ggtitle("PCoA con metodo Bray-Curtis")+
  theme_bw(base_rect_size = 1)+
  theme( text = element_text(size = 20))+
  stat_ellipse(aes(group=SampleType)) 

#PERMANOVA
samples <- data.frame(sample_data(Global_data.rarefied))
adonis2(bray_matrix ~ SampleType, data = samples)

#PCoA con Bray-Curtis per fase di fermentazione
plot_ordination(Global_data.rarefied, ord, color = "Phase") + 
  geom_text(label = samples_df$label[1:24], size = 5, nudge_y = 0.045)+ 
  geom_point(size = 4)+
  stat_ellipse(aes(group=Phase)) +
  ggtitle("PCoA con metodo Bray-Curtis")+
  theme_bw()+
  theme(text = element_text(size = 20))

#PERMANOVA
adonis2(bray_matrix ~ Phase, data = samples)


#PCoA con Bray-Curtis per tempi
plot_ordination(Global_data.rarefied, ord, color = "Time") + 
  geom_text(label = samples_df$label[1:24], size = 5, nudge_y = 0.035)+
  geom_point(size=4) + 
  stat_ellipse(aes(group=Time)) +
  ggtitle("PCoA con metodo Bray-Curtis: T0 vs T6")+
  theme_bw()+
  theme(text = element_text(size = 20))

# PERMANOVA
adonis2(bray_matrix ~ Time, data = samples)





#METODO UNWEIGHTED UNIFRAC

options(getClass.msg=FALSE)
#creo l'albero filogenetico come oggetto phyloseq
random_tree = rtree(ntaxa(Global_data.rarefied), rooted=TRUE, tip.label=taxa_names(Global_data.rarefied))
#plot(random_tree)

#Creo variabile con metodo Unweighted-Unifrac
Global_data.rarefied <- merge_phyloseq(Global_data.rarefied,random_tree)
ord = ordinate(Global_data.rarefied, method="PCoA", distance = "unifrac")
# ordinate e una funzione di phyloseq che ordina i dati in base al metodo che scegliamo
# ci dobbiamo creare l'ord che contiene i nostri 24 campioni.
# In questo caso uso "unifrac" che mi ordinate secondo il metodo Unweighted-Unifrac

#PCoA con Unweighted-Unifrac: CTRL vs LR vs LRH
plot_ordination(Global_data.rarefied, ord, color = "SampleType") + 
  geom_text(label = samples_df$label[1:24], size = 6, nudge_y = 0.03)+
  stat_ellipse(aes(group=SampleType)) +
  ggtitle("PCoA con metodo unweighted-UniFrac")+
  theme_bw(base_rect_size = 1)+
  theme(text = element_text(size = 20))

uni_matrix <- phyloseq::distance(Global_data.rarefied, method = "UniFrac")
uni_matrix <- as.matrix(uni_matrix)

#PERMANOVA
adonis2(uni_matrix ~ SampleType, data = samples)

# PCoA per "Phase"
sampletype <- samples_df$SampleType
sampletype <- sampletype[-18]
sampletype <- sampletype[-18]

#PCoA con Unweighted-Unifrac per fase di fermentazione
plot_ordination(Global_data.rarefied, ord, color = "Phase") + 
  geom_text(label = samples_df$label[1:24], size = 5, nudge_y = 0.04)+
  geom_point(size=4) + 
  stat_ellipse(aes(group=Phase)) +
  ggtitle("PCoA con metodo unweighted-UniFrac")+
  theme_bw()+
  theme(text = element_text(size = 20))

#PERMANOVA
adonis2(uni_matrix ~ Phase, data = samples)


#PCoA con Unweighted-Unifrac per tempi
plot_ordination(Global_data.rarefied, ord, color = "Time") + 
  geom_text(label = samples_df$label[1:24], size = 5, nudge_y = 0.035)+ 
  geom_point(size = 4)+
  stat_ellipse(aes(group=Time)) +
  ggtitle("PCoA con metodo unweighted-UniFrac: T0 vs T6")+
  theme_bw()+
  theme(text = element_text(size = 20))

#PERMANOVA
adonis2(uni_matrix ~ Time, data = samples)

## The Bray-Curtis dissimilarity is based on occurrence data (abundance), while
#the Jaccard distance is based on presence/absence data (does not include abundance information). 
#UniFrac distances take into account the occurrence table and the phylogeny diversity (sequence distance).

