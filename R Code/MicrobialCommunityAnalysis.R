#R code for the manuscript: Temporal variability of microbial response to crude oil exposure in the northern Gulf of Mexico

#Author: Melissa L Brock

#Code for microbial community analysis, community statistics, and figure building 

####Setup R Environment####
#Set working directory
getwd()
setwd("/path/to/directory")

#Load libraries
library(dada2)
library(seqinr)
library(phyloseq)
library(DESeq2)
library(hexbin)
library(fpc)
library(vegan)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dendextend)
library(RVAideMemoire)
library(car)
library(extrafont)
library(extrafontdb)
library(cowplot)
library(MicEco)
library(ggplotify)
library(betapart)

##Download fasta file from QIIME2 view using the Brock-2023-16S-rep-seqs.qzv file
#Convert Fasta file to csv
parsed <- read.fasta(file('sequences_365.fasta'), as.string = TRUE)
table <- data.frame(unlist(parsed), row.names = sapply(parsed, attr, 'Annot'))
write.csv(table, "sequences_365.csv")

#Import rep seq table
SS_Y1_ref_sequences <- read.csv("sequences_365.csv", header = T, sep = ",")

##Make ASV count table
#In terminal, export feature table from QIIME2 in BIOM format: qiime tools export --input-path Brock-16S-feature-table.qza --output-path 16S-table.biom
#Convert table from BIOM format to tsv format: biom convert -i 16S-table.biom/feature-table.biom -o asv_table.txt --to-tsv

##Import count table
SS_Y1_count_table <- read.csv("asv_table.csv", header = T, sep = ",")

#Make sure column names for OTU IDs match and that row numbers are the same
colnames(SS_Y1_ref_sequences)
colnames(SS_Y1_count_table)

dim(SS_Y1_ref_sequences)
#4597 2
dim(SS_Y1_count_table)
#4597 27

#Sort based on OTU IDs
SS_Y1_ref_sequences <- SS_Y1_ref_sequences[order(as.numeric(SS_Y1_ref_sequences$OTU.ID)),]
SS_Y1_count_table <- SS_Y1_count_table[order(as.numeric(SS_Y1_count_table$OTU.ID)),]

#Merge tables together and rearrange
SS_Y1_asv_table <- merge(SS_Y1_ref_sequences, SS_Y1_count_table, by = "OTU.ID")

rownames(SS_Y1_asv_table) <- SS_Y1_asv_table$Sequence
head(rownames(SS_Y1_asv_table))

colnames(SS_Y1_asv_table)
SS_Y1_asv_table <- SS_Y1_asv_table[ , -c(1:2)]
SS_Y1_asv_table_trans <- as.matrix(t(SS_Y1_asv_table))
colnames(SS_Y1_asv_table_trans) <- toupper(colnames(SS_Y1_asv_table_trans))
head(colnames(SS_Y1_asv_table_trans))


####Assign taxonomy using RDP's Naive Bayesian classifier and SILVA138 reference database####
SS_Y1_taxa_SILVA138 <- assignTaxonomy(SS_Y1_asv_table_trans, "silva_nr99_v138.1_train_set.fa.gz", verbose = TRUE, minBoot = 80, multithread = TRUE, outputBootstraps = TRUE)

####Manipulate taxa table and import metadata####
#Manipulate taxa chart
SS_Y1_taxa_SILVA138.df <- as.data.frame(SS_Y1_taxa_SILVA138)

colnames(SS_Y1_taxa_SILVA138.df)
SS_Y1_taxa_SILVA138.df2 <- SS_Y1_taxa_SILVA138.df[ , -c(7:12)]
colnames(SS_Y1_taxa_SILVA138.df2) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
SS_Y1_taxa_SILVA138.m <- as.matrix(SS_Y1_taxa_SILVA138.df2)

#Import metadata (Supplemental Table 1)
SS_Y1_metadata <- read.csv("SupplementalTable1.csv", header = T, row.names = 1, sep = ",")
View(SS_Y1_metadata)

#Add color to variables - Season
SS_Y1_metadata$Color_Season[SS_Y1_metadata$Season == "Spring"] <- "aquamarine4"
SS_Y1_metadata$Color_Season[SS_Y1_metadata$Season == "Summer"] <- "firebrick3"
SS_Y1_metadata$Color_Season[SS_Y1_metadata$Season == "Fall"] <- "darkorange2"
SS_Y1_metadata$Color_Season[SS_Y1_metadata$Season == "Winter"] <- "blue4"

####Determine proportion of taxa belonging to each Kingdom (L1)####
#Build phyloseq object
SS_Y1_metadata.phy <- sample_data(SS_Y1_metadata)
SS_Y1_asv_table.phy <- otu_table(SS_Y1_asv_table_trans, taxa_are_rows = F)
SS_Y1_taxa_Silva138.phy <- tax_table(SS_Y1_taxa_SILVA138.m)
SS_Y1138.physeq <- phyloseq(SS_Y1_taxa_Silva138.phy, SS_Y1_asv_table.phy, SS_Y1_metadata.phy)

#Make count table at the Kingdom level
SS_Y1_counts138 <- t(otu_table(tax_glom(SS_Y1138.physeq, taxrank = "Kingdom")))

#Agglomerate the taxa at the Kingdom level
SS_Y1_glom138 <- tax_glom(SS_Y1138.physeq, taxrank = "Kingdom")

#Store the Kingdom names in a vector
SS_Y1_vec138 <- as.data.frame(SS_Y1_glom138@tax_table@.Data)$Kingdom

#Assign the Kingdom names as row names in the count table
row.names(SS_Y1_counts138) <- as.vector(SS_Y1_vec138)
row.names(SS_Y1_counts138)

#Calculate the number of unannotated reads and add it to the count table
SS_Y1_unannot_counts138 <- colSums(t(SS_Y1_asv_table_trans)) - colSums(SS_Y1_counts138)
SS_Y1_kingdomandunannot138 <- rbind(SS_Y1_counts138, "Unannotated" = SS_Y1_unannot_counts138)

#Calculate the relative abundance of each Kingdom
SS_Y1_proportions138 <- apply(SS_Y1_kingdomandunannot138, 2, function(x) x/sum(x)*100)

#Get the dimensions of the table
dim(SS_Y1_proportions138)
#4 26

#View the table
View(SS_Y1_proportions138)
write.csv(SS_Y1_proportions138, "Seasonal Study 16S Kingdom Proportions 138.csv")

####Remove Eukaryotes, Archaea, and Unannotated####
View(SS_Y1_taxa_SILVA138.df2) #Nucleotide string is row name
View(SS_Y1_asv_table_trans) #Nucleotide string is column name

#Put nucleotide string of Eukaryotes, Archaea, and Unannotated into separate data frames
SS_Y1_Eukaryotes138.df <- data.frame(rownames(SS_Y1_taxa_SILVA138.df2)[which(grepl("Eukaryota", SS_Y1_taxa_SILVA138.df2$Kingdom))])
View(SS_Y1_Eukaryotes138.df)
colnames(SS_Y1_Eukaryotes138.df) <- "Sequence"

SS_Y1_Archaea138.df <- data.frame(rownames(SS_Y1_taxa_SILVA138.df2)[which(grepl("Archaea", SS_Y1_taxa_SILVA138.df2$Kingdom))])
View(SS_Y1_Archaea138.df)
colnames(SS_Y1_Archaea138.df) <- "Sequence"

SS_Y1_Unannotated.k.138.df <- data.frame(rownames(SS_Y1_taxa_SILVA138.df2)[which(is.na(SS_Y1_taxa_SILVA138.df2$Kingdom))])
View(SS_Y1_Unannotated.k.138.df)
colnames(SS_Y1_Unannotated.k.138.df) <- "Sequence"

#Put row numbers corresponding to Eukaryotes, Archaea, and Unannotated into separate data frames
SS_Y1_Eukaryotes138.rn <- data.frame(which(grepl("Eukaryota", SS_Y1_taxa_SILVA138.df2$Kingdom)))
View(SS_Y1_Eukaryotes138.rn)
colnames(SS_Y1_Eukaryotes138.rn) <- "RowNumber"

SS_Y1_Archaea138.rn <- data.frame(which(grepl("Archaea", SS_Y1_taxa_SILVA138.df2$Kingdom)))
View(SS_Y1_Archaea138.rn)
colnames(SS_Y1_Archaea138.rn) <- "RowNumber"

SS_Y1_Unannotated.k.138.rn <- data.frame(which(is.na(SS_Y1_taxa_SILVA138.df2$Kingdom)))
View(SS_Y1_Unannotated.k.138.rn)
colnames(SS_Y1_Unannotated.k.138.rn) <- "RowNumber"

#Identify the row numbers of the Eukaryotes, Archaea, and Unannotated in the ASV table and save them in separate data frames
SS_Y1_EuksinASVTable138 <- data.frame(which(colnames(SS_Y1_asv_table_trans) %in% SS_Y1_Eukaryotes138.df$Sequence))
View(SS_Y1_EuksinASVTable138)
colnames(SS_Y1_EuksinASVTable138) <- "Sequence"

SS_Y1_ArchinASVTable138 <- data.frame(which(colnames(SS_Y1_asv_table_trans) %in% SS_Y1_Archaea138.df$Sequence))
colnames(SS_Y1_ArchinASVTable138) <- "Sequence"

SS_Y1_UnannKinASVTable138 <- data.frame(which(colnames(SS_Y1_asv_table_trans) %in% SS_Y1_Unannotated.k.138.df$Sequence))
colnames(SS_Y1_UnannKinASVTable138) <- "Sequence"

#Combine dataframes
SS_Y1_ASV138.remove <- rbind(SS_Y1_EuksinASVTable138, SS_Y1_ArchinASVTable138, SS_Y1_UnannKinASVTable138)
SS_Y1_Taxa.remove138 <- rbind(SS_Y1_Eukaryotes138.rn, SS_Y1_Archaea138.rn, SS_Y1_Unannotated.k.138.rn)

#Remove the Eukaryotes, Archaea, and Unannotated Kingdom from the ASV table
SS_Y1_Proks_asv_table138 <- SS_Y1_asv_table_trans[, -c(SS_Y1_ASV138.remove$Sequence)]
dim(SS_Y1_Proks_asv_table138)
dim(SS_Y1_ASV138.remove)
dim(SS_Y1_asv_table_trans)

#Remove the Eukaryotes, Archaea, and Unannotated Kingdom from the taxa table
SS_Y1_Proks_taxa138 <- SS_Y1_taxa_SILVA138.df2[-c(SS_Y1_Taxa.remove138$RowNumber), ]
dim(SS_Y1_taxa_SILVA138.df2)
dim(SS_Y1_Proks_taxa138)
dim(SS_Y1_Taxa.remove138)

SS_Y1_Proks_taxa138.m <- as.matrix(SS_Y1_Proks_taxa138)

####Build taxonomic table at the Kingdom level to check that it contains only Bacteria ####
#Build phyloseq object
SS_Y1_Proks_metadata.phy <- sample_data(SS_Y1_metadata)
SS_Y1_Proks_asv_table138.phy <- otu_table(SS_Y1_Proks_asv_table138, taxa_are_rows = F)
SS_Y1_Proks_taxa_Silva138.phy <- tax_table(SS_Y1_Proks_taxa138.m)
SS_Y1_Proks138.physeq <- phyloseq(SS_Y1_Proks_taxa_Silva138.phy, SS_Y1_Proks_asv_table138.phy, SS_Y1_Proks_metadata.phy)

#Make count table at the Kingdom level
SS_Y1_Proks_counts138 <- t(otu_table(tax_glom(SS_Y1_Proks138.physeq, taxrank = "Kingdom")))

#Agglomerate the taxa at the Kingdom level
SS_Y1_Proks_glom138 <- tax_glom(SS_Y1_Proks138.physeq, taxrank = "Kingdom")

#Store the Kingdom names in a vector
SS_Y1_Proks_vec138 <- as.data.frame(SS_Y1_Proks_glom138@tax_table@.Data)$Kingdom

#Assign the Kingdom names as row names in the count table
row.names(SS_Y1_Proks_counts138) <- as.vector(SS_Y1_Proks_vec138)
row.names(SS_Y1_Proks_counts138)

#Calculate the number of unannotated reads and add it to the count table
SS_Y1_Proks_unannot_counts138 <- colSums(t(SS_Y1_Proks_asv_table138)) - colSums(SS_Y1_Proks_counts138)
SS_Y1_Proks_kingdomandunannot138 <- rbind(SS_Y1_Proks_counts138, "Unannotated" = SS_Y1_Proks_unannot_counts138)

#Calculate the relative abundance of each Kingdom
SS_Y1_Proks_proportions138 <- apply(SS_Y1_Proks_kingdomandunannot138, 2, function(x) x/sum(x)*100)

#Get the dimensions of the table
dim(SS_Y1_Proks_proportions138)
#2 26

#View the table
View(SS_Y1_Proks_proportions138)

####Remove Unannotated Phyla####
View(SS_Y1_Proks_taxa138) #Nucleotide string is row name
View(SS_Y1_Proks_asv_table138) #Nucleotide string is column name

#Put nucleotide string of Unannotated into a data frame
SS_Y1_Unannotated.p.138.df <- data.frame(rownames(SS_Y1_Proks_taxa138)[which(is.na(SS_Y1_Proks_taxa138$Phylum))])
colnames(SS_Y1_Unannotated.p.138.df) <- "Sequence"

#Put row numbers corresponding to Unannotated into a data frame
SS_Y1_Unannotated.p.138.rn <- data.frame(which(is.na(SS_Y1_Proks_taxa138$Phylum)))
colnames(SS_Y1_Unannotated.p.138.rn) <- "RowNumber"

#Identify the row numbers of the Unannotated Phyla in the ASV table and save them in a data frame
SS_Y1_UnannPinASVTable138 <- data.frame(which(colnames(SS_Y1_Proks_asv_table138) %in% SS_Y1_Unannotated.p.138.df$Sequence))
colnames(SS_Y1_UnannPinASVTable138) <- "Sequence"

#Remove the Unannotated Phyla from the ASV table
SS_Y1_Proks_asv_table2.138 <- SS_Y1_Proks_asv_table138[, -c(SS_Y1_UnannPinASVTable138$Sequence)]
dim(SS_Y1_Proks_asv_table2.138)
dim(SS_Y1_UnannPinASVTable138)
dim(SS_Y1_Proks_asv_table138)

#Remove the Unannotated Phyla from the taxa table
SS_Y1_Proks_taxa2.138 <- SS_Y1_Proks_taxa138[-c(SS_Y1_Unannotated.p.138.rn$RowNumber), ]
dim(SS_Y1_Unannotated.p.138.rn)
dim(SS_Y1_Proks_taxa2.138)
dim(SS_Y1_Proks_taxa138)

SS_Y1_Proks_taxa2.138.m <- as.matrix(SS_Y1_Proks_taxa2.138)

####Build taxonomic table at the Phylum level to check that Unannotated have been removed####
#Build phyloseq object
SS_Y1_Proks_asv_table2.138.phy <- otu_table(SS_Y1_Proks_asv_table2.138, taxa_are_rows = F)
SS_Y1_Proks_taxa2_Silva138.phy <- tax_table(SS_Y1_Proks_taxa2.138.m)
SS_Y1_Proks2.138.physeq <- phyloseq(SS_Y1_Proks_taxa2_Silva138.phy, SS_Y1_Proks_asv_table2.138.phy, SS_Y1_Proks_metadata.phy)

#Make count table at the Phylum level
SS_Y1_Proks2_counts138 <- t(otu_table(tax_glom(SS_Y1_Proks2.138.physeq, taxrank = "Phylum")))

#Agglomerate the taxa at the Kingdom level
SS_Y1_Proks2_glom138 <- tax_glom(SS_Y1_Proks2.138.physeq, taxrank = "Phylum")

#Store the Kingdom names in a vector
SS_Y1_Proks2_vec138 <- as.data.frame(SS_Y1_Proks2_glom138@tax_table@.Data)$Phylum

#Assign the Kingdom names as row names in the count table
row.names(SS_Y1_Proks2_counts138) <- as.vector(SS_Y1_Proks2_vec138)
row.names(SS_Y1_Proks2_counts138)

#Calculate the number of unannotated reads and add it to the count table
SS_Y1_Proks2_unannot_counts138 <- colSums(t(SS_Y1_Proks_asv_table2.138)) - colSums(SS_Y1_Proks2_counts138)
SS_Y1_Proks2_kingdomandunannot138 <- rbind(SS_Y1_Proks2_counts138, "Unannotated" = SS_Y1_Proks2_unannot_counts138)

#Calculate the relative abundance of each Kingdom
SS_Y1_Proks2_proportions138 <- apply(SS_Y1_Proks2_kingdomandunannot138, 2, function(x) x/sum(x)*100)

dim(SS_Y1_Proks2_proportions138)
#37 26

#View the table
View(SS_Y1_Proks2_proportions138)

####Taxa Barplot and Pie Charts - Genera####
SS_Y1_genus_counts138 <- t(otu_table(tax_glom(SS_Y1_Proks2.138.physeq, taxrank = "Genus")))
SS_Y1_genus_glom138 <- tax_glom(SS_Y1_Proks2.138.physeq, taxrank = "Genus")
SS_Y1_genus_vec138 <- as.data.frame(SS_Y1_genus_glom138@tax_table@.Data)$Genus
row.names(SS_Y1_genus_counts138) <- as.vector(SS_Y1_genus_vec138)
row.names(SS_Y1_genus_counts138)
SS_Y1_genus_unannot_counts138 <- colSums(t(SS_Y1_Proks_asv_table2.138)) - colSums(SS_Y1_genus_counts138)
SS_Y1_genusandunannot138 <- rbind(SS_Y1_genus_counts138, "Unassigned" = SS_Y1_genus_unannot_counts138)
SS_Y1_genus_proportions138 <- apply(SS_Y1_genusandunannot138, 2, function(x) x/sum(x)*100)
dim(SS_Y1_genus_proportions138)
#346 26
range(rowSums(SS_Y1_genus_proportions138))
colSums(SS_Y1_genus_proportions138)

SS_Y1_genus_proportions_sorted138 <- as.data.frame(SS_Y1_genus_proportions138[(order(rowSums(SS_Y1_genus_proportions138), decreasing = T)), ])
SS_Y1_genus_proportions138_top50 <- SS_Y1_genus_proportions_sorted138[c(1:51), ]
dim(SS_Y1_genus_proportions138_top50)
rownames(SS_Y1_genus_proportions138_top50)

SS_Y1_genus_proportions138_top25 <- SS_Y1_genus_proportions_sorted138[c(1:26), ]
dim(SS_Y1_genus_proportions138_top)
rownames(SS_Y1_genus_proportions138_top)

SS_Y1_genus_proportions_filt138_top25 <- colSums(SS_Y1_genus_proportions138) - colSums(SS_Y1_genus_proportions138_top25)
SS_Y1_genus_proportions_filt_tab138_top25 <- rbind(SS_Y1_genus_proportions138_top25, "Other" = SS_Y1_genus_proportions_filt138_top25)
rownames(SS_Y1_genus_proportions_filt_tab138_top25)
SS_Y1_genus_proportions_filt_tab2.138_top25 <- SS_Y1_genus_proportions_filt_tab138_top25[-1, ]
colSums(SS_Y1_genus_proportions_filt_tab2.138_top25)

SS_Y1_genus_proportions_filt_tab4plot138_top25 <- SS_Y1_genus_proportions_filt_tab2.138_top25
SS_Y1_genus_proportions_filt_tab4plot138_top25$Genus <- row.names(SS_Y1_genus_proportions_filt_tab4plot138_top25)
SS_Y1_genus_proportions_filt_tab4plot138_top25.g <- gather(SS_Y1_genus_proportions_filt_tab4plot138_top25, Sample, Proportion, -Genus)

SS_Y1_metadata$CollectionDate <- as.Date(SS_Y1_metadata$CollectionDate, format='%m/%d/%y')
colnames(SS_Y1_metadata)
SS_Y1_metadata$SampleID <- row.names(SS_Y1_metadata)

SS_Y1_sample_info <- data.frame("Sample" = SS_Y1_metadata$SampleID, "Season" = SS_Y1_metadata$Season, "Date" = SS_Y1_metadata$CollectionDate, stringsAsFactors = F)

SS_Y1_genus_proportions_filt_tab4plot138_top25.g2 <- merge(SS_Y1_genus_proportions_filt_tab4plot138_top25.g, SS_Y1_sample_info)

#Figure 3a
sort(unique(SS_Y1_genus_proportions_filt_tab4plot138_top25.g2$Date))

SS_Y1_taxaplot138.2_L6_top25 <- ggplot(SS_Y1_genus_proportions_filt_tab4plot138_top25.g2, aes(x = Date, y = Proportion, fill = Genus)) +
  geom_bar(width = 6, stat = "identity") + 
  geom_vline(xintercept = c(sort(as.Date("2015-09-18", "%Y-%m-%d")), sort(as.Date("2015-11-27", "%Y-%m-%d")), sort(as.Date("2016-03-21", "%Y-%m-%d")), sort(as.Date("2016-06-24", "%Y-%m-%d"))), linetype = "solid", linewidth = 0.5) + 
  annotate("rect", xmin = sort(as.Date("2015-09-11", "%Y-%m-%d")), xmax = sort(as.Date("2016-09-30", "%Y-%m-%d")), ymin = 73, ymax = Inf, fill = "white") + 
  annotate("text", y = c(78, 78, 78, 78, 78), x = c(sort(as.Date("2015-09-04", "%Y-%m-%d")), sort(as.Date("2015-10-23", "%Y-%m-%d")), sort(as.Date("2016-01-22", "%Y-%m-%d")), sort(as.Date("2016-05-06", "%Y-%m-%d")), sort(as.Date("2016-08-05", "%Y-%m-%d"))), label = c(" ", "Fall", "Winter", "Spring", "Summer"), size = 4.5, hjust = 0.5, family = "Times New Roman") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_text(size = 14), legend.text = element_text(size = 12), legend.key.size = unit(3.5, "mm"), text = element_text(family = "Times New Roman"), legend.text.align = 0, plot.title = element_text(size = 16, hjust = 0.5)) +  
  labs(x = "Date", y = "Relative Abundance (%)") +
  expand_limits(y = c(0, 80)) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_y_continuous(limits = c(0, 80), breaks = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70)) + 
  scale_fill_manual(guide = guide_legend(ncol = 1), name = "Genus:", values = c("pink4", "darkorange3", "pink1", "darkorange1", "goldenrod1", "lightgoldenrod1", "darkseagreen4", "darkseagreen1", "deeppink4", "deeppink1", "orchid4", "royalblue4", "royalblue1", "plum3", "violet", "plum1", "black", "seashell4", "seashell1", "turquoise4", "palegreen4", "lightsalmon1", "springgreen3", "turquoise", "tan4", "wheat3"), breaks = c("Acholeplasma", "Ascidiaceihabitans", "Aurantivirga", "Balneola", "Blastopirellula", "Candidatus Actinomarina", "Candidatus Aquiluna", "Clade Ia", "Clade Ib", "Cyanobium PCC-6307", "Fluviicola", "Formosa", "HIMB11", "MB11C04 marine group", "ML602J-51", "NS2b marine group", "NS4 marine group", "NS5 marine group", "OM43 clade", "OM60(NOR5) clade", "Other", "Rubripirellula", "Sva0996 marine group", "Synechococcus CC9902", "Winogradskyella", "Woeseia"), labels = c("Acholeplasma", "Ascidiaceihabitans", "Aurantivirga", "Balneola", "Blastopirellula", "Candidatus Actinomarina", "Candidatus Aquiluna", "SAR11 Clade Ia", "SAR11 Clade Ib", "Cyanobium PCC-6307", "Fluviicola", "Formosa", "HIMB11", "MB11C04 marine group", "ML602J-51", "NS2b marine group", "NS4 marine group", "NS5 marine group", "OM43 clade", "OM60(NOR5) clade", "Other", "Rubripirellula", "Sva0996 marine group", "Synechococcus", "Winogradskyella", "Woeseia")) 


#Pie Charts
SS_Y1_genus_proportions_filt_tab4plot138 <- SS_Y1_genus_proportions_filt_tab2.138
SS_Y1_genus_proportions_filt_tab4plot138$Genus <- row.names(SS_Y1_genus_proportions_filt_tab4plot138)
SS_Y1_genus_proportions_filt_tab4plot138.g <- gather(SS_Y1_genus_proportions_filt_tab4plot138, Sample, Proportion, -Genus)
SS_Y1_genus_proportions_filt_tab4plot138.g2 <- merge(SS_Y1_genus_proportions_filt_tab4plot138.g, SS_Y1_sample_info)

unique(SS_Y1_genus_proportions_filt_tab4plot138.g2$Season)
SS_Y1_genus_proportions_summer138 <- SS_Y1_genus_proportions_filt_tab4plot138.g2[which(SS_Y1_genus_proportions_filt_tab4plot138.g2$Season == "Summer"), ]
SS_Y1_genus_proportions_fall138 <- SS_Y1_genus_proportions_filt_tab4plot138.g2[which(SS_Y1_genus_proportions_filt_tab4plot138.g2$Season == "Fall"), ]
SS_Y1_genus_proportions_winter138 <- SS_Y1_genus_proportions_filt_tab4plot138.g2[which(SS_Y1_genus_proportions_filt_tab4plot138.g2$Season == "Winter"), ]
SS_Y1_genus_proportions_spring138 <- SS_Y1_genus_proportions_filt_tab4plot138.g2[which(SS_Y1_genus_proportions_filt_tab4plot138.g2$Season == "Spring"), ]

mean_by_Season_summer138 <- SS_Y1_genus_proportions_summer138 %>% 
  group_by(Genus) %>%
  summarize(averaged.Prop = mean(Proportion))
View(mean_by_Season_summer138)
#Identify genera with average relative abundance >= 1%

mean_by_Season_spring138 <- SS_Y1_genus_proportions_spring138 %>% 
  group_by(Genus) %>%
  summarize(averaged.Prop = mean(Proportion))
View(mean_by_Season_spring138)
#Identify genera with average relative abundance >= 1%

mean_by_Season_fall138 <- SS_Y1_genus_proportions_fall138 %>% 
  group_by(Genus) %>%
  summarize(averaged.Prop = mean(Proportion))
View(mean_by_Season_fall138)
#Identify genera with average relative abundance >= 1%

mean_by_Season_winter138 <- SS_Y1_genus_proportions_winter138 %>% 
  group_by(Genus) %>%
  summarize(averaged.Prop = mean(Proportion))
View(mean_by_Season_winter138)
#Identify genera with average relative abundance >= 1%

mean_by_Season_spring138$Genus
mean_by_Season_summer138$Genus
mean_by_Season_fall138$Genus
mean_by_Season_winter138$Genus

row.names(mean_by_Season_spring138) #Pull out genera identified above
mean_by_Season_spring_sub138 <- mean_by_Season_spring138[c(4, 6, 7, 8, 9, 12, 13, 15, 17, 20, 27, 30, 31, 32, 34, 35, 47), ]
mean_by_Season_summer_sub138 <- mean_by_Season_summer138[c(4, 6, 7, 8, 9, 12, 13, 15, 17, 20, 27, 30, 31, 32, 34, 35, 47), ]
mean_by_Season_fall_sub138 <- mean_by_Season_fall138[c(4, 6, 7, 8, 9, 12, 13, 15, 17, 20, 27, 30, 31, 32, 34, 35, 47), ]
mean_by_Season_winter_sub138 <- mean_by_Season_winter138[c(4, 6, 7, 8, 9, 12, 13, 15, 17, 20, 27, 30, 31, 32, 34, 35, 47), ]

#Build plots
sum(mean_by_Season_spring_sub138$averaged.Prop)
sort(unique(mean_by_Season_spring_sub138$Genus))
spring.avg.genus138 <- ggplot(mean_by_Season_spring_sub138, aes(x = "", y = averaged.Prop, fill = Genus))+
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0, direction = 1) + 
  theme_minimal() + 
  labs(title = "Spring") + 
  theme(axis.title =element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12), legend.key.size = unit(3.5, "mm"), text = element_text(family = "Times New Roman"), legend.text.align = 0, plot.title = element_text(size = 14, hjust = 0.5)) + 
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 5)) +
  scale_fill_manual(guide = guide_legend(ncol = 5), name = "Genus:", values = c("darkorange3", "darkorange1", "goldenrod1", "lightgoldenrod1", "darkseagreen4", "darkseagreen1", "deeppink4", "deeppink1", "royalblue4", "royalblue1", "plum3", "plum1", "black", "seashell4", "seashell1", "turquoise4", "turquoise"), breaks = c("Ascidiaceihabitans", "Balneola", "Blastopirellula", "Candidatus Actinomarina", "Candidatus Aquiluna", "Clade Ia", "Clade Ib", "Cyanobium PCC-6307", "Formosa", "HIMB11", "MB11C04 marine group", "NS2b marine group", "NS4 marine group", "NS5 marine group", "OM43 clade", "OM60(NOR5) clade", "Synechococcus CC9902"), labels = c("Ascidiaceihabitans", "Balneola", "Blastopirellula", "Candidatus Actinomarina", "Candidatus Aquiluna", "SAR11 Clade Ia", "SAR11 Clade Ib", "Cyanobium PCC-6307", "Formosa", "HIMB11", "MB11C04 marine group", "NS2b marine group", "NS4 marine group", "NS5 marine group", "OM43 clade", "OM60(NOR5) clade", "Synechococcus"))
spring.avg.genus138

sum(mean_by_Season_summer_sub138$averaged.Prop)
summer.avg.genus138 <- ggplot(mean_by_Season_summer_sub138, aes(x = "", y = averaged.Prop, fill = Genus))+
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0, direction = 1) + 
  theme_minimal() +
  labs(title = "Summer") + 
  theme(axis.title =element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12), legend.key.size = unit(3.5, "mm"), text = element_text(family = "Times New Roman"), legend.text.align = 0, plot.title = element_text(size = 14, hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 5)) +
  scale_fill_manual(guide = guide_legend(ncol = 5), name = "Genus:", values = c("darkorange3", "darkorange1", "goldenrod1", "lightgoldenrod1", "darkseagreen4", "darkseagreen1", "deeppink4", "deeppink1", "royalblue4", "royalblue1", "plum3", "plum1", "black", "seashell4", "seashell1", "turquoise4", "turquoise"), breaks = c("Ascidiaceihabitans", "Balneola", "Blastopirellula", "Candidatus Actinomarina", "Candidatus Aquiluna", "Clade Ia", "Clade Ib", "Cyanobium PCC-6307", "Formosa", "HIMB11", "MB11C04 marine group", "NS2b marine group", "NS4 marine group", "NS5 marine group", "OM43 clade", "OM60(NOR5) clade", "Synechococcus CC9902"), labels = c("Ascidiaceihabitans", "Balneola", "Blastopirellula", "Candidatus Actinomarina", "Candidatus Aquiluna", "SAR11 Clade Ia", "SAR11 Clade Ib", "Cyanobium PCC-6307", "Formosa", "HIMB11", "MB11C04 marine group", "NS2b marine group", "NS4 marine group", "NS5 marine group", "OM43 clade", "OM60(NOR5) clade", "Synechococcus"))
summer.avg.genus138

sum(mean_by_Season_fall_sub138$averaged.Prop)
fall.avg.genus138 <- ggplot(mean_by_Season_fall_sub138, aes(x = "", y = averaged.Prop, fill = Genus))+
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0, direction = 1) + 
  theme_minimal() + 
  labs(title = "Fall") + 
  theme(axis.title =element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12), legend.key.size = unit(3.5, "mm"), text = element_text(family = "Times New Roman"), legend.text.align = 0, plot.title = element_text(size = 14, hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 5)) +
  scale_fill_manual(guide = guide_legend(ncol = 5), name = "Genus:", values = c("darkorange3", "darkorange1", "goldenrod1", "lightgoldenrod1", "darkseagreen4", "darkseagreen1", "deeppink4", "deeppink1", "royalblue4", "royalblue1", "plum3", "plum1", "black", "seashell4", "seashell1", "turquoise4", "turquoise"), breaks = c("Ascidiaceihabitans", "Balneola", "Blastopirellula", "Candidatus Actinomarina", "Candidatus Aquiluna", "Clade Ia", "Clade Ib", "Cyanobium PCC-6307", "Formosa", "HIMB11", "MB11C04 marine group", "NS2b marine group", "NS4 marine group", "NS5 marine group", "OM43 clade", "OM60(NOR5) clade", "Synechococcus CC9902"), labels = c("Ascidiaceihabitans", "Balneola", "Blastopirellula", "Candidatus Actinomarina", "Candidatus Aquiluna", "SAR11 Clade Ia", "SAR11 Clade Ib", "Cyanobium PCC-6307", "Formosa", "HIMB11", "MB11C04 marine group", "NS2b marine group", "NS4 marine group", "NS5 marine group", "OM43 clade", "OM60(NOR5) clade", "Synechococcus"))
fall.avg.genus138

sum(mean_by_Season_winter_sub138$averaged.Prop)
winter.avg.genus138 <- ggplot(mean_by_Season_winter_sub138, aes(x = "", y = averaged.Prop, fill = Genus))+
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start = 0, direction = 1) + 
  theme_minimal() + 
  labs(title = "Winter") + 
  theme(axis.title =element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12), legend.key.size = unit(3.5, "mm"), text = element_text(family = "Times New Roman"), legend.text.align = 0, plot.title = element_text(size = 14, hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 5)) + 
  scale_fill_manual(guide = guide_legend(ncol = 5), name = "Genus:", values = c("darkorange3", "darkorange1", "goldenrod1", "lightgoldenrod1", "darkseagreen4", "darkseagreen1", "deeppink4", "deeppink1", "royalblue4", "royalblue1", "plum3", "plum1", "black", "seashell4", "seashell1", "turquoise4", "turquoise"), breaks = c("Ascidiaceihabitans", "Balneola", "Blastopirellula", "Candidatus Actinomarina", "Candidatus Aquiluna", "Clade Ia", "Clade Ib", "Cyanobium PCC-6307", "Formosa", "HIMB11", "MB11C04 marine group", "NS2b marine group", "NS4 marine group", "NS5 marine group", "OM43 clade", "OM60(NOR5) clade", "Synechococcus CC9902"), labels = c("Ascidiaceihabitans", "Balneola", "Blastopirellula", "Candidatus Actinomarina", "Candidatus Aquiluna", "SAR11 Clade Ia", "SAR11 Clade Ib", "Cyanobium PCC-6307", "Formosa", "HIMB11", "MB11C04 marine group", "NS2b marine group", "NS4 marine group", "NS5 marine group", "OM43 clade", "OM60(NOR5) clade", "Synechococcus"))
winter.avg.genus138

#Figure 3b
season.avg.genus138 <- ggarrange(fall.avg.genus138, winter.avg.genus138, spring.avg.genus138, summer.avg.genus138, nrow = 1, ncol = 4, common.legend = T, legend = "bottom")
season.avg.genus138

#Combined Figure 3
set_null_device("png")

png("Frontiers-Figure3-138.png", units="in", width=10, height=10, res=600)

ggarrange(SS_Y1_taxaplot138.2_L6_top25, season.avg.genus138, nrow = 2, ncol = 1, labels = c("(A)", "(B)"), font.label = list(size = 16, color = "black", face = "bold", family = "Times New Roman"), heights = c(1, 0.5))

dev.off()

####Normalize Data - Variance Stabilizing Transformation####
#Construct DESeqDataSet
SS_Y1.138.dds <- DESeqDataSetFromMatrix(countData = t(SS_Y1_Proks_asv_table2.138), colData = SS_Y1_metadata, design = ~ 1)
SS_Y1.138.dds

#Pre-filter the dataset: remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(t(SS_Y1_Proks_asv_table2.138))
sum(rowSums(t(SS_Y1_Proks_asv_table2.138))<2)

SS_Y1_filt138 <- SS_Y1.138.dds[rowSums(counts(SS_Y1.138.dds)) >1,]
nrow(SS_Y1_filt138)

#Normalize for sampling depth through variance stabilizing transformation
SS_Y1.vst138 <- varianceStabilizingTransformation(SS_Y1_filt138)
SS_Y1.vst138$sizeFactor

#Plot the effect of the transformation
SS_Y1.trans.plot138 <- as_tibble(assay(SS_Y1.vst138)[, 1:2]) %>% mutate(transformation = "vst")
colnames(SS_Y1.trans.plot138)[1:2] <- c("x","y")
transformations.plot138 <- ggplot(SS_Y1.trans.plot138, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation) + theme_classic()
transformations.plot138

boxplot(SS_Y1.trans.plot138$x, SS_Y1.trans.plot138$y)

#Pull out transformed table
SS_Y1.trans138 <- assay(SS_Y1.vst138)

####Normalize Data - Rarefaction####
#Plot histogram
hist(x = rowSums(SS_Y1_Proks_asv_table2.138), xlab = "Number of Sequences", main = "Histogram of Sequence Counts")
range(rowSums(SS_Y1_Proks_asv_table2.138))
median(rowSums(SS_Y1_Proks_asv_table2.138))
sort(rowSums(SS_Y1_Proks_asv_table2.138))

#Rarefy
SS_Y1_rarefied138 <- rrarefy(SS_Y1_Proks_asv_table2.138, 20000)

#Remove 0s and singletons
length(which(colSums(SS_Y1_rarefied138) < 2))
SS_Y1_rarefied2.138 <- SS_Y1_rarefied138[, (colSums(SS_Y1_rarefied138) > 1)]
dim(SS_Y1_rarefied2.138)

####PCA####
#Environmental data PCA
compiled.data.pca <- SS_Y1_metadata
compiled.data.pca2 <- compiled.data.pca[, c(13:17)] #pull out environmental variables of interest (temperature, salinity, TKN, NO3+NO2, and Orthophosphate)
SS_Y1_pca_env <- rda(compiled.data.pca2, scale = TRUE)
SS_Y1_pca_env

loadings <- scores(SS_Y1_pca_env, display = 'species', scaling = 0)
loadings
sort(abs(loadings[,1]), decreasing = TRUE)
sort(abs(loadings[,2]), decreasing = TRUE)

biplot(SS_Y1_pca_env, display = 'species', scaling = 'species')
ev <- SS_Y1_pca_env$CA$eig
evplot(ev)

#Convert the environmental loadings to vectors
arrowmat <- vegan::scores(SS_Y1_pca_env, display = "species")
arrowmat
rownames(arrowmat)

#Add labels, make a data.frame
arrowdf <- data.frame(labels = c("Salinity", "Temperature", "Nitrate+Nitrite", "Orthophosphate", "TKN"), arrowmat)
arrowdf

#Define the arrow aesthetic mapping
arrow_map <- aes(xend = 10*PC1, yend = 10*PC2, x = 0, y = 0, shape = NULL, color = NULL, label = labels)
arrow_map
label_map <- aes(x = 13*PC1, y = 13*PC2, shape = NULL, color = NULL, label = arrowdf$labels)
label_map
arrowhead <- arrow(length = unit(0.03, "npc"))

#Microbial Community PCA
SS_Y1.phy.pca138 <- otu_table(SS_Y1.trans138, taxa_are_rows = T)
SS_des_phy.pca138 <- sample_data(SS_Y1_metadata)
SS_Y1.physeq.pca138 <- phyloseq(SS_Y1.phy.pca138, SS_des_phy.pca138)

SS_Y1_pca138 <- ordinate(SS_Y1.physeq.pca138, method = "RDA")
SS_Y1_eigen.vals138.pca <- SS_Y1_pca138$CA$eig

SS_Y1_pca138.plot <- plot_ordination(SS_Y1.physeq.pca138, SS_Y1_pca138) + 
  geom_point(mapping = aes(color = Season, shape = Season), size = 4.5) + 
  coord_fixed(sqrt(SS_Y1_eigen.vals138.pca[2]/SS_Y1_eigen.vals138.pca[1])) + 
  #  labs(x = "Axis 1 [16.9%]", y = "Axis 2 [9.5%]") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 18, family = "Times New Roman")) +
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) +
  expand_limits(y = c(-20, 20), x = c(-20, 20)) + 
  geom_segment(arrow_map, size = 0.5, data = arrowdf, color = "grey40", arrow = arrowhead) +
  geom_text(label_map, size = 4, data = arrowdf, family = "Times New Roman")
SS_Y1_pca138.plot

SS_Y1_pca138.plot + stat_ellipse(mapping = aes(color = Season), type = "t", linetype = 1, level = 0.95)

#Figure 2b
SS_Y1_pca138.plotellipse <- SS_Y1_pca138.plot + stat_ellipse(mapping = aes(color = Season), type = "t", linetype = 1, level = 0.95)

#Extract PC1 and PC2
SS_Y1_cord.vals138.pca <- as.data.frame(SS_Y1_pca138$CA$u[,1:2])
SS_Y1_cord.vals138.pca$SampleID <- row.names(SS_Y1_cord.vals138.pca)
PCA.cor.df <- merge(SS_Y1_cord.vals138.pca, compiled.data, by = "SampleID")

#Plot
plot(y = PCA.cor.df$Temperature, x = PCA.cor.df$PC1, ylab = "Temperature", xlab = "PC1", pch = 16, cex = 1.5)
plot(y = PCA.cor.df$Temperature, x = PCA.cor.df$PC2, ylab = "Temperature", xlab = "PC2", pch = 16, cex = 1.5)

plot(y = PCA.cor.df$Total.Kjeldahl.Nitrogen, x = PCA.cor.df$PC1, ylab = "TKN", xlab = "PC1", pch = 16, cex = 1.5)
plot(y = PCA.cor.df$Total.Kjeldahl.Nitrogen, x = PCA.cor.df$PC2, ylab = "TKN", xlab = "PC2", pch = 16, cex = 1.5)

#Correlation tests
cor.test(y = PCA.cor.df$Temperature, x = PCA.cor.df$PC1, method = "pearson")
#p-value = 7.784e-06; cor 0.7564634
cor.test(y = PCA.cor.df$Temperature, x = PCA.cor.df$PC2, method = "pearson")
#p-value = 0.4176
cor.test(y = PCA.cor.df$Total.Kjeldahl.Nitrogen, x = PCA.cor.df$PC1, method = "pearson")
#p-value = 0.02524; cor -0.4379568
cor.test(y = PCA.cor.df$Total.Kjeldahl.Nitrogen, x = PCA.cor.df$PC2, method = "pearson")
#p-value = 0.07534

###
evplot <- function(ev)
{
  # Broken stick model (MacArthur 1957)
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}
#

####Permutational ANOVA####
#Calculate Euclidean distance matrix
SS_Y1.dist138 <- dist(t(SS_Y1.trans138))

#Check betadispersion
anova(betadisper(SS_Y1.dist138, SS_Y1_metadata$Season))
#Result: p-value = 0.2412
plot(betadisper(SS_Y1.dist138, SS_Y1_metadata$Season))

#Perform PERMANOVA
SS_Y1.adonis138 <- adonis(SS_Y1.dist138~SS_Y1_metadata$Season)
SS_Y1.adonis138$aov.tab
#Result: p-value = 0.001

#Perform post-hoc test
SS_Y1.adonis.ph138 <- pairwise.perm.manova(SS_Y1.dist138, SS_Y1_metadata$Season, p.method = "hochberg")
SS_Y1.adonis.ph138$p.value
#          Fall Spring Summer
#Spring   0.009     NA     NA
#Summer   0.012  0.009     NA
#Winter   0.009  0.034  0.006

####Alpha Diversity Plots, ANOVAs, and Correlation Analysis####
#Calculate Shannon Index
SS_Y1_shannon138 <- vegan::diversity(SS_Y1_rarefied2.138, index = "shannon")
SS_Y1_shannon138.df <- as.data.frame(SS_Y1_shannon138)

#Calculate richness
SS_Y1_richness <- microbiome::richness(t(SS_Y1_rarefied2.138))
SS_Y1_richness.df <- as.data.frame(SS_Y1_richness)
View(SS_Y1_richness.df)

#Calculate Pielou's evenness 
#Pielou’s evenness J = H′/log(S)
SS_Y1_evenness <- vegan::diversity(SS_Y1_rarefied2.138, index = "shannon")/log(specnumber(SS_Y1_rarefied2.138))
SS_Y1_evenness.df <- as.data.frame(SS_Y1_evenness)
View(SS_Y1_evenness.df)

#Merge with metadata
SS_Y1_shannon138.df$SampleID <- rownames(SS_Y1_shannon138.df)
colnames(SS_Y1_shannon138.df) <- c("Shannon", "SampleID")
SS_Y1_metadata_shannon138 <- merge(SS_Y1_metadata, SS_Y1_shannon138.df, by = "SampleID")
range(SS_Y1_metadata_shannon138$Shannon)

SS_Y1_richness.df$SampleID <- rownames(SS_Y1_richness.df)
SS_Y1_metadata_richness <- merge(SS_Y1_metadata, SS_Y1_richness.df, by = "SampleID")
range(SS_Y1_metadata_richness$observed)

SS_Y1_evenness.df$SampleID <- rownames(SS_Y1_evenness.df)
SS_Y1_metadata_evenness <- merge(SS_Y1_metadata, SS_Y1_evenness.df, by = "SampleID")
range(SS_Y1_metadata_evenness$SS_Y1_evenness)

#Plot
a.div138 <- ggplot(SS_Y1_metadata_shannon138, aes(x = CollectionDate, y = Shannon)) + 
  theme_classic() + 
  labs(y = "Shannon Diversity", x = " ") + 
  geom_point(mapping = aes(shape = Season, color = Season), size = 3) + 
  geom_smooth(method = "loess", se = F, span = 0.2, color = "black", linewidth = .5, linetype = "dashed") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "Times New Roman", size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 18), axis.text.x = element_blank()) + 
  expand_limits(y = c(3.5, 6.25)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) +
  scale_y_continuous(breaks = c(3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5), labels = c("3.5", "4.0", "4.5", "5.0", "5.5", "6.0", "6.5")) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month")
a.div138

asv.rich <- ggplot(SS_Y1_metadata_richness, aes(x = CollectionDate, y = observed)) + 
  theme_classic() + 
  labs(y = "ASV Richness", x = " ") + 
  geom_point(mapping = aes(shape = Season, color = Season), size = 3) + 
  geom_smooth(method = "loess", se = F, span = 0.2, color = "black", linewidth = .5, linetype = "dashed") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "Times New Roman", size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 18), axis.text.x = element_blank()) + 
  expand_limits(y = c(100, 900)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) +
  scale_y_continuous(breaks = seq(150, 900, 150)) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month")
asv.rich

#Figure 2a
asv.evenness <- ggplot(SS_Y1_metadata_evenness, aes(x = CollectionDate, y = SS_Y1_evenness)) + 
  theme_classic() + 
  labs(y = "Pielou's Evenness", x = "Date") + 
  geom_point(mapping = aes(shape = Season, color = Season), size = 3) + 
  geom_smooth(method = "loess", se = F, span = 0.2, color = "black", linewidth = .5, linetype = "dashed") + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "Times New Roman", size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  expand_limits(y = c(0.75, 0.9)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_y_continuous(breaks = seq(0.75, 0.9, 0.05)) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month")
asv.evenness
write.csv(SS_Y1_metadata_evenness_compiled, "ASV evenness.csv")

#Supplemental Figure 1
png("Frontiers-SuppFig1-138.png", units="in", width=10, height=12, res=600)

ggarrange(a.div138.2, asv.rich, asv.evenness, nrow = 3, ncol = 1, align = "hv", labels = c("(A)", "(B)", "(C)"), font.label = list(size = 16, color = "black", face = "bold", family = "Times New Roman"), common.legend = T, legend = "top", label.x = -0.0125)

dev.off()

####ANOVAs
#Test of normality
shapiro.test(SS_Y1_metadata_shannon138$Shannon)
#Result: p-value = 0.2306
shapiro.test(SS_Y1_metadata_richness$observed)
#Result: p-value = 0.01743
shapiro.test(SS_Y1_metadata_evenness$SS_Y1_evenness)
#Result: p-value = 0.405

#Test of HOV
leveneTest(y = SS_Y1_metadata_shannon138$Shannon, group = SS_Y1_metadata_shannon138$Season)
#Result: p-value = 0.2893
leveneTest(y = SS_Y1_metadata_richness$observed, group = SS_Y1_metadata_richness$Season)
#Result: p-value = 0.4093
leveneTest(y = SS_Y1_metadata_evenness$SS_Y1_evenness, group = SS_Y1_metadata_evenness$Season)
#Result: p-value = 0.601

#Build linear model
lm.adiv138 <- lm(SS_Y1_metadata_shannon138$Shannon ~ as.factor(SS_Y1_metadata_shannon138$Season))
lm.rich <- lm(SS_Y1_metadata_richness$observed ~ as.factor(SS_Y1_metadata_richness$Season))
lm.even <- lm(SS_Y1_metadata_evenness$SS_Y1_evenness ~ as.factor(SS_Y1_metadata_evenness$Season))

#Run anova
anova(lm.adiv138)
#Result: p-value = 0.2036

anova(lm.rich)
#Result: p-value = 0.4569

anova(lm.even)
#Result: p-value = 0.03328

#Post-hoc tests
require(graphics)

adiv.aov138 <- aov(SS_Y1_metadata_shannon138$Shannon ~ SS_Y1_metadata_shannon138$Season)
adiv.tukey138 <- TukeyHSD(adiv.aov138, conf.level = 0.95)
adiv.tukey138

#                     diff        lwr       upr     p adj
#Spring-Fall    0.03096989 -0.5495632 0.6115030 0.9988043
#Summer-Fall   -0.17162280 -0.7521559 0.4089103 0.8439588
#Winter-Fall    0.25285592 -0.3276772 0.8333890 0.6275323
#Summer-Spring -0.20259269 -0.7325445 0.3273591 0.7157875
#Winter-Spring  0.22188603 -0.3080658 0.7518378 0.6558442
#Winter-Summer  0.42447872 -0.1054731 0.9544305 0.1478293

rich.aov <- aov(SS_Y1_metadata_richness$observed ~ SS_Y1_metadata_richness$Season)
rich.tukey <- TukeyHSD(rich.aov, conf.level = 0.95)
rich.tukey

#                    diff        lwr      upr     p adj
#Spring-Fall    -3.142857 -238.55307 232.2674 0.9999811
#Summer-Fall   -43.142857 -278.55307 192.2674 0.9561002
#Winter-Fall    80.857143 -154.55307 316.2674 0.7764649
#Summer-Spring -40.000000 -254.89914 174.8991 0.9541618
#Winter-Spring  84.000000 -130.89914 298.8991 0.7018359
#Winter-Summer 124.000000  -90.89914 338.8991 0.3978991

even.aov <- aov(SS_Y1_metadata_evenness$SS_Y1_evenness ~ SS_Y1_metadata_evenness$Season)
even.tukey <- TukeyHSD(even.aov, conf.level = 0.95)
even.tukey

#                     diff          lwr         upr     p adj
#Spring-Fall    0.01002338 -0.026564094 0.046610856 0.8710969
#Summer-Fall   -0.01585581 -0.052443287 0.020731663 0.6312202
#Winter-Fall    0.02167608 -0.014911394 0.058263556 0.3753251
#Summer-Spring -0.02587919 -0.059278835 0.007520449 0.1683298
#Winter-Spring  0.01165270 -0.021746942 0.045052342 0.7682703
#Winter-Summer  0.03753189  0.004132251 0.070931535 0.0238538


#Correlation analysis
#Shannon Index
cor.test(x = SS_Y1_metadata_shannon138$Temperature, y = SS_Y1_metadata_shannon138$Shannon, method = "pearson")
#p-value = 0.0281; cor -0.4305838
cor.test(x = SS_Y1_metadata_shannon138$Total.Kjeldahl.Nitrogen, y = SS_Y1_metadata_shannon138$Shannon, method = "pearson")
#p-value = 0.04755; cor 0.392146


#Richness
cor.test(x = SS_Y1_metadata_richness$Temperature, y = SS_Y1_metadata_richness$observed, method = "pearson")
#p-value = 0.07337; cor -0.3570365
cor.test(x = SS_Y1_metadata_richness$Total.Kjeldahl.Nitrogen, y = SS_Y1_metadata_richness$observed, method = "pearson")
#p-value = 0.03506; cor 0.4149181


#Evenness
cor.test(x = SS_Y1_metadata_evenness$Temperature, y = SS_Y1_metadata_evenness$SS_Y1_evenness, method = "pearson")
#p-value = 0.00969; cor -0.4976156
cor.test(x = SS_Y1_metadata_evenness$Total.Kjeldahl.Nitrogen, y = SS_Y1_metadata_evenness$SS_Y1_evenness, method = "pearson")
#p-value =  0.145; cor 0.2939388 


####Euler Diagram of ASVs by Season####
#Data needs to be count data (cannot be vst table)
head(row.names(SS_Y1_rarefied138))
row.names(SS_Y1_rarefied138)
SS_Y1.phy.rare138 <- otu_table(SS_Y1_rarefied138, taxa_are_rows = F)
SS_des_phy.rare138 <- sample_data(SS_Y1_metadata)
SS_Y1_rare.physeq138 <- phyloseq(SS_Y1.phy.rare138, SS_des_phy.rare138)

#Euler diagram - Figure 2c
ps_euler_rarefied138 <- ps_euler(SS_Y1_rare.physeq138, group = "Season", shape = "circle", weight = FALSE, relative = TRUE, plot = TRUE, labels = list(fontfamily = "serif", cex = 1.5), fills = c("orange", "darkseagreen3", "firebrick3", "lightblue2"), quantities = list(fontfamily = "serif", type = c("counts", "percent"), cex = 1.25))
ps_euler_rarefied138

#Compiled Figure 2
set_null_device("png")

png("Frontiers-Figure2-138.png", units="in", width=12, height=10.5, res=600)

ggarrange(asv.evenness, ggarrange(SS_Y1_pca138.plotellipse, as.ggplot(ps_euler_rarefied138), nrow = 1, ncol = 2, labels = c("(B)", "(C)"), font.label = list(size = 16, color = "black", face = "bold", family = "Times New Roman"), label.x = -0.0125, align = "h"), labels = c("(A)"), nrow = 2, font.label = list(size = 16, color = "black", face = "bold", family = "Times New Roman"), label.x = -0.0125)

dev.off()

####Turnover and nestedness components of temporal change####

#At the ASV Level
#uses count data with site as row and taxa as columns
#Separate rarefied asv table into samples from Fall, Winter, Spring, and Summer 2016
#Take average of ASV counts across samples within each dataframe
#Convert to presence/absence format where 1 = present and 0 = absent
row.names(SS_Y1_rarefied2.138)

SS_Y1_asv_counts138_fall <- as.data.frame(SS_Y1_rarefied2.138[c(22:26), ])
row.names(SS_Y1_asv_counts138_fall)
SS_Y1_asv_counts138_fall.avg <- rbind(SS_Y1_asv_counts138_fall, Average = colMeans(SS_Y1_asv_counts138_fall))
Fall_asv_counts138_avg <- SS_Y1_asv_counts138_fall.avg["Average", , drop = F]
Fall_asv_counts138_presabs <- Fall_asv_counts138_avg
Fall_asv_counts138_presabs[Fall_asv_counts138_presabs > 0] <- 1 
View(Fall_asv_counts138_presabs)

SS_Y1_asv_counts138_winter <- as.data.frame(SS_Y1_rarefied2.138[c(28, 2:7), ])
row.names(SS_Y1_asv_counts138_winter)
SS_Y1_asv_counts138_winter.avg <- rbind(SS_Y1_asv_counts138_winter, Average = colMeans(SS_Y1_asv_counts138_winter))
Winter_asv_counts138_avg <- SS_Y1_asv_counts138_winter.avg["Average", , drop = F]
Winter_asv_counts138_presabs <- Winter_asv_counts138_avg
Winter_asv_counts138_presabs[Winter_asv_counts138_presabs > 0] <- 1 

SS_Y1_asv_counts138_spring <- as.data.frame(SS_Y1_rarefied2.138[c(8:11, 13:15), ])
row.names(SS_Y1_asv_counts138_spring)
SS_Y1_asv_counts138_spring.avg <- rbind(SS_Y1_asv_counts138_spring, Average = colMeans(SS_Y1_asv_counts138_spring))
Spring_asv_counts138_avg <- SS_Y1_asv_counts138_spring.avg["Average", , drop = F]
Spring_asv_counts138_presabs <- Spring_asv_counts138_avg
Spring_asv_counts138_presabs[Spring_asv_counts138_presabs > 0] <- 1 

SS_Y1_asv_counts138_summer <- as.data.frame(SS_Y1_rarefied2.138[c(16:21), ])
row.names(SS_Y1_asv_counts138_summer)
SS_Y1_asv_counts138_summer.avg <- rbind(SS_Y1_asv_counts138_summer, Average = colMeans(SS_Y1_asv_counts138_summer))
Summer_asv_counts138_avg <- SS_Y1_asv_counts138_summer.avg["Average", , drop = F]
Summer_asv_counts138_presabs <- Summer_asv_counts138_avg
Summer_asv_counts138_presabs[Summer_asv_counts138_presabs > 0] <- 1 

#Calculate nestedness and turnover components
fall_winter_asv_beta <- beta.temp(Fall_asv_counts138_presabs, Winter_asv_counts138_presabs, index.family="sor")
fall_winter_asv_beta
#beta.sim = turnover; beta.sne = nestedness; beta.sor = overall dissimilarity 
#         beta.sim  beta.sne beta.sor
#Average 0.5204403 0.1006778 0.621118

winter_spring_asv_beta <- beta.temp(Winter_asv_counts138_presabs, Spring_asv_counts138_presabs, index.family="sor")
winter_spring_asv_beta
#beta.sim = turnover; beta.sne = nestedness; beta.sor = overall dissimilarity 
#         beta.sim   beta.sne  beta.sor
#Average 0.4895397 0.03876735 0.5283071

spring_summer_asv_beta <- beta.temp(Spring_asv_counts138_presabs, Summer_asv_counts138_presabs, index.family="sor")
spring_summer_asv_beta
#beta.sim = turnover; beta.sne = nestedness; beta.sor = overall dissimilarity 
#         beta.sim   beta.sne  beta.sor
#Average 0.5322835 0.06404681 0.5963303

####Mantel Test####
#Bray-Curtis 
SS_Y1_rarefied138.t <- t(SS_Y1_rarefied138)
SS_Y1_rarefied138_bc <- vegdist(t(SS_Y1_rarefied138), method = "bray")

#Plot community distance matrix to look at the distance decay of similarity
plot(envdist, SS_Y1_rarefied138_bc, pch = 16, cex = 0.5, col = "black", xlab = "Spatial Distance", ylab = "Species composition dissimilarity")

#Generate a Euclidian distance matrix for bacterial inhibition data
colnames(SS_Y1_metadata)
bact <- SS_Y1_metadata[, c(1, 3)]
bact$Date <- NULL
bactdist <- vegdist(bact, method = "euclidean")

#Run Mantel test
mantel.inhib <- ecodist::mantel(bactdist ~ SS_Y1_rarefied138_bc, nperm = 10000, nboot = 1000, pboot = 0.95, cboot = 0.95)
mantel.inhib
#Mantel statistic r: 0.2810173
#pval = 0.0016000

####Figure Accessibility####
#Obtaining RGB codes
#Figure 2a and 2b
col2rgb(c("aquamarine4", "firebrick3", "darkorange2", "blue4"), alpha = FALSE)

#Figure 2c
col2rgb(c("orange", "darkseagreen3", "firebrick3", "lightblue2"), alpha = FALSE)

#Figure 3a
col2rgb(c("pink4", "darkorange3", "pink1", "darkorange1", "goldenrod1", "lightgoldenrod1", "darkseagreen4", "darkseagreen1", "deeppink4", "deeppink1", "orchid4", "royalblue4", "royalblue1", "plum3", "violet", "plum1", "black", "seashell4", "seashell1", "turquoise4", "palegreen4", "lightsalmon1", "springgreen3", "turquoise", "tan4", "wheat3"), alpha = FALSE)

#Figure 3b
col2rgb(c("darkorange3", "darkorange1", "goldenrod1", "lightgoldenrod1", "darkseagreen4", "darkseagreen1", "deeppink4", "deeppink1", "royalblue4", "royalblue1", "plum3", "plum1", "black", "seashell4", "seashell1", "turquoise4", "turquoise"), alpha = FALSE)
col2rgb("white", alpha = FALSE)

#Checked compatibility at: https://color.adobe.com/create/color-accessibility
#No conflicts found. Swatches are color blind safe.
