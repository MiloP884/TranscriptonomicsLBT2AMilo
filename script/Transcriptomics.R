setwd("/Users/milopostma/Documents/Transcriptomics")
unzip("Data_RA_raw.zip", exdir = "Data_raw")

install.packages('BiocManager')
BiocManager::install('Rsubread')
library(Rsubread)

browseVignettes('Rsubread')

unzip("ncbi_dataset-genome-2.zip", exdir = "Data_genome")

#Indexeren van genoom
buildindex(
  basename = 'ref_ecoli',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)

#Eerste data tegen genoom van de mens

align.controle1 <- align(index = "ref_ecoli", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2 = "SRR4785819_2_subset40k.fastq", output_file = "controle1.BAM")

align.controle2 <- align(index = "ref_ecoli", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2 = "SRR4785820_2_subset40k.fastq", output_file = "controle2.BAM")

align.controle3 <- align(index = "ref_ecoli", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2 = "SRR4785828_2_subset40k.fastq", output_file = "controle3.BAM")

align.controle4 <- align(index = "ref_ecoli", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2 = "SRR4785831_2_subset40k.fastq", output_file = "controle4.BAM")

align.RA1 <- align(index = "ref_ecoli", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2 = "SRR4785979_2_subset40k.fastq", output_file = "RA1.BAM")

align.RA2 <- align(index = "ref_ecoli", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2 = "SRR4785980_2_subset40k.fastq", output_file = "RA2.BAM")

align.RA3 <- align(index = "ref_ecoli", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2 = "SRR4785986_2_subset40k.fastq", output_file = "RA3.BAM")

align.RA4 <- align(index = "ref_ecoli", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2 = "SRR4785988_2_subset40k.fastq", output_file = "RA4.BAM")


#BAM bestanden maken

# Laad Rsamtools voor sorteren en indexeren
library(Rsamtools)

# Bestandsnamen van de monsters
samples <- c('controle1', 'controle2', 'controle3', 'controle4', 'RA1', 'RA2', 'RA3', 'RA4')

# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})

#.bam.bai bestanden genereren met Rsamtools

indexBam("controle1.sorted.bam")
indexBam("controle2.sorted.bam")
indexBam("controle3.sorted.bam")
indexBam("controle4.sorted.bam")
indexBam("RA1.sorted.bam")
indexBam("RA2.sorted.bam")
indexBam("RA3.sorted.bam")
indexBam("RA4.sorted.bam")
indexFa("GCF_000001405.40_GRCh38.p14_genomic.fna")

#Werkcollege 2
#Count matrix

library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)

#featurecounts uitvoeren

# Je definieert een vector met namen van BAM-bestanden. Elke BAM bevat reads van een RNA-seq-experiment (bijv. behandeld vs. controle).

allsamples <- c("controle1.BAM", "controle2.BAM", "controle3.BAM", "controle4.BAM", "RA1.BAM", "RA2.BAM", "RA3.BAM", "RA4.BAM")

unzip("ncbi_dataset-human-GTF.zip", exdir = "GTF_genome")

count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

head(count_matrix$annotation)
head(count_matrix$counts)

# Bekijk eerst de structuur van het object
str(count_matrix)

# Haal alleen de matrix met tellingen eruit
counts <- count_matrix$counts

colnames(counts) <- c("controle1", "controle2", "controle3", "controle4", "RA1", "RA2", "RA3", "RA4")

#opslaan van de matrix

write.csv(counts, "bewerkt_countmatrix.csv")

head(counts)

#Dag 3: statistiek en analyse

counts <- read.table("count_matrix.txt", header = TRUE, row.names = 1)


treatment <- c("controle", "controle", "controle", "controle", "RA", "RA", "RA", "RA")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c('controle1', 'controle2', 'controle3', 'controle4', 'RA1', 'RA2', 'RA3', 'RA4')


library(DESeq2)
library(KEGGREST)

# Maak DESeqDataSet aan
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)

# Voer analyse uit
dds <- DESeq(dds)
resultaten <- results(dds)

# Resultaten opslaan in een bestand
#Bij het opslaan van je tabel kan je opnieuw je pad instellen met `setwd()` of het gehele pad waar je de tabel wilt opslaan opgeven in de code.

write.table(resultaten, file = 'ResultatenWC3.csv', row.names = TRUE, col.names = TRUE)

#Welke genen zijn significant verdubbeld in expressie of gehalveerd?

sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

#Meest opvallende genen

hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]

head(laagste_p_waarde)

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')

# Alternatieve plot zonder p-waarde cutoff (alle genen zichtbaar)

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)

dev.copy(png, 'VolcanoplotWC.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

#Downloaden en installeren van de package pathview

if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}
library(pathview)

#Uitvoeren van een GO-analyse
install.packages("dplyr")
install.packages("tidyverse")
BiocManager::install('goseq')
BiocManager::install('geneLenDataBase')
BiocManager::install("org.Dm.eg.db")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("AnnotationDbi")
library("goseq")
library("geneLenDataBase")
library("org.Dm.eg.db")
library('org.Hs.eg.db')
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library('AnnotationDbi')
library(tidyverse)
library(dplyr)

#All bestand maken van de data
res_all <- results(dds, alpha = 0.05)
res_all_df <- as.data.frame(res_all)
write.csv(res_all_df, "ALL_genes.csv")

DEG <- read.table("ResultatenWC3.csv", header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)
ALL <- read.table("ALL_genes.csv", header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)

#Omzetten naar vector
res_df <- as.data.frame(resultaten)
filter_res <- filter(res_df, res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1 )

class(DEG)
DEG.vector <- rownames(filter_res)
ALL.vector <- rownames(resultaten)

gene.vector=as.integer(ALL.vector%in%DEG.vector)
names(gene.vector)=ALL.vector

pwf=nullp(gene.vector, 'hg19', 'geneSymbol')
GO.wall=goseq(pwf,"hg19", 'geneSymbol')

goResults <- goseq(pwf, "hg19", 'geneSymbol', test.cats = c("GO:BP"))

#Visualiseren van de data

library(pathview)
goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

library(GO.db)
GOTERM[[goResults$category[1]]]

#KEGG-pathway visualiseren

resultaten[1] <- NULL
resultaten[2:5] <- NULL
gene_vector <- resultaten$log2FoldChange
names(gene_vector) <- rownames(resultaten)
pathview(
  gene.data = gene_vector,
  pathway.id = "hsa04660",  # KEGG ID voor T-cell receptor signaling – Hompo Sapiens
  species = "hsa",          # 'hsa' = Humaan in KEGG
  gene.idtype = "SYMBOL",     # Geef aan dat het KEGG-ID's zijn
  limit = list(gene = 5)    # Kleurbereik voor log2FC van -5 tot +5
)

BiocManager::install('clusterProfiler')
BiocManager::install('enrichplot')
BiocManager::install('pheatmap')
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(pheatmap)

