# could be required for ChIPseeker installation:
# sudo chmod a+wr -R /usr/lib/R/library
# sudo apt-get install libxml2-dev
# sudo apt-get install libcurl4-openssl-dev libssl-dev
# install.packages("devtools", dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("ChIPseeker", quietly = TRUE))
  BiocManager::install("ChIPseeker")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

setwd("~/chipseq/workdir/downstream")
peaks <- readPeakFile("../macs2/GSM1102797_CD14_H3K4me3_hg19.chr15_broad0.1_peaks.broadPeak")
peaks

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(
  "../macs2/GSM1102797_CD14_H3K4me3_hg19.chr15_broad0.1_peaks.broadPeak",
  tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
