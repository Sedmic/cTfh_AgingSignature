library("genefilter")
library("ggplot2")
library("grid")
library("ggrepel")
library("RColorBrewer")
library("DESeq2")
library("BiocParallel")
library("WGCNA")
library("Rtsne")
library("pheatmap")
library("tools")
library("viridis")
library("ggpubr")
register(SnowParam(workers=7))
sessionInfo()

bestDataLog <- read.csv("bestDataLog.csv"); rownames(bestDataLog) <- bestDataLog$X; bestDataLog$X <- NULL

## ***************************     GSEA leading edge analysis    **************


leadingEdge <- list()
files <- list.files("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_Hallmark.GseaPreranked.1539515189801/", pattern='.xls$', full=TRUE)
files <- files[grep("HALLMARK",files)]
for(i in 1:length(files)){  leadingEdge[[i]] <- read.csv(file = files[i],sep="\t")  }
saveNames <- list()
for(i in 1:length(leadingEdge)) {  saveNames[[i]] <- leadingEdge[[i]]$PROBE }

# now i need to merge the elements of the list into a large matrix of counts
geneOccurrences <- as.data.frame(table(unlist(saveNames)))

Hallmark_HivsLo_v1_Pos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_Hallmark.GseaPreranked.1539515189801/gsea_report_for_na_pos_1539515189801.xls", sep="\t")
Hallmark_HivsLo_v1_signif <- Hallmark_HivsLo_v1_Pos[which(Hallmark_HivsLo_v1_Pos$FDR.q.val < 0.05),]

files

