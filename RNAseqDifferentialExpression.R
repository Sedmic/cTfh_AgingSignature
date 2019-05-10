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
library("scales")
library("gridExtra")
register(SnowParam(workers=7))
sessionInfo()

dataMatrixCounts <- read.table("RawProcessing/2018July_PORTnormalization/SPREADSHEETS/FINAL_master_list_of_gene_counts_MIN.FluVacAging.txt", sep="", header=T, stringsAsFactors = F)
rownames(dataMatrixCounts) <- make.names(dataMatrixCounts$geneSymbol, unique=TRUE)
dataMatrixCounts$id <- dataMatrixCounts$geneCoordinate <- dataMatrixCounts$geneSymbol <- NULL 

dataMatrixCounts$X222006_HiHi_v1 <- NULL  # based on what was learned from the RNAseq_QC script

bestDataCounts <- dataMatrixCounts


## ***************************     Log transformation  **************************************

metaData <- data.frame(row.names=colnames(bestDataCounts));  metaData$condition <- "empty"
#  1:36 are young, then 37:83 are elderly

metaData$ageGroup <- c(rep("Y",36),rep("E",47))
metaData$subject <- substr(rownames(metaData), start=2,stop=7)
metaData$condition <- substr(rownames(metaData),start=9, stop=20)
metaData$subgroup <- paste0(metaData$ageGroup,sep="_", metaData$condition)

fullDataset <- DESeqDataSetFromMatrix(countData=bestDataCounts, colData=metaData,design= ~ subgroup)

# vsd <- varianceStabilizingTransformation(fullDataset)
# bestDataLog <- as.data.frame(assay(vsd))
# write.csv(bestDataLog, file="bestDataLog.csv")
bestDataLog <- read.csv("bestDataLog.csv"); rownames(bestDataLog) <- bestDataLog$X; bestDataLog$X <- NULL



## ***************************     tSNE analysis  **************************************

set.seed(42)
tsneMap <- Rtsne(t(bestDataLog),epoch=50,perplexity=6,k=2,theta=0,verbosity=TRUE,max_iter=1000)
tsneMap <- as.data.frame(tsneMap$Y); tsneMap$subgroup <- metaData$subgroup; rownames(tsneMap) <- colnames(bestDataLog)
tsneReorder <- rbind(
  tsneMap[grep("Y_HiHi_v1",tsneMap$subgroup),], tsneMap[grep("Y_HiHi_v2",tsneMap$subgroup),], 
  tsneMap[grep("Y_LoLo_v1",tsneMap$subgroup),], tsneMap[grep("Y_LoLo_v2",tsneMap$subgroup),], 
  tsneMap[grep("Y_Naive_v1",tsneMap$subgroup),], tsneMap[grep("Y_Naive_v2",tsneMap$subgroup),],
  tsneMap[grep("E_HiHi_v1",tsneMap$subgroup),], tsneMap[grep("E_HiHi_v2",tsneMap$subgroup),], 
  tsneMap[grep("E_LoLo_v1",tsneMap$subgroup),], tsneMap[grep("E_LoLo_v2",tsneMap$subgroup),], 
  tsneMap[grep("E_Naive_v1",tsneMap$subgroup),], tsneMap[grep("E_Naive_v2",tsneMap$subgroup),])
#tsneReorder <- rbind(tsneMap[48:53,],tsneMap[66:71,],tsneMap[1:7,], tsneMap[24:31,],tsneMap[54:59,],tsneMap[72:77,],tsneMap[8:15,],tsneMap[32:39,],tsneMap[60:65,],tsneMap[78:83,],tsneMap[16:23,],tsneMap[40:47,])
tsneReorder$Class <- c(rep(1,length(grep("Y_HiHi_v1",tsneReorder$subgroup))),rep(2,length(grep("Y_HiHi_v2",tsneReorder$subgroup))),
                       rep(3,length(grep("Y_LoLo_v1",tsneReorder$subgroup))),rep(4,length(grep("Y_LoLo_v2",tsneReorder$subgroup))),
                       rep(5,length(grep("Y_Naive_v1",tsneReorder$subgroup))),rep(6,length(grep("Y_Naive_v2",tsneReorder$subgroup))),
                       rep(7,length(grep("E_HiHi_v1",tsneReorder$subgroup))),rep(8,length(grep("E_HiHi_v2",tsneReorder$subgroup))),
                       rep(9,length(grep("E_LoLo_v1",tsneReorder$subgroup))),rep(10,length(grep("E_LoLo_v2",tsneReorder$subgroup))),
                       rep(11,length(grep("E_Naive_v1",tsneReorder$subgroup))),rep(12,length(grep("E_Naive_v2",tsneReorder$subgroup))))
customPalette <- colorRampPalette(brewer.pal(12,"Paired"))(12)
customPalette[11] <- "#DE966D"  # modify paired palette to eliminate yellow
temp <- customPalette[c(7,8,3,4,1,1,9,10,1,1,11,12)]
temp[5:6] <- c("#fff849","#d6cf44"); temp[9:10] <- c("#03681e","#003d10"); temp[11:12] <- c("#99930a","#686408")
customPalette <- temp

ggplot(tsneReorder, aes(x=V1, y=V2)) + 
  geom_point(size=10,pch=21,colour="black",aes(fill=factor(Class,labels=c("Young ICOS+CD38+ Day 0","Young ICOS+CD38+ Day 7",
                                                                          "Young ICOS-CD38- Day 0","Young ICOS-CD38- Day 7",
                                                                          "Young Naive Day 0","Young Naive Day 7",
                                                                          "Elderly ICOS+CD38+ Day 0","Elderly ICOS+CD38+ Day 7",
                                                                          "Elderly ICOS-CD38- Day 0","Elderly ICOS-CD38- Day 7",
                                                                          "Elderly Naive Day 0","Elderly Naive Day 7")))) +  
  xlab("") + ylab("") + 
  ggtitle("t-SNE All samples") +  theme_light(base_size=15)  +  scale_fill_manual(values=customPalette[unique(tsneReorder$Class)]) +
  theme(strip.background = element_blank()) + labs(fill = "Sample cohort") + theme(legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) #+ xlim(-60,60) + ylim(-55,70)
# ggsave(filename="Images/tSNE_all_samples.pdf", device="pdf",width=8,height=6)


subsetTSNE <- tsneReorder[c(grep("Y_HiHi",tsneReorder$subgroup),grep("E_HiHi",tsneReorder$subgroup)),]  ## only view ICOS+CD38+ 
ggplot(subsetTSNE, aes(x=V1, y=V2)) + 
  geom_point(size=12,pch=21,colour="black",aes(fill=factor(Class, labels=c("Young ICOS+CD38+ Day 0","Young ICOS+CD38+ Day 7",
                                                                           "Elderly ICOS+CD38+ Day 0","Elderly ICOS+CD38+ Day 7")))) + 
  xlab("") + ylab("") + # geom_text(label=subsetTSNE$subgroup) +
  ggtitle("t-SNE ICOS+CD38+") +  theme_light(base_size=20)  +  scale_fill_manual(values=customPalette[unique(subsetTSNE$Class)]) +
  theme(strip.background = element_blank()) + labs(fill = "Sample cohort") + theme(legend.key = element_blank()) + 
  guides(colour=guide_legend(override.aes=list(size=8))) + xlim(-50,20) + ylim(-35,5)
# ggsave(filename="Images/tSNE_ICOShiCD38hi.pdf", device="pdf",width=8,height=6)


pca1 = prcomp(t(bestDataLog[,grep("HiHi", colnames(bestDataLog), value=F)]), scale = FALSE)
scores = as.data.frame(pca1$x) 
scores$color <- c("#FDBF6F","#FF7F00","#FDBF6F","#FF7F00","#FDBF6F","#FF7F00","#FDBF6F","#FF7F00","#FDBF6F","#FF7F00","#FDBF6F","#FF7F00",
                  "#6A3D9A","#CAB2D6", "#6A3D9A","#CAB2D6", "#6A3D9A","#CAB2D6", "#6A3D9A","#CAB2D6", "#6A3D9A","#CAB2D6", "#6A3D9A","#CAB2D6", "#6A3D9A","#CAB2D6", "#6A3D9A")

# plot of observations
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
  geom_hline(yintercept = 0, colour = "gray65") + theme_bw() +  geom_vline(xintercept = 0, colour = "gray65") +  geom_point(colour=scores$color, alpha=1, size=6) + 
#  geom_text_repel(colour = "black", alpha = 1, size = 3) +
  ggtitle("PCA plot for all samples in final pool") + theme(title = element_text(size=20))
# ggsave(filename = "Images/PCAplot_allSamples.pdf", device="pdf")




metaDataHiHi <- metaData[grep("HiHi", colnames(bestDataLog)),]
tsneMapHiHi <- Rtsne(t(bestDataLog[,grep("HiHi",colnames(bestDataLog))]),epoch=50,perplexity=6,k=2,theta=0,verbosity=TRUE,max_iter=1000)
tsneMapHiHi <- as.data.frame(tsneMapHiHi$Y); tsneMapHiHi$subgroup <- metaDataHiHi$subgroup; rownames(tsneMapHiHi) <- rownames(metaDataHiHi)
tsneReorderHiHi <- rbind(
  tsneMapHiHi[grep("Y_HiHi_v1",tsneMapHiHi$subgroup),], tsneMapHiHi[grep("Y_HiHi_v2",tsneMapHiHi$subgroup),], 
  tsneMapHiHi[grep("E_HiHi_v1",tsneMapHiHi$subgroup),], tsneMapHiHi[grep("E_HiHi_v2",tsneMapHiHi$subgroup),])
tsneReorderHiHi$Class <- c(rep(1,length(grep("Y_HiHi_v1",tsneReorderHiHi$subgroup))),rep(2,length(grep("Y_HiHi_v2",tsneReorderHiHi$subgroup))),
                       rep(3,length(grep("Y_LoLo_v1",tsneReorderHiHi$subgroup))),rep(4,length(grep("Y_LoLo_v2",tsneReorderHiHi$subgroup))),
                       rep(5,length(grep("Y_Naive_v1",tsneReorderHiHi$subgroup))),rep(6,length(grep("Y_Naive_v2",tsneReorderHiHi$subgroup))),
                       rep(7,length(grep("E_HiHi_v1",tsneReorderHiHi$subgroup))),rep(8,length(grep("E_HiHi_v2",tsneReorderHiHi$subgroup))),
                       rep(9,length(grep("E_LoLo_v1",tsneReorderHiHi$subgroup))),rep(10,length(grep("E_LoLo_v2",tsneReorderHiHi$subgroup))),
                       rep(11,length(grep("E_Naive_v1",tsneReorderHiHi$subgroup))),rep(12,length(grep("E_Naive_v2",tsneReorder$subgroup))))
customPalette <- colorRampPalette(brewer.pal(12,"Paired"))(12)
customPalette[11] <- "#DE966D"  # modify paired palette to eliminate yellow
temp <- customPalette[c(7,8,3,4,1,1,9,10,1,1,11,12)]
temp[5:6] <- c("#fff849","#d6cf44"); temp[9:10] <- c("#03681e","#003d10"); temp[11:12] <- c("#99930a","#686408")
customPalette <- temp

ggplot(tsneReorderHiHi, aes(x=V1, y=V2)) + 
  geom_point(size=10,pch=21,colour="black",aes(fill=factor(Class,labels=c("Young ICOS+CD38+ Day 0","Young ICOS+CD38+ Day 7",
                                                                          "Elderly ICOS+CD38+ Day 0","Elderly ICOS+CD38+ Day 7")))) +  
  xlab("") + ylab("") + 
  ggtitle("t-SNE ICOS+CD38+ cTfh") +  theme_light(base_size=15)  +  scale_fill_manual(values=customPalette[unique(tsneReorderHiHi$Class)]) +
  theme(strip.background = element_blank()) + labs(fill = "Sample cohort") + theme(legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) #+ xlim(-60,60) + ylim(-55,70)
# ggsave(filename="Images/tSNE_ICOShiCD38hi.pdf", device="pdf",width=8,height=6)




## ***************************     heatmap of DiffExp of HiHi vs LoLo for selected genes AllAges **************************************

probeList <- c("CXCR5", "PRDM1", "BCL6", "IL10",  "IFNG","MAF","CCR6","CXCR3","GATA3",
               "BTLA","TNFRSF4", "CD38","TIGIT","SLAMF1","POU2AF1","MKI67","SH2D1A", "TOX2", "BIRC5", "TBK1")
probeGenes <- bestDataLog[probeList,grep("v1",colnames(bestDataLog))]
Hi_v_Lo_allAges <- cbind(
  probeGenes[,grep("HiHi_v1",colnames(probeGenes))], probeGenes[,grep("HiHi_v2",colnames(probeGenes))], 
  probeGenes[,grep("LoLo_v1",colnames(probeGenes))], probeGenes[,grep("LoLo_v2",colnames(probeGenes))], 
  probeGenes[,grep("Naive_v1",colnames(probeGenes))], probeGenes[,grep("Naive_v2",colnames(probeGenes))])

annotateHeatmap <- data.frame(row.names = colnames(Hi_v_Lo_allAges), subset = c(rep("ICOS+CD38+ cTfh", 13),rep("ICOS-CD38- cTfh", 14),rep("Naive CD4", 14)), 
                              ageGroup = c(rep("Young",6), rep("Elderly", 7),rep("Young",6), rep("Elderly", 8),rep("Young",6), rep("Elderly", 8)))
ann_colors = list(  subset = c("ICOS+CD38+ cTfh" ="#0D0887", "ICOS-CD38- cTfh" = "#E16462", "Naive CD4" ="#F0F921"), ageGroup = c("Young"="orange3", "Elderly" = "purple")  )
pheatmap(Hi_v_Lo_allAges, scale="row", cluster_col=F, annotation_col = annotateHeatmap, show_colnames=F, main="Selected Tfh genes",
         gaps_col = c(13,27), annotation_colors = ann_colors, fontsize_row = 18, color=inferno(100), cellheight=30, cutree_rows=3, border_color = F, 
         #, filename = "Images/SelectedGenesHeatmapAllAges.pdf"
)



# dev.off(); dev.off(); 

## ***************************     heatmap of DiffExp of HiHi vs LoLo for selected genes  YOUNG **************************************

probeList <- c("CXCR5", "PRDM1", "BCL6", "IL10",  "IFNG","MAF","CCR6","CXCR3","GATA3",
               "BTLA","TNFRSF4", "CD38","TIGIT","SLAMF1","POU2AF1","LEF1","MKI67","SH2D1A", "TOX2", "BIRC5", "TBK1", "IL32")
probeGenes <- bestDataLog[probeList,grep("v1",colnames(bestDataLog[,1:36]))]
Hi_v_Lo_y2 <- cbind(
  probeGenes[,grep("HiHi_v1",colnames(probeGenes))], probeGenes[,grep("HiHi_v2",colnames(probeGenes))], 
  probeGenes[,grep("LoLo_v1",colnames(probeGenes))], probeGenes[,grep("LoLo_v2",colnames(probeGenes))], 
  probeGenes[,grep("Naive_v1",colnames(probeGenes))], probeGenes[,grep("Naive_v2",colnames(probeGenes))])

annotateHeatmap <- data.frame(row.names = colnames(Hi_v_Lo_y2), subset = c(rep("ICOS+CD38+ cTfh", 6),rep("ICOS-CD38- cTfh", 6),rep("Naive CD4", 6)))
ann_colors = list(  subset = c("ICOS+CD38+ cTfh" ="orange", "ICOS-CD38- cTfh" = "green", "Naive CD4" ="yellow")  )
pheatmap(Hi_v_Lo_y2, scale="row", cluster_col=F, annotation_col = annotateHeatmap, show_colnames=F, main="Selected Tfh genes",
         gaps_col = c(6,12), annotation_colors = ann_colors, cutree_rows = 3, fontsize_row = 16, color=inferno(100)
         #, filename = "Images/SelectedGenesHeatmap.pdf"
)

# dev.off(); dev.off(); 

## ***************************     heatmap of DiffExp of HiHi vs LoLo for selected genes  ELDERLY **************************************

probeList <- c("CXCR5", "PRDM1", "BCL6", "IL10",  "IFNG","MAF","CCR6","CXCR3","GATA3",
               "BTLA","TNFRSF4", "CD38","TIGIT","SLAMF1","POU2AF1","LEF1","MKI67","SH2D1A", "TOX2", "BIRC5", "TBK1", "IL32")
probeGenes <- bestDataLog[probeList,grep("v1",colnames(bestDataLog))]
probeGenes <- probeGenes[,grep("X222", colnames(probeGenes))]
Hi_v_Lo_y2 <- cbind(
  probeGenes[,grep("HiHi_v1",colnames(probeGenes))], probeGenes[,grep("HiHi_v2",colnames(probeGenes))], 
  probeGenes[,grep("LoLo_v1",colnames(probeGenes))], probeGenes[,grep("LoLo_v2",colnames(probeGenes))], 
  probeGenes[,grep("Naive_v1",colnames(probeGenes))], probeGenes[,grep("Naive_v2",colnames(probeGenes))])

annotateHeatmap <- data.frame(row.names = colnames(Hi_v_Lo_y2), subset = c(rep("ICOS+CD38+ cTfh", 7),rep("ICOS-CD38- cTfh", 8),rep("Naive CD4", 8)))
ann_colors = list(  subset = c("ICOS+CD38+ cTfh" ="purple", "ICOS-CD38- cTfh" = "darkgreen", "Naive CD4" ="gold")  )
pheatmap(Hi_v_Lo_y2, scale="row", cluster_col=F, annotation_col = annotateHeatmap, show_colnames=F, main="Selected Tfh genes - ELDERLY",
         gaps_col = c(7,15), annotation_colors = ann_colors, cutree_rows = 3, fontsize_row = 16, color=inferno(100)
        #, filename = "Images/SelectedGenesHeatmap_ELDERLY.pdf"
)

# dev.off(); dev.off(); 


## ******************************** Differential expression by subgroup and age category *********************************************
metaData
fullDataset <- DESeqDataSetFromMatrix(countData=bestDataCounts, colData=metaData, design= ~ subgroup)
dds <- estimateSizeFactors(fullDataset)
temp <- counts(dds,normalized=TRUE)
idx <- rowSums( counts(dds, normalized=TRUE) >= 20 ) >= 20  # filter for genes with at least 20 counts in 25% of remaining samples
# write.csv(bestDataCounts[idx,],file="FilteredBestDataCounts.csv")
# write.csv(bestDataLog[idx,],file="FilteredBestDataLog.csv")

fullDataset <- fullDataset[idx,]

DESdata_fullDataset <- DESeq(fullDataset, parallel=TRUE)

diffExpr <- list( Y_HiHivsLoLo_v1 = results(DESdata_fullDataset, contrast=c("subgroup","Y_HiHi_v1","Y_LoLo_v1"), parallel = T),
                  Y_HiHivsLoLo_v2 = results(DESdata_fullDataset, contrast=c("subgroup","Y_HiHi_v2","Y_LoLo_v2"), parallel = T), 
                  Y_HiHivsnaive_v1 = results(DESdata_fullDataset, contrast=c("subgroup","Y_HiHi_v1","Y_Naive_v1"), parallel = T), 
                  Y_HiHivsnaive_v2 = results(DESdata_fullDataset, contrast=c("subgroup","Y_HiHi_v2","Y_Naive_v2"), parallel = T), 
                  Y_LoLovsnaive_v1 = results(DESdata_fullDataset, contrast=c("subgroup","Y_LoLo_v1","Y_Naive_v1"), parallel = T), 
                  Y_LoLovsnaive_v2 = results(DESdata_fullDataset, contrast=c("subgroup","Y_LoLo_v2","Y_Naive_v2"), parallel = T), 
                  E_HiHivsLoLo_v1 = results(DESdata_fullDataset, contrast=c("subgroup","E_HiHi_v1","E_LoLo_v1"), parallel = T), 
                  E_HiHivsLoLo_v2 = results(DESdata_fullDataset, contrast=c("subgroup","E_HiHi_v2","E_LoLo_v2"), parallel = T), 
                  E_HiHivsnaive_v1 = results(DESdata_fullDataset, contrast=c("subgroup","E_HiHi_v1","E_Naive_v1"), parallel = T), 
                  E_HiHivsnaive_v2 = results(DESdata_fullDataset, contrast=c("subgroup","E_HiHi_v2","E_Naive_v2"), parallel = T), 
                  E_LoLovsnaive_v1 = results(DESdata_fullDataset, contrast=c("subgroup","E_LoLo_v1","E_Naive_v1"), parallel = T), 
                  E_LoLovsnaive_v2 = results(DESdata_fullDataset, contrast=c("subgroup","E_LoLo_v2","E_Naive_v2"), parallel = T), 
                  
                  Y_HiHi_v12 = results(DESdata_fullDataset, contrast=c("subgroup","Y_HiHi_v2","Y_HiHi_v1"), parallel = T), 
                  Y_LoLo_v12 = results(DESdata_fullDataset, contrast=c("subgroup","Y_LoLo_v2","Y_LoLo_v1"), parallel = T), 
                  Y_naivevnaive_v12 = results(DESdata_fullDataset, contrast=c("subgroup","Y_Naive_v2","Y_Naive_v1"), parallel = T), 
                  E_HiHi_v12 = results(DESdata_fullDataset, contrast=c("subgroup","E_HiHi_v2","E_HiHi_v1"), parallel = T), 
                  E_LoLo_v12 = results(DESdata_fullDataset, contrast=c("subgroup","E_LoLo_v2","E_LoLo_v1"), parallel = T), 
                  E_naivevnaive_v12 = results(DESdata_fullDataset, contrast=c("subgroup","E_Naive_v2","E_Naive_v1"), parallel = T), 
                  
                  YvE_HiHi_v1 = results(DESdata_fullDataset, contrast=c("subgroup","E_HiHi_v1","Y_HiHi_v1"), parallel = T), 
                  YvE_HiHi_v2 = results(DESdata_fullDataset, contrast=c("subgroup","E_HiHi_v2","Y_HiHi_v2"), parallel = T), 
                  YvE_LoLo_v1 = results(DESdata_fullDataset, contrast=c("subgroup","E_LoLo_v1","Y_LoLo_v1"), parallel = T), 
                  YvE_LoLo_v2 = results(DESdata_fullDataset, contrast=c("subgroup","E_LoLo_v2","Y_LoLo_v2"), parallel = T), 
                  YvE_naive_v1 = results(DESdata_fullDataset, contrast=c("subgroup","E_Naive_v1","Y_Naive_v1"), parallel = T), 
                  YvE_naive_v2 = results(DESdata_fullDataset, contrast=c("subgroup","E_Naive_v2","Y_Naive_v2"), parallel = T)
)

# for (i in 1:length(diffExpr)) { write.csv(diffExpr[[i]], file=paste0("DifferentialExpression/",names(diffExpr[i]),".csv")) }


## ******************************** Differential expression by subgroup only, ages mixed  *********************************************
metaData
fullDataset <- DESeqDataSetFromMatrix(countData=bestDataCounts, colData=metaData, design= ~ condition + subject)   # mixed ages but control for subject
dds <- estimateSizeFactors(fullDataset)
temp <- counts(dds,normalized=TRUE)
idx <- rowSums( counts(dds, normalized=TRUE) >= 20 ) >= 20  # filter for genes with at least 20 counts in 25% of remaining samples
# write.csv(bestDataCounts[idx,],file="FilteredBestDataCounts.csv")
# write.csv(bestDataLog[idx,],file="FilteredBestDataLog.csv")

fullDataset <- fullDataset[idx,]

DESdata_mixedAges <- DESeq(fullDataset, parallel=TRUE)

diffExprAll <- list( HiHivsLoLo_v1 = results(DESdata_mixedAges, contrast=c("condition","HiHi_v1","LoLo_v1"), parallel = T),
                  HiHivsLoLo_v2 = results(DESdata_mixedAges, contrast=c("condition","HiHi_v2","LoLo_v2"), parallel = T), 
                  HiHivsnaive_v1 = results(DESdata_mixedAges, contrast=c("condition","HiHi_v1","Naive_v1"), parallel = T), 
                  HiHivsnaive_v2 = results(DESdata_mixedAges, contrast=c("condition","HiHi_v2","Naive_v2"), parallel = T), 
                  LoLovsnaive_v1 = results(DESdata_mixedAges, contrast=c("condition","LoLo_v1","Naive_v1"), parallel = T), 
                  LoLovsnaive_v2 = results(DESdata_mixedAges, contrast=c("condition","LoLo_v2","Naive_v2"), parallel = T) 
)

# for (i in 1:length(diffExprAll)) { write.csv(diffExprAll[[i]], file=paste0("DifferentialExpression/",names(diffExprAll[i]),"_allAges.csv")) }


#  ************************* Differential expression for YOUNG controlling for subject *********************************
fullDataset_YoungOnly_controlSubj <- DESeqDataSetFromMatrix(countData=bestDataCounts[,1:36], colData=metaData[1:36,], design= ~ subgroup + subject)
dds <- estimateSizeFactors(fullDataset_YoungOnly_controlSubj)
temp <- counts(dds,normalized=TRUE)
idx <- rowSums( counts(dds, normalized=TRUE) >= 20 ) >= 9  # filter for genes with at least 10 counts in 25% of the remaining samples
fullDataset_YoungOnly_controlSubj <- fullDataset_YoungOnly_controlSubj[idx,]

DESdata_fullDataset_YoungOnly_controlSubj <- DESeq(fullDataset_YoungOnly_controlSubj, parallel=TRUE)
# diffExprSingleYoung <- list( Y_HiHivsLoLo_v1 = results(DESdata_fullDataset_YoungOnly_controlSubj, contrast=c("subgroup","Y_HiHi_v1","Y_LoLo_v1"), parallel = T))

diffExpr_YoungOnly_controlSubj <- list(
  Y_HiHivsLoLo_v1 = results(DESdata_fullDataset_YoungOnly_controlSubj, contrast=c("subgroup","Y_HiHi_v1","Y_LoLo_v1"), parallel = T), 
  Y_HiHivsLoLo_v2 = results(DESdata_fullDataset_YoungOnly_controlSubj, contrast=c("subgroup","Y_HiHi_v2","Y_LoLo_v2"), parallel = T), 
  Y_HiHivsnaive_v1 = results(DESdata_fullDataset_YoungOnly_controlSubj, contrast=c("subgroup","Y_HiHi_v1","Y_Naive_v1"), parallel = T), 
  Y_HiHivsnaive_v2 = results(DESdata_fullDataset_YoungOnly_controlSubj, contrast=c("subgroup","Y_HiHi_v2","Y_Naive_v2"), parallel = T), 
  Y_LoLovsnaive_v1 = results(DESdata_fullDataset_YoungOnly_controlSubj, contrast=c("subgroup","Y_LoLo_v1","Y_Naive_v1"), parallel = T), 
  Y_LoLovsnaive_v2 = results(DESdata_fullDataset_YoungOnly_controlSubj, contrast=c("subgroup","Y_LoLo_v2","Y_Naive_v2"), parallel = T), 
  Y_HiHi_v12 = results(DESdata_fullDataset_YoungOnly_controlSubj, contrast=c("subgroup","Y_HiHi_v2","Y_HiHi_v1"), parallel = T), 
  Y_LoLo_v12 = results(DESdata_fullDataset_YoungOnly_controlSubj, contrast=c("subgroup","Y_LoLo_v2","Y_LoLo_v1"), parallel = T), 
  Y_naivevnaive_v12 = results(DESdata_fullDataset_YoungOnly_controlSubj, contrast=c("subgroup","Y_Naive_v2","Y_Naive_v1"), parallel = T) )


# for (i in 1:length(diffExpr_YoungOnly_controlSubj)) { write.csv(diffExpr_YoungOnly_controlSubj[[i]], file=paste0("DifferentialExpression/controlSubject_",names(diffExpr_YoungOnly_controlSubj[i]),".csv")) }



#  ************************* Differential expression for ELDERLY controlling for subject *********************************
fullDataset_ElderlyOnly_controlSubj <- DESeqDataSetFromMatrix(countData=bestDataCounts[,37:83], colData=metaData[37:83,], design= ~ subgroup + subject)
dds <- estimateSizeFactors(fullDataset_ElderlyOnly_controlSubj)
temp <- counts(dds,normalized=TRUE)
idx <- rowSums( counts(dds, normalized=TRUE) >= 20 ) >= 11  # filter for genes with at least 10 counts in 25% of the remaining samples
fullDataset_ElderlyOnly_controlSubj <- fullDataset_ElderlyOnly_controlSubj[idx,]

DESdata_fullDataset_ElderlyOnly_controlSubj <- DESeq(fullDataset_ElderlyOnly_controlSubj, parallel=TRUE)

diffExpr_ElderlyOnly_controlSubj <- list(
  E_HiHivsLoLo_v1 = results(DESdata_fullDataset_ElderlyOnly_controlSubj, contrast=c("subgroup","E_HiHi_v1","E_LoLo_v1"), parallel = T), 
  E_HiHivsLoLo_v2 = results(DESdata_fullDataset_ElderlyOnly_controlSubj, contrast=c("subgroup","E_HiHi_v2","E_LoLo_v2"), parallel = T), 
  E_HiHivsnaive_v1 = results(DESdata_fullDataset_ElderlyOnly_controlSubj, contrast=c("subgroup","E_HiHi_v1","E_Naive_v1"), parallel = T), 
  E_HiHivsnaive_v2 = results(DESdata_fullDataset_ElderlyOnly_controlSubj, contrast=c("subgroup","E_HiHi_v2","E_Naive_v2"), parallel = T), 
  E_LoLovsnaive_v1 = results(DESdata_fullDataset_ElderlyOnly_controlSubj, contrast=c("subgroup","E_LoLo_v1","E_Naive_v1"), parallel = T), 
  E_LoLovsnaive_v2 = results(DESdata_fullDataset_ElderlyOnly_controlSubj, contrast=c("subgroup","E_LoLo_v2","E_Naive_v2"), parallel = T), 
  E_HiHi_v12 = results(DESdata_fullDataset_ElderlyOnly_controlSubj, contrast=c("subgroup","E_HiHi_v2","E_HiHi_v1"), parallel = T), 
  E_LoLo_v12 = results(DESdata_fullDataset_ElderlyOnly_controlSubj, contrast=c("subgroup","E_LoLo_v2","E_LoLo_v1"), parallel = T), 
  E_naivevnaive_v12 = results(DESdata_fullDataset_ElderlyOnly_controlSubj, contrast=c("subgroup","E_Naive_v2","E_Naive_v1"), parallel = T) )


# for (i in 1:length(diffExpr_YoungOnly_controlSubj)) { write.csv(diffExpr_ElderlyOnly_controlSubj[[i]], file=paste0("DifferentialExpression/controlSubject_",names(diffExpr_ElderlyOnly_controlSubj[i]),".csv")) }

# *********************************************************************************************************************


# save(diffExpr, diffExpr_YoungOnly_controlSubj, diffExpr_ElderlyOnly_controlSubj, file="Rsavedimages/20180706_DiffExp.Rdata")
# load(file="Rsavedimages/20180706_DiffExp.Rdata")

## *****************************     volcano plots for comparisons  ***************************************
data <- as.data.frame(diffExpr_YoungOnly_controlSubj[["Y_HiHivsLoLo_v1"]])
data <- data[which(data$lfcSE<2),]     # filter out very high SE which may be outliers
ggplot(data, aes(x=log2FoldChange, y=-log10(padj), label=row.names(data))) +
  geom_point(alpha=0.2, size=3, colour="black") + theme_bw() +  theme(legend.position="none") + 
#  geom_text(aes(label=row.names(data)), cex=4) +  
  geom_text_repel(data=subset(data, padj<0.001), aes(label=row.names(subset(data, padj<0.001)))) + 
  labs(title="HiHi vs LoLo visit 1") + xlab("log2 fold change") + ylab("-log10 p-adj") + ylim(0,13) + xlim(-11,11)
# ggsave(filename="DifferentialExpression/Images/Young_HiHi_vs_LoLo_v1.png", device="png")
plotMA(diffExpr_YoungOnly_controlSubj[["Y_HiHivsLoLo_v1"]],ylim=c(-11,11))

data <- as.data.frame(diffExpr_ElderlyOnly_controlSubj[["E_HiHivsLoLo_v1"]])
data <- data[which(data$lfcSE<2),]     # filter out very high SE which may be outliers
ggplot(data, aes(x=log2FoldChange, y=-log10(padj), label=row.names(data))) +
  geom_point(alpha=0.2, size=3, colour="black") + theme_bw() +  theme(legend.position="none") + 
  #  geom_text(aes(label=row.names(data)), cex=4) +  
  geom_text_repel(data=subset(data, padj<0.00005), aes(label=row.names(subset(data, padj<0.00005)))) + 
  labs(title="HiHi vs LoLo visit 1 - ELDERLY") + xlab("log2 fold change") + ylab("-log10 p-adj") + ylim(0,13) + xlim(-11,11)
# ggsave(filename="DifferentialExpression/Images/Elderly_HiHi_vs_LoLo_v1.png", device="png")
plotMA(diffExpr_ElderlyOnly_controlSubj[["E_HiHivsLoLo_v1"]],ylim=c(-11,11))




# now look at YOUNG day 7 vs day 0
  data <- as.data.frame(diffExpr_YoungOnly_controlSubj[["Y_HiHi_v12"]])
ggplot(data, aes(x=log2FoldChange, y=-log10(padj), label=row.names(data))) +
  geom_point(alpha=0.2, size=3, colour="black") + theme_bw() +  theme(legend.position="none") + 
  #  geom_text(aes(label=row.names(data)), cex=4) +  
  geom_text_repel(data=subset(data, padj<0.01), aes(label=row.names(subset(data, padj<0.01)))) + 
  labs(title="Young ICOS+CD38+ cTfh:  day 7 vs day 0") + xlab("log2 fold change") + ylab("-log10 p-adj") + ylim(0,5) + xlim(-11,11)
# ggsave(filename="DifferentialExpression/Images/Young_HiHi_v1-v2.png", device="png")
plotMA(diffExpr_YoungOnly_controlSubj[["Y_HiHi_v12"]],ylim=c(-11,11))

probeList <- as.data.frame(head(data[order(data$stat,decreasing=T),],20)); probeList <- rbind(probeList, head(data[order(data$stat,decreasing=F),],20))
probeList <- probeList[which(probeList$lfcSE < 2),]
probeGenes <- bestDataLog[rownames(probeList),grep("HiHi",colnames(bestDataLog[1:36]))]
probeGenes <- cbind(probeGenes[,grep("v1",colnames(probeGenes))], probeGenes[,grep("v2",colnames(probeGenes))])

annotateHeatmap <- data.frame(row.names = colnames(probeGenes), Day = c(rep("day 0", 6),rep("day 7", 6)))
ann_colors = list(  Day = c("day 0" ="orange", "day 7" = "orange4")  )
pheatmap(probeGenes, scale="row", cluster_col=F, annotation_col = annotateHeatmap, show_colnames=F, main="Differentially expressed - Young",
         gaps_col = c(6), annotation_colors = ann_colors, cutree_rows = 2, fontsize_row = 16, color=inferno(100)
#        , filename = "DifferentialExpression/Images/Y_hihi_V1vsV2_SelectedGenesHeatmap.pdf"
)

# now look at ELDERLY day 7 vs day 0
data <- as.data.frame(diffExpr_ElderlyOnly_controlSubj[["E_HiHi_v12"]])
data <- data[which(data$lfcSE<2),]     # filter out very high SE which may be outliers
ggplot(data, aes(x=log2FoldChange, y=-log10(padj), label=row.names(data))) +
  geom_point(alpha=0.2, size=3, colour="black") + theme_bw() +  theme(legend.position="none") + 
  #  geom_text(aes(label=row.names(data)), cex=4) +  
  geom_text_repel(data=subset(data, padj<0.05), aes(label=row.names(subset(data, padj<0.05)))) + 
  labs(title="Elderly ICOS+CD38+ cTfh:  day 7 vs day 0") + xlab("log2 fold change") + ylab("-log10 p-adj") + ylim(0,5) + xlim(-11,11)
# ggsave(filename="DifferentialExpression/Images/Elderly_HiHi_v1-v2.png", device="png")
plotMA(diffExpr_YoungOnly_controlSubj[["Y_HiHi_v12"]],ylim=c(-11,11))

probeList <- as.data.frame(head(data[order(data$stat,decreasing=T),],20)); probeList <- rbind(probeList, head(data[order(data$stat,decreasing=F),],20))
probeList <- probeList[which(probeList$lfcSE < 2),]
probeGenes <- bestDataLog[rownames(probeList),grep("HiHi",colnames(bestDataLog))]
probeGenes <- probeGenes[,grep("X2", colnames(probeGenes))]
probeGenes <- cbind(probeGenes[,grep("v1",colnames(probeGenes))], probeGenes[,grep("v2",colnames(probeGenes))])

annotateHeatmap <- data.frame(row.names = colnames(probeGenes), Day = c(rep("day 0", 7),rep("day 7", 8)))
ann_colors = list(  Day = c("day 0" ="mediumpurple1", "day 7" = "purple4")  )
pheatmap(probeGenes, scale="row", cluster_col=F, annotation_col = annotateHeatmap, show_colnames=F, main="Differentially expressed - Elderly",
         gaps_col = c(7), annotation_colors = ann_colors, cutree_rows = 2, fontsize_row = 12, color=inferno(100)
#         , filename = "DifferentialExpression/Images/E_hihi_V1vsV2_SelectedGenesHeatmap.pdf"
)

# now compare YOUNG and ELDERLY directly for ICOS+CD38+ cTfh at day 7
data <- as.data.frame(diffExpr[["YvE_HiHi_v2"]])
data <- data[which(data$lfcSE<1.5),]     # filter out very high SE which may be outliers
ggplot(data, aes(x=log2FoldChange, y=-log10(padj), label=row.names(data))) +
  geom_point(alpha=0.2, size=3, colour="black") + theme_bw() +  theme(legend.position="none") + 
#  geom_text(aes(label=row.names(data)), cex=4) +  
  geom_text_repel(data=subset(data, padj<0.05), aes(label=row.names(subset(data, padj<0.05)))) + 
  labs(title="Young vs Elderly ICOS+CD38+ cTfh at day 7") + xlab("log2 fold change") + ylab("-log10 p-adj") + ylim(0,3.5) + xlim(-11,11)
# ggsave(filename="DifferentialExpression/Images/YvE_HiHi_v2.png")
plotMA(diffExpr[["YvE_HiHi_v2"]],ylim=c(-11,11))


    
probeList <- as.data.frame(head(data[order(data$stat,decreasing=T),],20)); probeList <- rbind(probeList, head(data[order(data$stat,decreasing=F),],20))
probeList <- probeList[which(probeList$lfcSE < 1.25),]
probeGenes <- bestDataLog[rownames(probeList),grep("HiHi_v2",colnames(bestDataLog))]

annotateHeatmap <- data.frame(row.names = colnames(probeGenes), ageGroup = c(rep("Young", 6),rep("Elderly", 8)))
ann_colors = list(  ageGroup = c("Young" ="orange", "Elderly" = "purple")  )
pheatmap(probeGenes, scale="row", cluster_col=F, annotation_col = annotateHeatmap, show_colnames=F, main="Selected Tfh genes",
         gaps_col = c(6), annotation_colors = ann_colors, cutree_rows = 2, fontsize_row = 16, color=inferno(100)
         #, filename = "Images/YvE_hihi_v2_SelectedGenesHeatmap.pdf"
)

# dev.off(); dev.off(); dev.off(); dev.off();


## *****************************       Gene Ontology replot    ***********************************************
# take analyses from Metascape.org and replot them using pheatmap for uniformity
# metascape performed using all padj < 0.05 genes

geneOntology <- read.csv("DifferentialExpression/GeneOntology/20181113_AllAges_Hihi-v-Lolo_v1_Metascape_Comparative/Enrichment_heatmap/HeatmapSelectedGO.csv", stringsAsFactors = F)
rownames(geneOntology) <- geneOntology$Description  #paste0(geneOntology$GO, ": ", geneOntology$Description)
geneOntology <- geneOntology[-grep("hsa", geneOntology$GO),]; geneOntology <- geneOntology[-grep("HSA", geneOntology$GO),]
geneOntology$X_LogP_Enriched_HiHi <- -geneOntology$X_LogP_Enriched_HiHi;  geneOntology$X_LogP_Enriched_LoLo <- -geneOntology$X_LogP_Enriched_LoLo  # now legend is -log P value
geneOntology$GO <- geneOntology$Description <- NULL
annotateHeatmap <- data.frame(row.names = colnames(geneOntology), Subset = c(rep("ICOS+CD38+", 1),rep("ICOS-CD38-", 1)))
ann_colors = list(  Subset = c("ICOS+CD38+" ="#0D0887", "ICOS-CD38-" = "#E16462")  )
inferno_100 <- inferno(n = 100); infernoBias <- colorRampPalette(inferno_100,bias=1.65)(100)
pheatmap(geneOntology, scale="none", cluster_col=F, cluster_row=F,show_colnames=F, main="cTfh subsets",
         annotation_colors = ann_colors, annotation_col = annotateHeatmap, fontsize_row = 11, 
         color=infernoBias
         #, filename = "DifferentialExpression/GeneOntology/20181113_AllAges_Hihi-v-Lolo_v1_Metascape_Comparative/AllAges_HiHi-v-Lolo_v1_geneOntology.pdf"
)


## *****************************       GSEA results: Hallmark Genesets    ***********************************************
# all together at baseline
Hallmark_HivsLo_v1_Pos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_Hallmark.GseaPreranked.1539515189801/gsea_report_for_na_pos_1539515189801.xls", sep="\t")
Hallmark_HivsLo_v1_Neg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_Hallmark.GseaPreranked.1539515189801/gsea_report_for_na_neg_1539515189801.xls", sep="\t")
Hallmark_HivsLo_v1 <- rbind(Hallmark_HivsLo_v1_Pos, Hallmark_HivsLo_v1_Neg)

Hallmark_HivsLo_v1$NAME <- toTitleCase(tolower(substr(Hallmark_HivsLo_v1$NAME,start =10, stop=50)))
Hallmark_HivsLo_v1$NAME <- factor(Hallmark_HivsLo_v1$NAME, levels = Hallmark_HivsLo_v1$NAME[order(Hallmark_HivsLo_v1$NES, decreasing = T)])
ggplot(data=subset(Hallmark_HivsLo_v1, `FDR.q.val` < 0.05), aes(x=`NAME`, y=`NES`)) + geom_point(size=3) + 
  geom_bar( aes(x=`NAME`, y=`NES`) , stat="Identity", width=0.01, color="#0D0887", size=1) +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets\nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12))
# ggsave(file="DifferentialExpression/GSEA/Images/AllAges_hihi-v-lolo_v1_Hallmark.pdf", device="pdf", height=7, width=5.5)

# trim down to gene sets that are likely to be important for Tfh biology
genesetExclusion <- c(grep("Spermatogenesis",Hallmark_HivsLo_v1$NAME), grep("Bile_acid",Hallmark_HivsLo_v1$NAME), grep("Myogenesis",Hallmark_HivsLo_v1$NAME), 
                      grep("Coagulation",Hallmark_HivsLo_v1$NAME),grep("Angiogenesis",Hallmark_HivsLo_v1$NAME),grep("Xenob",Hallmark_HivsLo_v1$NAME), 
                      grep("Myc_targets_v1",Hallmark_HivsLo_v1$NAME), grep("Heme",Hallmark_HivsLo_v1$NAME),grep("Uv_response",Hallmark_HivsLo_v1$NAME),
                      grep("Complement",Hallmark_HivsLo_v1$NAME),grep("Apical",Hallmark_HivsLo_v1$NAME),grep("Estrogen",Hallmark_HivsLo_v1$NAME),
                      grep("Epithelial",Hallmark_HivsLo_v1$NAME),grep("Allograft",Hallmark_HivsLo_v1$NAME))
Hallmark_HivsLo_v1 <- Hallmark_HivsLo_v1[-genesetExclusion,]
ggplot(data=subset(Hallmark_HivsLo_v1, `FDR.q.val` < 0.05), aes(x=`NAME`, y=`NES`)) + geom_point(size=5) + 
  geom_bar( aes(x=`NAME`, y=`NES`) , stat="Identity", width=0.01, color="#0D0887", size=1.5) +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets\nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )

# ggsave(file="DifferentialExpression/GSEA/Images/AllAges_hihi-v-lolo_v1_Hallmark_limited.pdf", device="pdf", height=7, width=5.5)


# Oxidative Phosphorylation
AllAgesHiHivLoLov1 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_Hallmark.GseaPreranked.1539515189801/HALLMARK_OXIDATIVE_PHOSPHORYLATION.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(Hallmark_HivsLo_v1$NES[grep("OXIDATIVE",Hallmark_HivsLo_v1$NAME)],2), "\n", "FDR: ", formatC(Hallmark_HivsLo_v1$FDR.q.val[grep("OXIDATIVE",Hallmark_HivsLo_v1$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.75,  y=0.82, hjust=0, gp=gpar(col="black", fontsize=12)))
ggplot(data=AllAgesHiHivLoLov1, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="black", size=1) + 
  geom_rug(sides="b", size=0.75, alpha=0.5, color="black") +   
  ggtitle("Oxidative phosphorylation - Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0))+
  annotation_custom(my_grob1) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/AllAges_hihi-vs-lolo_d0_OxPhos.pdf", device="pdf", height=3.5, width=5)





# young at baseline
Hallmark_HivsLo_v1_Pos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_HALLMARK.GseaPreranked.1530926563560/gsea_report_for_na_pos_1530926563560.xls", sep="\t")
Hallmark_HivsLo_v1_Neg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_HALLMARK.GseaPreranked.1530926563560/gsea_report_for_na_neg_1530926563560.xls", sep="\t")
Hallmark_HivsLo_v1 <- rbind(Hallmark_HivsLo_v1_Pos, Hallmark_HivsLo_v1_Neg)

Hallmark_HivsLo_v1$NAME <- toTitleCase(tolower(substr(Hallmark_HivsLo_v1$NAME,start =10, stop=50)))
Hallmark_HivsLo_v1$NAME <- factor(Hallmark_HivsLo_v1$NAME, levels = Hallmark_HivsLo_v1$NAME[order(Hallmark_HivsLo_v1$NES, decreasing = T)])
ggplot(data=subset(Hallmark_HivsLo_v1, `FDR.q.val` < 0.05), aes(x=`NAME`, y=`NES`)) + geom_point(size=3) + 
  geom_bar( aes(x=`NAME`, y=`NES`) , stat="Identity", width=0.01, color="orange3", size=1) +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets YOUNG\nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12))
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_Hallmark.pdf", device="pdf", height=7, width=5.5, dpi=300)

# trim down to gene sets that are likely to be important for Tfh biology
genesetExclusion <- c(grep("Spermatogenesis",Hallmark_HivsLo_v1$NAME), grep("Bile_acid",Hallmark_HivsLo_v1$NAME), grep("Myogenesis",Hallmark_HivsLo_v1$NAME), 
                      grep("Coagulation",Hallmark_HivsLo_v1$NAME),grep("Angiogenesis",Hallmark_HivsLo_v1$NAME),grep("Xenob",Hallmark_HivsLo_v1$NAME), 
                      grep("Myc_targets_v1",Hallmark_HivsLo_v1$NAME), grep("Heme",Hallmark_HivsLo_v1$NAME),grep("Uv_response",Hallmark_HivsLo_v1$NAME),
                      grep("Complement",Hallmark_HivsLo_v1$NAME),grep("Apical",Hallmark_HivsLo_v1$NAME),grep("Estrogen",Hallmark_HivsLo_v1$NAME),
                      grep("Epithelial",Hallmark_HivsLo_v1$NAME),grep("Allograft",Hallmark_HivsLo_v1$NAME))
Hallmark_HivsLo_v1 <- Hallmark_HivsLo_v1[-genesetExclusion,]
ggplot(data=subset(Hallmark_HivsLo_v1, `FDR.q.val` < 0.05), aes(x=`NAME`, y=`NES`)) + geom_point(size=5) + 
  geom_bar( aes(x=`NAME`, y=`NES`) , stat="Identity", width=0.01, color="orange3", size=1.5) +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets YOUNG\nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )

# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_Hallmark_limited.pdf", device="pdf", height=7, width=5.5, dpi=300)


# elderly at baseline 
Hallmark_HivsLo_v1_Pos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_HALLMARK.GseaPreranked.1531837435104/gsea_report_for_na_pos_1531837435104.xls", sep="\t")
Hallmark_HivsLo_v1_Neg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_HALLMARK.GseaPreranked.1531837435104/gsea_report_for_na_neg_1531837435104.xls", sep="\t")
Hallmark_HivsLo_v1 <- rbind(Hallmark_HivsLo_v1_Pos, Hallmark_HivsLo_v1_Neg)

Hallmark_HivsLo_v1$NAME <- toTitleCase(tolower(substr(Hallmark_HivsLo_v1$NAME,start =10, stop=50)))
Hallmark_HivsLo_v1$NAME <- factor(Hallmark_HivsLo_v1$NAME, levels = Hallmark_HivsLo_v1$NAME[order(Hallmark_HivsLo_v1$NES, decreasing = T)])
ggplot(data=subset(Hallmark_HivsLo_v1, `FDR.q.val` < 0.05), aes(x=`NAME`, y=`NES`)) + geom_point(size=3) + 
  geom_bar( aes(x=`NAME`, y=`NES`) , stat="Identity", width=0.01, color="purple", size=1) +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets ELDERLY \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12))
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_Hallmark.pdf", device="pdf", height=7, width=5.5)

# trim down to gene sets that are likely to be important for Tfh biology
genesetExclusion <- c(grep("Spermatogenesis",Hallmark_HivsLo_v1$NAME), grep("Bile_acid",Hallmark_HivsLo_v1$NAME), grep("Myogenesis",Hallmark_HivsLo_v1$NAME), 
                      grep("Coagulation",Hallmark_HivsLo_v1$NAME),grep("Angiogenesis",Hallmark_HivsLo_v1$NAME),grep("Xenob",Hallmark_HivsLo_v1$NAME), 
                      grep("Myc_targets_v1",Hallmark_HivsLo_v1$NAME), grep("Heme",Hallmark_HivsLo_v1$NAME),grep("Uv_response",Hallmark_HivsLo_v1$NAME),
                      grep("Complement",Hallmark_HivsLo_v1$NAME),grep("Apical",Hallmark_HivsLo_v1$NAME),grep("Estrogen",Hallmark_HivsLo_v1$NAME),
                      grep("Epithelial",Hallmark_HivsLo_v1$NAME),grep("Allograft",Hallmark_HivsLo_v1$NAME))
Hallmark_HivsLo_v1 <- Hallmark_HivsLo_v1[-genesetExclusion,]
ggplot(data=subset(Hallmark_HivsLo_v1, `FDR.q.val` < 0.05), aes(x=`NAME`, y=`NES`)) + geom_point(size=5) + 
  geom_bar( aes(x=`NAME`, y=`NES`) , stat="Identity", width=0.01, color="purple", size=1.5) +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets ELDERLY \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )

# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_Hallmark_limited.pdf", device="pdf", height=7, width=5.5)




## *****************************       GSEA results: external Genesets    ***********************************************

# all together at baseline
ExternalGenesetsPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_ExternalGenesets.GseaPreranked.1539515274110/gsea_report_for_na_pos_1539515274110.xls", sep="\t")
ExternalGenesetsNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_ExternalGenesets.GseaPreranked.1539515274110/gsea_report_for_na_neg_1539515274110.xls", sep="\t")
ExternalGenesets <- rbind(ExternalGenesetsPos, ExternalGenesetsNeg)

LocciGSE50392 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_ExternalGenesets.GseaPreranked.1539515274110/GSE50392_GCTFH-VS-NAV.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(ExternalGenesets$NES[grep("GCTFH-VS-NAV",ExternalGenesets$NAME)],2), "\n", "FDR: ", formatC(ExternalGenesets$FDR.q.val[grep("GCTFH-VS-NAV",ExternalGenesets$NAME)], format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=LocciGSE50392, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE50392 HUMAN Tonsil GC-Tfh vs Naive") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/AllAges_hihi-v-lolo_v1_LocciGCTfh.pdf", device="pdf", height=3.5, width=5)
MOUSEGCTFHGSE16697 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_ExternalGenesets.GseaPreranked.1539515274110/GSE16697_MOUSEGCTFH.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(ExternalGenesets$NES[grep("16697",ExternalGenesets$NAME)],2), "\n", "FDR: ", formatC(ExternalGenesets$FDR.q.val[grep("16697",ExternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=MOUSEGCTFHGSE16697, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE16697 MOUSE Tfh vs non-Tfh CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/AllAges_hihi-v-lolo_v1_MouseGCTfh.pdf", device="pdf", height=3.5, width=5)
MOUSEGCTFHGSE32596 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_ExternalGenesets.GseaPreranked.1539515274110/GSE32596_MOUSEGCTFH.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(ExternalGenesets$NES[grep("32596",ExternalGenesets$NAME)],2), "\n", "FDR: ", formatC(ExternalGenesets$FDR.q.val[grep("32596",ExternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=MOUSEGCTFHGSE32596, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE32596 MOUSE Tfh vs Naive CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/AllAges_hihi-v-lolo_v1_MouseGCTfh2.pdf", device="pdf", height=3.5, width=5)


YexternalGenesetsPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/gsea_report_for_na_pos_1531001576554.xls", sep="\t")
YexternalGenesetsNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/gsea_report_for_na_neg_1531001576554.xls", sep="\t")
YexternalGenesets <- rbind(YexternalGenesetsPos, YexternalGenesetsNeg)
EexternalGenesetsPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_ExtraGeneSets.GseaPreranked.1531837454581/gsea_report_for_na_pos_1531837454581.xls", sep="\t")
EexternalGenesetsNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_ExtraGeneSets.GseaPreranked.1531837454581/gsea_report_for_na_neg_1531837454581.xls", sep="\t")
EexternalGenesets <- rbind(EexternalGenesetsPos, EexternalGenesetsNeg)

LocciGSE50392 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/GCTFH-VS-NAV.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(YexternalGenesets$NES[grep("GCTFH-VS-NAV",YexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(YexternalGenesets$FDR.q.val[grep("GCTFH-VS-NAV",YexternalGenesets$NAME)], format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=LocciGSE50392, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE50392 HUMAN Tonsil GC-Tfh vs Naive") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_LocciGCTfh.pdf", device="pdf", height=3.5, width=5)

LocciGSE50392 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_ExtraGeneSets.GseaPreranked.1531837454581/GSE50392_GCTFH-VS-NAV.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(EexternalGenesets$NES[grep("GCTFH-VS-NAV",EexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(EexternalGenesets$FDR.q.val[grep("GCTFH-VS-NAV",EexternalGenesets$NAME)], format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=LocciGSE50392, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE50392 HUMAN Tonsil GC-Tfh vs Naive") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_LocciGCTfh.pdf", device="pdf", height=3.5, width=5)

MOUSEGCTFHGSE16697 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/GSE16697_MOUSEGCTFH.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(YexternalGenesets$NES[grep("16697",YexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(YexternalGenesets$FDR.q.val[grep("16697",YexternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=MOUSEGCTFHGSE16697, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE16697 MOUSE Tfh vs non-Tfh CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_MouseGCTfh.pdf", device="pdf", height=3.5, width=5)

MOUSEGCTFHGSE16697 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_ExtraGeneSets.GseaPreranked.1531837454581/GSE16697_MOUSEGCTFH.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(EexternalGenesets$NES[grep("16697",EexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(EexternalGenesets$FDR.q.val[grep("16697",EexternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=MOUSEGCTFHGSE16697, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE16697 MOUSE Tfh vs non-Tfh CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_MouseGCTfh.pdf", device="pdf", height=3.5, width=5)

MOUSEGCTFHGSE32596 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/GSE32596_MOUSEGCTFH.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(YexternalGenesets$NES[grep("32596",YexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(YexternalGenesets$FDR.q.val[grep("32596",YexternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=MOUSEGCTFHGSE32596, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE32596 MOUSE Tfh vs Naive CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_MouseGCTfh2.pdf", device="pdf", height=3.5, width=5)

MOUSEGCTFHGSE32596 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_ExtraGeneSets.GseaPreranked.1531837454581/GSE32596_MOUSEGCTFH.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(EexternalGenesets$NES[grep("32596",EexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(EexternalGenesets$FDR.q.val[grep("32596",EexternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=MOUSEGCTFHGSE32596, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE32596 MOUSE Tfh vs Naive CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_MouseGCTfh2.pdf", device="pdf", height=3.5, width=5)



TCF1KOGSE65693 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/GSE65693_TCF1KO.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(YexternalGenesets$NES[grep("65693",YexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(YexternalGenesets$FDR.q.val[grep("65693",YexternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=TCF1KOGSE65693, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE65693 TCF1 knockout") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_TCF1ko.png", device="png", height=3.5, width=5)


## *****************************       GSEA results: Y and E D0-to-D7 for Hallmark    ***********************************************
Yd0tod7HallmarkPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/gsea_report_for_na_pos_1530926605818.xls", sep="\t")
Yd0tod7HallmarkNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/gsea_report_for_na_neg_1530926605818.xls", sep="\t")
Yd0tod7Hallmark <- rbind(Yd0tod7HallmarkPos, Yd0tod7HallmarkNeg)
Ed0tod7HallmarkPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/gsea_report_for_na_pos_1530926588879.xls", sep="\t")
Ed0tod7HallmarkNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/gsea_report_for_na_neg_1530926588879.xls", sep="\t")
Ed0tod7Hallmark <- rbind(Ed0tod7HallmarkPos, Ed0tod7HallmarkNeg)


# TGF-b Targets 
YoungHiHiv2TGF <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/HALLMARK_TGF_BETA_SIGNALING.xls", sep="\t")
ElderlyHiHiv2TGF <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/HALLMARK_TGF_BETA_SIGNALING.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(Yd0tod7Hallmark$NES[grep("TGF",Yd0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Yd0tod7Hallmark$FDR.q.val[grep("TGF",Yd0tod7Hallmark$NAME)], format = "e", digits = 1))
annotationInfo2 <- paste0("\nNES: ", round(Ed0tod7Hallmark$NES[grep("TGF",Ed0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Ed0tod7Hallmark$FDR.q.val[grep("TGF",Ed0tod7Hallmark$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.32, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.18, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=YoungHiHiv2TGF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="orange3") + 
  geom_line(data=ElderlyHiHiv2TGF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="purple")+
  geom_rug(sides="t", size=0.75, alpha=0.5, color="orange3") +  geom_rug(data=ElderlyHiHiv2TGF, sides="b", size=0.75, alpha=0.5, color="purple") +  
  ggtitle("TGF-b Signaling - Day 7 vs Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_TGFbeta.pdf", device="pdf", height=3.5, width=5)

# E2F Targets
YoungHiHiv2E2F <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/HALLMARK_E2F_TARGETS.xls", sep="\t")
ElderlyHiHiv2E2F <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/HALLMARK_E2F_TARGETS.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(Yd0tod7Hallmark$NES[grep("E2F_TARGETS",Yd0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Yd0tod7Hallmark$FDR.q.val[grep("E2F_TARGETS",Yd0tod7Hallmark$NAME)], format = "e", digits = 1))
annotationInfo2 <- paste0("\nNES: ", round(Ed0tod7Hallmark$NES[grep("E2F_TARGETS",Ed0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Ed0tod7Hallmark$FDR.q.val[grep("E2F_TARGETS",Ed0tod7Hallmark$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.88, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.7,  y=0.75, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=YoungHiHiv2E2F, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="orange3") + 
  geom_line(data=ElderlyHiHiv2E2F, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="purple")+
  geom_rug(sides="t", size=0.75, alpha=0.5, color="orange3") +  geom_rug(data=ElderlyHiHiv2E2F, sides="b", size=0.75, alpha=0.5, color="purple") +  
  ggtitle("E2F Targets - Day 7 vs Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_E2FTargets.pdf", device="pdf", height=3.5, width=5)

# IL2 signaling
YoungHiHiv2IL2 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/HALLMARK_IL2_STAT5_SIGNALING.xls", sep="\t")
ElderlyHiHiv2IL2 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/HALLMARK_IL2_STAT5_SIGNALING.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(Yd0tod7Hallmark$NES[grep("IL2",Yd0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Yd0tod7Hallmark$FDR.q.val[grep("IL2",Yd0tod7Hallmark$NAME)], format = "e", digits = 1))
annotationInfo2 <- paste0("\nNES: ", round(Ed0tod7Hallmark$NES[grep("IL2",Ed0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Ed0tod7Hallmark$FDR.q.val[grep("IL2",Ed0tod7Hallmark$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.88, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.7,  y=0.75, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=YoungHiHiv2IL2, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="orange3") + 
  geom_line(data=ElderlyHiHiv2IL2, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="purple")+
  geom_rug(sides="b", size=0.75, alpha=0.5, color="orange3") +  geom_rug(data=ElderlyHiHiv2IL2, sides="t", size=0.75, alpha=0.5, color="purple") +  
  ggtitle("IL2 signaling - Day 7 vs Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_IL2-STAT5.pdf", device="pdf", height=3.5, width=5)

# TNF - NFkb
YoungHiHiv2TNF <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
ElderlyHiHiv2TNF <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(Yd0tod7Hallmark$NES[grep("TNF",Yd0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Yd0tod7Hallmark$FDR.q.val[grep("TNF",Yd0tod7Hallmark$NAME)], format = "e", digits = 1))
annotationInfo2 <- paste0("\nNES: ", round(Ed0tod7Hallmark$NES[grep("TNF",Ed0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Ed0tod7Hallmark$FDR.q.val[grep("TNF",Ed0tod7Hallmark$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.88, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.7,  y=0.75, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=YoungHiHiv2TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="orange3") + 
  geom_line(data=ElderlyHiHiv2TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="purple")+
  geom_rug(sides="b", size=0.75, alpha=0.5, color="orange3") +  geom_rug(data=ElderlyHiHiv2TNF, sides="t", size=0.75, alpha=0.5, color="purple") +  
  ggtitle("TNF signaling - Day 7 vs Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_TNF-NFkb.pdf", device="pdf", height=3.5, width=5)

# Apoptosis
YoungHiHiv2Apoptosis <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/HALLMARK_APOPTOSIS.xls", sep="\t")
ElderlyHiHiv2Apoptosis <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/HALLMARK_APOPTOSIS.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(Yd0tod7Hallmark$NES[grep("APOPTOSIS",Yd0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Yd0tod7Hallmark$FDR.q.val[grep("APOPTOSIS",Yd0tod7Hallmark$NAME)], format = "e", digits = 1))
annotationInfo2 <- paste0("\nNES: ", round(Ed0tod7Hallmark$NES[grep("APOPTOSIS",Ed0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Ed0tod7Hallmark$FDR.q.val[grep("APOPTOSIS",Ed0tod7Hallmark$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.88, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.7,  y=0.75, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=YoungHiHiv2Apoptosis, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="orange3") + 
  geom_line(data=ElderlyHiHiv2Apoptosis, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="purple")+
  geom_rug(sides="b", size=0.75, alpha=0.5, color="orange3") +  geom_rug(data=ElderlyHiHiv2Apoptosis, sides="t", size=0.75, alpha=0.5, color="purple") +  
  ggtitle("Apoptosis - Day 7 vs Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_Apoptosis.pdf", device="pdf", height=3.5, width=5)


## *****************************       GSEA results: YvE parallel d7 vs d0 comparison - Hallmark    ***********************************************


Youngv1v2Pos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/gsea_report_for_na_pos_1530926605818.xls", sep="\t", stringsAsFactors = F)
Youngv1v2Neg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/gsea_report_for_na_neg_1530926605818.xls", sep="\t", stringsAsFactors = F)
Youngv1v2 <- rbind(Youngv1v2Pos, Youngv1v2Neg)

Elderlyv1v2Pos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/gsea_report_for_na_pos_1530926588879.xls", sep="\t", stringsAsFactors = F)
Elderlyv1v2Neg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/gsea_report_for_na_neg_1530926588879.xls", sep="\t", stringsAsFactors = F)
Elderlyv1v2 <- rbind(Elderlyv1v2Pos, Elderlyv1v2Neg)

Youngv1v2$NAME <- toTitleCase(tolower(substr(Youngv1v2$NAME,start =10, stop=50)))
Youngv1v2$NAME <- factor(Youngv1v2$NAME, levels = Youngv1v2$NAME[order(Youngv1v2$NES, decreasing = T)])

ggplot(data=subset(Youngv1v2, `FDR.q.val` < 0.95), aes(x=`NAME`, y=`NES`)) + geom_point(size=3) + 
  geom_segment( aes(x=`NAME`, xend=`NAME`, y=0, yend=`NES`) , size=1, color="orange", linetype="dashed" ) +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets \nICOS+CD38+ cTfh for Day 0 vs Day 7") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14))
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v1-vs-v2_Hallmark.png", device="png", height=5, width=5.5, dpi=300)

Elderlyv1v2$NAME <- toTitleCase(tolower(substr(Elderlyv1v2$NAME,start =10, stop=50)))
Elderlyv1v2$NAME <- factor(Elderlyv1v2$NAME, levels = Elderlyv1v2$NAME[order(Elderlyv1v2$NES, decreasing = T)])

ggplot(data=subset(Elderlyv1v2, `FDR.q.val` < 0.95), aes(x=`NAME`, y=`NES`)) + geom_point(size=3) + 
  geom_segment( aes(x=`NAME`, xend=`NAME`, y=0, yend=`NES`) , size=1, color="purple", linetype="dashed" ) +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets \nICOS+CD38+ cTfh for Day 0 vs Day 7") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14))
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v1-vs-v2_Hallmark.png", device="png", height=5, width=5.5, dpi=300)


# merge the two datasets into one
tempYoungv1v2 <- rbind(Youngv1v2Pos, Youngv1v2Neg); tempYoungv1v2 <- cbind(tempYoungv1v2, ageGroup = "Young"); tempYoungv1v2$GS.br..follow.link.to.MSigDB <- NULL
tempYoungv1v2 <- tempYoungv1v2[order(tempYoungv1v2$NES, decreasing=T),]
tempElderlyv1v2 <- rbind(Elderlyv1v2Pos, Elderlyv1v2Neg); tempElderlyv1v2 <- cbind(tempElderlyv1v2, ageGroup = "Elderly"); tempElderlyv1v2$GS.br..follow.link.to.MSigDB <- NULL
hallmarkMerged <- rbind(tempYoungv1v2[order(tempYoungv1v2$NAME),], tempElderlyv1v2[order(tempElderlyv1v2$NAME),])

hallmarkMerged$NAME <- toTitleCase(tolower(substr(hallmarkMerged$NAME,start =10, stop=50)))
hallmarkMerged <- hallmarkMerged[order(hallmarkMerged$NES, decreasing = T),]

hallmarkMergedWide<- reshape(hallmarkMerged, idvar = "NAME", timevar="ageGroup", direction="wide")
hallmarkMergedWide$EminusY <- hallmarkMergedWide$NES.Elderly - hallmarkMergedWide$NES.Young
hallmarkMergedWide <- hallmarkMergedWide[order(hallmarkMergedWide$EminusY, decreasing=T),]

hallmarkMergedWide <- hallmarkMergedWide[order(hallmarkMergedWide$NES.Young, decreasing=T),]
hallmarkFactored <- within(hallmarkMergedWide, NAME<-factor(NAME, levels=hallmarkMergedWide$NAME))
ggplot(data=hallmarkFactored) + 
  geom_point(aes(x=`NAME`, y=`NES.Young`),size=3, color="orange3") + 
  geom_bar(aes(x=`NAME`, y=`NES.Young`), stat="Identity", width=0.2, color="orange") +
  geom_point(aes(x=`NAME`, y=`NES.Elderly`),size=3, color="purple") + 
  geom_bar(aes(x=`NAME`, y=`NES.Elderly`), stat="Identity", width=0.2, color="purple") +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets \nICOS+CD38+ cTfh, Day 7 vs Day 0") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=16), plot.title = element_text(size=12)) + geom_hline(aes(yintercept = 0))
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_v2_hallmark_full.pdf", device="pdf", height=9, width=5)


# exclude genesets that do not have known role in Tfh biology

genesetExclusion <- c(grep("Spermatogenesis",hallmarkFactored$NAME), grep("Bile_acid",hallmarkFactored$NAME), grep("Myogenesis",hallmarkFactored$NAME), 
                      grep("Coagulation",hallmarkFactored$NAME),grep("Angiogenesis",hallmarkFactored$NAME),grep("Xenob",hallmarkFactored$NAME), 
                      grep("Myc_targets_v1",hallmarkFactored$NAME), grep("Heme",hallmarkFactored$NAME),grep("Uv_response",hallmarkFactored$NAME),
                      grep("Complement",hallmarkFactored$NAME),grep("Apical",hallmarkFactored$NAME),grep("Estrogen",hallmarkFactored$NAME),
                      grep("Epithelial",hallmarkFactored$NAME),grep("Allograft",hallmarkFactored$NAME))
hallmarkFactored <- hallmarkFactored[-genesetExclusion,]
ggplot(data=hallmarkFactored) + 
  geom_point(aes(x=`NAME`, y=`NES.Young`),size=4, color="orange3") + 
  geom_bar(aes(x=`NAME`, y=`NES.Young`), stat="Identity", width=0.01, color="orange3", size=2) +
  geom_point(aes(x=`NAME`, y=`NES.Elderly`),size=4, color="purple") + 
  geom_bar(aes(x=`NAME`, y=`NES.Elderly`), stat="Identity", width=0.01, color="purple", size=2) +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets \nICOS+CD38+ cTfh, Day 7 vs Day 0") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title = element_text(size=16), axis.text = element_text(size=14), plot.title = element_text(size=16)) + geom_hline(aes(yintercept = 0))
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_v2_hallmark_limited.pdf", device="pdf", height=12, width=7)



## *****************************       GSEA results: YvE direct comparison at d7 - Hallmark    ***********************************************


hallmarkv2Pos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/gsea_report_for_na_pos_1531000206530.xls", sep="\t", stringsAsFactors = F)
hallmarkv2Neg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/gsea_report_for_na_neg_1531000206530.xls", sep="\t", stringsAsFactors = F)
hallmarkv2 <- rbind(hallmarkv2Pos, hallmarkv2Neg)

hallmarkv2$NAME <- toTitleCase(tolower(substr(hallmarkv2$NAME,start =10, stop=50)))
hallmarkv2$NAME <- factor(hallmarkv2$NAME, levels = hallmarkv2$NAME[order(hallmarkv2$NES, decreasing = T)])

ggplot(data=subset(hallmarkv2, `FDR.q.val` < 0.05), aes(x=`NAME`, y=`NES`)) + geom_point(size=3) + 
  geom_bar( aes(x=`NAME`, y=`NES`) , stat="Identity", width=0.01, color="#0D0887", size=1) +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets FDR<0.05\nICOS+CD38+\nYoung vs Elderly at day 7") + ylab("Normalized Enrichment Score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12))
# ggsave(file="DifferentialExpression/GSEA/Images/YvEhihi-v2_Hallmark.pdf", device="pdf", width=6, height=5)


# Apoptosis
YvEHiHiv2Apoptosis <- read.csv("DifferentialExpression/GSEA/GSEA_Results/YvE_HiHi_v2_HALLMARK.GseaPreranked.1540118672676/HALLMARK_APOPTOSIS.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(hallmarkv2$NES[grep("Apoptosis",hallmarkv2$NAME)],2), "\n", "FDR: ", formatC(hallmarkv2$FDR.q.val[grep("Apoptosis",hallmarkv2$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.63,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=YvEHiHiv2Apoptosis, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Apoptosis ICOS+CD38+ cTfh at day 7") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YvE_HiHi_d7_Hallmark_Apoptosis.pdf", device="pdf", height=3.5, width=5)

# TNF
YvEHiHiv2TNF<- read.csv("DifferentialExpression/GSEA/GSEA_Results/YvE_HiHi_v2_HALLMARK.GseaPreranked.1540118672676/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(hallmarkv2$NES[grep("Tnf",hallmarkv2$NAME)],2), "\n", "FDR: ", formatC(hallmarkv2$FDR.q.val[grep("Tnf",hallmarkv2$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.63,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=YvEHiHiv2TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
 ggtitle("TNF-NFkB ICOS+CD38+ cTfh at day 7") + ylab("Enrichment score") + xlab("Rank in gene list") + 
 # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
 theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
 annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YvE_HiHi_d7_Hallmark_TNF-NFkB.pdf", device="pdf", height=3.5, width=5)
 
 
# TNF and Apoptosis on one plot
annotationInfo <- paste0("NES: ", round(hallmarkv2$NES[grep("Apoptosis",hallmarkv2$NAME)],2), "\n", "FDR: ", formatC(hallmarkv2$FDR.q.val[grep("Apoptosis",hallmarkv2$NAME)], format = "e", digits = 1))
annotationInfo2 <- paste0("\nNES: ", round(hallmarkv2$NES[grep("Tnf",hallmarkv2$NAME)],2), "\n", "FDR: ", formatC(hallmarkv2$FDR.q.val[grep("Tnf",hallmarkv2$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.88, hjust=0, gp=gpar(col="#a52c60", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.7,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=12)))
ggplot(data=YvEHiHiv2Apoptosis, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="#a52c60") + 
  geom_line(data=YvEHiHiv2TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="black")+
  geom_rug(sides="b", size=0.75, alpha=0.5, color="#a52c60") +  geom_rug(data=YvEHiHiv2TNF, sides="t", size=0.75, alpha=0.5, color="black") +  
  ggtitle("ICOS+CD38+ cTfh at day 7") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YvE_HiHi_d7_TNFandApoptosis.pdf", device="pdf", height=3.5, width=5)


## *****************************       GSEA results: YvE direct comparison - HexagonPlots    ***********************************************

# viridis color scheme for discrete
# library(scales)
# show_col(viridis_pal(option="B")(10))
# show_col(colorRampPalette(c("#CF4446","white", "#ED6925"))(50))

hallmarkYvE_v1_Hi <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v1_HALLMARK.GseaPreranked.1531000213854/gsea_report_for_na_pos_1531000213854.xls", sep="\t", stringsAsFactors = F),
                           read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v1_HALLMARK.GseaPreranked.1531000213854/gsea_report_for_na_neg_1531000213854.xls", sep="\t", stringsAsFactors = F)
                           )
hallmarkYvE_v1_Hi$NAME <- toTitleCase(tolower(substr(hallmarkYvE_v1_Hi$NAME,start =10, stop=50)))
hallmarkYvE_v1_Hi$NAME <- factor(hallmarkYvE_v1_Hi$NAME, levels = hallmarkYvE_v1_Hi$NAME[order(hallmarkYvE_v1_Hi$NES, decreasing = T)])
rownames(hallmarkYvE_v1_Hi) <- hallmarkYvE_v1_Hi$NAME

hallmarkYvE_v2_Hi <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/gsea_report_for_na_pos_1531000206530.xls", sep="\t", stringsAsFactors = F),
                           read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/gsea_report_for_na_neg_1531000206530.xls", sep="\t", stringsAsFactors = F)
                           )
hallmarkYvE_v2_Hi$NAME <- toTitleCase(tolower(substr(hallmarkYvE_v2_Hi$NAME,start =10, stop=50)))
hallmarkYvE_v2_Hi$NAME <- factor(hallmarkYvE_v2_Hi$NAME, levels = hallmarkYvE_v2_Hi$NAME[order(hallmarkYvE_v2_Hi$NES, decreasing = T)])
rownames(hallmarkYvE_v2_Hi) <- hallmarkYvE_v2_Hi$NAME

hallmarkYvE_v1_Lo <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_LoLo_vs_E_LoLo_v1_Hallmark.GseaPreranked.1540914053735/gsea_report_for_na_pos_1540914053735.xls", sep="\t", stringsAsFactors = F),
                           read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_LoLo_vs_E_LoLo_v1_Hallmark.GseaPreranked.1540914053735/gsea_report_for_na_neg_1540914053735.xls", sep="\t", stringsAsFactors = F)
                           )
hallmarkYvE_v1_Lo$NAME <- toTitleCase(tolower(substr(hallmarkYvE_v1_Lo$NAME,start =10, stop=50)))
hallmarkYvE_v1_Lo$NAME <- factor(hallmarkYvE_v1_Lo$NAME, levels = hallmarkYvE_v1_Lo$NAME[order(hallmarkYvE_v1_Lo$NES, decreasing = T)])
rownames(hallmarkYvE_v1_Lo) <- hallmarkYvE_v1_Lo$NAME

hallmarkYvE_v2_Lo <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_LoLo_vs_E_LoLo_v2_Hallmark.GseaPreranked.1540914059884/gsea_report_for_na_pos_1540914059884.xls", sep="\t", stringsAsFactors = F),
                           read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_LoLo_vs_E_LoLo_v2_Hallmark.GseaPreranked.1540914059884/gsea_report_for_na_neg_1540914059884.xls", sep="\t", stringsAsFactors = F)
                           )
hallmarkYvE_v2_Lo$NAME <- toTitleCase(tolower(substr(hallmarkYvE_v2_Lo$NAME,start =10, stop=50)))
hallmarkYvE_v2_Lo$NAME <- factor(hallmarkYvE_v2_Lo$NAME, levels = hallmarkYvE_v2_Lo$NAME[order(hallmarkYvE_v2_Lo$NES, decreasing = T)])
rownames(hallmarkYvE_v2_Lo) <- hallmarkYvE_v2_Lo$NAME

hallmarkYvE_v1_Nav <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v1_Hallmark.GseaPreranked.1540914080384/gsea_report_for_na_pos_1540914080384.xls", sep="\t", stringsAsFactors = F),
                            read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v1_Hallmark.GseaPreranked.1540914080384/gsea_report_for_na_neg_1540914080384.xls", sep="\t", stringsAsFactors = F)
                            )
hallmarkYvE_v1_Nav$NAME <- toTitleCase(tolower(substr(hallmarkYvE_v1_Nav$NAME,start =10, stop=50)))
hallmarkYvE_v1_Nav$NAME <- factor(hallmarkYvE_v1_Nav$NAME, levels = hallmarkYvE_v1_Nav$NAME[order(hallmarkYvE_v1_Nav$NES, decreasing = T)])
rownames(hallmarkYvE_v1_Nav) <- hallmarkYvE_v1_Nav$NAME

hallmarkYvE_v2_Nav <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v2_Hallmark.GseaPreranked.1540914108109/gsea_report_for_na_pos_1540914108109.xls", sep="\t", stringsAsFactors = F), 
                            read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v2_Hallmark.GseaPreranked.1540914108109/gsea_report_for_na_neg_1540914108109.xls", sep="\t", stringsAsFactors = F)
                            )
hallmarkYvE_v2_Nav$NAME <- toTitleCase(tolower(substr(hallmarkYvE_v2_Nav$NAME,start =10, stop=50)))
hallmarkYvE_v2_Nav$NAME <- factor(hallmarkYvE_v2_Nav$NAME, levels = hallmarkYvE_v2_Nav$NAME[order(hallmarkYvE_v2_Nav$NES, decreasing = T)])
rownames(hallmarkYvE_v2_Nav) <- hallmarkYvE_v2_Nav$NAME

hallmark_subsets <- data.frame(row.names=hallmarkYvE_v1_Hi$NAME[order(hallmarkYvE_v1_Hi$NAME)],
                               v1Hi = hallmarkYvE_v1_Hi$NES[order(hallmarkYvE_v1_Hi$NAME)],
                               v2Hi = hallmarkYvE_v2_Hi$NES[order(hallmarkYvE_v2_Hi$NAME)],
                               v1Lo = hallmarkYvE_v1_Lo$NES[order(hallmarkYvE_v1_Lo$NAME)],
                               v2Lo = hallmarkYvE_v2_Lo$NES[order(hallmarkYvE_v2_Lo$NAME)],
                               v1Nav = hallmarkYvE_v1_Nav$NES[order(hallmarkYvE_v1_Nav$NAME)],
                               v2Nav = hallmarkYvE_v2_Nav$NES[order(hallmarkYvE_v2_Nav$NAME)]
                               )

hallmark_subsetsREVISED <- merge(hallmarkYvE_v1_Hi[,c("NAME","NES")], hallmarkYvE_v2_Hi[,c("NAME","NES")], by='NAME'); 
hallmark_subsetsREVISED <- merge(hallmark_subsetsREVISED, hallmarkYvE_v1_Lo[,c("NAME","NES")], by='NAME'); names(hallmark_subsetsREVISED) <- c("NAME","Hi_v1", "Hi_v2", "Lo_v1")
hallmark_subsetsREVISED <- merge(hallmark_subsetsREVISED, hallmarkYvE_v2_Lo[,c("NAME","NES")], by='NAME')
hallmark_subsetsREVISED <- merge(hallmark_subsetsREVISED, hallmarkYvE_v1_Nav[,c("NAME","NES")], by='NAME')
hallmark_subsetsREVISED <- merge(hallmark_subsetsREVISED, hallmarkYvE_v2_Nav[,c("NAME","NES")], by='NAME'); names(hallmark_subsetsREVISED) <- c("NAME","Hi_v1", "Hi_v2", "Lo_v1", "Lo_v2", "Nav_v1", "Nav_v2")

hallmark_subsets <- hallmark_subsetsREVISED;  rownames(hallmark_subsets) <- hallmark_subsets$NAME;  hallmark_subsets$NAME <- NULL; rm(hallmark_subsetsREVISED)
hallmark_NFkB <- hallmark_subsets[44,]
Xpositions <- c(2.75, 3.25, 2.5, 3.5, 2.75, 3.25); Ypositions <- c(2.15, 2.15, 2, 2, 1.85, 1.85)
miniMatrix <- data.frame(NES = as.vector(t(hallmark_NFkB)), Xpos=Xpositions, Ypos = Ypositions, display=T)
miniMatrix[c("refLow","refHigh","refZero"),] <- matrix(c(-2.5,2.5,0,3,3,3,2,2,2,F,F,F),nrow=3)
b<-ggplot(data=miniMatrix, aes(x=`Xpos`, y=`Ypos`,size=abs(`NES`))) + #geom_point() + annotation_custom(grob(label="PlaceHolder", x=2,  y=3))
  geom_point(data=subset(miniMatrix, display==1),aes(color=`NES`)) + 
  scale_size(trans="exp", range=c(5,60))+
  scale_color_gradient2(low= "#67a9cf", high="#ef8a62") + 
  geom_vline(xintercept=3, linetype="dashed", color="grey", size=1) +
  theme_bw() + ggtitle(paste0(rownames(hallmark_NFkB))) + ylab(NULL) + xlab(NULL) + ylim(c(1.5,2.5)) + xlim(c(2,4)) +
  theme(plot.title = element_text(hjust = 0.5, size=36)) + 
  theme(legend.position="right", legend.text = element_text(size=14), legend.title=element_text(size=18)) + 
#  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border = element_blank(),axis.text = element_blank(), axis.ticks = element_blank()) + 
  geom_blank() # no idea why this is required but it doesn't work if I don't include it
b
# ggsave(filename = "DifferentialExpression/GSEA/Images/HexagonPlots/prototype.pdf", width=10, height=12)

savePlots <- list(); 
for (i in 1:nrow(hallmark_subsets)) 
{
  miniMatrix <- data.frame(NES = as.vector(t(hallmark_subsets[i,])), Xpos=Xpositions, Ypos = Ypositions, display=T)
  miniMatrix[c("refLow","refHigh","refZero"),] <- matrix(c(-2.5,2.5,0,3,3,3,2,2,2,F,F,F),nrow=3)
  savePlots[[i]] <- ggplot(data=miniMatrix, aes(x=`Xpos`, y=`Ypos`,color=`NES`, size=abs(`NES`))) + 
    geom_vline(xintercept=3, linetype="dashed", color="grey", size=1) +
    geom_point(data=subset(miniMatrix, display==1)) + 
    scale_size(trans="exp", range=c(5,60))+
    scale_color_gradient2(low= "#67a9cf", high="#ef8a62") + 
    theme_bw() + ggtitle(paste0(rownames(hallmark_subsets)[i])) + ylab(NULL) + xlab(NULL) + ylim(c(1.7,2.3)) + xlim(c(2,4)) +
    theme(plot.title = element_text(hjust = 0.5, size=24)) + 
    theme(legend.position="none", legend.text = element_text(size=24)) + 
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border = element_blank(),axis.text = element_blank(), axis.ticks = element_blank()) + 
    geom_blank()  # no idea why this is required but it doesn't work if I don't include it
}
names(savePlots) <- rownames(hallmark_subsets)
 
a <- ggarrange(savePlots[["Tnfa_signaling_via_nfkb"]], savePlots[["Inflammatory_response"]],savePlots[["Il2_stat5_signaling"]] ,
               savePlots[["Myc_targets_v2"]], savePlots[["Mitotic_spindle"]], savePlots[["Oxidative_phosphorylation"]],
               savePlots[["Hypoxia"]],savePlots[["Tgf_beta_signaling"]], savePlots[["Apoptosis"]],
               common.legend = TRUE, legend="right")
a
# ggexport(a, filename = "DifferentialExpression/GSEA/Images/HexagonPlots/extremesByV2Hi.pdf", width = 20, height=14)



#extract legend    https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(b)
blankPlot <- ggplot(miniMatrix) + geom_blank() + theme_bw() + theme(panel.border=element_blank())
savePlots[[50]] <- blankPlot; savePlots[[51]] <- mylegend
ml <- marrangeGrob(savePlots, nrow=13, ncol=4); ml
# ggexport(ml,  filename = "DifferentialExpression/GSEA/Images/HexagonPlots/HexagonPlots_allgeneSets.pdf", width = 20, height=70)
 
 

# hallmark_subsets <- data.frame(row.names=hallmarkYvE_v1_Hi$NAME[order(hallmarkYvE_v1_Hi$NAME)],
#                                v1Hi = hallmarkYvE_v1_Hi$NES[order(hallmarkYvE_v1_Hi$NAME)],
#                                v2Hi = hallmarkYvE_v2_Hi$NES[order(hallmarkYvE_v2_Hi$NAME)],
#                                v1Lo = hallmarkYvE_v1_Lo$NES[order(hallmarkYvE_v1_Lo$NAME)],
#                                v2Lo = hallmarkYvE_v2_Lo$NES[order(hallmarkYvE_v2_Lo$NAME)],
#                                v1Nav = hallmarkYvE_v1_Nav$NES[order(hallmarkYvE_v1_Nav$NAME)],
#                                v2Nav = hallmarkYvE_v2_Nav$NES[order(hallmarkYvE_v2_Nav$NAME)]
# )

# ICOS+CD38+ cTfh at day 7 for Young vs Elderly
YvEHiHiv2TNF<- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(hallmarkYvE_v2_Hi$NES[grep("Tnf",hallmarkYvE_v2_Hi$NAME)],2), "\n", "FDR: ", formatC(hallmarkYvE_v2_Hi$FDR.q.val[grep("Tnf",hallmarkYvE_v2_Hi$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.63,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=YvEHiHiv2TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("TNF-NFkB ICOS+CD38+ cTfh at day 7") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YvE_HiHi_d7_Hallmark_TNF-NFkB.pdf", device="pdf", height=3.5, width=5)

# ICOS-CD38- cTfh at day 7 for Young vs Elderly
YvELoLov2TNF<- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_LoLo_vs_E_LoLo_v2_Hallmark.GseaPreranked.1540914059884/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(hallmarkYvE_v2_Lo$NES[grep("Tnf",hallmarkYvE_v2_Lo$NAME)],2), "\n", "FDR: ", formatC(hallmarkYvE_v2_Lo$FDR.q.val[grep("Tnf",hallmarkYvE_v2_Lo$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.63,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=YvELoLov2TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("TNF-NFkB ICOS-CD38- cTfh at day 7") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YvE_LoLo_d7_Hallmark_TNF-NFkB.pdf", device="pdf", height=3.5, width=5)

# Naive CD4 at day 7 for Young vs Elderly
YvENaivev2TNF<- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v2_Hallmark.GseaPreranked.1540914108109/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(hallmarkYvE_v2_Nav$NES[grep("Tnf",hallmarkYvE_v2_Nav$NAME)],2), "\n", "FDR: ", formatC(hallmarkYvE_v2_Nav$FDR.q.val[grep("Tnf",hallmarkYvE_v2_Nav$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.63,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=YvENaivev2TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("TNF-NFkB Naive CD4 at day 7") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YvE_Nav_d7_Hallmark_TNF-NFkB.pdf", device="pdf", height=3.5, width=5)

# ICOS+CD38+ cTfh at day 0 for Young vs Elderly
YvEHiHiv1TNF<- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v1_HALLMARK.GseaPreranked.1531000213854/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(hallmarkYvE_v1_Hi$NES[grep("Tnf",hallmarkYvE_v1_Hi$NAME)],2), "\n", "FDR: ", formatC(hallmarkYvE_v1_Hi$FDR.q.val[grep("Tnf",hallmarkYvE_v1_Hi$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.63,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=YvEHiHiv1TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("TNF-NFkB ICOS+CD38+ cTfh at day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YvE_HiHi_d0_Hallmark_TNF-NFkB.pdf", device="pdf", height=3.5, width=5)

# ICOS-CD38- cTfh at day 0 for Young vs Elderly
YvELoLov1TNF<- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_LoLo_vs_E_LoLo_v1_Hallmark.GseaPreranked.1540914053735/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(hallmarkYvE_v1_Lo$NES[grep("Tnf",hallmarkYvE_v1_Lo$NAME)],2), "\n", "FDR: ", formatC(hallmarkYvE_v1_Lo$FDR.q.val[grep("Tnf",hallmarkYvE_v1_Lo$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.63,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=YvELoLov1TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("TNF-NFkB ICOS-CD38- cTfh at day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YvE_LoLo_d0_Hallmark_TNF-NFkB.pdf", device="pdf", height=3.5, width=5)

# Naive CD4 at day 0 for Young vs Elderly
YvENaivev1TNF<- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v1_Hallmark.GseaPreranked.1540914080384/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(hallmarkYvE_v1_Nav$NES[grep("Tnf",hallmarkYvE_v1_Nav$NAME)],2), "\n", "FDR: ", formatC(hallmarkYvE_v1_Nav$FDR.q.val[grep("Tnf",hallmarkYvE_v1_Nav$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.63,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=YvENaivev1TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("TNF-NFkB Naive CD4 at day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YvE_Nav_d0_Hallmark_TNF-NFkB.pdf", device="pdf", height=3.5, width=5)


