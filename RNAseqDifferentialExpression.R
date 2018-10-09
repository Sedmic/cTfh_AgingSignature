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
register(SnowParam(workers=7))
sessionInfo()


# cust <- colorRampPalette(c("#7bb2b2","gray15","orange"))(400) ##F1a340  rev(brewer.pal(9,"YlGn")))(400)

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
                                                                          "Young ICOS-CD38+- Day 0","Young ICOS-CD38+- Day 7",
                                                                          "Young Naive Day 0","Young Naive Day 7",
                                                                          "Elderly ICOS+CD38+ Day 0","Elderly ICOS+CD38+ Day 7",
                                                                          "Elderly ICOS-CD38- Day 0","Elderly ICOS-CD38- Day 7",
                                                                          "Elderly Naive Day 0","Elderly Naive Day 7")))) +  
  xlab("") + ylab("") + 
  ggtitle("t-SNE All samples") +  theme_light(base_size=15)  +  scale_fill_manual(values=customPalette[unique(tsneReorder$Class)]) +
  theme(strip.background = element_blank()) + labs(fill = "Sample cohort") + theme(legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) #+ xlim(-60,60) + ylim(-55,70)
# ggsave(filename="Images/tSNE_all_samples.png", width=8,height=6)
# ggsave(filename="Images/tSNE_all_samples.eps", device="eps",width=8,height=6)

subsetTSNE <- tsneReorder[grep("Y_HiHi",tsneReorder$subgroup),]  ## only consider ICOS+CD38+ from young 
ggplot(subsetTSNE, aes(x=V1, y=V2)) + 
  geom_point(size=10,pch=21,colour="black",aes(fill=factor(Class, labels=c("YngD0+/+","YngD7+/+")))) + 
  xlab("") + ylab("") + # geom_text(label=subsetTSNE$subgroup) +
  ggtitle("t-SNE ICOS+CD38+ Yng") +  theme_light(base_size=20)  +  scale_fill_manual(values=customPalette[unique(subsetTSNE$Class)]) +
  theme(strip.background = element_blank()) + labs(fill = "Sample cohort") + theme(legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) + xlim(-50,50) + ylim(-55,70)
# ggsave(filename="Images/tSNE_ICOShiCD38hi_young.png", width=8,height=6)

subsetTSNE <- tsneReorder[c(grep("Y_HiHi",tsneReorder$subgroup),grep("E_HiHi",tsneReorder$subgroup)),]  ## only consider ICOS+CD38+ from elderly
ggplot(subsetTSNE, aes(x=V1, y=V2)) + 
  geom_point(size=12,pch=21,colour="black",aes(fill=factor(Class, labels=c("Young ICOS+CD38+ Day 0","Young ICOS+CD38+ Day 7",
                                                                           "Elderly ICOS+CD38+ Day 0","Elderly ICOS+CD38+ Day 7")))) + 
  xlab("") + ylab("") + # geom_text(label=subsetTSNE$subgroup) +
  ggtitle("t-SNE ICOS+CD38+ Eld") +  theme_light(base_size=20)  +  scale_fill_manual(values=customPalette[unique(subsetTSNE$Class)]) +
  theme(strip.background = element_blank()) + labs(fill = "Sample cohort") + theme(legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) + xlim(-60,60) + ylim(-55,70)
# ggsave(filename="Images/tSNE_ICOShiCD38hi_elderly.png", width=8,height=6)
# ggsave(filename="Images/tSNE_ICOShiCD38hi_elderly.eps", device="eps",width=8,height=6)

subsetTSNE <- tsneReorder[grep("LoLo",tsneReorder$subgroup),]  ## only consider ICOS-CD38- from young and elderly
ggplot(subsetTSNE, aes(x=V1, y=V2)) + 
  geom_point(size=12,pch=21,colour="black",aes(fill=factor(Class, labels=c("YngD0-/-","YngD7-/-","EldD0-/-","EldD7-/-"))))+
  xlab("") + ylab("") +
  ggtitle("t-SNE ICOS-CD38- ") +  theme_light(base_size=20)  +  scale_fill_manual(values=customPalette[unique(subsetTSNE$Class)]) +
  theme(strip.background = element_blank()) + labs(fill = "Sample cohort") + theme(legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) + xlim(-50,50) + ylim(-55,70)
# ggsave(filename="Images/tSNE_ICOSloCD38lo.png", width=8,height=6)

subsetTSNE <- tsneReorder[grep("Naive",tsneReorder$subgroup),]  ## only consider Naive from young and elderly
ggplot(subsetTSNE, aes(x=V1, y=V2)) + 
  geom_point(size=12,pch=21,colour="black",aes(fill=factor(Class, labels=c("YngD0Nav","YngD7Nav","EldD0Nav","EldD7Nav"))))+
  xlab("") + ylab("") +
  ggtitle("t-SNE Naive") +  theme_light(base_size=20)  +  scale_fill_manual(values=customPalette[unique(subsetTSNE$Class)]) +
  theme(strip.background = element_blank()) + labs(fill = "Sample cohort") + theme(legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) + xlim(-50,50) + ylim(-55,70)
# ggsave(filename="Images/tSNE_Naive.png", width=8,height=6)

subsetTSNE <- tsneReorder[c(grep("Y_HiHi_v1",tsneReorder$subgroup), grep("Y_LoLo_v1",tsneReorder$subgroup), grep("Y_Naive_v1",tsneReorder$subgroup)),]  ## only consider young samples day 0
ggplot(subsetTSNE, aes(x=V1, y=V2)) +    
  geom_point(size=12, pch=21, colour="black", aes(fill=factor(Class, labels=c("ICOShiCD38hi cTfh","ICOSloCD38lo cTfh","Naive CD4")))) + 
  xlab("") + ylab("") + 
  ggtitle("t-SNE Young Pre-Post") +  theme_light(base_size=20)  +  scale_fill_manual(values=customPalette[unique(subsetTSNE$Class)]) +
  theme(strip.background = element_blank()) + labs(fill = "Subset") + theme(legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) + xlim(-50,50) + ylim(-55,70)
# ggsave(filename="Images/tSNE_Yng_CD4subsets.png",width=8.5,height=6)


subsetTSNE <- tsneReorder[c(grep("Y_HiHi_",tsneReorder$subgroup), grep("Y_LoLo_",tsneReorder$subgroup), grep("Y_Naive_",tsneReorder$subgroup)),]  ## only consider young samples day 0
ggplot(subsetTSNE, aes(x=V1, y=V2)) +    
  geom_point(size=12, pch=21, colour="black", aes(fill=factor(Class, labels=c("ICOShiCD38hi D0","ICOShiCD38hi D7","ICOSloCD38lo D0","ICOSloCD38lo D7", "Naive CD4 D0", "Naive CD4 D7")))) + 
  xlab("") + ylab("") + 
  ggtitle("t-SNE Young Pre-Post") +  theme_light(base_size=20)  +  scale_fill_manual(values=customPalette[unique(subsetTSNE$Class)]) +
  theme(strip.background = element_blank()) + labs(fill = "Subset") + theme(legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) + xlim(-50,50) + ylim(-55,70)
# ggsave(filename="Images/tSNE_Yng_prepost.png",width=8,height=6)





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
         gaps_col = c(6,12), annotation_colors = ann_colors, cutree_rows = 3, fontsize_row = 16, color=viridis(100)
         , filename = "Images/SelectedGenesHeatmap.pdf"
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
         gaps_col = c(7,15), annotation_colors = ann_colors, cutree_rows = 3, fontsize_row = 16, color=viridis(100)
        , filename = "Images/SelectedGenesHeatmap_ELDERLY.pdf"
)

# dev.off(); dev.off(); 


## ******************************** Differential expression *********************************************
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

for (i in 1:length(diffExpr)) { write.csv(diffExpr[[i]], file=paste0("DifferentialExpression/",names(diffExpr[i]),".csv")) }


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
  geom_point(alpha=0.2, size=5, colour="black") + theme_bw() +  theme(legend.position="none") + 
#  geom_text(aes(label=row.names(data)), cex=4) +  
  geom_text_repel(data=subset(data, padj<0.001), aes(label=row.names(subset(data, padj<0.001)))) + 
  labs(title="HiHi vs LoLo visit 1") + xlab("log2 fold change") + ylab("-log10 p-adj") + ylim(0,13) + xlim(-11,11)
# ggsave(filename="DifferentialExpression/Images/Young_HiHi_vs_LoLo_v1.png", device="png")
plotMA(diffExpr_YoungOnly_controlSubj[["Y_HiHivsLoLo_v1"]],ylim=c(-11,11))

data <- as.data.frame(diffExpr_ElderlyOnly_controlSubj[["E_HiHivsLoLo_v1"]])
data <- data[which(data$lfcSE<2),]     # filter out very high SE which may be outliers
ggplot(data, aes(x=log2FoldChange, y=-log10(padj), label=row.names(data))) +
  geom_point(alpha=0.2, size=5, colour="black") + theme_bw() +  theme(legend.position="none") + 
  #  geom_text(aes(label=row.names(data)), cex=4) +  
  geom_text_repel(data=subset(data, padj<0.00005), aes(label=row.names(subset(data, padj<0.00005)))) + 
  labs(title="HiHi vs LoLo visit 1 - ELDERLY") + xlab("log2 fold change") + ylab("-log10 p-adj") + ylim(0,13) + xlim(-11,11)
# ggsave(filename="DifferentialExpression/Images/Young_HiHi_vs_LoLo_v1_ELDERLY.png", device="png")
plotMA(diffExpr_ElderlyOnly_controlSubj[["E_HiHivsLoLo_v1"]],ylim=c(-11,11))




# now look at YOUNG day 7 vs day 0
data <- as.data.frame(diffExpr_YoungOnly_controlSubj[["Y_HiHi_v12"]])
ggplot(data, aes(x=log2FoldChange, y=-log10(padj), label=row.names(data))) +
  geom_point(alpha=0.2, size=5, colour="black") + theme_bw() +  theme(legend.position="none") + 
  #  geom_text(aes(label=row.names(data)), cex=4) +  
  geom_text_repel(data=subset(data, padj<0.01), aes(label=row.names(subset(data, padj<0.01)))) + 
  labs(title="Young ICOS+CD38+ cTfh:  day 7 vs day 0") + xlab("log2 fold change") + ylab("-log10 p-adj") + ylim(0,5) + xlim(-11,11)
# ggsave(filename="DifferentialExpression/Images/Young_HiHi_v1-v2.png", device="png")
plotMA(diffExpr_YoungOnly_controlSubj[["Y_HiHi_v12"]],ylim=c(-11,11))

probeList <- as.data.frame(head(data[order(data$stat,decreasing=T),],20)); probeList <- rbind(probeList, head(data[order(data$stat,decreasing=F),],20))
probeList <- probeList[which(probeList$lfcSE < 2),]
probeGenes <- bestDataLog[rownames(probeList),grep("HiHi",colnames(bestDataLog[1:36]))]
probeGenes <- cbind(probeGenes[,grep("v1",colnames(probeGenes))], probeGenes[,grep("v2",colnames(probeGenes))])

annotateHeatmap <- data.frame(row.names = colnames(probeGenes), ageGroup = c(rep("day 0", 6),rep("day 7", 6)))
ann_colors = list(  ageGroup = c("day 0" ="orange", "day 7" = "orange4")  )
pheatmap(probeGenes, scale="row", cluster_col=F, annotation_col = annotateHeatmap, show_colnames=F, main="Differentially expressed - Young",
         gaps_col = c(6), annotation_colors = ann_colors, cutree_rows = 2, fontsize_row = 16, color=viridis(100)
        #, filename = "DifferentialExpression/Images/Y_hihi_V1vsV2_SelectedGenesHeatmap.png"
)

# now look at ELDERLY day 7 vs day 0
data <- as.data.frame(diffExpr_ElderlyOnly_controlSubj[["E_HiHi_v12"]])
data <- data[which(data$lfcSE<2),]     # filter out very high SE which may be outliers
ggplot(data, aes(x=log2FoldChange, y=-log10(padj), label=row.names(data))) +
  geom_point(alpha=0.2, size=5, colour="black") + theme_bw() +  theme(legend.position="none") + 
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

annotateHeatmap <- data.frame(row.names = colnames(probeGenes), ageGroup = c(rep("day 0", 7),rep("day 7", 8)))
ann_colors = list(  ageGroup = c("day 0" ="mediumpurple1", "day 7" = "purple4")  )
pheatmap(probeGenes, scale="row", cluster_col=F, annotation_col = annotateHeatmap, show_colnames=F, main="Differentially expressed - Elderly",
         gaps_col = c(7), annotation_colors = ann_colors, cutree_rows = 2, fontsize_row = 12, color=viridis(100)
         #, filename = "DifferentialExpression/Images/E_hihi_V1vsV2_SelectedGenesHeatmap.png"
)

# now compare YOUNG and ELDERLY directly for ICOS+CD38+ cTfh at day 7
data <- as.data.frame(diffExpr[["YvE_HiHi_v2"]])
data <- data[which(data$lfcSE<1.5),]     # filter out very high SE which may be outliers
ggplot(data, aes(x=log2FoldChange, y=-log10(padj), label=row.names(data))) +
  geom_point(alpha=0.2, size=5, colour="black") + theme_bw() +  theme(legend.position="none") + 
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
         gaps_col = c(6), annotation_colors = ann_colors, cutree_rows = 2, fontsize_row = 16, color=viridis(100)
         #, filename = "Images/YvE_hihi_v2_SelectedGenesHeatmap.png"
)

# dev.off(); dev.off(); dev.off(); dev.off();


## *****************************       GSEA results: Hallmark Genesets    ***********************************************
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

YexternalGenesetsPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/gsea_report_for_na_pos_1531001576554.xls", sep="\t")
YexternalGenesetsNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/gsea_report_for_na_neg_1531001576554.xls", sep="\t")
YexternalGenesets <- rbind(YexternalGenesetsPos, YexternalGenesetsNeg)
EexternalGenesetsPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_ExtraGeneSets.GseaPreranked.1531837454581/gsea_report_for_na_pos_1531837454581.xls", sep="\t")
EexternalGenesetsNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_ExtraGeneSets.GseaPreranked.1531837454581/gsea_report_for_na_neg_1531837454581.xls", sep="\t")
EexternalGenesets <- rbind(EexternalGenesetsPos, EexternalGenesetsNeg)

LocciGSE50392 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/GCTFH-VS-NAV.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(YexternalGenesets$NES[grep("GCTFH-VS-NAV",YexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(YexternalGenesets$FDR.q.val[grep("GCTFH-VS-NAV",YexternalGenesets$NAME)], format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=LocciGSE50392, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="orange", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE50392 HUMAN Tonsil GC-Tfh vs Naive") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_LocciGCTfh.png", device="png", height=3.5, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_LocciGCTfh.pdf", device="pdf", height=3.5, width=5, dpi=300)

LocciGSE50392 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_ExtraGeneSets.GseaPreranked.1531837454581/GSE50392_GCTFH-VS-NAV.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(EexternalGenesets$NES[grep("GCTFH-VS-NAV",EexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(EexternalGenesets$FDR.q.val[grep("GCTFH-VS-NAV",EexternalGenesets$NAME)], format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=LocciGSE50392, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="purple", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE50392 HUMAN Tonsil GC-Tfh vs Naive") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_LocciGCTfh.png", device="png", height=3.5, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_LocciGCTfh.pdf", device="pdf", height=3.5, width=5, dpi=300)

MOUSEGCTFHGSE16697 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/GSE16697_MOUSEGCTFH.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(YexternalGenesets$NES[grep("16697",YexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(YexternalGenesets$FDR.q.val[grep("16697",YexternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=MOUSEGCTFHGSE16697, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="orange", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE16697 MOUSE Tfh vs non-Tfh CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_MouseGCTfh.png", device="png", height=3.5, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_MouseGCTfh.pdf", device="pdf", height=3.5, width=5, dpi=300)

MOUSEGCTFHGSE16697 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_ExtraGeneSets.GseaPreranked.1531837454581/GSE16697_MOUSEGCTFH.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(EexternalGenesets$NES[grep("16697",EexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(EexternalGenesets$FDR.q.val[grep("16697",EexternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=MOUSEGCTFHGSE16697, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="purple", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE16697 MOUSE Tfh vs non-Tfh CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_MouseGCTfh.png", device="png", height=3.5, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_MouseGCTfh.pdf", device="pdf", height=3.5, width=5, dpi=300)

MOUSEGCTFHGSE32596 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/GSE32596_MOUSEGCTFH.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(YexternalGenesets$NES[grep("32596",YexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(YexternalGenesets$FDR.q.val[grep("32596",YexternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=MOUSEGCTFHGSE32596, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="orange", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE32596 MOUSE Tfh vs Naive CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_MouseGCTfh2.png", device="png", height=3.5, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_MouseGCTfh2.pdf", device="pdf", height=3.5, width=5, dpi=300)

MOUSEGCTFHGSE32596 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_ExtraGeneSets.GseaPreranked.1531837454581/GSE32596_MOUSEGCTFH.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(EexternalGenesets$NES[grep("32596",EexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(EexternalGenesets$FDR.q.val[grep("32596",EexternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=MOUSEGCTFHGSE32596, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="purple", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE32596 MOUSE Tfh vs Naive CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_MouseGCTfh2.png", device="png", height=3.5, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/Ehihi-v-lolo_v1_MouseGCTfh2.pdf", device="pdf", height=3.5, width=5, dpi=300)



TCF1KOGSE65693 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_ExternalGeneSets.GseaPreranked.1531001576554/GSE65693_TCF1KO.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(YexternalGenesets$NES[grep("65693",YexternalGenesets$NAME)],2), "\n", "FDR: ", formatC(YexternalGenesets$FDR.q.val[grep("65693",YexternalGenesets$NAME)], format = "e", digits = 1))
my_grob = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=TCF1KOGSE65693, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="orange", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE65693 TCF1 knockout") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/Yhihi-v-lolo_v1_TCF1ko.png", device="png", height=3.5, width=5, dpi=300)


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
  geom_line(color="orange3", size=1) + 
  geom_line(data=ElderlyHiHiv2TGF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="purple")+
  geom_rug(sides="t", size=0.75, alpha=0.5, color="orange3") +  geom_rug(data=ElderlyHiHiv2TGF, sides="b", size=0.75, alpha=0.5, color="purple") +  
  ggtitle("TGF-b Signaling - Day 7 vs Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_TGFbeta.png", device="png", height=3.5, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_TGFbeta.pdf", device="pdf", height=3.5, width=5, dpi=300)

# E2F Targets
YoungHiHiv2E2F <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/HALLMARK_E2F_TARGETS.xls", sep="\t")
ElderlyHiHiv2E2F <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/HALLMARK_E2F_TARGETS.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(Yd0tod7Hallmark$NES[grep("E2F_TARGETS",Yd0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Yd0tod7Hallmark$FDR.q.val[grep("E2F_TARGETS",Yd0tod7Hallmark$NAME)], format = "e", digits = 1))
annotationInfo2 <- paste0("\nNES: ", round(Ed0tod7Hallmark$NES[grep("E2F_TARGETS",Ed0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Ed0tod7Hallmark$FDR.q.val[grep("E2F_TARGETS",Ed0tod7Hallmark$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.88, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.7,  y=0.75, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=YoungHiHiv2E2F, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="orange3", size=1) + 
  geom_line(data=ElderlyHiHiv2E2F, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="purple")+
  geom_rug(sides="t", size=0.75, alpha=0.5, color="orange3") +  geom_rug(data=ElderlyHiHiv2E2F, sides="b", size=0.75, alpha=0.5, color="purple") +  
  ggtitle("E2F Targets - Day 7 vs Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_E2FTargets.png", device="png", height=3.5, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_E2FTargets.pdf", device="pdf", height=3.5, width=5, dpi=300)

# IL2 signaling
YoungHiHiv2IL2 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/HALLMARK_IL2_STAT5_SIGNALING.xls", sep="\t")
ElderlyHiHiv2IL2 <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/HALLMARK_IL2_STAT5_SIGNALING.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(Yd0tod7Hallmark$NES[grep("IL2",Yd0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Yd0tod7Hallmark$FDR.q.val[grep("IL2",Yd0tod7Hallmark$NAME)], format = "e", digits = 1))
annotationInfo2 <- paste0("\nNES: ", round(Ed0tod7Hallmark$NES[grep("IL2",Ed0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Ed0tod7Hallmark$FDR.q.val[grep("IL2",Ed0tod7Hallmark$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.88, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.7,  y=0.75, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=YoungHiHiv2IL2, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="orange3", size=1) + 
  geom_line(data=ElderlyHiHiv2IL2, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="purple")+
  geom_rug(sides="b", size=0.75, alpha=0.5, color="orange3") +  geom_rug(data=ElderlyHiHiv2IL2, sides="t", size=0.75, alpha=0.5, color="purple") +  
  ggtitle("IL2 signaling - Day 7 vs Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_IL2-STAT5.png", device="png", height=3.5, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_IL2-STAT5.pdf", device="pdf", height=3.5, width=5, dpi=300)

# TNF - NFkb
YoungHiHiv2TNF <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926605818/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
ElderlyHiHiv2TNF <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v1_vs_v2_HALLMARK.GseaPreranked.1530926588879/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(Yd0tod7Hallmark$NES[grep("TNF",Yd0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Yd0tod7Hallmark$FDR.q.val[grep("TNF",Yd0tod7Hallmark$NAME)], format = "e", digits = 1))
annotationInfo2 <- paste0("\nNES: ", round(Ed0tod7Hallmark$NES[grep("TNF",Ed0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Ed0tod7Hallmark$FDR.q.val[grep("TNF",Ed0tod7Hallmark$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.88, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.7,  y=0.75, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=YoungHiHiv2TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="orange3", size=1) + 
  geom_line(data=ElderlyHiHiv2TNF, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="purple")+
  geom_rug(sides="b", size=0.75, alpha=0.5, color="orange3") +  geom_rug(data=ElderlyHiHiv2TNF, sides="t", size=0.75, alpha=0.5, color="purple") +  
  ggtitle("TNF signaling - Day 7 vs Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_TNF-NFkb.png", device="png", height=3.5, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_d7vsd0_TNF-NFkb.pdf", device="pdf", height=3.5, width=5, dpi=300)

## *****************************       GSEA results: YvE direct comparison - Hallmark    ***********************************************


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
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_v2_hallmark_full.png", device="png", height=7, width=5, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_v2_hallmark_full.pdf", device="pdf", height=7, width=5, dpi=300)


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
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_v2_hallmark_limited.png", device="png", height=10, width=7, dpi=300)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi_v2_hallmark_limited.pdf", device="pdf", height=10, width=7, dpi=300)

# what is the most differentially altered between the two groups? need to do a subtraction

ggplot(data=hallmarkMergedWide, aes(x=reorder(`NAME`,-`EminusY`), y=`EminusY`)) + geom_point(size=3) + 
  geom_bar( stat="Identity", width=0.2, color="blue") +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets \nICOS+CD38+ cTfh for Day 7 vs Day 0") + ylab("E-Y Difference in NES") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=16)) + geom_hline(aes(yintercept = 0))
# ggsave(file="DifferentialExpression/GSEA/Images/EminusY_hihi_v2_hallmark_full.png", device="png", height=7, width=5, dpi=300)


ggplot(data=hallmarkFactored, aes(x=reorder(`NAME`,-`EminusY`), y=`EminusY`)) + geom_point(size=3) + 
  geom_bar( stat="Identity", width=0.2, color="blue") +
  coord_flip() + theme_bw() + ggtitle("HALLMARK Genesets \nICOS+CD38+ cTfh for Day 7 vs Day 0") + ylab("E-Y Difference in NES") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=16)) + geom_hline(aes(yintercept = 0))
# ggsave(file="DifferentialExpression/GSEA/Images/EminusY_hihi_v2_hallmark_limited.png", device="png", height=7, width=5, dpi=300)


