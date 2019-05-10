library("genefilter")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("DESeq2")
library("BiocParallel")
library(WGCNA)
register(MulticoreParam(6))
sessionInfo()


bestDataLog <- read.csv(file="bestDataLog.csv", stringsAsFactors = FALSE)




##  ************************* Heatmaps of few genes for verification ***************


# lmat = rbind(c(2,3),c(4,1)); lwid = c(0.5,4); lhei = c(1.5,4); layout(mat = lmat, widths = lwid, heights = lhei)

probeList <- c("ICOS","IRF4","TBK1","MAFF","SOCS3","SMAD2","SMAD3","SKI","IFNGR2","TGFB1","TNFSF10","ANXA4","IL2RB","HMOX1","CBL","CD36","ITGB3")
probeGenes <- bestDataLog[probeList,]

#geneMatrix <- probeGenes[,c(16:23,40:47,8:15,32:39,1:7,24:31)]   # naive LoLo HiHi in elderly 
geneMatrix <- probeGenes[,c(66:71,24:31)] 
gsg <- goodSamplesGenes(t(geneMatrix), verbose = 6)
geneMatrix <- geneMatrix[gsg[[1]],]

# png(file = "Images/Elderly_v1-v-v2_Heatmap.png", height=600,width=1500, res=100)
heatmap.2(as.matrix(t(geneMatrix)),trace="none", col=RdYlBu, cexRow=0.5, cexCol=1,key.title=NA, scale="column", 
          main="Elderly adults",key.ylab=NA, key.xlab=NA, keysize=1, density.info="none", Colv = TRUE,
          margins=c(8,5), Rowv=FALSE )# lmat = lmat, lwid = lwid, lhei = lhei)
# dev.off()

geneMatrix <- probeGenes[,c(60:65,78:83,54:59,72:77,48:53,66:71)] 
gsg <- goodSamplesGenes(t(geneMatrix), verbose = 6)
geneMatrix <- geneMatrix[gsg[[1]],]


png(file = "Images/Young_v1-v-v2_Heatmap.png", height=600,width=1500, res=100)
heatmap.2(as.matrix(t(geneMatrix)),trace="none", col=RdYlBu, cexRow=0.5, cexCol=1,key.title=NA, scale="column", 
          main="Young adults",key.ylab=NA, key.xlab=NA, keysize=1, density.info="none", Colv = FALSE,
          margins=c(8,5), Rowv=FALSE )# lmat = lmat, lwid = lwid, lhei = lhei)
dev.off()





##  TRANSCRIPTION FACTORS
probeGenes <- read.csv(file="TransFac.csv",stringsAsFactors = FALSE)
row.names(probeGenes) <- probeGenes$Gene
probeGenes <- bestDataLog[row.names(probeGenes),]

YvE_v2 <- probeGenes[,c(24:31,66:71)] #bestDataLog[row.names(bestDataLog) %in% row.names(highLow),c(66:83)]
library(WGCNA)
gsg <- goodSamplesGenes(t(YvE_v2), verbose = 6)
YvE_v2 <- YvE_v2[gsg[[1]],]

lmat = rbind(c(2,3),c(4,1))
lwid = c(0.5,4)
lhei = c(1.5,4)
layout(mat = lmat, widths = lwid, heights = lhei)

# png(file = "Images/YvE_HiHi_transFactor.png", height=7000,width=10000, res=100)
heatmap.2(as.matrix(t(YvE_v2)),trace="none", col=RdYlBu, cexRow=0.5, cexCol=0.5,key.title=NA, scale="column", 
          main="Selected Genes",key.ylab=NA, key.xlab=NA, keysize=1, density.info="none", 
          margins=c(8,5), Rowv=FALSE )# lmat = lmat, lwid = lwid, lhei = lhei)
# dev.off()





##  ********************* heatmap for SelectedGenes in Figure 3A *****************

# prove d7 vs d0 in Y and E by heatmap --  TNFa leading edge

GenesFiltered <- read.csv(file="DifferentialExpression/GSEA/Inflammation/TNFa_GeneSet_filteredList.csv",stringsAsFactors = FALSE)

GeneList <- c("F2RL1","CCL20","KYNU","PMEPA1","PPAP2B","BMP2","CXCR7","EFNA1","LDLR","IL23A","NR4A1","CSF1")
probeGenes <- bestDataLog[GeneList,]
probeSet <- probeGenes[,c(1:7,24:31,48:53,66:71)] 

heatmap.2(as.matrix(probeSet),trace="none", col=RdYlBu, cexRow=1, cexCol=1,key.title=NA, scale="row", 
          main="TNFa geneset",key.ylab=NA, key.xlab=NA, keysize=1, density.info="none", 
          margins=c(8,5), Colv=FALSE, Rowv=TRUE ) # lmat = lmat, lwid = lwid, lhei = lhei)


probeSet <- probeGenes[,c(2,25,5,28,7,31)] 
heatmap.2(as.matrix(probeSet),trace="none", col=RdYlBu, cexRow=1, cexCol=1,key.title=NA, scale="row", 
          main="TNFa geneset",key.ylab=NA, key.xlab=NA, keysize=1, density.info="none", 
          margins=c(8,5), Colv=FALSE, Rowv=FALSE ) # lmat = lmat, lwid = lwid, lhei = lhei)

# subtract d7 from d0 in elderly topThird and then d0-d7 in botThird to get summary of change
AlphabetizedBestDataLog <- bestDataLog[,c(1:7,24:29,31,48:53,66:71)]
AlphabetizedBestDataLog <- AlphabetizedBestDataLog[,order(names(AlphabetizedBestDataLog))]
diffMatrix <- data.frame("222039"=rep(0,22257),"222049"=rep(0,22257),"222060"=rep(0,22257),"222062"=rep(0,22257),"222065"=rep(0,22257), 
                         "222067"=rep(0,22257),"222068"=rep(0,22257),"111004"=rep(0,22257),"111005"=rep(0,22257),"111039"=rep(0,22257),
                         "111043"=rep(0,22257),"111049"=rep(0,22257),"111055"=rep(0,22257),stringsAsFactors = FALSE)
for(i in 1:13) { diffMatrix[,i] <- AlphabetizedBestDataLog[,i*2]-AlphabetizedBestDataLog[,i*2-1]; }
rownames(diffMatrix) <- rownames(bestDataLog)

probeTNF <- diffMatrix[c(GenesFiltered$TopThird),]
heatmap.2(as.matrix(probeTNF),trace="none", col=RdYlBu, cexRow=1, cexCol=1,key.title=NA, scale="row", 
          main="TNFa major",key.ylab=NA, key.xlab=NA, keysize=1, density.info="none", 
          margins=c(8,5), Colv=TRUE, Rowv=TRUE )


