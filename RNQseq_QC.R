library("gplots")
library("ggplot2")
library("DESeq2")
library("viridis")
sessionInfo()


dataMatrix <- read.table("RawProcessing/2018July_PORTnormalization/SPREADSHEETS/FINAL_master_list_of_gene_counts_MIN.FluVacAging.txt", sep="", header=T, stringsAsFactors = F)
rownames(dataMatrix) <- make.names(dataMatrix$geneSymbol, unique=TRUE)
dataMatrix$id <- dataMatrix$geneCoordinate <- dataMatrix$geneSymbol <- NULL 

bestData <- dataMatrix 
colSums2(as.matrix(dataMatrix))

sampleDists <- dist(t(bestData))
#sampleDists <- dist(t(bestData[nrow(bestData),4:ncol(bestData)]))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- inferno(255) 
hc <- hclust(sampleDists)
# png(filename = "Images/AllSamples_bestTranscripts_ClusterDendrogram.png", width=1500, res = 75)
  plot(hc, main="Cluster dendrogram of best transcripts")
# dev.off()
# pdf(file = "Images/AllSamples_SampleDistanceMatrix_heatmap.pdf")
heatmap.2(sampleDistMatrix, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=colors, margins=c(2,10), labCol=FALSE,cexRow = 0.4)
# dev.off()



pca1 = prcomp(t(bestData), scale = FALSE)
# create data frame with scores
scores = as.data.frame(pca1$x)

# plot of observations
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
  geom_hline(yintercept = 0, colour = "gray65") + theme_bw() +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(colour="black", alpha=1, size=3) + 
#  geom_text(colour = "tomato", alpha = 1, size = 3) +
  ggtitle("PCA plot for all samples in final pool") + theme(title = element_text(size=20))
# ggsave(filename = "Images/PCAplot_allSamples.pdf", device="pdf")






