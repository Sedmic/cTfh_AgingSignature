library(WGCNA)
library(DESeq2)

#  gene filtering approach
dataMatrixCounts <- read.table("RawProcessing/2018July_PORTnormalization/SPREADSHEETS/FINAL_master_list_of_gene_counts_MIN.FluVacAging.txt", sep="", header=T, stringsAsFactors = F)
rownames(dataMatrixCounts) <- make.names(dataMatrixCounts$geneSymbol, unique=TRUE)
dataMatrixCounts$id <- dataMatrixCounts$geneCoordinate <- dataMatrixCounts$geneSymbol <- NULL 
dataMatrixCounts$X222006_HiHi_v1 <- NULL  # based on what was learned from the RNAseq_QC script
bestDataCounts <- dataMatrixCounts
metaData <- data.frame(row.names=colnames(bestDataCounts));  metaData$condition <- "empty"
#  1:36 are young, then 37:83 are elderly
metaData$ageGroup <- c(rep("Y",36),rep("E",47))
metaData$subject <- substr(rownames(metaData), start=2,stop=7)
metaData$condition <- substr(rownames(metaData),start=9, stop=20)
metaData$subgroup <- paste0(metaData$ageGroup,sep="_", metaData$condition)

bestDataLog <- read.csv(file = "bestDataLog.csv", stringsAsFactors = F); row.names(bestDataLog) <- bestDataLog$X; bestDataLog$X <- NULL
fullDataset <- DESeqDataSetFromMatrix(countData=bestDataCounts, colData=metaData, design= ~ subgroup)
dds <- estimateSizeFactors(fullDataset)
temp <- counts(dds,normalized=TRUE)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 10  # filter for genes with at least 5 counts in 10% of all samples
bestDataLog <- bestDataLog[idx,]


NetworkData <- data.frame(bestDataLog[,grep("HiHi_v2", colnames(bestDataLog))])  
gsg <- goodSamplesGenes(t(NetworkData), verbose = 6)
NetworkData <- NetworkData[gsg[[1]],]  # only keep genes with variability and nonzero reads
NetworkData <- NetworkData[,grep("X222", colnames(NetworkData))]   #elderly v1-v2 hihi 


sampleTree <- hclust(dist(t(NetworkData)),method="average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 8000, col = "red");
clust <- cutreeStatic(sampleTree,cutHeight=300,minSize=1)
table(clust)
keepSamples <- (clust==1)
DataMatrix <- t(NetworkData[,keepSamples])
nGenes <- ncol(DataMatrix)
nSamples <- nrow(DataMatrix)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(DataMatrix, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



#  *****   MANUAL METHOD *******   
 
adjacency <- adjacency(DataMatrix,power=20)
TOM = TOMsimilarity(adjacency); 
# save.image(file="Network\2018912_E_hihi_v2.Rdata")
dissTOM <- 1-TOM; rm(TOM)
geneTree <- hclust(as.dist(dissTOM), method = "average");

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);

#Identification of modules from the dendogram
minModuleSize = 100;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = minModuleSize);

dynamicColors = labels2colors(dynamicMods)
collectGarbage()
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, 
                    guideHang = 0.05,main = "Gene dendrogram and module colors")


#Identification of module eigen genes
MEList = moduleEigengenes(DataMatrix, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average");

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
abline(h=0, col = "red")

#Merging of modules
MEDissThres = 0.0
merge <- mergeCloseModules(DataMatrix, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors;
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)

#pdf(file="GeneDendroGram_Modules_merged_Ehihiv2.pdf")
plotDendroAndColors(geneTree, mergedColors, main="Modules after merging", c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
#dev.off()


moduleColors = mergedColors
colorOrder = c(standardColors(50), "grey");
moduleLabels = match(moduleColors, colorOrder)-1;

sort(table(dynamicColors))

#Reading clinical trait data and matching rownames with expression data
traitData<-read.csv("Network/ClinicalTraits.csv")
Samples <- rownames(DataMatrix)
traitRows <- match(Samples, traitData$Subject)
datTraits <- traitData[traitRows, -1]
rownames(datTraits) <- traitData[traitRows, 1]

#Creating colors for various traits
traitColors <- labels2colors(datTraits, zeroIsGrey = TRUE, colorSeq = NULL, naColor = "grey", commonColorCode = TRUE)

#Plotting dendogram with the trait heatmap
# plotDendroAndColors(sampleTree, traitColors,groupLabels = names(datTraits),main = "Sample dendrogram and trait heatmap")


#Quantifying module trait associations
MEs = mergedMEs
MEs0 = moduleEigengenes(DataMatrix, moduleColors)$eigengenes
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# modules to assay:  lightgreen, green, brown, pink, black, turquoise, yellow, blue
toMatch <- c("MElightgreen", "MEgreen$", "MEbrown","MEpink","MEblack","MEturquoise","MEyellow$","MEblue$")
keyModules <- grep(paste(toMatch,collapse="|"), rownames(moduleTraitCor), value=F)
moduleTraitCor <- moduleTraitCor[keyModules,]; moduleTraitPvalue <- moduleTraitPvalue[keyModules,]
MEs <- MEs[,keyModules]
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


#Plot for new module trait relationship

sizeGrWindow(10,6)
par(mar = c(8, 12, 3, 3))
# pdf(file = "Network/Eld_HiHiV2/Module-trait_relationships.pdf")
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), 
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, cex.lab.y=0.7, 
               yColorWidth=0.02, cex.lab.x=0.7,
               setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
# dev.off()


#Plotting TOM
plotTOM <- (dissTOM)^7
diag(plotTOM) = NA;
sizeGrWindow(9,9)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

#Calculation of gene relationships to traits,gene significance and module membership
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(DataMatrix, MEs, use = "p"))  #same as signedKME function
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Relate module to trait
trait <- as.data.frame(datTraits$ICOSCD38_freqTfh_V2)
names(trait)="ICOSCD38_freqTfh_V2"
geneTraitSignificance = as.data.frame(cor(DataMatrix, trait, use = "p"))   #calculate p values for MM
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))  #calculate p values for GS

names(geneTraitSignificance) = paste("GS.", names(trait), sep="") #pre-pend label for Y column
names(GSPvalue) = paste("p.GS.", names(trait), sep="")   #pre-pend label for Y column

module="turquoise"
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),  abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for blue",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'black')#module)

MET = orderMEs(cbind(MEs, trait))

#sizeGrWindow(5,7.5);
#par(cex = 0.9)
#plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 180)


##  ------------------------------- Investigate per-module gene list -------------------------------------------

DataMatrixGenes <- colnames(DataMatrix)
turquoise <- DataMatrixGenes[moduleColors=="turquoise"]
ranked <- geneModuleMembership[order(geneModuleMembership$MMturquoise, decreasing=TRUE),]
ranked <- ranked[,grep(paste(c("MMgreen","MMlightgreen","MMbrown","MMturquoise","MMyellow","MMblue", "MMpink","MMblack"), collapse="$|^"), colnames(ranked), value=T)]
# write.csv(ranked,"Network/Eld_HiHiV2/20180912_RankedGenes_moduleMembership.csv")
# save(geneTree, merge, MEList, METree, file="Network/Eld_HiHiV2/20180912_diffModuleComparison_E_hihi_v2.Rdata")
# save.image(file="Network/Eld_HiHiV2/20180913_FinalModules.Rdata")


FilterGenes <- abs(geneTraitSignificance)> .9 & abs(geneModuleMembership$MMblue)>.98
FilterGenes <- subset(FilterGenes,FilterGenes[,1] == TRUE)


#Gene annotation and output of genes with statistics 
#annot = read.csv(file = "GeneAnnotation.csv");
#dim(annot)
#names(annot)
#probes = names(datExpr0)
#probes2annot = match(probes, annot$SYMBOL)

geneInfo0 = data.frame(SYMBOL = DataMatrixGenes, moduleColor = moduleColors, geneTraitSignificance, GSPvalue)


geneOrder = order(geneInfo0$moduleColor)#, -abs(geneInfo0$GS.treatment))
geneInfo = geneInfo0[geneOrder, ]

#write.csv(geneInfo, file = "geneInfo.csv")



# --------------------------------------------




# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);









save(DataMatrix, NetworkData, TOM, net, sft, modProbes, file = "wgcna_20150725.RData");
save(DataMatrix, TOM, sft, sampleTree, file = "wgcna_20150818.RData")
lnames <- load(file="wgcna_20150818.RData")

save(geneTree, merge, MEList, METree, file = "wgcna_20150915.RData")

