library("grid")
library("genefilter")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("DESeq2")
library("GSVA")
library("GSEABase")
library("pheatmap")
library("viridis")
library("tools")
library("ggcorrplot")
library("ggrepel")
sessionInfo()

bestDataLog <- read.csv(file="bestDataLog.csv",stringsAsFactors = FALSE); rownames(bestDataLog) <- bestDataLog$X; bestDataLog$X <- NULL
phenotypeMatrix <- read.csv(file="../../Analysis/BAA_yr4_AllMergedData_FC.csv"); phenotypeMatrix <- phenotypeMatrix[-c(62:73),]; rownames(phenotypeMatrix) <- phenotypeMatrix$Subject; 
phenotypeMatrix$Subject <- NULL
gsets <- getGmt("DifferentialExpression/GSVA/h.all.v6.1.symbols.gmt")
HiHiv2 <- as.matrix(bestDataLog[,grep("HiHi_v2",colnames(bestDataLog))])
GSVAhallmark <- gsva(HiHiv2, gsets, method="gsva")
# write.csv(GSVAhallmark, file="DifferentialExpression/GSVA/GSVA_HiHi_v2_hallmark.csv")

colnames(GSVAhallmark) <- substr(colnames(GSVAhallmark),start=2,stop=7)

rownames(GSVAhallmark) <- toTitleCase(tolower(substr(rownames(GSVAhallmark),start =10, stop=50)))

annotateHeatmap <- data.frame(row.names = colnames(GSVAhallmark), ageGroup = c(rep("Young", 6),rep("Elderly", 8)))
ann_colors = list(  ageGroup = c("Young" ="orange", "Elderly" = "purple")  )
pheatmap(GSVAhallmark, scale="none", cluster_col=T, annotation_col = annotateHeatmap, show_colnames=F, main="GSVA scores - Hallmark genesets",
         annotation_colors = ann_colors, cutree_cols = 3, cutree_rows=2, fontsize_row = 8, color=inferno(100)
#         , filename = "DifferentialExpression/GSVA/Images/GSVAhallmark_allGeneSets.pdf"
)


# ****************************************************** What genes are driving the clustering? ***************************************
genesHallmark <- read.csv(file = "DifferentialExpression/GSVA/h.all.v6.1.symbols.gmt", sep="\t", stringsAsFactors = F, header = F); rownames(genesHallmark) <- genesHallmark$V1
genesHallmark$V2 <- genesHallmark$V1 <- NULL

# now take away gene sets that don't contribute to the clustering
genesHallmark <- genesHallmark[-c(grep("MYC_TARGETS_V1", rownames(genesHallmark)), grep("OXIDATIVE",  rownames(genesHallmark)), 
                                  grep("DNA_REPAIR",rownames(genesHallmark)), grep("ADIPOGENESIS",  rownames(genesHallmark)), 
                                  grep("FATTY_ACID", rownames(genesHallmark)), grep("BILE_ACID",  rownames(genesHallmark)),
                                  grep("PEROXISOME",  rownames(genesHallmark)), grep("KRAS_SIGNALING",  rownames(genesHallmark)),
                                  grep("PANCREAS", rownames(genesHallmark))
                                  ),
                               ]

genesHallmark <- unique(unlist(genesHallmark))  # now a linear concatenation of all of the genes in the Hallmark dataset
result <- data.frame(table(unlist(genesHallmark)))
result <- result[order(result$Freq, decreasing = T),]

filteredBestDataLog <- read.csv(file="FilteredBestDataLog.csv",stringsAsFactors = F); rownames(filteredBestDataLog) <- filteredBestDataLog$X; filteredBestDataLog$X <- NULL
probeList <- as.character(result$Var1[2:100]); probeGenes <- filteredBestDataLog[probeList,grep("HiHi_v2",colnames(filteredBestDataLog))]; probeGenes <- na.omit(probeGenes)

# probeGenes <- probeGenes
annotateHeatmap <- data.frame(row.names = colnames(probeGenes), ageGroup = c(rep("Young", 6),rep("Elderly", 8)))
ann_colors = list(  ageGroup = c("Young" ="orange", "Elderly" = "purple")  )
pheatmap(probeGenes, scale="row", cluster_col=T, annotation_col = annotateHeatmap, show_colnames=T, main="Selected genes",
         annotation_colors = ann_colors #, cutree_rows=2, # gaps_col = c(0), fontsize_row = 16
         #, filename = "Images/YvE_hihi_v2_SelectedGenesHeatmap.png"
)
# this approach was not fruitful, did not cluster the subjects in the same way as the GSVA on all genesets


# ********************************************* correlate GSVA against serum and flow cytometry ***********************************************

#  external work: paste in GSVA scores for TNF into spreadsheet

phenotypeMatrix <- read.csv(file="../../Analysis/BAA_yr4_AllMergedData_FC.csv"); 
phenotypeMatrix <- phenotypeMatrix[-c(62:73),]; rownames(phenotypeMatrix) <- phenotypeMatrix$Subject; # not sure why there are NA rows introduced into the end
phenotypeMatrix$Subject <- NULL
phenotypeMatrixYoung <- subset(phenotypeMatrix,Identifier=="Young",stat="identity")
phenotypeMatrixElderly <- subset(phenotypeMatrix,Identifier=="Elderly",stat="identity")


#  TNF vs GSVA-TNFNFkB score
ggplot(data=phenotypeMatrix, aes(x=`TNFa`,y=`GSVAscoreTNF`, fill="black")) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1) + 
  geom_point(size=6, pch=21) +  ylim(-1,1) + xlim(0,30)+
  theme(axis.text = element_text(size=24,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))+
  ggtitle("GSVA for TNF-NFkB vs serum TNF") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("GSVA score for \nTNFa-NFkB geneset")  + xlab("TNF (pg/mL)")+
  scale_fill_manual(values=c('black','#E69F00')) + scale_color_manual(values=c('black', '#E69F00'))
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_serumTNF_all.pdf", device="pdf", width=8, height=6)
cor.test(phenotypeMatrix$TNFa, phenotypeMatrix$GSVAscoreTNF, use="complete")


# annotationInfo <- paste0("NES: ", round(Yd0tod7Hallmark$NES[grep("IL2",Yd0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Yd0tod7Hallmark$FDR.q.val[grep("IL2",Yd0tod7Hallmark$NAME)], format = "e", digits = 1))
# annotationInfo2 <- paste0("\nNES: ", round(Ed0tod7Hallmark$NES[grep("IL2",Ed0tod7Hallmark$NAME)],2), "\n", "FDR: ", formatC(Ed0tod7Hallmark$FDR.q.val[grep("IL2",Ed0tod7Hallmark$NAME)], format = "e", digits = 1))
# my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.85, hjust=0, gp=gpar(col="orange3", fontsize=12)))
# my_grob2 = grobTree(textGrob(annotationInfo2, x=0.7,  y=0.75, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`TNFa`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  ylim(-1,1) + xlim(0,30)+
  ggtitle("GSVA for TNFa-NFkB vs serum TNF") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA score for TNFa-NFkB geneset")  + xlab("TNF (pg/mL)")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_serumTNF_byAge.png", device="png")
cor.test(phenotypeMatrixYoung$TNFa, phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
cor.test(phenotypeMatrixElderly$TNFa, phenotypeMatrixElderly$GSVAscoreTNF, use="complete")


# ***************************************  Comparison to the B cell response ***********************************************
#  ASC vs GSVA-TNFNFkB score 
grep("ASC",colnames(phenotypeMatrix), value=T)
a <- cor.test(phenotypeMatrixYoung$ASC.freqLive, phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
b <- cor.test(phenotypeMatrixElderly$ASC.freqLive, phenotypeMatrixElderly$GSVAscoreTNF, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), "\n", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`ASC.freqLive`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  ylim(-1,1) + xlim(0,3)+
  ggtitle("GSVA TNF-NFkB vs ASC freq FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA for TNFa-NFkB geneset")  + xlab("ASC freq foldchange")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_ASC_byAge.png", device="png")



#  H3N2 nAb vs GSVA-TNFNFkB score 
grep("H3N2",colnames(phenotypeMatrix), value=T)
a <- cor.test(phenotypeMatrixYoung$H3N2.nAb.FCd28, phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
b <- cor.test(phenotypeMatrixElderly$H3N2.nAb.FCd28, phenotypeMatrixElderly$GSVAscoreTNF, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), "\n", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`H3N2.nAb.FCd28`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  ylim(-1,1) + xlim(0,20)+
  ggtitle("GSVA TNF-NFkB vs H3N2.nAb.FCd28") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA for TNFa-NFkB geneset")  + xlab("H3N2.nAb.FCd28")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_H3N2nAb_byAge.png", device="png")

#  H1N1 nAb vs GSVA-TNFNFkB score 
grep("H1N1",colnames(phenotypeMatrix), value=T)
a <- cor.test(phenotypeMatrixYoung$H1N1.nAb.FCd28, phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
b <- cor.test(phenotypeMatrixElderly$H1N1.nAb.FCd28, phenotypeMatrixElderly$GSVAscoreTNF, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), "\n", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`H1N1.nAb.FCd28`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  ylim(-1,1) + xlim(0,10)+
  ggtitle("GSVA TNF-NFkB vs H1N1.nAb.FCd28") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA for TNFa-NFkB geneset")  + xlab("H1N1.nAb.FCd28")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_H1N1nAb_byAge.png", device="png")



#  H3N2 IgG vs GSVA-TNFNFkB score 
grep("H3N2",colnames(phenotypeMatrix), value=T)
a <- cor.test(phenotypeMatrixYoung$H3N2.IgG.FCd28, phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
b <- cor.test(phenotypeMatrixElderly$H3N2.IgG.FCd28, phenotypeMatrixElderly$GSVAscoreTNF, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), "\n", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`H3N2.IgG.FCd28`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,20)+
  ggtitle("GSVA TNF-NFkB vs H3N2.IgG.FCd28") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA for TNFa-NFkB geneset")  + xlab("H3N2.IgG.FCd28")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_H3N2IgG_byAge.png", device="png")

#  H1N1 IgG vs GSVA-TNFNFkB score 
grep("H1N1",colnames(phenotypeMatrix), value=T)
a <- cor.test(phenotypeMatrixYoung$H1N1.IgG.FCd28, phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
b <- cor.test(phenotypeMatrixElderly$H1N1.IgG.FCd28, phenotypeMatrixElderly$GSVAscoreTNF, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), "\n", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`H1N1.IgG.FCd28`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  ylim(-1,1) + xlim(1,1.25)+
  ggtitle("GSVA TNF-NFkB vs H1N1.IgG.FCd28") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA for TNFa-NFkB geneset")  + xlab("H1N1.IgG.FCd28")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_H1N1IgG_byAge.png", device="png")





# ***************************************  Comparison to the Tfh response ***********************************************

#  PD-1 mfi on cTfh vs GSVA-TNFNFkB score 
 
a <- cor.test(phenotypeMatrix$CD4hiNonnaive.CXCR5hi.PD1hi.....MFI_PD1., phenotypeMatrix$GSVAscoreTNF, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`CD4hiNonnaive.CXCR5hi.PD1hi.....MFI_PD1.`,y=`GSVAscoreTNF`, fill="black")) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1) + 
  geom_point(size=6, pch=21) +  #ylim(-1,1) + xlim(0,30)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=18,hjust = 0.5))+
  ggtitle("GSVA for TNFa-NFkB vs cTfh MFI PD-1 FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + 
  ylab("GSVA score for TNFa-NFkB geneset")  + xlab("MFI PD-1 foldchange")+ annotation_custom(my_grob1)+
  scale_fill_manual(values=c('black','#E69F00')) + scale_color_manual(values=c('black', '#E69F00'))
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_MFI-PD1-cTfh.png", device="png")



a <- cor.test(phenotypeMatrixYoung$CD4hiNonnaive.CXCR5hi.PD1hi.....MFI_PD1., phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
b <- cor.test(phenotypeMatrixElderly$CD4hiNonnaive.CXCR5hi.PD1hi.....MFI_PD1., phenotypeMatrixElderly$GSVAscoreTNF, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), "\n", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`CD4hiNonnaive.CXCR5hi.PD1hi.....MFI_PD1.`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,30)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=18,hjust = 0.5))+
  ggtitle("GSVA TNF-NFkB vs cTfh MFI PD-1 FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA for TNFa-NFkB geneset")  + xlab("MFI PD-1 foldchange")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_MFI-PD1-cTfh_byAge.png", device="png")



#  HiHi freq in cTfh vs GSVA-TNFNFkB score 
grep("CD38hi",colnames(phenotypeMatrix), value=T)
a <- cor.test(phenotypeMatrixYoung$Tfh_ICOShi.CD38hi...Freq..of...., phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
b <- cor.test(phenotypeMatrixElderly$Tfh_ICOShi.CD38hi...Freq..of...., phenotypeMatrixElderly$GSVAscoreTNF, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), "\n", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`Tfh_ICOShi.CD38hi...Freq..of....`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,30)+
  ggtitle("GSVA TNF-NFkB vs +/+ frequency FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA for TNFa-NFkB geneset")  + xlab("ICOS+CD38+ Freq foldchange")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_FreqHiHi_byAge.png", device="png")


#  CD27 mfi on ICOS+CD38+ cTfh vs GSVA-TNFNFkB score 
grep("MFI_PD1",colnames(phenotypeMatrix), value=T)
a <- cor.test(phenotypeMatrixYoung$Tfh_ICOShi.CD38hi.....MFI_CD27., phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
b <- cor.test(phenotypeMatrixElderly$Tfh_ICOShi.CD38hi.....MFI_CD27., phenotypeMatrixElderly$GSVAscoreTNF, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), "\n", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`Tfh_ICOShi.CD38hi.....MFI_CD27.`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,30)+
  ggtitle("GSVA TNF-NFkB vs +/+ MFI CD27 FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA for TNFa-NFkB geneset")  + xlab("MFI CD27 foldchange")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_MFI-CD27-HiHi_byAge.png", device="png")


#  PD-1 mfi on ICOS+CD38+ cTfh vs GSVA-TNFNFkB score 
grep("MFI_PD1",colnames(phenotypeMatrix), value=T)
a <- cor.test(phenotypeMatrixYoung$Tfh_ICOShi.CD38hi.....MFI_PD1., phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
b <- cor.test(phenotypeMatrixElderly$Tfh_ICOShi.CD38hi.....MFI_PD1., phenotypeMatrixElderly$GSVAscoreTNF, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18)))
ggplot(data=phenotypeMatrix, aes(x=`Tfh_ICOShi.CD38hi.....MFI_PD1.`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,30)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=18,hjust = 0.5))+
  ggtitle("GSVA TNF-NFkB vs +/+ MFI PD-1 FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA for TNFa-NFkB geneset")  + xlab("MFI PD-1 foldchange")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_MFI-PD1-HiHi_byAge.png", device="png")




#  PD-1 mfi on ICOS+CD38+ cTfh vs ASC
grep("MFI_PD1",colnames(phenotypeMatrix), value=T)
a <- cor.test(phenotypeMatrixYoung$Tfh_ICOShi.CD38hi.....MFI_PD1., phenotypeMatrixYoung$ASC.freqLive, use="complete")
b <- cor.test(phenotypeMatrixElderly$Tfh_ICOShi.CD38hi.....MFI_PD1., phenotypeMatrixElderly$ASC.freqLive, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), "\n", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`Tfh_ICOShi.CD38hi.....MFI_PD1.`,y=`ASC.freqLive`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,30)+
  ggtitle("ASC vs +/+ MFI PD-1 FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("ASC freq foldchange")  + xlab("MFI PD-1 foldchange")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/ASC_vs_MFI-PD1-HiHi_byAge.png", device="png")



#  Ki67 mfi on ICOS+CD38+ cTfh vs GSVA-TNFNFkB score 
grep("Ki67hi...FreqParent",colnames(phenotypeMatrix), value=T)
a <- cor.test(phenotypeMatrixYoung$Tfh_ICOShi.CD38hi.Ki67hi...FreqParent, phenotypeMatrixYoung$GSVAscoreTNF, use="complete")
b <- cor.test(phenotypeMatrixElderly$Tfh_ICOShi.CD38hi.Ki67hi...FreqParent, phenotypeMatrixElderly$GSVAscoreTNF, use="complete")

annotationInfo <- paste0("r= ", round(a$estimate,2), "\n", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), "\n", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=phenotypeMatrix, aes(x=`Tfh_ICOShi.CD38hi.....MFI_PD1.`,y=`GSVAscoreTNF`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,30)+
  ggtitle("GSVA TNF-NFkB vs +/+ Ki67 Freq FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA for TNFa-NFkB geneset")  + xlab("Freq Ki67 foldchange")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_Ki67-HiHi_byAge.png", device="png")



# ***************************************  GSVA vs clinical parameters ***********************************************

clinicalCharacteristics <- read.csv(file = "DifferentialExpression/GSVA/ClinicalCharacteristics_plusGSVA.csv", stringsAsFactors = F)
clinicalCharacteristics$full.id <- paste0(substr(clinicalCharacteristics$full.id,1,3), substr(clinicalCharacteristics$full.id,5,8))
rownames(clinicalCharacteristics) <- clinicalCharacteristics$full.id; clinicalCharacteristics$full.id <- NULL

summary(lm(data = clinicalCharacteristics, HALLMARK_CHOLESTEROL_HOMEOSTASIS ~ sex + age))
summary(lm(data = clinicalCharacteristics, HALLMARK_IL6_JAK_STAT3_SIGNALING ~ sex + age))
summary(lm(data = clinicalCharacteristics, HALLMARK_HEDGEHOG_SIGNALING ~ age))
summary(lm(data = clinicalCharacteristics, HALLMARK_INFLAMMATORY_RESPONSE ~ BMI))






# ****************************************************** Correlate GSVA scores against phenotype ***************************************

bestDataLog <- read.csv(file="bestDataLog.csv",stringsAsFactors = FALSE); rownames(bestDataLog) <- bestDataLog$X; bestDataLog$X <- NULL

gsets <- getGmt("DifferentialExpression/GSVA/h.all.v6.1.symbols.gmt")
HiHiv2 <- as.matrix(bestDataLog[,grep("HiHi_v2",colnames(bestDataLog))])
GSVAhallmark <- gsva(HiHiv2, gsets, method="gsva")
colnames(GSVAhallmark) <- substr(colnames(GSVAhallmark),start=2,stop=7)
rownames(GSVAhallmark) <- toTitleCase(tolower(substr(rownames(GSVAhallmark),start =10, stop=50)))

phenotypeMatrix <- read.csv(file="../../Analysis/BAA_yr4_AllMergedData_FC.csv"); phenotypeMatrix <- phenotypeMatrix[-c(62:73),]; rownames(phenotypeMatrix) <- phenotypeMatrix$Subject; 
phenotypeMatrix$Subject <- NULL

phenotypeGSVA <- merge(phenotypeMatrix, t(GSVAhallmark), by=0)
rownames(phenotypeGSVA) <- phenotypeGSVA$Row.names
phenotypeGSVA$Row.names <- phenotypeGSVA$Visit <- phenotypeGSVA$GSVA_TNF_category <- NULL
# write.csv(phenotypeGSVA, file = "phenotypeGSVA.csv")


phenotypeGSVAYoung <- phenotypeGSVA[1:6,]
phenotypeGSVAElderly <- phenotypeGSVA[7:14,]


corrAll <- cor(phenotypeGSVA[,c(2:23, 530:579)], use="complete")
p.matAll <- cor_pmat(phenotypeGSVA[,c(2:23, 530:579)], use="complete")
colorScheme <- colorRampPalette(c("#67a9cf","white", "#ef8a62"))(3)

a<- ggcorrplot(corrAll, hc.order=FALSE, type="lower",outline.col="white",lab=FALSE, p.mat=p.matAll, insig="blank", tl.cex=7, colors=colorScheme, title = "All GSVA scores")
a

corrY <- cor(phenotypeGSVAYoung[, c(2:23, 530:579)], use="complete")
p.matY <- cor_pmat(phenotypeGSVAYoung[, c(2:23, 530:579)], use="complete")
b<- ggcorrplot(corrY, hc.order=FALSE, type="lower",outline.col="white",lab=FALSE, p.mat=p.matY, insig="blank", tl.cex=5, colors=colorScheme, title = "Young GSVA scores" )
b

corrE <- cor(phenotypeGSVAElderly[, c(2:23, 530:579)], use="complete")
p.matE <- cor_pmat(phenotypeGSVAElderly[, c(2:23, 530:579)], use="complete")
c<- ggcorrplot(corrE, hc.order=FALSE, type="lower",outline.col="white",lab=FALSE, p.mat=p.matE, insig="blank", tl.cex=5, colors=colorScheme,title = "Elderly GSVA scores" )
c
# ggsave(b,filename="DifferentialExpression/GSVA/Images/interCategorycorrelations_Young.svg")
# ggsave(c,filename="DifferentialExpression/GSVA/Images/interCategorycorrelations_Elderly.svg")
# ggsave(a,filename="DifferentialExpression/GSVA/Images/interCategorycorrelations_All.svg")



#  GSVA-TNF vs GSVA-Apoptosis response 
a <- cor.test(phenotypeGSVAYoung$Tnfa_signaling_via_nfkb, phenotypeGSVAYoung$Apoptosis, use="complete")
b <- cor.test(phenotypeGSVAElderly$Tnfa_signaling_via_nfkb, phenotypeGSVAElderly$Apoptosis, use="complete")
annotationInfo <- paste0("r = ", round(a$estimate,2), ";   ", "P = ",  formatC(a$p.value, format = "e", digits = 1))
annotationInfo2 <- paste0("r = ", round(b$estimate,2), ";   ", "P = ",  formatC(b$p.value, format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.52,  y=0.12, hjust=0, gp=gpar(col="orange3", fontsize=24)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.52,  y=0.05, hjust=0, gp=gpar(col="purple", fontsize=24)))
ggplot(data=phenotypeGSVA, aes(x=`Tnfa_signaling_via_nfkb`,y=`Apoptosis`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeGSVA, Identifier=="Elderly", stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeGSVA, Identifier=="Young", stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("GSVA scores: Apoptosis vs TNF-NFkB") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("Apoptosis GSVA score")  + xlab("TNF-NFkB GSVA score")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/Apoptosis_vs_TNF_GSVAscores.pdf", device="pdf", width=8, height=6)


#  Plasmablast vs GSVA-Hedgehog response 
a <- cor.test(phenotypeGSVAYoung$Hedgehog_signaling, phenotypeGSVAYoung$`ASC.freqLive`, use="complete")
b <- cor.test(phenotypeGSVAElderly$Hedgehog_signaling, phenotypeGSVAElderly$`ASC.freqLive`, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", round(b$p.value,3))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.10, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.05, hjust=0, gp=gpar(col="purple", fontsize=18)))
ggplot(data=phenotypeGSVA, aes(x=`Hedgehog_signaling`,y=`ASC.freqLive`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeGSVA, Identifier=="Elderly", stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeGSVA, Identifier=="Young", stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("GSVA: Plasmablast resp vs Hedgehog") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("Plasmablast FC")  + xlab("Hedgehog GSVA score")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/Plasmablast_vs_Hedgehog_GSVAscores.pdf", device="pdf", width=8, height=6)


#  Tfh resp vs GSVA-Hedgehog_signaling response      
a <- cor.test(phenotypeGSVAYoung$Hedgehog_signaling, phenotypeGSVAYoung$`Tfh_ICOShi.CD38hi...Freq..of....`, use="complete")
b <- cor.test(phenotypeGSVAElderly$Hedgehog_signaling, phenotypeGSVAElderly$`Tfh_ICOShi.CD38hi...Freq..of....`, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", round(b$p.value,3))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.10, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.05, hjust=0, gp=gpar(col="purple", fontsize=18)))
ggplot(data=phenotypeGSVA, aes(x=`Hedgehog_signaling`,y=`Tfh_ICOShi.CD38hi...Freq..of....`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeGSVA, Identifier=="Elderly", stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeGSVA, Identifier=="Young", stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("GSVA: Tfh resp vs Hedgehog") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("Tfh response FC")  + xlab("Hedgehog GSVA score")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/HiHi_vs_Hedgehog_GSVAscores.pdf", device="pdf", width=8, height=6)


summary(lm(data = phenotypeGSVA, Inflammatory_response ~ Apoptosis))
summary(lm(data = phenotypeGSVA, Hedgehog_signaling ~ TNFa))



# ****************************************************** Correlate GSVA scores against transcripts ***************************************
phenotypeGSVA <- merge(phenotypeMatrix, t(GSVAhallmark), by=0); rownames(phenotypeGSVA) <- phenotypeGSVA$Row.names
phenotypeGSVA$Row.names <- phenotypeGSVA$Visit <- phenotypeGSVA$GSVA_TNF_category <- NULL

geneList <- c("BIRC5", "BCL2", "BCL2L2","BCL2A1", "BCL2L4", "BID", "BAD", "BCL2L11")
geneOfInterest <- bestDataLog[geneList,grep("HiHi_v2",colnames(bestDataLog))]
colnames(geneOfInterest) <- substr(colnames(geneOfInterest),start=2,stop=7)
phenotypeGSVA <- merge(phenotypeGSVA, t(geneOfInterest), by=0); phenotypeGSVAYoung <- phenotypeGSVA[1:6,]; phenotypeGSVAElderly <- phenotypeGSVA[7:14,]

a <- cor.test(phenotypeGSVAYoung$Tnfa_signaling_via_nfkb, phenotypeGSVAYoung$`BCL2A1`, use="complete")
b <- cor.test(phenotypeGSVAElderly$Tnfa_signaling_via_nfkb, phenotypeGSVAElderly$`BCL2A1`, use="complete")
annotationInfo <- paste0("r = ", round(a$estimate,2), ";   ", "P = ", formatC(a$p.value, format = "e", digits = 1))
annotationInfo2 <- paste0("r = ", round(b$estimate,2), ";   ", "P = ", formatC(b$p.value,format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.4,  y=0.12, hjust=0, gp=gpar(col="orange3", fontsize=24)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.4,  y=0.05, hjust=0, gp=gpar(col="purple", fontsize=24)))
ggplot(data=phenotypeGSVA, aes(x=`Tnfa_signaling_via_nfkb`,y=`BCL2A1`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeGSVA, Identifier=="Elderly", stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeGSVA, Identifier=="Young", stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("BCL2A1 vs GSVA:TNF-NFkB") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("log2 counts of BCL2A1")  + xlab("TNF-NFkB GSVA score")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "DifferentialExpression/GSVA/Images/BCL2A1_vs_TNF_GSVAscores.pdf", device="pdf")


proSurvivalGenes <- c("BCL2","BCL2L1","MCL1","BCL2A1","BCL2L12","BCL2L2","BCL2L10")
proSurvival <- bestDataLog[proSurvivalGenes,grep("HiHi_v2",colnames(bestDataLog))]
colnames(proSurvival) <- substr(colnames(proSurvival),start=2,stop=7)
a <- cor(phenotypeGSVAYoung$`Tnfa_signaling_via_nfkb`, t(proSurvival[,1:6]))
b <- cor(phenotypeGSVAElderly$`Tnfa_signaling_via_nfkb`, t(proSurvival[,7:14]))
proSurvivalYvE <- data.frame(Young=t(a),Elderly=t(b),Category="ProSurvival")

ApoptoticGenes <- c("BAK1","BCLAF1","BCL2L13","BCL2L11","BBC3","PMAIP1","BCL2L14", "BAX", "BID")
Apoptotic <- bestDataLog[ApoptoticGenes,grep("HiHi_v2",colnames(bestDataLog))]
colnames(Apoptotic) <- substr(colnames(Apoptotic),start=2,stop=7)
a<- cor(phenotypeGSVAYoung$`Tnfa_signaling_via_nfkb`, t(Apoptotic[,1:6]))
b<- cor(phenotypeGSVAElderly$`Tnfa_signaling_via_nfkb`, t(Apoptotic[,7:14]))
ApoptoticYvE <-  data.frame(Young=t(a),Elderly=t(b), Category="Apoptotic")

cellLifeDeath <- rbind(proSurvivalYvE, ApoptoticYvE)

ggplot(data=cellLifeDeath, aes(x=`Young`,y=`Elderly`, color=`Category`)) + theme_bw() +  theme(legend.position = c(0.2,0.1), legend.text = element_text(size=20)) +
  geom_point(size=4) +
  geom_text_repel(size=7, point.padding = 0.5, show.legend = F, aes(label=row.names(cellLifeDeath))) + 
  ggtitle("Correlation to GSVA TNF-NFkB score") + theme(plot.title = element_text(size=24,hjust = 0.5)) + xlim(-1,1) + ylim(-1,1)+
  ylab("Pearson r, Elderly")  + xlab("Pearson r, Young") +
  geom_vline(xintercept=0, linetype="dashed", color="grey", size=1) + geom_hline(yintercept=0, linetype="dashed", color="grey", size=1) +
  theme(axis.text = element_text(size=16,hjust = 0.5), axis.title = element_text(size=24,hjust = 0.5), legend.title = element_blank()) +
  scale_color_viridis_d(begin=0, end=0.65)
# ggsave(filename="DifferentialExpression/GSVA/Images/Bcl_family_GSVATNF_hihi_v2.pdf", device="pdf")




