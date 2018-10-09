library("grid")
library("genefilter")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("DESeq2")
library("GSVA")
library("GSEABase")
library("pheatmap")

bestDataLog <- read.csv(file="bestDataLog.csv",stringsAsFactors = FALSE); rownames(bestDataLog) <- bestDataLog$X; bestDataLog$X <- NULL

gsets <- getGmt("DifferentialExpression/GSVA/h.all.v6.1.symbols.gmt")
HiHiv2 <- as.matrix(bestDataLog[,grep("HiHi_v2",colnames(bestDataLog))])
GSVAhallmark <- gsva(HiHiv2, gsets, method="gsva")
# write.csv(GSVAhallmark, file="DifferentialExpression/GSVA/GSVA_HiHi_v2_hallmark.csv")

colnames(GSVAhallmark) <- substr(colnames(GSVAhallmark),start=2,stop=7)

rownames(GSVAhallmark) <- toTitleCase(tolower(substr(rownames(GSVAhallmark),start =10, stop=50)))

annotateHeatmap <- data.frame(row.names = colnames(GSVAhallmark), ageGroup = c(rep("Young", 6),rep("Elderly", 8)))
ann_colors = list(  ageGroup = c("Young" ="orange", "Elderly" = "purple")  )
pheatmap(GSVAhallmark, scale="none", cluster_col=T, annotation_col = annotateHeatmap, show_colnames=F, main="GSVA scores - Hallmark genesets",
         annotation_colors = ann_colors, cutree_cols = 3, cutree_rows=2, fontsize_row = 8, color=viridis(100)
        #, filename = "DifferentialExpression/GSVA/Images/GSVAhallmark_allGeneSets.png"
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
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=18,hjust = 0.5))+
  ggtitle("GSVA for TNFa-NFkB vs serum TNF") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("GSVA score for TNFa-NFkB geneset")  + xlab("TNF (pg/mL)")+
  scale_fill_manual(values=c('black','#E69F00')) + scale_color_manual(values=c('black', '#E69F00'))
# ggsave(filename = "DifferentialExpression/GSVA/Images/GSVA_TNF-NFkB_vs_serumTNF_all.png", device="png")
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


