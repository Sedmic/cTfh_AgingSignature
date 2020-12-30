library("grid")
library("genefilter")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("DESeq2")
library("GSVA")
library("GSEABase")
library("pheatmap")
library("reshape2")
library("ggcorrplot")
library("viridis")
sessionInfo()

bestDataLog <- read.csv(file="bestDataLog.csv", row.names=1)

phenotypeMatrix <- read.csv(file="20180828_phenotypeMatrix.csv", row.names=1); #rownames(bestDataLog) <- bestDataLog$X; bestDataLog$X <- NULL
phenotypeMatrixYoung <- subset(phenotypeMatrix,Identifier=="Young",stat="identity")
phenotypeMatrixElderly <- subset(phenotypeMatrix,Identifier=="Elderly",stat="identity")

# ***************************************  TNFRSF1A and B from the bestDataLog  ***********************************************

# pull TNFR1 and TNFR2 data and then combine with relevant lines from phenotypeMatrix
# ICOS+CD38+ cTfh
bestDataLog_hihi_v2 <- bestDataLog[,grep("HiHi_v2", colnames(bestDataLog))]
bestDataLog_hihi_v2 <- data.frame(t(bestDataLog_hihi_v2))
bestDataLog_hihi_v2$subject <- substr(rownames(bestDataLog_hihi_v2), start=2, stop=7)
extractGenes <- bestDataLog_hihi_v2[,c("TNFRSF1A", "TNFRSF1B", "subject")]
phenotypeMatrix$TNFRSF1Ahihid7 <- extractGenes$TNFRSF1A[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSF1Bhihid7 <- extractGenes$TNFRSF1B[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSFratiohihid7 <- phenotypeMatrix$TNFRSF1Ahihid7/phenotypeMatrix$TNFRSF1Bhihid7

bestDataLog_hihi_v1 <- bestDataLog[,grep("HiHi_v1", colnames(bestDataLog))]
bestDataLog_hihi_v1 <- data.frame(t(bestDataLog_hihi_v1))
bestDataLog_hihi_v1$subject <- substr(rownames(bestDataLog_hihi_v1), start=2, stop=7)
extractGenes <- bestDataLog_hihi_v1[,c("TNFRSF1A", "TNFRSF1B", "subject")]
phenotypeMatrix$TNFRSF1Ahihid0 <- extractGenes$TNFRSF1A[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSF1Bhihid0 <- extractGenes$TNFRSF1B[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSFratiohihid0 <- phenotypeMatrix$TNFRSF1Ahihid0/phenotypeMatrix$TNFRSF1Bhihid0
# write.csv(phenotypeMatrix, file="20180820_phenotypeMatrix.csv")   

# ICOS-CD38- cTfh
bestDataLog_lolo_v2 <- bestDataLog[,grep("LoLo_v2", colnames(bestDataLog))]
bestDataLog_lolo_v2 <- data.frame(t(bestDataLog_lolo_v2))
bestDataLog_lolo_v2$subject <- substr(rownames(bestDataLog_lolo_v2), start=2, stop=7)
extractGenes <- bestDataLog_lolo_v2[,c("TNFRSF1A", "TNFRSF1B", "subject")]
phenotypeMatrix$TNFRSF1Alolod7 <- extractGenes$TNFRSF1A[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSF1Blolod7 <- extractGenes$TNFRSF1B[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSFratiololod7 <- phenotypeMatrix$TNFRSF1Alolod7/phenotypeMatrix$TNFRSF1Blolod7

bestDataLog_lolo_v1 <- bestDataLog[,grep("LoLo_v1", colnames(bestDataLog))]
bestDataLog_lolo_v1 <- data.frame(t(bestDataLog_lolo_v1))
bestDataLog_lolo_v1$subject <- substr(rownames(bestDataLog_lolo_v1), start=2, stop=7)
extractGenes <- bestDataLog_lolo_v1[,c("TNFRSF1A", "TNFRSF1B", "subject")]
phenotypeMatrix$TNFRSF1Alolod0 <- extractGenes$TNFRSF1A[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSF1Blolod0 <- extractGenes$TNFRSF1B[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSFratiololod0 <- phenotypeMatrix$TNFRSF1Alolod0/phenotypeMatrix$TNFRSF1Blolod0
# write.csv(phenotypeMatrix, file="20180820_phenotypeMatrix.csv")   

# Naive
bestDataLog_naive_v2 <- bestDataLog[,grep("Naive_v2", colnames(bestDataLog))]
bestDataLog_naive_v2 <- data.frame(t(bestDataLog_naive_v2))
bestDataLog_naive_v2$subject <- substr(rownames(bestDataLog_naive_v2), start=2, stop=7)
extractGenes <- bestDataLog_naive_v2[,c("TNFRSF1A", "TNFRSF1B", "subject")]
phenotypeMatrix$TNFRSF1Anaived7 <- extractGenes$TNFRSF1A[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSF1Bnaived7 <- extractGenes$TNFRSF1B[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSFrationaived7 <- phenotypeMatrix$TNFRSF1Anaived7/phenotypeMatrix$TNFRSF1Bnaived7

bestDataLog_naive_v1 <- bestDataLog[,grep("HiHi_v1", colnames(bestDataLog))]
bestDataLog_naive_v1 <- data.frame(t(bestDataLog_naive_v1))
bestDataLog_naive_v1$subject <- substr(rownames(bestDataLog_naive_v1), start=2, stop=7)
extractGenes <- bestDataLog_naive_v1[,c("TNFRSF1A", "TNFRSF1B", "subject")]
phenotypeMatrix$TNFRSF1Anaived0 <- extractGenes$TNFRSF1A[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSF1Bnaived0 <- extractGenes$TNFRSF1B[match(rownames(phenotypeMatrix), extractGenes$subject)]  # match the values of TNFR1 over to phenotypeMatrix
phenotypeMatrix$TNFRSFrationaived0 <- phenotypeMatrix$TNFRSF1Anaived0/phenotypeMatrix$TNFRSF1Bnaived0
# write.csv(phenotypeMatrix, file="20180828_phenotypeMatrix.csv")   

phenotypeMatrixYoung <- subset(phenotypeMatrix,Identifier=="Young",stat="identity")
phenotypeMatrixElderly <- subset(phenotypeMatrix,Identifier=="Elderly",stat="identity")



# ********************************  compare TNFR1 and 2 by subset at day 0  *****************************************

probeList <- c("TNFRSF1A", "TNFRSF1B")
probeGenes <- bestDataLog[probeList,grep("v1",colnames(bestDataLog))]
logDataSelectedv1 <- cbind(
  probeGenes[,grep("HiHi_v1",colnames(probeGenes))], probeGenes[,grep("HiHi_v2",colnames(probeGenes))], 
  probeGenes[,grep("LoLo_v1",colnames(probeGenes))], probeGenes[,grep("LoLo_v2",colnames(probeGenes))], 
  probeGenes[,grep("Naive_v1",colnames(probeGenes))], probeGenes[,grep("Naive_v2",colnames(probeGenes))])


logDataSelectedv1 <- data.frame(t(logDataSelectedv1))
logDataSelectedv1$ageGroup <- logDataSelectedv1$subset <- "x"
logDataSelectedv1$ageGroup[grep("X111", rownames(logDataSelectedv1))] <- "Young"
logDataSelectedv1$ageGroup[grep("X222", rownames(logDataSelectedv1))] <- "Elderly"
logDataSelectedv1$subset[grep("HiHi", rownames(logDataSelectedv1))] <- "ICOS+CD38+ cTfh"
logDataSelectedv1$subset[grep("LoLo", rownames(logDataSelectedv1))] <- "ICOS-CD38- cTfh"
logDataSelectedv1$subset[grep("Naive", rownames(logDataSelectedv1))] <- "Naive CD4"
logDataSelectedv1$subset <- factor(logDataSelectedv1$subset, levels=c("ICOS+CD38+ cTfh", "ICOS-CD38- cTfh", "Naive CD4"))
logDataSelectedv1$labels <- paste0(logDataSelectedv1$ageGroup,"_",logDataSelectedv1$subset)

logDataSelectedv1$labels <- factor(logDataSelectedv1$labels, levels = c("Young_ICOS+CD38+ cTfh", "Young_ICOS-CD38- cTfh", "Young_Naive CD4", 
                                                                    "Elderly_ICOS+CD38+ cTfh", "Elderly_ICOS-CD38- cTfh", "Elderly_Naive CD4"  ))

customPalette <- c("#FDBF6F", "#FF7F00", "#B2DF8A", "#33A02C", "#fff849", "#d6cf44", "#CAB2D6", "#6A3D9A", "#03681e", "#003d10", "#99930a", "#686408")
colors <- customPalette[c(1,3,5,7,9,11)]
# just TNFR1 and 2 in the HiHi at day 0
ggplot(data=logDataSelectedv1, aes(x=subset,y=TNFRSF1A, fill=`labels`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=labels)) + 
  geom_boxplot(fill="white", outlier.shape = NA) + geom_jitter(size=6, pch=21, width=0.1) + 
  ggtitle("TNFRSF1A at day 0") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("Log-transformed counts")  + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+ theme(axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank())+
  scale_fill_manual(values=colors) + ylim(3.5, 11) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename = "Images/TNFR1_allSubsets_v1.pdf", device="pdf")

ggplot(data=logDataSelectedv1, aes(x=subset,y=TNFRSF1B, fill=`labels`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=labels)) + 
  geom_boxplot(fill="white", outlier.shape = NA) + geom_jitter(size=6, pch=21, width=0.1) + 
  ggtitle("TNFRSF1B at day 0") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("Log-transformed counts")  + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+ theme(axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank())+
  scale_fill_manual(values=colors) + ylim(7.5,16) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename = "Images/TNFR2_allSubsets_v1.pdf", device="pdf")


TukeyHSD(aov(data = logDataSelectedv1, TNFRSF1A ~ subset))
TukeyHSD(aov(data = logDataSelectedv1, TNFRSF1B ~ subset))


# ********************************  compare TNFR1 and 2 to vaccination response *****************************************

# what about the day 0 to day 7 change in TNFR1 and 2? 
# just ICOS+CD38+ cTfh
probeList <- c("TNFRSF1A", "TNFRSF1B", "BCL2A1", "MCL1", "BID")
probeGenes <- bestDataLog[probeList,grep("HiHi",colnames(bestDataLog))]  
logDataSelected <- cbind(
  probeGenes[,grep("HiHi_v1",colnames(probeGenes))], probeGenes[,grep("HiHi_v2",colnames(probeGenes))], 
  probeGenes[,grep("LoLo_v1",colnames(probeGenes))], probeGenes[,grep("LoLo_v2",colnames(probeGenes))], 
  probeGenes[,grep("Naive_v1",colnames(probeGenes))], probeGenes[,grep("Naive_v2",colnames(probeGenes))])
logDataSelected <- data.frame(t(logDataSelected))
logDataSelected$ageGroup <- c(rep("Young", 6), rep("Elderly",7),rep("Young", 6), rep("Elderly",8) )
logDataSelected$visit <- substr(rownames(logDataSelected), start=15,stop=16)
logDataSelected$subject <- substr(rownames(logDataSelected), start=2,stop=7)

# TNFR1 pre-post vaccine
TNFRSF1A_beforeAfter <- dcast(logDataSelected, subject~visit, value.var = "TNFRSF1A" )  # give me a before-after view of the matrix
TNFRSF1A_beforeAfter <- na.omit(TNFRSF1A_beforeAfter); colnames(TNFRSF1A_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=TNFRSF1A_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=TNFRSF1A_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=TNFRSF1A_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of TNFRSF1A") + ggtitle("TNFRSF1A Pre-Post vaccine") + 
  theme(axis.text = element_text(size=28,hjust = 0.5))+ theme(axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=36,hjust = 0.5)) 
# ggsave(filename="Images/TNFR1_HiHi_beforeAfter.pdf", device="pdf")

t.test(TNFRSF1A_beforeAfter$visit1[1:6],TNFRSF1A_beforeAfter$visit2[1:6],paired = T)

# TNFR2 pre-post vaccine
TNFRSF1B_beforeAfter <- dcast(logDataSelected, subject~visit, value.var = "TNFRSF1B" )  # give me a before-after view of the matrix
TNFRSF1B_beforeAfter <- na.omit(TNFRSF1B_beforeAfter); colnames(TNFRSF1B_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=TNFRSF1B_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=TNFRSF1B_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=TNFRSF1B_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of TNFRSF1B") + ggtitle("TNFRSF1B Pre-Post vaccine") + 
  theme(axis.text = element_text(size=24,hjust = 0.5))+ theme(axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=30,hjust = 0.5)) 
# ggsave(filename="Images/TNFR2_HiHi_beforeAfter.pdf", device="pdf")


# ratio of TNFR1 to TNFR2 pre-post vaccine
temp <- logDataSelected
temp$ratio <- temp$TNFRSF1A / temp$TNFRSF1B
TNFRSFratio_beforeAfter <- dcast(temp, subject~visit, value.var = "ratio" )  # give me a before-after view of the matrix
TNFRSFratio_beforeAfter <- na.omit(TNFRSFratio_beforeAfter); colnames(TNFRSFratio_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=TNFRSFratio_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=TNFRSFratio_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=TNFRSFratio_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of TNFRSF1A/TNFRSF1B") + ggtitle("Ratio Pre-Post vaccine") + 
  theme(axis.text = element_text(size=24,hjust = 0.5))+ theme(axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=30,hjust = 0.5)) 
# ggsave(filename="Images/TNFRratio_HiHi_beforeAfter.pdf", device="pdf")



# BCL2A1 pre-post vaccine
BCL2A1_beforeAfter <- dcast(logDataSelected, subject~visit, value.var = "BCL2A1" )  # give me a before-after view of the matrix
BCL2A1_beforeAfter <- na.omit(BCL2A1_beforeAfter); colnames(BCL2A1_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=BCL2A1_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=BCL2A1_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=BCL2A1_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of BCL2A1") + ggtitle("BCL2A1 Pre-Post vaccine") + 
  theme(axis.text = element_text(size=24,hjust = 0.5))+ theme(axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=30,hjust = 0.5)) 
# ggsave(filename="Images/BCL2A1_HiHi_beforeAfter.pdf", device="pdf")


# MCL1 pre-post vaccine
MCL1_beforeAfter <- dcast(logDataSelected, subject~visit, value.var = "MCL1" )  # give me a before-after view of the matrix
MCL1_beforeAfter <- na.omit(MCL1_beforeAfter); colnames(MCL1_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=MCL1_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=MCL1_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=MCL1_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of MCL1") + ggtitle("MCL1 Pre-Post vaccine") + 
  theme(axis.text = element_text(size=24,hjust = 0.5))+ theme(axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=30,hjust = 0.5)) 
# ggsave(filename="Images/MCL1_HiHi_beforeAfter.pdf", device="pdf")


# BID pre-post vaccine
BID_beforeAfter <- dcast(logDataSelected, subject~visit, value.var = "BID" )  # give me a before-after view of the matrix
BID_beforeAfter <- na.omit(BID_beforeAfter); colnames(BID_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=BID_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=BID_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=BID_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of BID") + ggtitle("BID Pre-Post vaccine") + 
  theme(axis.text = element_text(size=24,hjust = 0.5))+ theme(axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=30,hjust = 0.5)) 
# ggsave(filename="Images/BID_HiHi_beforeAfter.pdf", device="pdf")


#  ****************************** ICOS-CD38- cTfh subset as pre-post  *************************************

# just ICOS-CD38- cTfh
probeList <- c("TNFRSF1A", "TNFRSF1B")
probeGenes <- bestDataLog[probeList,grep("LoLo",colnames(bestDataLog))]  
logDataSelected <- cbind(
  probeGenes[,grep("HiHi_v1",colnames(probeGenes))], probeGenes[,grep("HiHi_v2",colnames(probeGenes))], 
  probeGenes[,grep("LoLo_v1",colnames(probeGenes))], probeGenes[,grep("LoLo_v2",colnames(probeGenes))], 
  probeGenes[,grep("Naive_v1",colnames(probeGenes))], probeGenes[,grep("Naive_v2",colnames(probeGenes))])

# TNFR1 pre-post vaccine
logDataSelected <- data.frame(t(logDataSelected))
logDataSelected$ageGroup <- c(rep("Young", 6), rep("Elderly",8),rep("Young", 6), rep("Elderly",8) )
logDataSelected$visit <- substr(rownames(logDataSelected), start=15,stop=16)
logDataSelected$subject <- substr(rownames(logDataSelected), start=2,stop=7)
TNFRSF1A_beforeAfter <- dcast(logDataSelected, subject~visit, value.var = "TNFRSF1A" )  # give me a before-after view of the matrix
TNFRSF1A_beforeAfter <- na.omit(TNFRSF1A_beforeAfter); colnames(TNFRSF1A_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=TNFRSF1A_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=TNFRSF1A_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=TNFRSF1A_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of TNFRSF1A") + ggtitle("TNFRSF1A Pre-Post - ICOS-CD38- cTfh") + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+ theme(axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=28,hjust = 0.5)) 
# ggsave(filename="Images/TNFR1_LoLo_beforeAfter.png", device="png")


# TNFR2 pre-post vaccine
TNFRSF1B_beforeAfter <- dcast(logDataSelected, subject~visit, value.var = "TNFRSF1B" )  # give me a before-after view of the matrix
TNFRSF1B_beforeAfter <- na.omit(TNFRSF1B_beforeAfter); colnames(TNFRSF1B_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=TNFRSF1B_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=TNFRSF1B_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=TNFRSF1B_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of TNFRSF1A") + ggtitle("TNFRSF1B Pre-Post - ICOS-CD38- cTfh") + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+ theme(axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=28,hjust = 0.5)) 
# ggsave(filename="Images/TNFR2_LoLo_beforeAfter.png", device="png")


# ratio of TNFR1 to TNFR2 pre-post vaccine
temp <- logDataSelected
temp$ratio <- temp$TNFRSF1A / temp$TNFRSF1B
TNFRSFratio_beforeAfter <- dcast(temp, subject~visit, value.var = "ratio" )  # give me a before-after view of the matrix
TNFRSFratio_beforeAfter <- na.omit(TNFRSFratio_beforeAfter); colnames(TNFRSFratio_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=TNFRSFratio_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=TNFRSFratio_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=TNFRSFratio_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of TNFRSF1A/TNFRSF1B") + ggtitle("Ratio Pre-Post - ICOS-CD38- cTfh") + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+ theme(axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=28,hjust = 0.5)) 
# ggsave(filename="Images/TNFRratio_LoLo_beforeAfter.png", device="png")

#  ******************************************************************************


#  ****************************** Naive CD4 subset as pre-post  *************************************

# just Naive CD4
probeList <- c("TNFRSF1A", "TNFRSF1B")
probeGenes <- bestDataLog[probeList,grep("Naive",colnames(bestDataLog))]  
logDataSelected <- cbind(
  probeGenes[,grep("HiHi_v1",colnames(probeGenes))], probeGenes[,grep("HiHi_v2",colnames(probeGenes))], 
  probeGenes[,grep("LoLo_v1",colnames(probeGenes))], probeGenes[,grep("LoLo_v2",colnames(probeGenes))], 
  probeGenes[,grep("Naive_v1",colnames(probeGenes))], probeGenes[,grep("Naive_v2",colnames(probeGenes))])

# TNFR1 pre-post vaccine
logDataSelected <- data.frame(t(logDataSelected))
logDataSelected$ageGroup <- c(rep("Young", 6), rep("Elderly",8),rep("Young", 6), rep("Elderly",8) )
logDataSelected$visit <- substr(rownames(logDataSelected), start=15,stop=16)
logDataSelected$subject <- substr(rownames(logDataSelected), start=2,stop=7)
TNFRSF1A_beforeAfter <- dcast(logDataSelected, subject~visit, value.var = "TNFRSF1A" )  # give me a before-after view of the matrix
TNFRSF1A_beforeAfter <- na.omit(TNFRSF1A_beforeAfter); colnames(TNFRSF1A_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=TNFRSF1A_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=TNFRSF1A_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=TNFRSF1A_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of TNFRSF1A") + ggtitle("TNFRSF1A Pre-Post - Naive CD4") + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+ theme(axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=28,hjust = 0.5)) 
# ggsave(filename="Images/TNFR1_Naive_beforeAfter.png", device="png")


# TNFR2 pre-post vaccine
TNFRSF1B_beforeAfter <- dcast(logDataSelected, subject~visit, value.var = "TNFRSF1B" )  # give me a before-after view of the matrix
TNFRSF1B_beforeAfter <- na.omit(TNFRSF1B_beforeAfter); colnames(TNFRSF1B_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=TNFRSF1B_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=TNFRSF1B_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=TNFRSF1B_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of TNFRSF1A") + ggtitle("TNFRSF1B Pre-Post - Naive CD4") + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+ theme(axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=28,hjust = 0.5)) 
# ggsave(filename="Images/TNFR2_Naive_beforeAfter.png", device="png")


# ratio of TNFR1 to TNFR2 pre-post vaccine
temp <- logDataSelected
temp$ratio <- temp$TNFRSF1A / temp$TNFRSF1B
TNFRSFratio_beforeAfter <- dcast(temp, subject~visit, value.var = "ratio" )  # give me a before-after view of the matrix
TNFRSFratio_beforeAfter <- na.omit(TNFRSFratio_beforeAfter); colnames(TNFRSFratio_beforeAfter) <- c("subject","visit1","visit2")
ggplot(data=TNFRSFratio_beforeAfter, aes()) + theme_bw() + 
  geom_segment(data=TNFRSFratio_beforeAfter[1:6,], aes(x=" Day0", xend=" Day7", y=visit1, yend=visit2), color="orange",size=1) + 
  geom_segment(data=TNFRSFratio_beforeAfter[7:13,], aes(x="Day0", xend="Day7", y=visit1, yend=visit2), colour="purple",size=1) +
  ylab("Log Counts of TNFRSF1A/TNFRSF1B") + ggtitle("Ratio Pre-Post - Naive CD4") + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+ theme(axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank())+
  theme(plot.title = element_text(size=28,hjust = 0.5)) 
# ggsave(filename="Images/TNFRratio_Naive_beforeAfter.png", device="png")



#  *********************************** t-test comparison plots - ICOS+CD38+ cTfh ********************************

# now just day 7 for ICOS+CD38+ cTfh
probeList <- c("TNFRSF1A", "TNFRSF1B")
probeGenes <- bestDataLog[probeList,grep("HiHi_v2",colnames(bestDataLog))]
logDataSelected <- cbind(
  probeGenes[,grep("HiHi_v1",colnames(probeGenes))], probeGenes[,grep("HiHi_v2",colnames(probeGenes))], 
  probeGenes[,grep("LoLo_v1",colnames(probeGenes))], probeGenes[,grep("LoLo_v2",colnames(probeGenes))], 
  probeGenes[,grep("Naive_v1",colnames(probeGenes))], probeGenes[,grep("Naive_v2",colnames(probeGenes))])


logDataSelected <- data.frame(t(logDataSelected))
logDataSelected$ageGroup <- c(rep("Young", 6), rep("Elderly",8))

a <- t.test(logDataSelected$TNFRSF1A[1:6], logDataSelected$TNFRSF1A[7:14])
annotationInfo <- paste0("p= ", round(a$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=32)))
ggplot(data=logDataSelected, aes(x=ageGroup,y=TNFRSF1A, fill=`ageGroup`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=ageGroup)) + 
  geom_boxplot(fill="white") + geom_point(size=6, pch=21) + 
  ggtitle("TNFRSF1A at day 7") + theme(plot.title = element_text(size=36,hjust = 0.5)) + ylab("Log counts TNFRSF1A")  + 
  theme(axis.text = element_text(size=32,hjust = 0.5))+theme(axis.title = element_text(size=36,hjust = 0.5))+theme(axis.title.x = element_blank())+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1)
# ggsave(filename = "Images/TNFR1_HiHi_v2.pdf", device="pdf")

a <- t.test(logDataSelected$TNFRSF1B[1:6], logDataSelected$TNFRSF1B[7:14])
annotationInfo <- paste0("p= ", round(a$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.7,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=32)))
ggplot(data=logDataSelected, aes(x=ageGroup,y=TNFRSF1B, fill=`ageGroup`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=ageGroup)) + 
  geom_boxplot(fill="white") + geom_point(size=6, pch=21) + 
  ggtitle("TNFRSF1B at day 7") + theme(plot.title = element_text(size=36,hjust = 0.5)) + ylab("Log counts TNFRSF1B")  + 
  theme(axis.text = element_text(size=36,hjust = 0.5))+theme(axis.title = element_text(size=36,hjust = 0.5))+theme(axis.title.x = element_blank())+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1)
# ggsave(filename = "Images/TNFR2_HiHi_v2.pdf", device="pdf")

logDataSelected$ratio <- logDataSelected$TNFRSF1A/logDataSelected$TNFRSF1B
a <- t.test(logDataSelected$ratio[1:6], logDataSelected$ratio[7:14])
annotationInfo <- paste0("p= ", round(a$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.8,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=24)))
ggplot(data=logDataSelected, aes(x=ageGroup,y=TNFRSF1A/TNFRSF1B, fill=`ageGroup`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=ageGroup)) + 
  geom_boxplot(fill="white") + geom_point(size=6, pch=21) + 
  ggtitle("TNFR1-to-TNFR2 Ratio") + theme(plot.title = element_text(size=36,hjust = 0.5)) + ylab("TNFRSF1A/TNFRSF1B")  + 
  theme(axis.text = element_text(size=28,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))+theme(axis.title.x = element_blank())+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1)
# ggsave(filename = "Images/TNFR1-to-TNFR2_HiHi_v2.pdf", device="pdf")



#  *********************************** t-test comparison plots - ICOS-CD38- cTfh ********************************

# now just day 7 for ICOS-CD38- cTfh
probeList <- c("TNFRSF1A", "TNFRSF1B")
probeGenes <- bestDataLog[probeList,grep("LoLo_v2",colnames(bestDataLog))]
logDataSelected <- cbind(
  probeGenes[,grep("HiHi_v1",colnames(probeGenes))], probeGenes[,grep("HiHi_v2",colnames(probeGenes))], 
  probeGenes[,grep("LoLo_v1",colnames(probeGenes))], probeGenes[,grep("LoLo_v2",colnames(probeGenes))], 
  probeGenes[,grep("Naive_v1",colnames(probeGenes))], probeGenes[,grep("Naive_v2",colnames(probeGenes))])


logDataSelected <- data.frame(t(logDataSelected))
logDataSelected$ageGroup <- c(rep("Young", 6), rep("Elderly",8))

a <- t.test(logDataSelected$TNFRSF1A[1:6], logDataSelected$TNFRSF1A[7:14])
annotationInfo <- paste0("p= ", round(a$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.85,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=logDataSelected, aes(x=ageGroup,y=TNFRSF1A, fill=`ageGroup`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=ageGroup)) + 
  geom_boxplot(fill="white") + geom_point(size=6, pch=21) + 
  ggtitle("TNFRSF1A at day 7 ICOS-CD38- cTfh") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("Log counts TNFRSF1A")  + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+theme(axis.title.x = element_blank())+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1)
# ggsave(filename = "Images/TNFR1_LoLo_v2.eps", device="eps")

a <- t.test(logDataSelected$TNFRSF1B[1:6], logDataSelected$TNFRSF1B[7:14])
annotationInfo <- paste0("p= ", round(a$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.85,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=logDataSelected, aes(x=ageGroup,y=TNFRSF1B, fill=`ageGroup`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=ageGroup)) + 
  geom_boxplot(fill="white") + geom_point(size=6, pch=21) + 
  ggtitle("TNFRSF1B at day 7 ICOS-CD38- cTfh") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("Log counts TNFRSF1B")  + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+theme(axis.title.x = element_blank())+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1)
# ggsave(filename = "Images/TNFR2_LoLo_v2.png", device="png")

logDataSelected$ratio <- logDataSelected$TNFRSF1A/logDataSelected$TNFRSF1B
a <- t.test(logDataSelected$ratio[1:6], logDataSelected$ratio[7:14])
annotationInfo <- paste0("p= ", round(a$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.85,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=14)))
ggplot(data=logDataSelected, aes(x=ageGroup,y=TNFRSF1A/TNFRSF1B, fill=`ageGroup`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=ageGroup)) + 
  geom_boxplot(fill="white") + geom_point(size=6, pch=21) + 
  ggtitle("TNFR1-to-TNFR2 Ratio ICOS-CD38- cTfh") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("TNFRSF1A/TNFRSF1B")  + 
  theme(axis.text = element_text(size=22,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+theme(axis.title.x = element_blank())+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1)
# ggsave(filename = "Images/TNFR1-to-TNFR2_LoLo_v2.png", device="png")




# ********************************  compare TNFR1 and 2 to plasmablast response *****************************************

a <- cor.test(phenotypeMatrixYoung$ASC.freqLive, phenotypeMatrixYoung$TNFRSF1Ahihid7, use="complete")
b <- cor.test(phenotypeMatrixElderly$ASC.freqLive, phenotypeMatrixElderly$TNFRSF1Ahihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,3), ";   ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.55,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=24)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.55,  y=0.87, hjust=0, gp=gpar(col="purple", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`ASC.freqLive`,y=`TNFRSF1Ahihid7`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  
  ggtitle("TNFRSF1A vs Plasmablasts") + theme(plot.title = element_text(size=36,hjust = 0.5)) + xlim(0,2.5) +
  ylab("TNFRSF1A at day 7")  + xlab("Plasmablast freq foldchange")+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/TNFRSF1A_vs_ASC_byAge.pdf", device="pdf")


a <- cor.test(phenotypeMatrixYoung$ASC.freqLive, phenotypeMatrixYoung$TNFRSF1Bhihid7, use="complete")
b <- cor.test(phenotypeMatrixElderly$ASC.freqLive, phenotypeMatrixElderly$TNFRSF1Bhihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=24)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.87, hjust=0, gp=gpar(col="purple", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`ASC.freqLive`,y=`TNFRSF1Bhihid7`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  
  ggtitle("TNFRSF1B vs Plasmablasts") + theme(plot.title = element_text(size=36,hjust = 0.5)) + xlim(0,2.5) +
  ylab("TNFRSF1B at day 7")  + xlab("Plasmablast freq foldchange")+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/TNFRSF1B_vs_ASC_byAge.pdf", device="pdf")


# TNFR ratio vs plasmablast responses
a <- cor.test(phenotypeMatrixYoung$ASC.freqLive, phenotypeMatrixYoung$TNFRSFratiohihid7, use="complete")
b <- cor.test(phenotypeMatrixElderly$ASC.freqLive, phenotypeMatrixElderly$TNFRSFratiohihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";     ", "p = ", formatC(a$p.value, format="e", digits=1))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.49,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=24)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.49,  y=0.87, hjust=0, gp=gpar(col="purple", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`ASC.freqLive`,y=`TNFRSFratiohihid7`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  
  ggtitle("TNFR ratio vs Plasmablasts") + theme(plot.title = element_text(size=36,hjust = 0.5)) + xlim(0,2.5) +
  ylab("TNFRSF1A/TNFRSF1B, day 7")  + xlab("Plasmablast foldchange")+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/TNFRSFratio_vs_ASC_byAge.pdf", device="pdf")



a <- cor.test(phenotypeMatrix$ASC.freqLive, phenotypeMatrix$TNFRSF1Ahihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.58,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`ASC.freqLive`,y=`TNFRSF1Ahihid7`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1) + 
  geom_point(size=6, pch=21, fill="black") + 
  ggtitle("ICOS+CD38+ cTfh d7 vs Plasmablast FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + xlim(0,2.5) +
  ylab("TNFRSF1A at day 7")  + xlab("Plasmablast freq foldchange")+
  #scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1)+ 
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))
# ggsave(filename = "Images/TNFRSF1A_vs_ASC_allTogether.pdf", device="pdf")

a <- cor.test(phenotypeMatrix$ASC.freqLive, phenotypeMatrix$TNFRSF1Bhihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.58,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`ASC.freqLive`,y=`TNFRSF1Bhihid7`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1) + 
  geom_point(size=6, pch=21, fill="black") + 
  ggtitle("ICOS+CD38+ cTfh d7 vs Plasmablast FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + xlim(0,2.5) +
  ylab("TNFRSF1B at day 7")  + xlab("Plasmablast freq foldchange")+
  #scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1)+ 
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))
# ggsave(filename = "Images/TNFRSF1B_vs_ASC_allTogether.pdf", device="pdf")

a <- cor.test(phenotypeMatrix$ASC.freqLive, phenotypeMatrix$TNFRSFratiohihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.48,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`ASC.freqLive`,y=`TNFRSFratiohihid7`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1) + 
  geom_point(size=6, pch=21, fill="black") + 
  ggtitle("ICOS+CD38 cTfh d7 vs Plasmablast FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + xlim(0,2.5) +
  ylab("TNFRSF1A/TNFRSF1B, day 7")  + xlab("Plasmablast freq foldchange")+
  #scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1)+ 
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))
# ggsave(filename = "Images/TNFRSFratio_vs_ASC_allTogether.pdf", device="pdf")


## ****************************************** correlation matrix & heatmap view ********************************************

keyColumns <- grep("^ASC.freqLive", colnames(phenotypeMatrix))  # only match the if ASC is at the beginning of the string
keyColumns <- c(keyColumns, grep("TNFRSF", colnames(phenotypeMatrix)))
keyColumns <- c(keyColumns, grep("nAb", colnames(phenotypeMatrix)))
day0Columns <- keyColumns[grep("d0$",colnames(phenotypeMatrix[,keyColumns]))] #d0 match at the end of the string
keyColumnsd7 <- setdiff(keyColumns,day0Columns)

correlMatrix <- round(cor(phenotypeMatrix[,keyColumnsd7],use="complete.obs", method="pearson"),2)
pMatrix <- cor_pmat(phenotypeMatrix[,keyColumnsd7],use="complete.obs", method="pearson")
ggcorrplot(correlMatrix, hc.order = T, type="lower", p.mat = pMatrix, insig = "blank")



simplifyCorrel <- correlMatrix[c(11,12),]  # column 1 would have been for plasmablasts
colsToRemove <- grep("ASC",colnames(simplifyCorrel))
colsToRemove <- c(colsToRemove, grep("nAb",colnames(simplifyCorrel)))
simplifyCorrel <- simplifyCorrel[,-colsToRemove]
rownames(simplifyCorrel) <- c("H3N2.nAb.FCd28","H1N1.nAb.FCd28")  # removed first element of this vector which was Plasmablast FC
annotateHeatmap <- data.frame(row.names = colnames(simplifyCorrel), 
                              CD4subset = c(rep("ICOS+CD38+ cTfh", 3), rep("ICOS-CD38- cTfh", 3),rep("Naive CD4", 3)),
                              TNFR= c("TNFRSF1A","TNFRSF1B","TNFRSF1A/TNFRSF1B","TNFRSF1A","TNFRSF1B","TNFRSF1A/TNFRSF1B",
                                      "TNFRSF1A","TNFRSF1B","TNFRSF1A/TNFRSF1B")
)
ann_colors = list(  CD4subset = c("ICOS+CD38+ cTfh" ="orange", "ICOS-CD38- cTfh" = "green", "Naive CD4" ="yellow"),
                    TNFR = c("TNFRSF1A"="cyan","TNFRSF1B"="beige","TNFRSF1A/TNFRSF1B"="black"))
breaksList <- seq(from=-1,to=1, by=0.02)
pheatmap(simplifyCorrel, cluster_row=T, cluster_col=F, show_rownames=T, breaks=breaksList, show_colnames=F,
         colors=colorRampPalette(rev(brewer.pal(n=3,name="RdYlBu")))(length(breaksList)), gaps_col = c(3,6)
         , annotation_col = annotateHeatmap, annotation_colors = ann_colors,height=2.2, color=inferno(100)
         #, filename="Images/TNFRSF_correlations.pdf"
)




## ****************************************** Multivariate model ********************************************

summary(lm(data=phenotypeMatrix, `ASC.freqLive`~ H1N1.nAb.FCd28)) # no statistical significance for the interaction
plot(data=phenotypeMatrix, `ASC.freqLive`~ H1N1.nAb.FCd28)
summary(lm(data=phenotypeMatrix, `ASC.freqLive`~ H3N2.nAb.FCd28)) # no statistical significance for the interaction
plot(data=phenotypeMatrix, `ASC.freqLive`~ H3N2.nAb.FCd28)

summary(lm(data=phenotypeMatrix, `ASC.freqLive`~ TNFa*TNFRSF1Ahihid7)) # no statistical significance for the interaction
summary(lm(data=phenotypeMatrix, `ASC.freqLive`~ TNFa*TNFRSF1Bhihid7))   # no statistical significance for the interaction 
summary(lm(data=phenotypeMatrix, `ASC.freqLive`~ TNFa*TNFRSFratiohihid7)) # significant!
summary(lm(data=phenotypeMatrix, `ASC.freqLive`~ TNFRSF1Ahihid7)) # significant!
summary(lm(data=phenotypeMatrix, `ASC.freqLive`~ TNFRSF1Bhihid7))   # no statistical significance 
summary(lm(data=phenotypeMatrix, `ASC.freqLive`~ TNFRSFratiohihid7))  # significant!

antibodyColumns <- grep("H1N1",colnames(phenotypeMatrix), value = F)
antibodyColumns <- c(antibodyColumns, grep("H3N2",colnames(phenotypeMatrix), value = F))
outputLinearModels <- list()

# look at linear models to compare antibody readouts vs TNFR1 
for(i in 1:length(antibodyColumns)) {   outputLinearModels[[i]] <- lm(data=phenotypeMatrix, phenotypeMatrix[,antibodyColumns[i]] ~ TNFRSF1Ahihid7) }
temp <- sapply(outputLinearModels, function(x) summary(x)$coefficients[,4])
significantHits <- which(temp[2,] < 0.05)
colnames(phenotypeMatrix)[antibodyColumns[significantHits]]
# look at linear models to compare antibody readouts vs TNFR2 
for(i in 1:length(antibodyColumns)) {   outputLinearModels[[i]] <- lm(data=phenotypeMatrix, phenotypeMatrix[,antibodyColumns[i]] ~ TNFRSF1Bhihid7) }
temp <- sapply(outputLinearModels, function(x) summary(x)$coefficients[,4])
significantHits <- which(temp[2,] < 0.05)
colnames(phenotypeMatrix)[antibodyColumns[significantHits]]
# look at linear models to compare antibody readouts vs TNFRratio 
for(i in 1:length(antibodyColumns)) {   outputLinearModels[[i]] <- lm(data=phenotypeMatrix, phenotypeMatrix[,antibodyColumns[i]] ~ TNFRSFratiohihid7) }
temp <- sapply(outputLinearModels, function(x) summary(x)$coefficients[,4])
significantHits <- which(temp[2,] < 0.05)
colnames(phenotypeMatrix)[antibodyColumns[significantHits]]
summary(lm(data=phenotypeMatrix, `H1N1.nAb.FCd28`~ TNFRSFratiohihid7))  # validation
summary(lm(data=phenotypeMatrix, `H3N2.nAb.FCd28`~ TNFRSFratiohihid7))
summary(lm(data=phenotypeMatrix, `H3N2.IgG.FCd7`~ TNFRSFratiohihid7))
cor(phenotypeMatrix$H3N2.IgG.FCd7, phenotypeMatrix$TNFRSFratiohihid7, use="complete.obs")
cor(phenotypeMatrix$H1N1.IgG.FCd7, phenotypeMatrix$TNFRSFratiohihid7, use="complete.obs")

plot(data=phenotypeMatrix, `H3N2.IgG.FCd7`~ TNFRSFratiohihid7)
plot(data=phenotypeMatrix, `H1N1.IgG.FCd7`~ TNFRSFratiohihid7)


summary(lm(data=phenotypeMatrix, `H3N2.nAb.FCd28`~ CXCL11*TNFRSF1Ahihid7)) # no statistical significance for the interaction
summary(lm(data=phenotypeMatrix, `ASC.freqLive`~ TNFRSF1Bhihid7 ))   # no statistical significance for the interaction 

summary(lm(data=phenotypeMatrix, `H3N2.nAb.FCd28`~ TNFRSFratiohihid7))
summary(lm(data=phenotypeMatrix, `H3N2.nAb.FCd28`~ TNFa*TNFRSFratiohihid7))
summary(lm(data=phenotypeMatrix, `H1N1.nAb.FCd28`~ TNFRSFratiohihid7))
summary(lm(data=phenotypeMatrix, `H1N1.nAb.FCd28`~ TNFa*TNFRSFratiohihid7))

a <- lm(data=phenotypeMatrix, `H1N1.nAb.FCd28`~ TNFa*TNFRSFratiohihid7)
b <- lm(data=phenotypeMatrix, `H3N2.nAb.FCd28`~ TNFa*TNFRSFratiohihid7)
c <- lm(data=phenotypeMatrix, `H1N1.nAb.FCd28`~ TNFRSFratiohihid7)
d <- lm(data=phenotypeMatrix, `H3N2.nAb.FCd28`~ TNFRSFratiohihid7)
e <- lm(data=phenotypeMatrix, `H1N1.nAb.FCd28`~ TNFa)
f <- lm(data=phenotypeMatrix, `H3N2.nAb.FCd28`~ TNFa)


a <- cor.test(phenotypeMatrix$TNFa, phenotypeMatrix$TNFRSFratiohihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.55,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=phenotypeMatrix, aes(x=`TNFa`,y=`TNFRSFratiohihid7`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1) + 
  geom_point(size=6, pch=21, fill="black") + 
  ggtitle("ICOS+CD38 cTfh d7 vs Plasmablast FC") + theme(plot.title = element_text(size=24,hjust = 0.5)) + #xlim(0,2.5) +
  ylab("TNFRSF1A/TNFRSF1B, day 7")  + xlab("TNFa serum (pg/mL)")+
  #scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1)+ 
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))
# ggsave(filename = "Images/TNFRSFratio_vs_ASC_allTogether.pdf", device="pdf")




## ****************************************** nAb responses and correlations ********************************************


# TNFR ratio vs nAb responses
a <- cor.test(phenotypeMatrixYoung$`H1N1.nAb.FCd28`, phenotypeMatrixYoung$TNFRSFratiohihid7, use="complete")
b <- cor.test(phenotypeMatrixElderly$`H1N1.nAb.FCd28`, phenotypeMatrixElderly$TNFRSFratiohihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=24)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.87, hjust=0, gp=gpar(col="purple", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`H1N1.nAb.FCd28`,y=`TNFRSFratiohihid7`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  
  ggtitle("TNFR1-to-2 vs H1N1 HAI Ab") + theme(plot.title = element_text(size=36,hjust = 0.5)) + xlim(0,10) +
  ylab("TNFRSF1A/TNFRSF1B, day 7")  + xlab("H1N1.nAb.FCd28")+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/TNFRSFratio_vs_H1N1nAb_byAge.pdf", device="pdf")


a <- cor.test(phenotypeMatrix$`H1N1.nAb.FCd28`, phenotypeMatrix$TNFRSFratiohihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`H1N1.nAb.FCd28`,y=`TNFRSFratiohihid7`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1) + 
  geom_point(size=6, pch=21, fill="black") + 
  ggtitle("TNFR1-to-2 vs H1N1 HAI Ab") + theme(plot.title = element_text(size=24,hjust = 0.5)) + xlim(0,10) +
  ylab("TNFRSF1A/TNFRSF1B, day 7")  + xlab("H1N1.nAb.FCd28 foldchange")+
  #scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1)+ 
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))
# ggsave(filename = "Images/TNFRSFratio_vs_H1N1nAb_allTogether.pdf", device="pdf")


a <- cor.test(phenotypeMatrixYoung$`H3N2.nAb.FCd28`, phenotypeMatrixYoung$TNFRSFratiohihid7, use="complete")
b <- cor.test(phenotypeMatrixElderly$`H3N2.nAb.FCd28`, phenotypeMatrixElderly$TNFRSFratiohihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=24)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.87, hjust=0, gp=gpar(col="purple", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`H3N2.nAb.FCd28`,y=`TNFRSFratiohihid7`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  
  ggtitle("TNFR1-to-2 vs H3N2 HAI Ab") + theme(plot.title = element_text(size=36,hjust = 0.5)) + xlim(0,33) + ylim(0.5,0.9)+
  ylab("TNFRSF1A/TNFRSF1B, day 7")  + xlab("H3N2.nAb.FCd28")+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/TNFRSFratio_vs_H3N2nAb_byAge.pdf", device="pdf")


a <- cor.test(phenotypeMatrix$`H3N2.nAb.FCd28`, phenotypeMatrix$TNFRSFratiohihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`H3N2.nAb.FCd28`,y=`TNFRSFratiohihid7`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1) + 
  geom_point(size=6, pch=21, fill="black") + 
  ggtitle("TNFR1-to-2 vs H3N2 HAI Ab") + theme(plot.title = element_text(size=24,hjust = 0.5)) + xlim(0,35) +
  ylab("TNFRSF1A/TNFRSF1B, day 7")  + xlab("H3N2.nAb.FCd28 foldchange")+
  #scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1)+ 
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))
# ggsave(filename = "Images/TNFRSFratio_vs_H3N2nAb_allTogether.pdf", device="pdf")



# TNFR1 vs nAb responses
a <- cor.test(phenotypeMatrixYoung$`H1N1.nAb.FCd28`, phenotypeMatrixYoung$TNFRSF1Ahihid7, use="complete")
b <- cor.test(phenotypeMatrixElderly$`H1N1.nAb.FCd28`, phenotypeMatrixElderly$TNFRSF1Ahihid7, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.4,  y=0.11, hjust=0, gp=gpar(col="orange3", fontsize=24)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.4,  y=0.05, hjust=0, gp=gpar(col="purple", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`H1N1.nAb.FCd28`,y=`TNFRSF1Ahihid7`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  
  ggtitle("TNFRSF1A vs H1N1 HAI Ab") + theme(plot.title = element_text(size=32,hjust = 0.5)) + xlim(0,10) +
  ylab("TNFRSF1A, day 7")  + xlab("H1N1.nAb.FCd28")+
  theme(axis.text = element_text(size=24,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/TNFRSF1A_vs_H1N1nAb_byAge.pdf", device="pdf")


a <- cor.test(phenotypeMatrix$`H1N1.nAb.FCd28`, phenotypeMatrix$TNFRSF1Ahihid7, use="complete")
annotationInfo <- paste0("r = ", round(a$estimate,3), "   \n", "P = ", formatC(a$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.6,  y=0.15, hjust=0, gp=gpar(col="black", fontsize=32)))
ggplot(data=phenotypeMatrix, aes(x=`H1N1.nAb.FCd28`,y=`TNFRSF1Ahihid7`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1) + 
  geom_point(size=6, pch=21, fill="black") + 
  ggtitle("TNFRSF1A vs H1N1 HAI Ab") + theme(plot.title = element_text(size=24,hjust = 0.5)) + xlim(0,10) +
  ylab("Log counts of TNFRSF1A, day 7")  + xlab("H1N1.nAb.FCd28 foldchange")+
  #scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1)+ 
  theme(axis.text = element_text(size=24,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))
# ggsave(filename = "Images/TNFRSF1A_vs_H1N1nAb_allTogether.pdf", device="pdf", width=7, height=7)


a <- cor.test(phenotypeMatrix$`H1N1.nAb.FCd28`, phenotypeMatrix$TNFRSF1Ahihid7, use="complete")
annotationInfo <- paste0("r = ", round(a$estimate,3), "   \n", "P = ", formatC(a$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.6,  y=0.15, hjust=0, gp=gpar(col="black", fontsize=32)))
ggplot(data=phenotypeMatrix, aes(x=`H1N1.nAb.FCd28`,y=`TNFRSF1Ahihid7`,fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(fill=`Visit`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=8, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=8, pch=21) +  
  ggtitle("TNFRSF1A vs H1N1 HAI Ab") + theme(plot.title = element_text(size=24,hjust = 0.5)) + xlim(0,10) +
  ylab("Log counts of TNFRSF1A, day 7")  + xlab("H1N1.nAb.FCd28 foldchange")+
  scale_fill_manual(values=c('purple','black','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1)+ 
  theme(axis.text = element_text(size=24,hjust = 0.5), axis.title = element_text(size=28,hjust = 0.5))
# ggsave(filename = "Images/TNFRSF1A_vs_H1N1nAb_allTogether-colorCohorts.pdf", device="pdf", width=7, height=7)




a <- cor.test(phenotypeMatrix$`H3N2.nAb.FCd28`, phenotypeMatrix$TNFRSF1Ahihid7, use="complete")
annotationInfo <- paste0("r = ", round(a$estimate,3), "   \n", "P = ", formatC(a$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.6,  y=0.15, hjust=0, gp=gpar(col="black", fontsize=32)))
ggplot(data=phenotypeMatrix, aes(x=`H3N2.nAb.FCd28`,y=`TNFRSF1Ahihid7`,fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(fill=`Visit`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=8, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=8, pch=21) +  
  ggtitle("TNFRSF1A vs H3N2 HAI Ab") + theme(plot.title = element_text(size=24,hjust = 0.5)) + xlim(0,35) +
  ylab("Log counts of TNFRSF1A, day 7")  + xlab("H3N2.nAb.FCd28 foldchange")+
  scale_fill_manual(values=c('purple','black','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + 
  annotation_custom(my_grob1)+ 
  theme(axis.text = element_text(size=24,hjust = 0.5), axis.title = element_text(size=28,hjust = 0.5))
# ggsave(filename = "Images/TNFRSF1A_vs_H3N2nAb_allTogether-colorCohorts.pdf", device="pdf", width=7, height=7)

