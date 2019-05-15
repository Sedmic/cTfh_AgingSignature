library("grid")
library("genefilter")
library("ggplot2")
library("gplots")
library("viridis")
library("DESeq2")
library("reshape2")
library("ggpubr")
library("ggcorrplot")
library("gridExtra")
library("RColorBrewer")
library("scales")



##  bring in the Luminex data

luminex <- read.csv(file="../../../Luminex/Year4_Luminex/LuminexSummary.csv",stringsAsFactors = FALSE)
luminexYoung <- subset(luminex, AgeCategory=="Young" & Day=="0")
luminexElderly <- subset(luminex, AgeCategory=="Elderly"& Day=="0")

qCXCL11 <-   ggplot(data=luminex, aes(x=Day,y=CXCL11)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) +  scale_fill_brewer(palette="Dark2") + #scale_fill_manual(values=c("#FDBF6F", "#FF7F00")) +
  ggtitle("CXCL11") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("CXCL11 (pg/mL)") + ylim(c(0,3000))
qFractaline <- ggplot(data=luminex, aes(x=Day,y=Fractaline)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("Fractaline") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("Fractaline (pg/mL)") + ylim(c(0,500))
qMIP3a <- ggplot(data=luminex, aes(x=Day,y=MIP3a)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("MIP3a") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("MIP3a (pg/mL)") + ylim(c(0,100))
qIL7 <- ggplot(data=luminex, aes(x=Day,y=IL7)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("IL7") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("IL7 (pg/mL)")  + ylim(c(0,80))
qIL8 <- ggplot(data=luminex, aes(x=Day,y=IL8)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("IL8") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("IL8 (pg/mL)") + ylim(c(0,260))
qMIP1b <- ggplot(data=luminex, aes(x=Day,y=MIP1b)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("MIP1b") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("MIP1b (pg/mL)") + ylim(c(0,35))
qTNF <- ggplot(data=luminex, aes(x=Day,y=TNFa)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("TNF") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("TNF (pg/mL)") + ylim(c(0,27))
qLeptin <- ggplot(data=luminex, aes(x=Day,y=Leptin)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("Leptin") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("Leptin (pg/mL)") + ylim(c(0,8200))
qHGF <- ggplot(data=luminex, aes(x=Day,y=HGF)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("HGF") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("HGF (pg/mL)") + ylim(c(0,160))
qAdiponectin <- ggplot(data=luminex, aes(x=Day,y=Adiponectin)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("Adiponectin") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("Adiponectin (pg/mL)") + ylim(c(0,7000000)) 
qMCP1 <- ggplot(data=luminex, aes(x=Day,y=MCP1)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("MCP1") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("MCP1 (pg/mL)") + ylim(c(0,180))
qResistin <- ggplot(data=luminex, aes(x=Day,y=Resistin)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Young",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("Resistin") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("Resistin (pg/mL)") + ylim(c(0,9000))


a <- ggarrange(qCXCL11, qFractaline, qMIP3a, qIL7, qIL8, qMIP1b, qTNF, qLeptin, qHGF, qAdiponectin, qMCP1, qResistin, common.legend = TRUE, legend="right")
ggexport(a, filename = "../../../Luminex/Year4_Luminex/Images/AllChemokinesByDay_Young.pdf")


qCXCL11 <- ggplot(data=luminex, aes(x=Day,y=CXCL11)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("CXCL11") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("CXCL11 (pg/mL)") + ylim(c(0,3000))
qFractaline <- ggplot(data=luminex, aes(x=Day,y=Fractaline)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("Fractaline") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("Fractaline (pg/mL)") + ylim(c(0,500))
qMIP3a <- ggplot(data=luminex, aes(x=Day,y=MIP3a)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("MIP3a") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("MIP3a (pg/mL)") + ylim(c(0,100))
qIL7 <- ggplot(data=luminex, aes(x=Day,y=IL7)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("IL7") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("IL7 (pg/mL)") + ylim(c(0,80))
qIL8 <- ggplot(data=luminex, aes(x=Day,y=IL8)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("IL8") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("IL8 (pg/mL)") + ylim(c(0,260))
qMIP1b <- ggplot(data=luminex, aes(x=Day,y=MIP1b)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("MIP1b") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("MIP1b (pg/mL)") + ylim(c(0,35))
qTNF <- ggplot(data=luminex, aes(x=Day,y=TNFa)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("TNF") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("TNF (pg/mL)") + ylim(c(0,27))
qLeptin <- ggplot(data=luminex, aes(x=Day,y=Leptin)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("Leptin") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("Leptin (pg/mL)") + ylim(c(0,8200))
qHGF <- ggplot(data=luminex, aes(x=Day,y=HGF)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("HGF") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("HGF (pg/mL)") + ylim(c(0,160))
qAdiponectin <- ggplot(data=luminex, aes(x=Day,y=Adiponectin)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("Adiponectin") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("Adiponectin (pg/mL)") + ylim(c(0,7000000)) 
qMCP1 <- ggplot(data=luminex, aes(x=Day,y=MCP1)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("MCP1") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("MCP1 (pg/mL)") + ylim(c(0,180))
qResistin <- ggplot(data=luminex, aes(x=Day,y=Resistin)) + geom_violin(data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity"), aes(fill=Day)) + theme_bw() +
  geom_boxplot(width=0.1, data=subset(luminex,Day==c(0,7) & AgeCategory=="Elderly",stat="identity")) + scale_fill_brewer(palette="Dark2") +
  ggtitle("Resistin") + theme(plot.title = element_text(size=16,hjust = 0.5)) + ylab("Resistin (pg/mL)") + ylim(c(0,9000))

a <- ggarrange(qCXCL11, qFractaline, qMIP3a, qIL7, qIL8, qMIP1b, qTNF, qLeptin, qHGF, qAdiponectin, qMCP1, qResistin, common.legend = TRUE, legend="right")
ggexport(a, filename = "../../../Luminex/Year4_Luminex/Images/AllChemokinesByDay_Elderly.pdf")


# CXCL11 vs TNF
a <- cor.test(luminexYoung$CXCL11, luminexYoung$TNFa, use="complete")
b <- cor.test(luminexElderly$CXCL11, luminexElderly$TNFa, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value, 2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=20)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=20)))
ggplot(data=subset(luminex, Day==0), aes(x=`CXCL11`,y=TNFa, fill=`AgeCategory`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`AgeCategory`)) + 
  geom_point(data=subset(luminex,AgeCategory=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(luminex,AgeCategory=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("TNF vs CXCL11 at day 0") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("TNF (pg/mL)")  + xlab("CXCL11 (pg/mL)")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/TNFvsCXCL11_YandE_scatter.png", device="png")
# ggsave(filename = "../../Analysis/Antibody data/Images/TNFvsCXCL11_YandE_scatter.pdf", device="pdf")

# MIP1b vs TNF
a <- cor.test(luminexYoung$MIP1b, luminexYoung$TNFa, use="complete")
b <- cor.test(luminexElderly$MIP1b, luminexElderly$TNFa, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value, 2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18)))
ggplot(data=subset(luminex, Day==0), aes(x=`MIP1b`,y=TNFa, fill=`AgeCategory`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`AgeCategory`)) + 
  geom_point(data=subset(luminex,AgeCategory=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(luminex,AgeCategory=="Young",stat="identity"), size=6, pch=21) +  
  ggtitle("TNF vs MIP1b at day 0") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("TNF (pg/mL)")  + xlab("MIP1b (pg/mL)")+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/TNFvsMIP1b_YandE_scatter.png", device="png")
# ggsave(filename = "../../Analysis/Antibody data/Images/TNFvsMIP1b_YandE_scatter.pdf", device="pdf")


# **************************** Antibody plus Luminex data ****************************************

phenotypeMatrix <- read.csv(file="../../Analysis/BAA_yr4_AllMergedData_FC.csv"); phenotypeMatrix <- phenotypeMatrix[-c(62:73),]; rownames(phenotypeMatrix) <- phenotypeMatrix$Subject; 
phenotypeMatrix$Subject <- NULL
phenotypeMatrixYoung <- subset(phenotypeMatrix,Identifier=="Young",stat="identity")
phenotypeMatrixElderly <- subset(phenotypeMatrix,Identifier=="Elderly",stat="identity")


#  ASC vs TNF 
a <- cor.test(phenotypeMatrixYoung$ASC.freqLive, phenotypeMatrixYoung$TNFa, use="complete")
b <- cor.test(phenotypeMatrixElderly$ASC.freqLive, phenotypeMatrixElderly$TNFa, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.9, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.63,  y=0.85, hjust=0, gp=gpar(col="purple", fontsize=18)))
ggplot(data=phenotypeMatrix, aes(x=`TNFa`,y=`ASC.freqLive`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("Plasmablast response vs TNF d0") + theme(plot.title = element_text(size=28,hjust = 0.5)) + xlab("TNF (ng/mL)")  + ylab("Plasmablast frequency foldchange")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/TNF_vs_ASC_byAge.pdf", device="pdf")

summary(lm(data = phenotypeMatrixYoung, ASC.freqLive ~ cTfh_ICOShi.CD38hi...Freq..of + TNFa))
summary(lm(data = phenotypeMatrixElderly, ASC.freqLive ~ cTfh_ICOShi.CD38hi...Freq..of + TNFa))

#  ASC vs Tfh response 
a <- cor.test(phenotypeMatrixYoung$ASC.freqLive, phenotypeMatrixYoung$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
b <- cor.test(phenotypeMatrixElderly$ASC.freqLive, phenotypeMatrixElderly$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value,format="e", digits=1))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.55,  y=0.12, hjust=0, gp=gpar(col="orange3", fontsize=24)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.55,  y=0.05, hjust=0, gp=gpar(col="purple", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`ASC.freqLive`,y=`cTfh_ICOShi.CD38hi...Freq..of`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("ICOS+CD38+ cTfh vs Plasmablast") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("ICOShiCD38hi cTfh foldchange")  + xlab("Plasmablast frequency foldchange")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/HiHiresponse_vs_ASC_byAge.pdf", device="pdf", width=8, height=6)

summary(lm(data = phenotypeMatrixYoung, ASC.freqLive ~ cTfh_ICOShi.CD38hi...Freq..of + TNFa))
summary(lm(data = phenotypeMatrixElderly, ASC.freqLive ~ cTfh_ICOShi.CD38hi...Freq..of + TNFa))

#  ASC vs Tfh response - OUTLIERS REMOVED
phenotypeMatrix <- phenotypeMatrix[-which(phenotypeMatrix$ASC.freqLive == max(phenotypeMatrix$ASC.freqLive, na.rm=T)),]
phenotypeMatrix <- phenotypeMatrix[-which(phenotypeMatrix$ASC.freqLive == max(phenotypeMatrix$ASC.freqLive, na.rm=T)),]
phenotypeMatrixYoung <- phenotypeMatrixYoung[-which(phenotypeMatrixYoung$ASC.freqLive==max(phenotypeMatrixYoung$ASC.freqLive,na.rm=T)),]
phenotypeMatrixElderly <- phenotypeMatrixElderly[-which(phenotypeMatrixElderly$ASC.freqLive==max(phenotypeMatrixElderly$ASC.freqLive,na.rm=T)),]
a <- cor.test(phenotypeMatrixYoung$ASC.freqLive, phenotypeMatrixYoung$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
b <- cor.test(phenotypeMatrixElderly$ASC.freqLive, phenotypeMatrixElderly$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.55,  y=0.12, hjust=0, gp=gpar(col="orange3", fontsize=24)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.55,  y=0.05, hjust=0, gp=gpar(col="purple", fontsize=24)))
ggplot(data=phenotypeMatrix, aes(x=`ASC.freqLive`,y=`cTfh_ICOShi.CD38hi...Freq..of`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("ICOS+CD38+ cTfh vs Plasmablasts") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("ICOShiCD38hi cTfh foldchange")  + xlab("Plasmablast frequency foldchange")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/HiHiresponse_vs_ASC_byAge_OUTLIERREMOVED.pdf", device="pdf", width=8, height=6)
phenotypeMatrix <- read.csv(file="../../Analysis/BAA_yr4_AllMergedData_FC.csv"); phenotypeMatrix <- phenotypeMatrix[-c(62:73),]; rownames(phenotypeMatrix) <- phenotypeMatrix$Subject; 
phenotypeMatrix$Subject <- NULL
phenotypeMatrixYoung <- subset(phenotypeMatrix,Identifier=="Young",stat="identity")
phenotypeMatrixElderly <- subset(phenotypeMatrix,Identifier=="Elderly",stat="identity")


summary(lm(data = phenotypeMatrixYoung, ASC.freqLive ~ cTfh_ICOShi.CD38hi...Freq..of + TNFa))
summary(lm(data = phenotypeMatrixElderly, ASC.freqLive ~ cTfh_ICOShi.CD38hi...Freq..of + TNFa))


#  nAb vs Tfh response
a <- cor.test(phenotypeMatrixYoung$H1N1.nAb.FCd28, phenotypeMatrixYoung$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
b <- cor.test(phenotypeMatrixElderly$H1N1.nAb.FCd28, phenotypeMatrixElderly$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.6,  y=0.10, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.6,  y=0.05, hjust=0, gp=gpar(col="purple", fontsize=18)))
ggplot(data=phenotypeMatrix, aes(x=`H1N1.nAb.FCd28`,y=`cTfh_ICOShi.CD38hi...Freq..of`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))+
  ggtitle("ICOS+CD38+ cTfh vs H1N1 HAI") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("ICOShiCD38hi cTfh foldchange")  + xlab("H1N1.nAb.FCd28")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/HiHiresponse_vs_H1N1_nAb_byAge.pdf", device="pdf")

a <- cor.test(phenotypeMatrixYoung$H3N2.nAb.FCd28, phenotypeMatrixYoung$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
b <- cor.test(phenotypeMatrixElderly$H3N2.nAb.FCd28, phenotypeMatrixElderly$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
annotationInfo2 <- paste0("r= ", round(b$estimate,3), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.55,  y=0.10, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.55,  y=0.05, hjust=0, gp=gpar(col="purple", fontsize=18)))
ggplot(data=phenotypeMatrix, aes(x=`H3N2.nAb.FCd28`,y=`cTfh_ICOShi.CD38hi...Freq..of`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5))+
  ggtitle("ICOS+CD38+ cTfh vs H3N2 HAI") + theme(plot.title = element_text(size=24,hjust = 0.5)) + ylab("ICOShiCD38hi cTfh foldchange")  + xlab("H3N2.nAb.FCd28")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/HiHiresponse_vs_H3N2_nAb_byAge.pdf", device="pdf")



#  ASC vs Tfh response - YOUNG ONLY  not for manuscript 
phenotypeMatrix <- phenotypeMatrix[-which(phenotypeMatrix$ASC.freqLive == max(phenotypeMatrix$ASC.freqLive, na.rm=T)),]
phenotypeMatrix <- phenotypeMatrix[-which(phenotypeMatrix$ASC.freqLive == max(phenotypeMatrix$ASC.freqLive, na.rm=T)),]
phenotypeMatrixYoung <- phenotypeMatrixYoung[-which(phenotypeMatrixYoung$ASC.freqLive==max(phenotypeMatrixYoung$ASC.freqLive,na.rm=T)),]
phenotypeMatrixElderly <- phenotypeMatrixElderly[-which(phenotypeMatrixElderly$ASC.freqLive==max(phenotypeMatrixElderly$ASC.freqLive,na.rm=T)),]
a <- cor.test(phenotypeMatrixYoung$ASC.freqLive, phenotypeMatrixYoung$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
b <- cor.test(phenotypeMatrixElderly$ASC.freqLive, phenotypeMatrixElderly$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", formatC(a$p.value, format="e", digits=1))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.9, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.62,  y=0.05, hjust=0, gp=gpar(col="purple", fontsize=18)))
phenotypeMatrixYoungReduced <- phenotypeMatrixYoung[sample.int(nrow(phenotypeMatrixYoung), 12, replace=F),]
ggplot(data=phenotypeMatrixYoungReduced, aes(x=`ASC.freqLive`,y=`cTfh_ICOShi.CD38hi...Freq..of`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
#  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrixYoungReduced,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("ICOS+CD38+ cTfh vs Plasmablasts") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("ICOShiCD38hi cTfh foldchange")  + xlab("Plasmablast frequency foldchange")+
  scale_fill_manual(values=c('#E69F00','#E69F00')) + scale_color_manual(values=c('#E69F00', '#E69F00')) + annotation_custom(my_grob1) #+ annotation_custom(my_grob2)
# ggsave(filename = "Images/HiHiresponse_vs_ASC_YoungOnly_OUTLIERREMOVED.png", device="png")
phenotypeMatrix <- read.csv(file="../../Analysis/BAA_yr4_AllMergedData_FC.csv"); phenotypeMatrix <- phenotypeMatrix[-c(62:73),]; rownames(phenotypeMatrix) <- phenotypeMatrix$Subject; 
phenotypeMatrix$Subject <- NULL
phenotypeMatrixYoung <- subset(phenotypeMatrix,Identifier=="Young",stat="identity")
phenotypeMatrixElderly <- subset(phenotypeMatrix,Identifier=="Elderly",stat="identity")


#  TNF vs Tfh response 
a <- cor.test(phenotypeMatrixYoung$TNFa, phenotypeMatrixYoung$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
b <- cor.test(phenotypeMatrixElderly$TNFa, phenotypeMatrixElderly$cTfh_ICOShi.CD38hi...Freq..of, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";   ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";   ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18)))
ggplot(data=phenotypeMatrix, aes(x=`TNFa`,y=`cTfh_ICOShi.CD38hi...Freq..of`, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  #ylim(-1,1) + xlim(0,3)+
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("ICOS+CD38+ cTfh FC vs TNFa d0") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("ICOShiCD38hi cTfh foldchange")  + xlab("TNF (ng/mL)")+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00')) + annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "Images/HiHiresponse_vs_TNF_byAge.pdf", device="pdf")



# H1N1 IgM d28 vs TNF
a <- cor.test(phenotypeMatrixYoung$H1N1.IgM.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H1N1.IgM.FCd28, phenotypeMatrixElderly$TNFa, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2),";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18)))
ggplot(data=phenotypeMatrix[-1,], aes(x=`H1N1.IgM.FCd28`,y=TNFa, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`Identifier`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +  
  ggtitle("TNF vs H1N1 IgM FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("TNF (pg/mL)") + xlab("Fold change of day 28 vs 0") + 
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H1N1_IgM_FCd28_YandE_scatter.png", device="png")
# ggsave(filename = "../../Analysis/Antibody data/Images/H1N1_IgM_FCd28_YandE_scatter.pdf", device="pdf")
cor.test(phenotypeMatrixYoung$H1N1.IgM.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
cor.test(phenotypeMatrixElderly$H1N1.IgM.FCd28, phenotypeMatrixElderly$TNFa, use="complete")



# H3N2 IgM d28 vs TNF
a <- cor.test(phenotypeMatrixYoung$H3N2.IgM.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H3N2.IgM.FCd28, phenotypeMatrixElderly$TNFa, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", formatC(b$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18)))
ggplot(data=phenotypeMatrix[-1,], aes(x=`H3N2.IgM.FCd28`,y=TNFa, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Young",stat="identity"), size=6, pch=21) +   
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("TNF vs H3N2 IgM FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("TNF (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H3N2_IgM_FCd28_YandE_scatter.png", device="png")
# ggsave(filename = "../../Analysis/Antibody data/Images/H3N2_IgM_FCd28_YandE_scatter.pdf", device="pdf")
cor.test(phenotypeMatrixYoung$H3N2.IgM.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
cor.test(phenotypeMatrixElderly$H3N2.IgM.FCd28, phenotypeMatrixElderly$TNFa, use="complete")




# H1N1 IgG d28 vs TNF
a <- cor.test(phenotypeMatrixYoung$H1N1.IgG.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H1N1.IgG.FCd28, phenotypeMatrixElderly$TNFa, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2),  ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H1N1.IgG.FCd28`,y=TNFa, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`Identifier`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("TNF vs H1N1 IgG FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("TNF (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H1N1_IgG_FCd28_YandE_scatter.pdf", device="pdf")
cor.test(phenotypeMatrixYoung$H1N1.IgG.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
cor.test(phenotypeMatrixElderly$H1N1.IgG.FCd28, phenotypeMatrixElderly$TNFa, use="complete")


# H1N1 IgG d28 vs MIP1b
a <- cor.test(phenotypeMatrixYoung$H1N1.IgG.FCd28[-1], phenotypeMatrixYoung$MIP1b[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H1N1.IgG.FCd28, phenotypeMatrixElderly$MIP1b, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2),  ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H1N1.IgG.FCd28`,y=MIP1b, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`Identifier`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("MIP1b vs H1N1 IgG FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("MIP1b (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H1N1_IgG_FCd28_Mip1b_YandE_scatter.pdf", device="pdf")
cor.test(phenotypeMatrixYoung$H1N1.IgG.FCd28[-1], phenotypeMatrixYoung$MIP1b[-1], use="complete")
cor.test(phenotypeMatrixElderly$H1N1.IgG.FCd28, phenotypeMatrixElderly$MIP1b, use="complete")


# H1N1 IgG d28 vs CXCL11
a <- cor.test(phenotypeMatrixYoung$H1N1.IgG.FCd28[-1], phenotypeMatrixYoung$CXCL11[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H1N1.IgG.FCd28, phenotypeMatrixElderly$CXCL11, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2),  ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H1N1.IgG.FCd28`,y=CXCL11, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`Identifier`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("CXCL11 vs H1N1 IgG FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("CXCL11 (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H1N1_IgG_FCd28_CXCL11_YandE_scatter.pdf", device="pdf")
cor.test(phenotypeMatrixYoung$H1N1.IgG.FCd28[-1], phenotypeMatrixYoung$CXCL11[-1], use="complete")
cor.test(phenotypeMatrixElderly$H1N1.IgG.FCd28, phenotypeMatrixElderly$CXCL11, use="complete")



# H3N2 IgG d28 vs TNF
a <- cor.test(phenotypeMatrixYoung$H3N2.IgG.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H3N2.IgG.FCd28, phenotypeMatrixElderly$TNFa, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";  ","p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H3N2.IgG.FCd28`,y=TNFa, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +   
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("TNF vs H3N2 IgG FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("TNF (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H3N2_IgG_FCd28_YandE_scatter.png", device="png")
# ggsave(filename = "../../Analysis/Antibody data/Images/H3N2_IgG_FCd28_YandE_scatter.pdf", device="pdf")
cor.test(phenotypeMatrixYoung$H3N2.IgG.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
cor.test(phenotypeMatrixElderly$H3N2.IgG.FCd28, phenotypeMatrixElderly$TNFa, use="complete")



# H3N2 IgG d28 vs MIP1b
a <- cor.test(phenotypeMatrixYoung$H3N2.IgG.FCd28[-1], phenotypeMatrixYoung$MIP1b[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H3N2.IgG.FCd28, phenotypeMatrixElderly$MIP1b, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2),  ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H3N2.IgG.FCd28`,y=MIP1b, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`Identifier`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("MIP1b vs H3N2 IgG FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("MIP1b (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H3N2_IgG_FCd28_Mip1b_YandE_scatter.pdf", device="pdf")
cor.test(phenotypeMatrixYoung$H3N2.IgG.FCd28[-1], phenotypeMatrixYoung$MIP1b[-1], use="complete")
cor.test(phenotypeMatrixElderly$H3N2.IgG.FCd28, phenotypeMatrixElderly$MIP1b, use="complete")


# H3N2 IgG d28 vs CXCL11
a <- cor.test(phenotypeMatrixYoung$H3N2.IgG.FCd28[-1], phenotypeMatrixYoung$CXCL11[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H3N2.IgG.FCd28, phenotypeMatrixElderly$CXCL11, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2),  ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H3N2.IgG.FCd28`,y=CXCL11, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`Identifier`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("CXCL11 vs H3N2 IgG FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("CXCL11 (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H3N2_IgG_FCd28_CXCL11_YandE_scatter.pdf", device="pdf")
cor.test(phenotypeMatrixYoung$H3N2.IgG.FCd28[-1], phenotypeMatrixYoung$CXCL11[-1], use="complete")
cor.test(phenotypeMatrixElderly$H3N2.IgG.FCd28, phenotypeMatrixElderly$CXCL11, use="complete")




# H1N1 nAb d28 vs TNF
a <- cor.test(phenotypeMatrixYoung$H1N1.nAb.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H1N1.nAb.FCd28, phenotypeMatrixElderly$TNFa, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H1N1.nAb.FCd28`,y=TNFa, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`Identifier`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("TNF vs H1N1 nAb FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("TNF (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H1N1_nAb_FCd28_YandE_scatter.png", device="png")
# ggsave(filename = "../../Analysis/Antibody data/Images/H1N1_nAb_FCd28_YandE_scatter.pdf", device="pdf")
cor(phenotypeMatrixYoung$H1N1.nAb.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
cor(phenotypeMatrixElderly$H1N1.nAb.FCd28, phenotypeMatrixElderly$TNFa, use="complete")

# H3N2 nAb d28 vs TNF
a <- cor.test(phenotypeMatrixYoung$H3N2.nAb.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H3N2.nAb.FCd28, phenotypeMatrixElderly$TNFa, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H3N2.nAb.FCd28`,y=TNFa, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("TNF vs H3N2 nAb FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("TNF (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H3N2_nAb_FCd28_YandE_scatter.pdf", device="pdf")
cor(phenotypeMatrixYoung$H3N2.nAb.FCd28[-1], phenotypeMatrixYoung$TNFa[-1], use="complete")
cor(phenotypeMatrixElderly$H3N2.nAb.FCd28, phenotypeMatrixElderly$TNFa, use="complete")


# H1N1 nAb d28 vs MIP1b
a <- cor.test(phenotypeMatrixYoung$H1N1.nAb.FCd28[-1], phenotypeMatrixYoung$MIP1b[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H1N1.nAb.FCd28, phenotypeMatrixElderly$MIP1b, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H1N1.nAb.FCd28`,y=MIP1b, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`Identifier`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("Mip1b vs H1N1 nAb FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("Mip1b (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H1N1_nAb_FCd28_mip1b_YandE_scatter.pdf", device="pdf")
cor(phenotypeMatrixYoung$H1N1.nAb.FCd28[-1], phenotypeMatrixYoung$MIP1b[-1], use="complete")
cor(phenotypeMatrixElderly$H1N1.nAb.FCd28, phenotypeMatrixElderly$MIP1b, use="complete")

# H3N2 nAb d28 vs MIP1b
a <- cor.test(phenotypeMatrixYoung$H3N2.nAb.FCd28[-1], phenotypeMatrixYoung$MIP1b[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H3N2.nAb.FCd28, phenotypeMatrixElderly$MIP1b, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H3N2.nAb.FCd28`,y=MIP1b, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("Mip1b vs H3N2 nAb FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("Mip1b (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H3N2_nAb_FCd28_mip1b_YandE_scatter.pdf", device="pdf")
cor(phenotypeMatrixYoung$H3N2.nAb.FCd28[-1], phenotypeMatrixYoung$MIP1b[-1], use="complete")
cor(phenotypeMatrixElderly$H3N2.nAb.FCd28, phenotypeMatrixElderly$MIP1b, use="complete")



# H1N1 nAb d28 vs CXCL11
a <- cor.test(phenotypeMatrixYoung$H1N1.nAb.FCd28[-1], phenotypeMatrixYoung$CXCL11[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H1N1.nAb.FCd28, phenotypeMatrixElderly$CXCL11, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H1N1.nAb.FCd28`,y=CXCL11, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +  # exclude row 1 because d28 sample not available for that person
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=`Identifier`)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("CXCL11 vs H1N1 nAb FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("CXCL11 (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H1N1_nAb_FCd28_CXCL11_YandE_scatter.pdf", device="pdf")
cor(phenotypeMatrixYoung$H1N1.nAb.FCd28[-1], phenotypeMatrixYoung$CXCL11[-1], use="complete")
cor(phenotypeMatrixElderly$H1N1.nAb.FCd28, phenotypeMatrixElderly$CXCL11, use="complete")

# H3N2 nAb d28 vs CXCL11
a <- cor.test(phenotypeMatrixYoung$H3N2.nAb.FCd28[-1], phenotypeMatrixYoung$CXCL11[-1], use="complete")
b <- cor.test(phenotypeMatrixElderly$H3N2.nAb.FCd28, phenotypeMatrixElderly$CXCL11, use="complete")
annotationInfo <- paste0("r= ", round(a$estimate,2), ";  ", "p = ", round(a$p.value,2))
annotationInfo2 <- paste0("r= ", round(b$estimate,2), ";  ", "p = ", round(b$p.value,2))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="orange3", fontsize=18)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.65,  y=0.88, hjust=0, gp=gpar(col="purple", fontsize=18))) 
ggplot(data=phenotypeMatrix[-1,], aes(x=`H3N2.nAb.FCd28`,y=CXCL11, fill=`Identifier`)) + theme_bw() +  theme(legend.position = "none") +
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1, aes(color=Identifier)) + 
  geom_point(data=subset(phenotypeMatrix,Identifier=="Elderly",stat="identity"), size=6, pch=21) + 
  geom_point(data=subset(phenotypeMatrix[-1,],Identifier=="Young",stat="identity"), size=6, pch=21) +  
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=24,hjust = 0.5))+
  ggtitle("CXCL11 vs H3N2 nAb FCd28") + theme(plot.title = element_text(size=28,hjust = 0.5)) + ylab("CXCL11 (pg/mL)")  + xlab("Fold change of day 28 vs 0") + 
  scale_fill_manual(values=c('purple','#E69F00')) + scale_color_manual(values=c('purple', '#E69F00'))+ annotation_custom(my_grob1) + annotation_custom(my_grob2)
# ggsave(filename = "../../Analysis/Antibody data/Images/H3N2_nAb_FCd28_CXCL11_YandE_scatter.pdf", device="pdf")
cor(phenotypeMatrixYoung$H3N2.nAb.FCd28[-1], phenotypeMatrixYoung$CXCL11[-1], use="complete")
cor(phenotypeMatrixElderly$H3N2.nAb.FCd28, phenotypeMatrixElderly$CXCL11, use="complete")



corrY <- cor(phenotypeMatrixYoung[,c(4:6, 10, 17:26)],use="complete") 
p.matY <- cor_pmat(phenotypeMatrixYoung[,c(4:6, 10, 17:26)],use="complete")
a<-ggcorrplot(corrY, hc.order=FALSE, type="full",outline.col="white",lab=FALSE, p.mat=p.matY, insig="blank", tl.cex=19, colors=inferno(3) )
a
# ggsave(plot=a, filename="../../Analysis/Antibody data/Images/Luminex_and_antibody_correlationsYoung.pdf", device="pdf")

corrE <- cor(phenotypeMatrixElderly[,c(4:6, 10, 17:26)],use="complete") 
p.matE <- cor_pmat(phenotypeMatrixElderly[,c(4:6, 10, 17:26)],use="complete")
b <- ggcorrplot(corrE, hc.order=FALSE, type="full",outline.col="white",lab=FALSE, p.mat=p.matE, insig="blank", tl.cex=19, colors=inferno(3) )
b
# ggsave(plot=b, filename="../../Analysis/Antibody data/Images/Luminex_and_antibody_correlationsElderly.pdf", device="pdf")

grid.arrange(a,b, nrow=1)


