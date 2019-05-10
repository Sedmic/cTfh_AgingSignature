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




## *****************************       GSEA results: hallmark at baseline    ***********************************************


Yd0HallmarkPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_HALLMARK.GseaPreranked.1530926563560/gsea_report_for_na_pos_1530926563560.xls", sep="\t")
Yd0HallmarkNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_HALLMARK.GseaPreranked.1530926563560/gsea_report_for_na_neg_1530926563560.xls", sep="\t")
Yd0Hallmark <- rbind(Yd0HallmarkPos, Yd0HallmarkNeg)
Ed0HallmarkPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_HALLMARK.GseaPreranked.1531837435104/gsea_report_for_na_pos_1531837435104.xls", sep="\t")
Ed0HallmarkNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_HALLMARK.GseaPreranked.1531837435104/gsea_report_for_na_neg_1531837435104.xls", sep="\t")
Ed0Hallmark <- rbind(Ed0HallmarkPos, Ed0HallmarkNeg)


# Ox-Phos  
YoungHiHiv1Oxphos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_lolo_v1_HALLMARK.GseaPreranked.1530926563560/HALLMARK_OXIDATIVE_PHOSPHORYLATION.xls", sep="\t")
ElderlyHiHiv1Oxphos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/E_hihi_v_LoLo_v1_HALLMARK.GseaPreranked.1531837435104//HALLMARK_OXIDATIVE_PHOSPHORYLATION.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(Yd0Hallmark$NES[grep("OXIDATIVE",Yd0Hallmark$NAME)],2), "\n", "FDR: ", formatC(Yd0Hallmark$FDR.q.val[grep("OXIDATIVE",Yd0Hallmark$NAME)], format = "e", digits = 1))
annotationInfo2 <- paste0("\nNES: ", round(Ed0Hallmark$NES[grep("OXIDATIVE",Ed0Hallmark$NAME)],2), "\n", "FDR: ", formatC(Ed0Hallmark$FDR.q.val[grep("OXIDATIVE",Ed0Hallmark$NAME)], format = "e", digits = 1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.75,  y=0.82, hjust=0, gp=gpar(col="orange3", fontsize=12)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.75,  y=0.65, hjust=0, gp=gpar(col="purple", fontsize=12)))
ggplot(data=YoungHiHiv1Oxphos, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) +  theme_bw() +
  geom_line(color="orange3", size=1) + 
  geom_line(data=ElderlyHiHiv1Oxphos, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`), color="purple")+
  geom_rug(sides="t", size=0.75, alpha=0.5, color="orange3") +  geom_rug(data=ElderlyHiHiv1Oxphos, sides="b", size=0.75, alpha=0.5, color="purple") +  
  ggtitle("Oxidative Phosphorylation - Day 0") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0))+
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/Images/YandE_hihi-v-lolo_d0_oxphos.pdf", device="pdf", height=3.5, width=5)




## *****************************       IPA: AllAges canonical pathways    ***********************************************


AllAgesHivLov1_can <- read.csv("IPA/IPAresults/20181115_AllAges_HiHi-v-LoLo_v1_IPAcanonicalPathways.csv", stringsAsFactors = F )
AllAgesHivLov1_can <- AllAgesHivLov1_can[-grep("#NUM!", AllAgesHivLov1_can$z.score),]  # suppress the rows with #NUM! in the z-score column
AllAgesHivLov1_can$z.score <- as.numeric(AllAgesHivLov1_can$z.score)                  # convert to numeric
AllAgesHivLov1_can <- AllAgesHivLov1_can[which(abs(AllAgesHivLov1_can$z.score)>=2),]  # retain only rows with z score >2 or <-2
AllAgesHivLov1_can <- AllAgesHivLov1_can[which(abs(AllAgesHivLov1_can$X.log.p.value.)>=1.3),]  # retain only rows with -log(pvalue) >1.3 

subsetMetab <- AllAgesHivLov1_can[grep("Glycolysis|Oxidative Phosphorylation|Gluconeogenesis|Pentose", AllAgesHivLov1_can$Ingenuity.Canonical.Pathways),]  # suppress the rows with #NUM! in the z-score column
subsetMetab$Ingenuity.Canonical.Pathways <- factor(subsetMetab$Ingenuity.Canonical.Pathways, levels = subsetMetab$Ingenuity.Canonical.Pathways[order(subsetMetab$z.score, decreasing = F)])
ggplot(data=subsetMetab, aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`)) + 
  geom_bar( aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`) , stat="Identity", width=0.01, color="#e16462", size=1.5) + geom_point(size=5) + 
  coord_flip() + theme_bw() + ggtitle("IPA Canonical Pathways \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("z score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )

# ggsave(file="IPA/Images/AllAges_HivLo_v1_canonical_Metab.pdf", device="pdf",height=2.5,width=5)



subsetTraffic <- AllAgesHivLov1_can[grep("Actin|Extravasation", AllAgesHivLov1_can$Ingenuity.Canonical.Pathways),]  # suppress the rows with #NUM! in the z-score column
subsetTraffic$Ingenuity.Canonical.Pathways <- factor(subsetTraffic$Ingenuity.Canonical.Pathways, levels = subsetTraffic$Ingenuity.Canonical.Pathways[order(subsetTraffic$z.score, decreasing = F)])
ggplot(data=subsetTraffic, aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`)) + 
  geom_bar( aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`) , stat="Identity", width=0.01, color="#e16462", size=1.5) + geom_point(size=5) + 
  coord_flip() + theme_bw() + ggtitle("IPA Canonical Pathways \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("z score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )

# ggsave(file="IPA/Images/AllAges_HivLo_v1_canonical_Trafficking.pdf", device="pdf",height=2.5,width=6)


subsetInflam <- AllAgesHivLov1_can[grep("NF-kB|Inflammasome|Interferon|Toll", AllAgesHivLov1_can$Ingenuity.Canonical.Pathways),]  # suppress the rows with #NUM! in the z-score column
subsetInflam$Ingenuity.Canonical.Pathways <- factor(subsetInflam$Ingenuity.Canonical.Pathways, levels = subsetInflam$Ingenuity.Canonical.Pathways[order(subsetInflam$z.score, decreasing = F)])
ggplot(data=subsetInflam, aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`)) + 
  geom_bar( aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`) , stat="Identity", width=0.01, color="#e16462", size=1.5) + geom_point(size=5) + 
  coord_flip() + theme_bw() + ggtitle("IPA Canonical Pathways \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("z score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )

# ggsave(file="IPA/Images/AllAges_HivLo_v1_canonical_Inflamm.pdf", device="pdf",height=2.5,width=5.5)



## *****************************       IPA: YoungOnly canonical pathways    ***********************************************


YoungHivLov1_can <- read.csv("IPA/IPAresults/20181116_Young_HiHi-v-LoLo_v1_IPAcanonicalPathways.csv", stringsAsFactors = F )
YoungHivLov1_can <- YoungHivLov1_can[-grep("#NUM!", YoungHivLov1_can$z.score),]  # suppress the rows with #NUM! in the z-score column
YoungHivLov1_can$z.score <- as.numeric(YoungHivLov1_can$z.score)                  # convert to numeric
YoungHivLov1_can <- YoungHivLov1_can[which(abs(YoungHivLov1_can$z.score)>=2),]  # retain only rows with z score >2 or <-2
YoungHivLov1_can <- YoungHivLov1_can[which(abs(YoungHivLov1_can$X.log.p.value.)>=1.3),]  # retain only rows with -log(pvalue) >1.3 

subsetMetab <- YoungHivLov1_can[grep("Glycolysis|Oxidative Phosphorylation|Gluconeogenesis|Pentose", YoungHivLov1_can$Ingenuity.Canonical.Pathways),]  # suppress the rows with #NUM! in the z-score column
subsetMetab$Ingenuity.Canonical.Pathways <- factor(subsetMetab$Ingenuity.Canonical.Pathways, levels = subsetMetab$Ingenuity.Canonical.Pathways[order(subsetMetab$z.score, decreasing = F)])
ggplot(data=subsetMetab, aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`)) + 
  geom_bar( aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`) , stat="Identity", width=0.01, color="#e16462", size=1.5) + geom_point(size=5) + 
  coord_flip() + theme_bw() + ggtitle("IPA Canonical Pathways \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("z score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )
# no plot because these genesets are not statistically significant

subsetTraffic <- YoungHivLov1_can[grep("Actin|Extravasation", YoungHivLov1_can$Ingenuity.Canonical.Pathways),]  # suppress the rows with #NUM! in the z-score column
subsetTraffic$Ingenuity.Canonical.Pathways <- factor(subsetTraffic$Ingenuity.Canonical.Pathways, levels = subsetTraffic$Ingenuity.Canonical.Pathways[order(subsetTraffic$z.score, decreasing = F)])
ggplot(data=subsetTraffic, aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`)) + 
  geom_bar( aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`) , stat="Identity", width=0.01, color="#e16462", size=1.5) + geom_point(size=5) + 
  coord_flip() + theme_bw() + ggtitle("IPA Canonical Pathways \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("z score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )

subsetInflam <- YoungHivLov1_can[grep("NF-kB|Inflammasome|Interferon|Toll", YoungHivLov1_can$Ingenuity.Canonical.Pathways),]  # suppress the rows with #NUM! in the z-score column
subsetInflam$Ingenuity.Canonical.Pathways <- factor(subsetInflam$Ingenuity.Canonical.Pathways, levels = subsetInflam$Ingenuity.Canonical.Pathways[order(subsetInflam$z.score, decreasing = F)])
ggplot(data=subsetInflam, aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`)) + 
  geom_bar( aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`) , stat="Identity", width=0.01, color="#e16462", size=1.5) + geom_point(size=5) + 
  coord_flip() + theme_bw() + ggtitle("IPA Canonical Pathways \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("z score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )


## *****************************       IPA: ElderlyOnly canonical pathways    ***********************************************


EldHivLov1_can <- read.csv("IPA/IPAresults/20181116_Elderly_HiHi-v-LoLo_v1_IPAcanonicalPathways.csv", stringsAsFactors = F )
EldHivLov1_can <- EldHivLov1_can[-grep("#NUM!", EldHivLov1_can$z.score),]  # suppress the rows with #NUM! in the z-score column
EldHivLov1_can$z.score <- as.numeric(EldHivLov1_can$z.score)                  # convert to numeric
EldHivLov1_can <- EldHivLov1_can[which(abs(EldHivLov1_can$z.score)>=2),]  # retain only rows with z score >2 or <-2
EldHivLov1_can <- EldHivLov1_can[which(abs(EldHivLov1_can$X.log.p.value.)>=1.3),]  # retain only rows with -log(pvalue) >1.3 

subsetMetab <- EldHivLov1_can[grep("Glycolysis|Oxidative Phosphorylation|Gluconeogenesis|Pentose", EldHivLov1_can$Ingenuity.Canonical.Pathways),]  # suppress the rows with #NUM! in the z-score column
subsetMetab$Ingenuity.Canonical.Pathways <- factor(subsetMetab$Ingenuity.Canonical.Pathways, levels = subsetMetab$Ingenuity.Canonical.Pathways[order(subsetMetab$z.score, decreasing = F)])
ggplot(data=subsetMetab, aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`)) + 
  geom_bar( aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`) , stat="Identity", width=0.01, color="#e16462", size=1.5) + geom_point(size=5) + 
  coord_flip() + theme_bw() + ggtitle("IPA Canonical Pathways \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("z score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )

subsetTraffic <- EldHivLov1_can[grep("Actin|Extravasation", EldHivLov1_can$Ingenuity.Canonical.Pathways),]  # suppress the rows with #NUM! in the z-score column
subsetTraffic$Ingenuity.Canonical.Pathways <- factor(subsetTraffic$Ingenuity.Canonical.Pathways, levels = subsetTraffic$Ingenuity.Canonical.Pathways[order(subsetTraffic$z.score, decreasing = F)])
ggplot(data=subsetTraffic, aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`)) + 
  geom_bar( aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`) , stat="Identity", width=0.01, color="#e16462", size=1.5) + geom_point(size=5) + 
  coord_flip() + theme_bw() + ggtitle("IPA Canonical Pathways \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("z score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )

subsetInflam <- EldHivLov1_can[grep("NF-kB|Inflammasome|Interferon|Toll", EldHivLov1_can$Ingenuity.Canonical.Pathways),]  # suppress the rows with #NUM! in the z-score column
subsetInflam$Ingenuity.Canonical.Pathways <- factor(subsetInflam$Ingenuity.Canonical.Pathways, levels = subsetInflam$Ingenuity.Canonical.Pathways[order(subsetInflam$z.score, decreasing = F)])
ggplot(data=subsetInflam, aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`)) + 
  geom_bar( aes(x=`Ingenuity.Canonical.Pathways`, y=`z.score`) , stat="Identity", width=0.01, color="#e16462", size=1.5) + geom_point(size=5) + 
  coord_flip() + theme_bw() + ggtitle("IPA Canonical Pathways \nICOS+CD38+ vs ICOS-CD38- cTfh") + ylab("z score") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=14), axis.text = element_text(size=12), plot.title=element_text(size=12) )



## *****************************       Survivin and caspases    ***********************************************
birc5 <- bestDataLog["BIRC5",grep("HiHi",colnames(bestDataLog))]
casp3 <- bestDataLog["CASP3",grep("HiHi",colnames(bestDataLog))]
apoptosis <- rbind(birc5,casp3)
apoptosis <- as.data.frame(t(apoptosis))
ggplot(apoptosis[grep("v1",rownames(apoptosis)),], aes(x=`BIRC5`,y=`CASP3`)) + geom_point()
a <- lm(data = apoptosis, BIRC5~CASP3)



