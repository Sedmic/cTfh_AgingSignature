library("genefilter")
library("ggplot2")
library("grid")
library("ggrepel")
library("RColorBrewer")
library("DESeq2")
library("BiocParallel")
library("Rtsne")
library("pheatmap")
library("tools")
library("viridis")
library("ggpubr")
library("scales")
library("gridExtra")
library("GSVA"); library("GSEABase")
library("limma"); # library("GEOquery"); library("Biobase")
library("edgeR")
sessionInfo()




bestDataLog <- read.csv("bestDataLog.csv"); rownames(bestDataLog) <- bestDataLog$X; bestDataLog$X <- NULL


hallmarkYvE_v1_Hi <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v1_HALLMARK.GseaPreranked.1531000213854/gsea_report_for_na_pos_1531000213854.xls", sep="\t", stringsAsFactors = F),
                           read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v1_HALLMARK.GseaPreranked.1531000213854/gsea_report_for_na_neg_1531000213854.xls", sep="\t", stringsAsFactors = F)
)
hallmarkYvE_v1_Hi$NAME <- toTitleCase(tolower(substr(hallmarkYvE_v1_Hi$NAME,start =10, stop=50)))
hallmarkYvE_v1_Hi$NAME <- factor(hallmarkYvE_v1_Hi$NAME, levels = hallmarkYvE_v1_Hi$NAME[order(hallmarkYvE_v1_Hi$NES, decreasing = T)])
rownames(hallmarkYvE_v1_Hi) <- hallmarkYvE_v1_Hi$NAME

hallmarkYvE_v2_Hi <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1562121168676/gsea_report_for_na_pos_1562121168676.xls", sep="\t", stringsAsFactors = F),
                           read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1562121168676/gsea_report_for_na_neg_1562121168676.xls", sep="\t", stringsAsFactors = F)
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

hallmarkYvE_v2_Nav <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v2_Hallmark.GseaPreranked.1562122339545/gsea_report_for_na_pos_1562122339545.xls", sep="\t", stringsAsFactors = F), 
                            read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v2_Hallmark.GseaPreranked.1562122339545/gsea_report_for_na_neg_1562122339545.xls", sep="\t", stringsAsFactors = F)
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


## *****************************         deltaNES          ***********************************************



deltaNES <- hallmark_subsets
deltaNES$Hi_delta <- abs(deltaNES$Hi_v2 - deltaNES$Hi_v1)
deltaNES$Lo_delta <- abs(deltaNES$Lo_v2 - deltaNES$Lo_v1)
deltaNES$Nav_delta <- abs(deltaNES$Nav_v2 - deltaNES$Nav_v1)

deltaNES <- deltaNES[,-c(1:6)]
deltaNES$geneset <- rownames(deltaNES)
deltaNES <- deltaNES[order(deltaNES$Hi_delta, decreasing = T),]
deltaNES$geneset <- factor(deltaNES$geneset, levels=deltaNES$geneset[order(deltaNES$Hi_delta, decreasing = T)], ordered=T)
deltaNESmelted <- reshape2::melt(deltaNES, id.vars = "geneset") 


# plot using violin plots
ggplot(deltaNESmelted, aes(x=variable, y=value)) + geom_violin(color="grey") + theme_bw() + #facet_wrap(~variable) +
  # geom_sina(data=to_lower_ascii(deltaNESmelted), aes(x=variable, y=value), alpha=0.4, scale=F, method="density", maxwidth = .6) +   # fails in ggforce 0.2.1, awaiting update
   stat_summary(fun.y=median, geom="point", color="#f66334", size=50, shape='-') + #geom_dotplot(stackdir='center', binaxis='y', dotsize=0.3, binwidth=0.1) +
  geom_jitter(width = 0.05, color="gray50") + 
  theme(axis.text = element_text(size=14,hjust = 0.5), axis.title = element_blank(), plot.title = element_text(size=18,hjust = 0.5), panel.border = element_blank()) +
  ggtitle("delta NES for hallmark genesets") + theme(axis.text.x = element_text(angle=45,hjust=1))
# ggsave(filename = "DifferentialExpression/GSEA/Images/deltaNES_violinplot.pdf")

# plot in waterfall-bar graph form
ggplot(deltaNESmelted, aes(x=geneset, y=value)) + geom_point() + facet_wrap(~variable) + theme_bw() + 
  theme(axis.text = element_text(size=12,hjust = 0.5)) + theme(axis.title = element_text(size=14,hjust = 0.5)) + theme(plot.title = element_text(size=18,hjust = 0.5))+
  ggtitle("delta NES for hallmark genesets") + theme(axis.text.x = element_text(angle=45,hjust=1))


# ****************************************************** milieu Interior Nanostring dataset   ***************************************

nanostringExpr <- read.csv(file="DifferentialExpression/ComparePublishedGeneSets/MilieuInterior/Nano_1000_NULL.csv", stringsAsFactors = F)
nanostringDemo <- read.csv(file="DifferentialExpression/ComparePublishedGeneSets/MilieuInterior/Nano_demographics.csv", stringsAsFactors = F)
nano <- merge(nanostringDemo, nanostringExpr, by="SUBJID")
colnames(nano) <- sub('\\_.*', '', colnames(nano)) # take away _NULL from each gene name

nano$SUBJID <- paste0("SUBJ_",nano$SUBJID)
nano$AgeCateg <- NA;  nano$AgeCateg[which(nano$AGE.V0 > 65)] <- "E";    nano$AgeCateg[which(nano$AGE.V0 < 40)] <- "Y"   
nanoAging <- nano[which(nano$AGE.V0 > 65 | nano$AGE.V0 < 40), -grep(paste("Batch","SEX","CMV","X", sep="|"),colnames(nano))];   rownames(nanoAging) <- nanoAging$SUBJID;  nanoAging$SUBJID <- NULL

metaData <- data.frame(row.names=rownames(nanoAging));       metaData$ageGroup <- nanoAging$AgeCateg;      labels <- metaData$ageGroup

nanoAgingExprs <- as.data.frame(t(nanoAging[,-grep("AgeCateg", colnames(nanoAging), value = F),]));   nanoAgingExprs <- nanoAgingExprs[-1,];
mm <- model.matrix(~0 + labels);        y <- voom(nanoAgingExprs, mm, plot = T);          fit <- lmFit(y, mm);   # head(coef(fit))
contr <- makeContrasts(labelsE - labelsY, levels = colnames(coef(fit)));        tmp <- contrasts.fit(fit, contr);       tmp <- eBayes(tmp)
top.tableMI <- topTable(tmp, sort.by = "T", n = Inf) ;      
head(top.tableMI, 20); 
top.tableMI <- top.tableMI[order(top.tableMI$t, decreasing=F),];     
# write.csv(top.tableMI, file = "DifferentialExpression/ComparePublishedGeneSets/MilieuInterior/diffExpr_MI_YvE.csv")
# write.table(cbind(rownames(top.tableMI),as.numeric(top.tableMI$t)), file = "DifferentialExpression/ComparePublishedGeneSets/MilieuInterior/diffExpr_MI_YvE.rnk.txt", col.names = F, row.names=F, sep="\t")




# ******************************************************  Wistar BAA datasets from Yrs 2, 3, 4, and 5  ***************************************

BAAyr2 <- read.csv(file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/data_BAA/Yr2_d0.csv", stringsAsFactors = F); rownames(BAAyr2) <- BAAyr2$X; BAAyr2$X <- NULL
BAAyr3 <- read.csv(file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/data_BAA/Yr3_d0.csv", stringsAsFactors = F); rownames(BAAyr3) <- BAAyr3$X; BAAyr3$X <- NULL
BAAyr4 <- read.csv(file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/data_BAA/Yr4_d0.csv", stringsAsFactors = F); rownames(BAAyr4) <- BAAyr4$X; BAAyr4$X <- NULL
BAAyr5 <- read.csv(file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/data_BAA/Yr5_d0.csv", stringsAsFactors = F); rownames(BAAyr5) <- BAAyr5$X; BAAyr5$X <- NULL

metaData <- data.frame(row.names=colnames(BAAyr4));       metaData$ageGroup <- substr(colnames(BAAyr4),1,1);      labels <- metaData$ageGroup
mm <- model.matrix(~0 + labels);        y <- voom(BAAyr4, mm, plot = T);          fit <- lmFit(y, mm); 
contr <- makeContrasts(labelsA - labelsY, levels = colnames(coef(fit)));        tmp <- contrasts.fit(fit, contr);       tmp <- eBayes(tmp)
top.table4 <- topTable(tmp, sort.by = "T", n = Inf);      head(top.table4, 20); tail(top.table4, 20)
# write.csv(top.table4, file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/diffExpr_Yr4_YvE.csv")]
# write.csv(rownames(top.table4[which(top.table4$t > 0 & top.table4$adj.P.Val < 0.2),]), file="DifferentialExpression/ComparePublishedGeneSets/WistarBAA/agingSignature.csv")      # ***** aging signature ******
# write.csv(rownames(top.table4[which(top.table4$t < 0 & top.table4$adj.P.Val < 0.2),]), file="DifferentialExpression/ComparePublishedGeneSets/WistarBAA/youthSignature.csv")      # ***** youth signature ******


metaData <- data.frame(row.names=colnames(BAAyr2));       metaData$ageGroup <- substr(colnames(BAAyr2),1,1);      labels <- metaData$ageGroup
mm <- model.matrix(~0 + labels);        y <- voom(BAAyr2, mm, plot = T);          fit <- lmFit(y, mm); 
contr <- makeContrasts(labelsA - labelsY, levels = colnames(coef(fit)));        tmp <- contrasts.fit(fit, contr);       tmp <- eBayes(tmp)
top.table2 <- topTable(tmp, sort.by = "T", n = Inf) ;      head(top.table2, 20); tail(top.table2, 20)
# write.csv(top.table2, file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/diffExpr_Yr2_YvE.csv")

metaData <- data.frame(row.names=colnames(BAAyr3));       metaData$ageGroup <- substr(colnames(BAAyr3),1,1);      labels <- metaData$ageGroup
mm <- model.matrix(~0 + labels);        y <- voom(BAAyr3, mm, plot = T);          fit <- lmFit(y, mm); 
contr <- makeContrasts(labelsA - labelsY, levels = colnames(coef(fit)));        tmp <- contrasts.fit(fit, contr);       tmp <- eBayes(tmp)
top.table3 <- topTable(tmp, sort.by = "T", n = Inf);      head(top.table3, 20); tail(top.table3, 20)
# write.csv(top.table3, file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/diffExpr_Yr3_YvE.csv")

metaData <- data.frame(row.names=colnames(BAAyr5));       metaData$ageGroup <- substr(colnames(BAAyr5),1,1);      labels <- metaData$ageGroup
mm <- model.matrix(~0 + labels);        y <- voom(BAAyr5, mm, plot = T);          fit <- lmFit(y, mm); 
contr <- makeContrasts(labelsA - labelsY, levels = colnames(coef(fit)));        tmp <- contrasts.fit(fit, contr);       tmp <- eBayes(tmp)
top.table5 <- topTable(tmp, sort.by = "T", n = Inf);      head(top.table5, 20); tail(top.table5, 20)
# write.csv(top.table5, file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/diffExpr_Yr5_YvE.csv")



# ******************************************************  Plot all of the individual GSEA results  ***************************************



# ********************************************************************
# *******************  WISTAR BAA YEAR 2  ****************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetBAAyr2_YvE.GseaPreranked.1557956545938/")
BAAyr2gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557956545938.xls"), sep="\t")
BAAyr2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t")
NES <- BAAyr2gseaFDRd$NES[grep("WISTAR",BAAyr2gseaFDRd$NAME)]
FDR <- BAAyr2gseaFDRd$FDR.q.val[grep("YOUTH",BAAyr2gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
  {   rpt <- read.delim(file=paste0(path, "revprobeBAAyr4_targetBAAyr2_YvE.GseaPreranked.1557956545938.rpt"),row.names=NULL, stringsAsFactors = F);  FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3]) }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=BAAyr2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY622 - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr2_YOUTH.pdf", device="pdf", height=3.5, width=5)

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetBAAyr2_YvE.GseaPreranked.1557855277477/")
BAAyr2gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855277477.xls"), sep="\t")
BAAyr2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t")
NES <- BAAyr2gseaFDR$NES[grep("WISTAR",BAAyr2gseaFDR$NAME)]
FDR <- BAAyr2gseaFDR$FDR.q.val[grep("AGING",BAAyr2gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
  {  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetBAAyr2_YvE.GseaPreranked.1557855277477.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3]) }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=BAAyr2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY622 - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr2_AGING.pdf", device="pdf", height=3.5, width=5)





# ********************************************************************
# *******************  WISTAR BAA YEAR 3  ****************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetBAAyr3_YvE.GseaPreranked.1558267709905/")
BAAyr3gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1558267709905.xls"), sep="\t")
BAAyr3gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1558267709905.xls"), sep="\t")); BAAyr3gsea <- rbind(BAAyr3gsea, c(0,0,0,0,numGenes,0,0,0,0))
zeroRankRow <- c(0,0,0,0,1,0,0,0,0) 
BAAyr3gsea <- rbind(zeroRankRow,BAAyr3gsea)
NES <- BAAyr3gseaFDRd$NES[grep("WISTAR",BAAyr3gseaFDRd$NAME)]
FDR <- BAAyr3gseaFDRd$FDR.q.val[grep("YOUTH",BAAyr3gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
  {  rpt <- read.delim(file=paste0(path, "revprobeBAAyr4_targetBAAyr3_YvE.GseaPreranked.1558267709905.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=BAAyr3gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY648 - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr3_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetBAAyr3_YvE.GseaPreranked.1558267679137/")
BAAyr3gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1558267679137.xls"), sep="\t")
BAAyr3gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1558267679137.xls"), sep="\t")); BAAyr3gsea <- rbind(BAAyr3gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- BAAyr3gseaFDR$NES[grep("WISTAR",BAAyr3gseaFDR$NAME)]
FDR <- BAAyr3gseaFDR$FDR.q.val[grep("AGING",BAAyr3gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
  {  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetBAAyr3_YvE.GseaPreranked.1558267679137.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=BAAyr3gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY648 - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr3_AGING.pdf", device="pdf", height=3.5, width=5)






# ********************************************************************
# *******************  WISTAR BAA YEAR 5  ****************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetBAAyr5_YvE.GseaPreranked.1557956604086/")
BAAyr5gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557956604086.xls"), sep="\t")
BAAyr5gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557956604086.xls"), sep="\t")); BAAyr5gsea <- rbind(BAAyr5gsea, c(0,0,0,0,numGenes,0,0,0,0))
zeroRankRow <- c(0,0,0,0,1,0,0,0,0) 
BAAyr5gsea <- rbind(zeroRankRow,BAAyr5gsea)
NES <- BAAyr5gseaFDRd$NES[grep("WISTAR",BAAyr5gseaFDRd$NAME)]
FDR <- BAAyr5gseaFDRd$FDR.q.val[grep("YOUTH",BAAyr5gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "revprobeBAAyr4_targetBAAyr5_YvE.GseaPreranked.1557956604086.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=BAAyr5gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY819 - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr5_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetBAAyr5_YvE.GseaPreranked.1557855086230/")
BAAyr5gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855086230.xls"), sep="\t")
BAAyr5gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855086230.xls"), sep="\t")); BAAyr5gsea <- rbind(BAAyr5gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- BAAyr5gseaFDR$NES[grep("WISTAR",BAAyr5gseaFDR$NAME)]
FDR <- BAAyr5gseaFDR$FDR.q.val[grep("AGING",BAAyr5gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetBAAyr5_YvE.GseaPreranked.1557855086230.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=BAAyr5gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY819 - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr5_AGING.pdf", device="pdf", height=3.5, width=5)








# ********************************************************************
# *******************  GSE 79396 visit 1  ****************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetGSE79396_v1.GseaPreranked.1557956302721/")
GSE79396gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557956302721.xls"), sep="\t")
GSE79396gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# GSE79396gsea <- rbind(zeroRankRow,GSE79396gsea)
NES <- GSE79396gseaFDRd$NES[grep("WISTAR",GSE79396gseaFDRd$NAME)]
FDR <- GSE79396gseaFDRd$FDR.q.val[grep("YOUTH",GSE79396gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "revprobeBAAyr4_targetGSE79396_v1.GseaPreranked.1557956302721.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=GSE79396gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE79396 ZostaVax day 0 - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE79396_d0_YOUTH.pdf", device="pdf", height=3.5, width=5)



path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetGSE79396_v1.GseaPreranked.1557854598532/")
GSE79396gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557854598532.xls"), sep="\t")
GSE79396gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557854598532.xls"), sep="\t")); GSE79396gsea <- rbind(GSE79396gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- GSE79396gseaFDR$NES[grep("WISTAR",GSE79396gseaFDR$NAME)]
FDR <- GSE79396gseaFDR$FDR.q.val[grep("AGING",GSE79396gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE79396_v1.GseaPreranked.1557854598532.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=GSE79396gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE79396 ZostaVax day 0 - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE79396_d0_AGING.pdf", device="pdf", height=3.5, width=5)






# ********************************************************************
# *******************  GSE 123696 2013  ******************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetGSE123696_2013.GseaPreranked.1557951377531/")
GSE123696gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557951377531.xls"), sep="\t")
GSE123696gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# GSE123696gsea <- rbind(zeroRankRow,GSE123696gsea)
NES <- GSE123696gseaFDRd$NES[grep("WISTAR",GSE123696gseaFDRd$NAME)]
FDR <- GSE123696gseaFDRd$FDR.q.val[grep("YOUTH",GSE123696gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "revprobeBAAyr4_targetGSE123696_2013.GseaPreranked.1557951377531.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=GSE123696gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123696 (2013) - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123696_2013_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetGSE123696_2013.GseaPreranked.1557855311404/")
GSE123696gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855311404.xls"), sep="\t")
GSE123696gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855311404.xls"), sep="\t")); GSE123696gsea <- rbind(GSE123696gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- GSE123696gseaFDR$NES[grep("WISTAR",GSE123696gseaFDR$NAME)]
FDR <- GSE123696gseaFDR$FDR.q.val[grep("AGING",GSE123696gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE123696_2013.GseaPreranked.1557855311404.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=GSE123696gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123696 (2013) - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123696_2013_AGING.pdf", device="pdf", height=3.5, width=5)







# ********************************************************************
# *******************  GSE 123697 2014  ******************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetGSE123697_2014.GseaPreranked.1558267819173/")
GSE123697gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1558267819173.xls"), sep="\t")
GSE123697gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# GSE123697gsea <- rbind(zeroRankRow,GSE123697gsea)
NES <- GSE123697gseaFDRd$NES[grep("WISTAR",GSE123697gseaFDRd$NAME)]
FDR <- GSE123697gseaFDRd$FDR.q.val[grep("YOUTH",GSE123697gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "revprobeBAAyr4_targetGSE123697_2014.GseaPreranked.1558267819173.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=GSE123697gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123697 (2014) - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123697_2014_YOUTH.pdf", device="pdf", height=3.5, width=5)



path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetGSE123697_2014.GseaPreranked.1558267851571/")
GSE123697gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1558267851571.xls"), sep="\t")
GSE123697gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1558267851571.xls"), sep="\t")); GSE123697gsea <- rbind(GSE123697gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- GSE123697gseaFDR$NES[grep("WISTAR",GSE123697gseaFDR$NAME)]
FDR <- GSE123697gseaFDR$FDR.q.val[grep("AGING",GSE123697gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE123697_2014.GseaPreranked.1558267851571.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=GSE123697gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123697 (2014) - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123697_2014_AGING.pdf", device="pdf", height=3.5, width=5)




# ********************************************************************
# *******************  GSE 123698 2015  ******************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetGSE123698_2015.GseaPreranked.1557951429483/")
GSE123698gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557951429483.xls"), sep="\t")
GSE123698gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# GSE123698gsea <- rbind(zeroRankRow,GSE123698gsea)
NES <- GSE123698gseaFDRd$NES[grep("WISTAR",GSE123698gseaFDRd$NAME)]
FDR <- GSE123698gseaFDRd$FDR.q.val[grep("YOUTH",GSE123698gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "revprobeBAAyr4_targetGSE123698_2015.GseaPreranked.1557951429483.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=GSE123698gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123698 (2015) - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123698_2015_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetGSE123698_2015.GseaPreranked.1557855351693/")
GSE123698gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855351693.xls"), sep="\t")
GSE123698gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855351693.xls"), sep="\t")); GSE123698gsea <- rbind(GSE123698gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- GSE123698gseaFDR$NES[grep("WISTAR",GSE123698gseaFDR$NAME)]
FDR <- GSE123698gseaFDR$FDR.q.val[grep("AGING",GSE123698gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE123698_2015.GseaPreranked.1557855351693.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=GSE123698gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123698 (2015) - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123698_2015_AGING.pdf", device="pdf", height=3.5, width=5)






# ********************************************************************
# *******************  milieu Interieur  ****************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetMilieauInterieur.GseaPreranked.1557956370920/")
GSEmilieuInterieurgseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557956370920.xls"), sep="\t")
GSEmilieuInterieurgsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557956370920.xls"), sep="\t")); GSEmilieuInterieurgsea <- rbind(GSEmilieuInterieurgsea, c(0,0,0,0,numGenes,0,0,0,0))
GSEmilieuInterieurgsea <- rbind(zeroRankRow,GSEmilieuInterieurgsea)
NES <- GSEmilieuInterieurgseaFDRd$NES[grep("WISTAR",GSEmilieuInterieurgseaFDRd$NAME)]
FDR <- GSEmilieuInterieurgseaFDRd$FDR.q.val[grep("YOUTH",GSEmilieuInterieurgseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "revprobeBAAyr4_targetMilieauInterieur.GseaPreranked.1557956370920.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=GSEmilieuInterieurgsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("MilieuInterieur - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSEMilieauInterieurgsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetMilieauInterieur.GseaPreranked.1557855965199/")
GSEmilieuInterieurgseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855965199.xls"), sep="\t")
GSEmilieuInterieurgsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855965199.xls"), sep="\t")); GSEmilieuInterieurgsea <- rbind(GSEmilieuInterieurgsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- GSEmilieuInterieurgseaFDR$NES[grep("WISTAR",GSEmilieuInterieurgseaFDR$NAME)]
FDR <- GSEmilieuInterieurgseaFDR$FDR.q.val[grep("AGING",GSEmilieuInterieurgseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetMilieauInterieur.GseaPreranked.1557855965199.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=GSEmilieuInterieurgsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("MilieuInterieur - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSEMilieauInterieurgsea_AGING.pdf", device="pdf", height=3.5, width=5)



# ********************************************************************
# *********************  +/+ YvE day 7  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetHiHi_v2_YvE.GseaPreranked.1557945703384/")
HiHiv2gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557945703384.xls"), sep="\t")
HiHiv2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
HiHiv2gsea <- rbind(zeroRankRow,HiHiv2gsea)
NES <- HiHiv2gseaFDRd$NES[grep("WISTAR",HiHiv2gseaFDRd$NAME)]
FDR <- HiHiv2gseaFDRd$FDR.q.val[grep("YOUTH",HiHiv2gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetHiHi_v2_YvE.GseaPreranked.1557945703384.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.61, y=0.8, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=HiHiv2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS+CD38+ cTfh day 7 - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetHiHiv2gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetHiHi_v2_YvE.GseaPreranked.1557855723962/")
HiHiv2gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855723962.xls"), sep="\t")
HiHiv2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855723962.xls"), sep="\t")); HiHiv2gsea <- rbind(HiHiv2gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- HiHiv2gseaFDR$NES[grep("WISTAR",HiHiv2gseaFDR$NAME)]
FDR <- HiHiv2gseaFDR$FDR.q.val[grep("AGING",HiHiv2gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetHiHi_v2_YvE.GseaPreranked.1557855723962.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.61, y=0.8, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=HiHiv2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS+CD38+ cTfh day 7 - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetHiHiv2gsea_AGING.pdf", device="pdf", height=3.5, width=5)





# ********************************************************************
# *********************  +/+ YvE day 0  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetHiHi_v1_YvE.GseaPreranked.1557943543782/")
HiHiv1gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557943543782.xls"), sep="\t"); HiHiv1gseaFDRd <- HiHiv1gseaFDRd[-grep("AGING", HiHiv1gseaFDRd$NAME, value = F),]
HiHiv1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# HiHiv1gsea <- rbind(zeroRankRow,HiHiv1gsea)
NES <- HiHiv1gseaFDRd$NES[grep("WISTAR",HiHiv1gseaFDRd$NAME)]
FDR <- HiHiv1gseaFDRd$FDR.q.val[grep("YOUTH",HiHiv1gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetHiHi_v1_YvE.GseaPreranked.1557943543782.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.61, y=0.8, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=HiHiv1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS+CD38+ cTfh day 0 - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetHiHiv1gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetHiHi_v1_YvE.GseaPreranked.1557855402299/")
HiHiv1gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855402299.xls"), sep="\t"); HiHiv1gseaFDR <- HiHiv1gseaFDR[-grep("YOUTH", HiHiv1gseaFDR$NAME, value = F),]
HiHiv1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855402299.xls"), sep="\t")); HiHiv1gsea <- rbind(HiHiv1gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- HiHiv1gseaFDR$NES[grep("WISTAR",HiHiv1gseaFDR$NAME)]
FDR <- HiHiv1gseaFDR$FDR.q.val[grep("AGING",HiHiv1gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetHiHi_v1_YvE.GseaPreranked.1557855402299.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.21, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=HiHiv1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS+CD38+ cTfh day 0 - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetHiHiv1gsea_AGING.pdf", device="pdf", height=3.5, width=5)




# ********************************************************************
# *********************  -/- YvE day 7  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetLoLo_v2_YvE.GseaPreranked.1557944196589/")
LoLov2gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557944196589.xls"), sep="\t")
LoLov2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# LoLov2gsea <- rbind(zeroRankRow,LoLov2gsea)
NES <- LoLov2gseaFDRd$NES[grep("WISTAR",LoLov2gseaFDRd$NAME)]
FDR <- LoLov2gseaFDRd$FDR.q.val[grep("YOUTH",LoLov2gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetLoLo_v2_YvE.GseaPreranked.1557944196589.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.61, y=0.8, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=LoLov2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS-CD38- cTfh day 7 - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetLoLov2gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetLoLo_v2_YvE.GseaPreranked.1557855890936/")
LoLov2gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855890936.xls"), sep="\t")
LoLov2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855890936.xls"), sep="\t")); LoLov2gsea <- rbind(LoLov2gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- LoLov2gseaFDR$NES[grep("WISTAR",LoLov2gseaFDR$NAME)]
FDR <- LoLov2gseaFDR$FDR.q.val[grep("AGING",LoLov2gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetLoLo_v2_YvE.GseaPreranked.1557855890936.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.61, y=0.8, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=LoLov2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS-CD38- cTfh day 7 - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetLoLov2gsea_AGING.pdf", device="pdf", height=3.5, width=5)

# ********************************************************************
# *********************  -/- YvE day 0  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetLoLo_v1_YvE.GseaPreranked.1557944216635/")
LoLov1gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557944216635.xls"), sep="\t"); LoLov1gseaFDRd <- LoLov1gseaFDRd[-grep("AGING", LoLov1gseaFDRd$NAME, value = F),]
LoLov1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# LoLov1gsea <- rbind(zeroRankRow,LoLov1gsea)
NES <- LoLov1gseaFDRd$NES[grep("WISTAR",LoLov1gseaFDRd$NAME)]
FDR <- LoLov1gseaFDRd$FDR.q.val[grep("YOUTH",LoLov1gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetLoLo_v1_YvE.GseaPreranked.1557944216635.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.61, y=0.8, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=LoLov1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS-CD38- cTfh day 0 - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetLoLov1gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetLoLo_v1_YvE.GseaPreranked.1557855843503/")
LoLov1gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855843503.xls"), sep="\t"); LoLov1gseaFDR <- LoLov1gseaFDR[-grep("YOUTH", LoLov1gseaFDR$NAME, value = F),]
LoLov1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855843503.xls"), sep="\t")); LoLov1gsea <- rbind(LoLov1gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- LoLov1gseaFDR$NES[grep("WISTAR",LoLov1gseaFDR$NAME)]
FDR <- LoLov1gseaFDR$FDR.q.val[grep("AGING",LoLov1gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetLoLo_v1_YvE.GseaPreranked.1557855843503.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.21, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=LoLov1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS-CD38- cTfh day 0 - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetLoLov1gsea_AGING.pdf", device="pdf", height=3.5, width=5)



# ********************************************************************
# *********************  Nav YvE day 7  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetNaive_v2_YvE.GseaPreranked.1557944257266/")
Naivev2gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557944257266.xls"), sep="\t")
Naivev2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# Naivev2gsea <- rbind(zeroRankRow,Naivev2gsea)
NES <- Naivev2gseaFDRd$NES[grep("WISTAR",Naivev2gseaFDRd$NAME)]
FDR <- Naivev2gseaFDRd$FDR.q.val[grep("YOUTH",Naivev2gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetNaive_v2_YvE.GseaPreranked.1557944257266.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.61, y=0.8, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=Naivev2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Naive CD4 day 7 - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetNaivev2gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetNaive_v2_YvE.GseaPreranked.1557855934146/")
Naivev2gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855934146.xls"), sep="\t")
Naivev2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855934146.xls"), sep="\t")); Naivev2gsea <- rbind(Naivev2gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- Naivev2gseaFDR$NES[grep("WISTAR",Naivev2gseaFDR$NAME)]
FDR <- Naivev2gseaFDR$FDR.q.val[grep("AGING",Naivev2gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetNaive_v2_YvE.GseaPreranked.1557855934146.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.61, y=0.8, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=Naivev2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Naive CD4 day 7 - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetNaivev2gsea_AGING.pdf", device="pdf", height=3.5, width=5)


# ********************************************************************
# *********************  Nav YvE day 0  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/revprobeBAAyr4_targetNaive_v1_YvE.GseaPreranked.1557944236383/")
Naivev1gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_pos_1557944236383.xls"), sep="\t")
Naivev1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# Naivev1gsea <- rbind(zeroRankRow,Naivev1gsea)
NES <- Naivev1gseaFDRd$NES[grep("WISTAR",Naivev1gseaFDRd$NAME)]
FDR <- Naivev1gseaFDRd$FDR.q.val[grep("YOUTH",Naivev1gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetNaive_v1_YvE.GseaPreranked.1557944236383.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.61, y=0.8, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=Naivev1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Naive CD4 day 0 - Youth") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetNaivev1gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetNaive_v1_YvE.GseaPreranked.1557855915997/")
Naivev1gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855915997.xls"), sep="\t")
Naivev1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855915997.xls"), sep="\t")); Naivev1gsea <- rbind(Naivev1gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- Naivev1gseaFDR$NES[grep("WISTAR",Naivev1gseaFDR$NAME)]
FDR <- Naivev1gseaFDR$FDR.q.val[grep("AGING",Naivev1gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the Agingper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetNaive_v1_YvE.GseaPreranked.1557855915997.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.61, y=0.8, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=Naivev1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Naive CD4 day 0 - Aging") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetNaivev1gsea_AGING.pdf", device="pdf", height=3.5, width=5)





# ****************************************************** milieu Interieur correlations  ***************************************


gsets <- getGmt("DifferentialExpression/GSEA/AgingSignature/WistarBAAyr4_GeneSets.gmt.txt")
nanoGSVA <- as.data.frame(GSVA::gsva(as.matrix(t(nano[,7:600])), gsets, method="gsva"))

nanoGSVA <- as.data.frame( t(nanoGSVA) ); nanoGSVA$SUBJID <- nano$SUBJID
nanostringDemo$SUBJID <- paste0("SUBJ_",nanostringDemo$SUBJID)
nanoGSVAdemo <- merge(nanostringDemo, nanoGSVA, by="SUBJID")


# plot(nanoGSVAdemo$AGE.V0, nanoGSVAdemo$WistarBAA_Year4_YOUTH);  cor.test(nanoGSVAdemo$AGE.V0, nanoGSVAdemo$WistarBAA_Year4_YOUTH)
fit <- cor( nanoGSVAdemo$WistarBAA_Year4_YOUTH, nanoGSVAdemo$AGE.V0, method="pearson", use="complete.obs")
fit.p <-   cor.test( nanoGSVAdemo$WistarBAA_Year4_YOUTH, nanoGSVAdemo$AGE.V0, method="pearson", use="complete.obs")
annotationInfo <- paste0("r = ", round(fit, 2), ", ", "P = ", formatC(fit.p$p.value, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.05, hjust=0, gp=gpar(col="black", fontsize=36)))
ggplot(nanoGSVAdemo, aes(x=AGE.V0, y=WistarBAA_Year4_YOUTH)) + geom_point(size=3, alpha=0.5) + theme_bw() + 
  theme(axis.text = element_text(size=30,hjust = 0.5))+theme(axis.title = element_text(size=30,hjust = 0.5))+theme(plot.title = element_text(size=36,hjust = 0.5))+
  ggtitle("Milieu Interieur - Youth signature") + scale_x_continuous(name = "Age", breaks=seq(20,80,10)) + scale_y_continuous(name="GSVA scores for Youth signature") + 
  stat_smooth(method="lm", col="red") + annotation_custom(my_grob)
# ggsave(filename = "DifferentialExpression/GSEA/AgingSignature/Images/MilieauInterieur_vs_YouthSignature.pdf", width=9,height=9)

# plot(nanoGSVAdemo$AGE.V0, nanoGSVAdemo$WistarBAA_Year4_AGING); cor.test(nanoGSVAdemo$AGE.V0, nanoGSVAdemo$WistarBAA_Year4_AGING)
fit <- cor( nanoGSVAdemo$WistarBAA_Year4_AGING, nanoGSVAdemo$AGE.V0, method="pearson", use="complete.obs")
fit.p <-   cor.test( nanoGSVAdemo$WistarBAA_Year4_AGING, nanoGSVAdemo$AGE.V0, method="pearson", use="complete.obs")
annotationInfo <- paste0("r = ", round(fit, 2), ", ", "P = ", formatC(fit.p$p.value, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.05, hjust=0, gp=gpar(col="black", fontsize=36)))
ggplot(nanoGSVAdemo, aes(x=AGE.V0, y=WistarBAA_Year4_AGING)) + geom_point(size=3, alpha=0.5) + theme_bw() + 
  theme(axis.text = element_text(size=30,hjust = 0.5))+theme(axis.title = element_text(size=30,hjust = 0.5))+theme(plot.title = element_text(size=36,hjust = 0.5))+
  ggtitle("Milieu Interieur - Aging signature") + scale_x_continuous(name = "Age", breaks=seq(20,80,10)) + scale_y_continuous(name="GSVA scores for Aging signature") + 
  stat_smooth(method="lm", col="red") + annotation_custom(my_grob)
# ggsave(filename = "DifferentialExpression/GSEA/AgingSignature/Images/MilieauInterieur_vs_AgingSignature.pdf", width=9,height=9)


summary(lm ( nanoGSVAdemo$WistarBAA_Year4_YOUTH ~ nanoGSVAdemo$CMV + nanoGSVAdemo$AGE.V0 + nanoGSVAdemo$SEX)) 
ggplot(nanoGSVAdemo, aes(x=SEX, y=WistarBAA_Year4_YOUTH)) + geom_violin() + geom_boxplot(width=0.1) + theme_bw() 
summary(lm ( nanoGSVAdemo$WistarBAA_Year4_AGING ~ nanoGSVAdemo$CMV + nanoGSVAdemo$AGE.V0 +  nanoGSVAdemo$SEX))  
ggplot(nanoGSVAdemo, aes(x=SEX, y=WistarBAA_Year4_AGING)) + geom_violin() + geom_boxplot(width=0.1) + theme_bw()



# ****************************************************** Summary plots of the GSEA results  ***************************************





extMicroarraysAging <- rbind(BAAyr2gseaFDR, BAAyr3gseaFDR, BAAyr5gseaFDR, GSE79396gseaFDR, GSE123697gseaFDR, GSE123696gseaFDR, GSE123698gseaFDR, GSEmilieuInterieurgseaFDR)
extMicroarraysAging$NAME <- c("Immport SDY622", "Immport SDY648", "Immport SDY819", "GSE79396 baseline", "GSE123697 (2014)", "GSE123696 (2013)", "GSE123698 (2015)", "Milieu Interieur")
extMicroarraysAging$NAME <- factor(extMicroarraysAging$NAME, levels = extMicroarraysAging$NAME[order(extMicroarraysAging$NES, decreasing = F)])

ggplot(extMicroarraysAging, aes(x=NAME, y=NES)) + geom_point(size=8) + geom_bar(stat="identity", width=0.1, fill="black") + theme_bw() + ylab("Normalized Enrichment \nScore") + 
  theme(axis.text = element_text(size=16,hjust = 0.5), axis.title = element_text(size=16,hjust = 0.5), plot.title = element_text(size=24,hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.y = element_blank()) + coord_flip() + scale_y_continuous(breaks = seq(-6,6,1)) + ggtitle("Aging signature")
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/externalStudies_wholeBloodArrays_Aging.pdf", width=5, height=8)


extMicroarraysYouth <- rbind(BAAyr2gseaFDRd, BAAyr3gseaFDRd, BAAyr5gseaFDRd, GSE79396gseaFDRd, GSE123697gseaFDRd, GSE123696gseaFDRd, GSE123698gseaFDRd, GSEmilieuInterieurgseaFDRd)
extMicroarraysYouth$NAME <- c("Immport SDY622", "Immport SDY648", "Immport SDY819", "GSE79396 baseline", "GSE123697 (2014)", "GSE123696 (2013)", "GSE123698 (2015)", "Milieu Interieur")
extMicroarraysYouth$NAME <- factor(extMicroarraysYouth$NAME, levels = extMicroarraysYouth$NAME[order(extMicroarraysYouth$NES, decreasing = F)])

ggplot(extMicroarraysYouth, aes(x=NAME, y=NES)) + geom_point(size=8) + geom_bar(stat="identity", width=0.1, fill="black") + theme_bw() + ylab("Normalized Enrichment \nScore") + 
  theme(axis.text = element_text(size=16,hjust = 0.5), axis.title = element_text(size=16,hjust = 0.5), plot.title = element_text(size=24,hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.y = element_blank()) + coord_flip() + scale_y_continuous(breaks = seq(-6,6,1)) + ggtitle("Youth signature")
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/externalStudies_wholeBloodArrays_Youth.pdf", width=5, height=8)




customPalette <- colorRampPalette(brewer.pal(12,"Paired"))(12)
customPalette[11] <- "#DE966D"  # modify paired palette to eliminate yellow
temp <- customPalette[c(7,8,3,4,1,1,9,10,1,1,11,12)]
temp[5:6] <- c("#fff849","#d6cf44"); temp[9:10] <- c("#03681e","#003d10"); temp[11:12] <- c("#99930a","#686408")
customPalette <- temp;  show_col(customPalette)
CD4subsetsAgingSig <- rbind(Naivev1gseaFDR, Naivev2gseaFDR, LoLov1gseaFDR, LoLov2gseaFDR, HiHiv1gseaFDR, HiHiv2gseaFDR)
CD4subsetsAgingSig$NAME <- c("Naive CD4 d0", "Naive CD4 d7", "ICOS-CD38- cTfh d0", "ICOS-CD38- cTfh d7", "ICOS+CD38+ cTfh d0", "ICOS+CD38+ cTfh d7" )
CD4subsetsAgingSig$DAY <- c(0,7,0,7,0,7)
CD4subsetsAgingSig$NAME <- factor(CD4subsetsAgingSig$NAME, levels = CD4subsetsAgingSig$NAME)
myColors <- customPalette[c(5,6,3,4,1,2)];   names(myColors) <- levels(CD4subsetsAgingSig$NAME)
CD4subsetsAgingSig$FDR.q.val <- -log10(CD4subsetsAgingSig$FDR.q.val)
ggplot(CD4subsetsAgingSig, aes(x=NES, y=FDR.q.val, fill=NAME)) + geom_point(shape=21, aes(size=5)) + theme_bw() + #geom_text_repel(label=CD4subsetsAgingSig$NAME, point.padding = 0.25, size=5) +   
  scale_fill_manual(name = "NAME", values=myColors) +   # and in the main ggplot call,          
  scale_size(range=c(3,10)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme(panel.border = element_blank()) + 
  theme(axis.text = element_text(size=16,hjust = 0.5), axis.title = element_text(size=16,hjust = 0.5), plot.title = element_text(size=36,hjust = 0.5), axis.ticks.y = element_blank()) + 
  xlab("Normalized Enrichment Score") + ylab("-log10 False Discovery Rate") + ggtitle("Aging signature") + 
  theme(axis.title.x = element_text(margin=margin(t=5,r=0,b=0,l=0)), legend.position = "none", plot.title = element_text(margin=margin(l=0,r=0,t=0,b=20))) + 
  scale_y_continuous(limits=c(0,4)) + scale_x_continuous(limits=c(-3,3), breaks=seq(-3,3,1)) + geom_rect(aes(xmin=-3,xmax=3,ymin=-log10(0.05),ymax=4), fill="yellow", alpha=0.02) + 
  geom_hline(yintercept = -log10(0.05),linetype="dashed", color="grey") + annotate("text", x=-2.3, y=1.42, label="FDR <0.05", size=5, color="grey")
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/AgingSignature_CD4subsets.pdf", width=5, height=6)




CD4subsetsYouthSig <- rbind(Naivev1gseaFDRd, Naivev2gseaFDRd, LoLov1gseaFDRd, LoLov2gseaFDRd, HiHiv1gseaFDRd, HiHiv2gseaFDRd)
CD4subsetsYouthSig$NAME <- c("Naive CD4 d0", "Naive CD4 d7", "ICOS-CD38- cTfh d0", "ICOS-CD38- cTfh d7", "ICOS+CD38+ cTfh d0", "ICOS+CD38+ cTfh d7"  )
CD4subsetsAgingSig$DAY <- c(0,7,0,7,0,7)
CD4subsetsYouthSig$NAME <- factor(CD4subsetsYouthSig$NAME, levels = CD4subsetsYouthSig$NAME)
myColors <- customPalette[c(5,6,3,4,1,2)];   names(myColors) <- levels(CD4subsetsAgingSig$NAME)
CD4subsetsYouthSig$FDR.q.val <- -log10(CD4subsetsYouthSig$FDR.q.val)
ggplot(CD4subsetsYouthSig, aes(x=NES, y=FDR.q.val, fill=NAME)) + geom_point(shape=21, size=10) + theme_bw() + #geom_text_repel(label=CD4subsetsYouthSig$NAME, point.padding = 0.25, size=5) +   
  scale_fill_manual(name = "NAME", values=myColors) + 
  scale_size(range=c(3,10)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme(panel.border = element_blank()) + 
  theme(axis.text = element_text(size=16,hjust = 0.5), axis.title = element_text(size=16,hjust = 0.5), plot.title = element_text(size=36,hjust = 0.5), axis.ticks.y = element_blank()) + 
  xlab("Normalized Enrichment Score") + ylab("-log10 False Discovery Rate") + ggtitle("Youth signature") + 
  theme(axis.title.x = element_text(margin=margin(t=5,r=0,b=0,l=0)), legend.position = "none", plot.title = element_text(margin=margin(l=0,r=0,t=0,b=20))) + 
  scale_y_continuous(limits=c(0,4)) + scale_x_continuous(limits=c(-3,3), breaks=seq(-3,3,1))+ geom_rect(aes(xmin=-3,xmax=3,ymin=-log10(0.05),ymax=4), fill="yellow", alpha=0.02) + 
  geom_hline(yintercept = -log10(0.05),linetype="dashed", color="grey") + annotate("text", x=-2.3, y=1.42, label="FDR <0.05", size=5, color="grey")
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/YouthSignature_CD4subsets.pdf", width=5, height=6)


# ******************************************************  What is the aging gene signature?   ***************************************


#  *************************     c7 signatures by Fisher's exact  *********************************   

agingSignature <- read.csv(file="DifferentialExpression/GSEA/AgingSignature/WistarBAAyr4_GeneSets.gmx.csv", stringsAsFactors = F); agingSignature <- agingSignature[-1,]
c7signatures <- read.csv(file="DifferentialExpression/GSEA/C7_msigdb_gmx.csv", stringsAsFactors = F)
c7signatures <- c7signatures[-1,]
hallmarksets <- read.csv(file="DifferentialExpression/GSEA/Hallmark_msigdb_gmx.csv", stringsAsFactors = F); hallmarksets <- hallmarksets[-1,]

# compare agingSignature to C7 signatures by column, calculate Fishers for overlap

makeContingency <- function (genelist1, genelist2)
{ #             in A   not in A
  # in B          a       b
  # not in B      c       d=23000
  a <- length(which(genelist1 %in% genelist2));   b <- length(genelist2) - a;   c <- length(genelist1) - a;   d <- 33000
  return( matrix(ncol=2,  c( a, b, c, d)  ) )   }

overlapResults <- sapply(c7signatures, function (x) {    return( fisher.test( makeContingency(agingSignature[,1], x) )$p.value)      })
overlapResults <- as.data.frame(overlapResults[order(overlapResults, decreasing = F)])
# write.csv(overlapResults, file = "DifferentialExpression/GSEA/AgingSignature/comparisonToMsigDB_C7_youth.csv")

overlapResults <- sapply(c7signatures, function (x) {    return( fisher.test( makeContingency(agingSignature[,2], x) )$p.value)      })
overlapResults <- as.data.frame(overlapResults[order(overlapResults, decreasing = F)])
# write.csv(overlapResults, file = "DifferentialExpression/GSEA/AgingSignature/comparisonToMsigDB_C7_aging.csv")



overlapResults <- sapply(hallmarksets, function (x) {    return( fisher.test( makeContingency(agingSignature[,2], x) )$p.value)      })  # overlap of signature of aging vs Hallmark sets
overlapResults <- as.data.frame(overlapResults[order(overlapResults, decreasing = F)] )
# write.csv(overlapResults, file = "DifferentialExpression/GSEA/AgingSignature/comparisonToMsigDB_hallmark_aging.csv")

overlapResults <- sapply(hallmarksets, function (x) {    return( fisher.test( makeContingency(agingSignature[,1], x) )$p.value)      })  # overlap of signature of youth vs Hallmark sets 
overlapResults <- as.data.frame(overlapResults[order(overlapResults, decreasing = F)] )
# write.csv(overlapResults, file = "DifferentialExpression/GSEA/AgingSignature/comparisonToMsigDB_hallmark_youth.csv")



overlapHiHi_leadingEdge <- sapply(c7signatures, function (x) {  return (fisher.test(makeContingency( HiHiv2gsea$PROBE[which(HiHiv2gsea$CORE.ENRICHMENT == "Yes")] , x) )$p.value)})
overlapHiHi_leadingEdge <- as.data.frame(overlapHiHi_leadingEdge[order(overlapHiHi_leadingEdge, decreasing = F)])
# write.csv(overlapHiHi_leadingEdge, file = "DifferentialExpression/GSEA/AgingSignature/HiHiv2_leadingEdge_comparisonToMsigDB_C7.csv")
























