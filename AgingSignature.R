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
library("limma"); library("GEOquery"); library("Biobase")
sessionInfo()




# ****************************************************** Milieau Interior Nanostring dataset   ***************************************

nanostringExpr <- read.csv(file="DifferentialExpression/ComparePublishedGeneSets/MilieuInterior/Nano_1000_NULL.csv", stringsAsFactors = F)
nanostringDemo <- read.csv(file="DifferentialExpression/ComparePublishedGeneSets/MilieuInterior/Nano_demographics.csv", stringsAsFactors = F)

nano <- merge(nanostringDemo, nanostringExpr, by="SUBJID")
colnames(nano) <- sub('\\_.*', '', colnames(nano)) # take away _NULL from each gene name

gsets <- getGmt("DifferentialExpression/GSEA/deltaNES/YvE_datasets_deltaNES.gmt.txt")
# gsets <- getGmt("DifferentialExpression/GSEA/AgingSignature/merge_aging_horiz.gmt.txt")
gsets <- getGmt("DifferentialExpression/GSEA/GSEA_Results/AllAges_HiHi_vs_LoLo_v1_Hallmark.GseaPreranked.1539515189801/edb/gene_sets.gmt")
nanoGSVA <- as.data.frame(GSVA::gsva(as.matrix(t(nano[,7:600])), gsets, method="gsva"))

nanoGSVA <- as.data.frame( t(nanoGSVA) ); nanoGSVA$SUBJID <- nano$SUBJID
nanoGSVAdemo <- merge(nanostringDemo, nanoGSVA, by="SUBJID")


plot(nanoGSVAdemo$AGE.V0, nanoGSVAdemo$deltaNES_HiHi_cut2.25)
cor.test(nanoGSVAdemo$AGE.V0, nanoGSVAdemo$deltaNES_LoLo_cut2.5)
summary(lm ( nanoGSVAdemo$d7_hihi ~ nanoGSVAdemo$CMV + nanoGSVAdemo$AGE.V0 + nanoGSVAdemo$CMV*nanoGSVAdemo$AGE.V0 + nanoGSVAdemo$SEX))
ggplot(nanoGSVAdemo, aes(x=SEX, y=d7_hihi)) + geom_violin()

pheatmap(nano[,7:600], scale="column", cluster_col=T, cluster_row=T, #annotation_col = annotateHeatmap, 
         show_colnames=F, main="deltaNES",
         # annotation_colors = ann_colors, cutree_rows = 7,
         fontsize_row = 9, color=inferno(100)
         #, filename = "DifferentialExpression/GSEA/Images/____.pdf"
)
nano$SUBJID <- paste0("SUBJ_",nano$SUBJID)
nano$AgeCateg <- NA;  nano$AgeCateg[which(nano$AGE.V0 > 65)] <- "E";    nano$AgeCateg[which(nano$AGE.V0 < 40)] <- "Y"   
nanoAging <- nano[which(nano$AGE.V0 > 65 | nano$AGE.V0 < 40), -grep(paste("Batch","SEX","CMV","X", sep="|"),colnames(nano))];   rownames(nanoAging) <- nanoAging$SUBJID;  nanoAging$SUBJID <- NULL

metaData <- data.frame(row.names=rownames(nanoAging));       metaData$ageGroup <- nanoAging$AgeCateg;      labels <- metaData$ageGroup

nanoAgingExprs <- as.data.frame(t(nanoAging[,-grep("AgeCateg", rownames(nanoAgingExprs), value = F),]));   nanoAgingExprs <- nanoAgingExprs[-1,];
mm <- model.matrix(~0 + labels);        y <- voom(nanoAgingExprs, mm, plot = T);          fit <- lmFit(y, mm);   # head(coef(fit))
contr <- makeContrasts(labelsE - labelsY, levels = colnames(coef(fit)));        tmp <- contrasts.fit(fit, contr);       tmp <- eBayes(tmp)
top.tableMI <- topTable(tmp, sort.by = "T", n = Inf) ;      
head(top.tableMI, 20); 
top.tableMI <- top.tableMI[order(top.tableMI$t, decreasing=F),];     # write.csv(top.tableMI, file = "DifferentialExpression/ComparePublishedGeneSets/MilieuInterior/diffExpr_MI_YvE.csv")
# write.table(cbind(rownames(top.tableMI),as.numeric(top.tableMI$t)), file = "DifferentialExpression/ComparePublishedGeneSets/MilieuInterior/diffExpr_MI_YvE.rnk.txt", col.names = F, row.names=F, sep="\t")




# ****************************************************** Correl GSVA of deltaNES genes against VZV IgG titers  ***************************************


# library(GEOquery)
# GSE79396 <- GEOquery::getGEO('GSE79396',GSEMatrix=TRUE)
# exp_design <- pData(phenoData(GSE79396[[1]]))
# exp_design$group <- exp_design$title
# exp_design$group <- gsub(" ","_",exp_design$group)
# design_file <- exp_design[,c("geo_accession","group")]
# colnames(design_file) <- c("sample","condition")
# expr_data <- data.frame(exprs(GSE79396[[1]]))
# write.csv(expr_data, file = "DifferentialExpression/ComparePublishedGeneSets/GSE79396_ZosterVaccineAging/expr_data.csv")
expr_data <- read.csv(file="DifferentialExpression/ComparePublishedGeneSets/GSE79396_ZosterVaccineAging/expr_data.csv", stringsAsFactors = F)

affyAnnotations <- read.csv(file = "DifferentialExpression/ComparePublishedGeneSets/GSE79396_ZosterVaccineAging/AffymetrixDefinitions/probe-symbol.csv", stringsAsFactors=F)

# now convert the expr_data probes to gene symbols based on the lookup table
expr_data <- expr_data[1:54613,]  # rows beyond this do not have gene symbols according to Affy file
setdiff(rownames(expr_data)[1:54613], affyAnnotations[,1]) 
identical(rownames(expr_data)[1:54613], affyAnnotations[,1])    # since true, then rewrite every Affy probe with the gene symbol

expr_data$symbol <- affyAnnotations[,2]
expr_data <- expr_data[!duplicated(expr_data$symbol),]
rownames(expr_data) <- expr_data$symbol; expr_data$symbol <- NULL

sampleNames <- read.csv(file="DifferentialExpression/ComparePublishedGeneSets/GSE79396_ZosterVaccineAging/SampleSet.csv", stringsAsFactors = F)  # copy/pasted from GEO into spreadsheet
sampleNames <- sampleNames[-c(289:296),]
setdiff(colnames(expr_data), sampleNames$Accession)
identical(colnames(expr_data), sampleNames$Accession)
# now only really care about visit 1 and visit 4 from GEO data

sampleNames$X <- substr(sampleNames$Title, 12, 19)
colnames(expr_data) <- sampleNames$X

gsets <- getGmt("DifferentialExpression/GSEA/deltaNES/YvE_datasets_literature.gmt.txt")
GSE79396day0 <- as.matrix(expr_data[grep("D0",colnames(expr_data)),])
GSVAagingSignature <- as.data.frame(GSVA::gsva(GSE79396day0, gsets, method="gsva"))
GSVAagingSignatureDay0 <- GSVAagingSignature[,grep("D0", colnames(GSVAagingSignature))]
colnames(GSVAagingSignatureDay0) <- substr(colnames(GSVAagingSignatureDay0), 1, 5)
GSVAagingSignatureDay0 <- t(GSVAagingSignatureDay0)
elisaSubject <- read.csv("DifferentialExpression/ComparePublishedGeneSets/GSE79396_ZosterVaccineAging/elisa-subject.csv", stringsAsFactors = F)
# go for d14 IgG data
d14IgG <- elisaSubject[which(elisaSubject$STUDY_TIME_COLLECTED == 14 & elisaSubject$ANALYTE_REPORTED == "IgG"),]
d14IgG <- merge(d14IgG, GSVAagingSignatureDay0, by.x="SubjectLabel", by.y =0)
ggplot(d14IgG, aes(x=VALUE_PREFERRED, y=deltaNES_HiHi_cut1.6)) + geom_point(size=4) + theme_bw()
ggplot(d14IgG, aes(x=VALUE_PREFERRED, y=deltaNES_LoLo_cut1.6)) + geom_point(size=4) + theme_bw()
cor(d14IgG$VALUE_PREFERRED, d14IgG$deltaNES_HiHi_cut1.6)

elderlyOnly <- sampleNames[which(sampleNames$Age.group == "elderly"),]
elderlyOnly <- elderlyOnly[grep("D0", elderlyOnly$Title),]
length(which(d14IgG$SubjectLabel %in% substr(elderlyOnly$X, 1, 5)))
elderlyOnlyd14IgG <- d14IgG[which(d14IgG$SubjectLabel %in% substr(elderlyOnly$X, 1, 5)),]
ggplot(elderlyOnlyd14IgG, aes(x=VALUE_PREFERRED, y=deltaNES_HiHi_cut1.6)) + geom_point(size=4) + theme_bw()
cor.test(elderlyOnlyd14IgG$VALUE_PREFERRED, elderlyOnlyd14IgG$deltaNES_HiHi_cut2.5)


set1 <- as.data.frame(gsets[[1]])




# ******************************************************  Wistar BAA datasets from Yrs 2, 3, 4, and 5  ***************************************
library("limma")
library("edgeR")
BAAyr2 <- read.csv(file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/data_BAA/Yr2_d0.csv", stringsAsFactors = F); rownames(BAAyr2) <- BAAyr2$X; BAAyr2$X <- NULL
BAAyr3 <- read.csv(file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/data_BAA/Yr3_d0.csv", stringsAsFactors = F); rownames(BAAyr3) <- BAAyr3$X; BAAyr3$X <- NULL
BAAyr4 <- read.csv(file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/data_BAA/Yr4_d0.csv", stringsAsFactors = F); rownames(BAAyr4) <- BAAyr4$X; BAAyr4$X <- NULL
BAAyr5 <- read.csv(file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/data_BAA/Yr5_d0.csv", stringsAsFactors = F); rownames(BAAyr5) <- BAAyr5$X; BAAyr5$X <- NULL


metaData <- data.frame(row.names=colnames(BAAyr2));       metaData$ageGroup <- substr(colnames(BAAyr2),1,1);      labels <- metaData$ageGroup
mm <- model.matrix(~0 + labels);        y <- voom(BAAyr2, mm, plot = T);          fit <- lmFit(y, mm); 
contr <- makeContrasts(labelsA - labelsY, levels = colnames(coef(fit)));        tmp <- contrasts.fit(fit, contr);       tmp <- eBayes(tmp)
top.table2 <- topTable(tmp, sort.by = "T", n = Inf) ;      # write.csv(top.table2, file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/diffExpr_Yr2_YvE.csv")
head(top.table2, 20); tail(top.table2, 20)

metaData <- data.frame(row.names=colnames(BAAyr3));       metaData$ageGroup <- substr(colnames(BAAyr3),1,1);      labels <- metaData$ageGroup
mm <- model.matrix(~0 + labels);        y <- voom(BAAyr3, mm, plot = T);          fit <- lmFit(y, mm); 
contr <- makeContrasts(labelsA - labelsY, levels = colnames(coef(fit)));        tmp <- contrasts.fit(fit, contr);       tmp <- eBayes(tmp)
top.table3 <- topTable(tmp, sort.by = "T", n = Inf);      # write.csv(top.table3, file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/diffExpr_Yr3_YvE.csv")
head(top.table3, 20); tail(top.table3, 20)

metaData <- data.frame(row.names=colnames(BAAyr4));       metaData$ageGroup <- substr(colnames(BAAyr4),1,1);      labels <- metaData$ageGroup
mm <- model.matrix(~0 + labels);        y <- voom(BAAyr4, mm, plot = T);          fit <- lmFit(y, mm); 
contr <- makeContrasts(labelsA - labelsY, levels = colnames(coef(fit)));        tmp <- contrasts.fit(fit, contr);       tmp <- eBayes(tmp)
top.table4 <- topTable(tmp, sort.by = "T", n = Inf);      # write.csv(top.table4, file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/diffExpr_Yr4_YvE.csv")
head(top.table4, 20); tail(top.table4, 20)

metaData <- data.frame(row.names=colnames(BAAyr5));       metaData$ageGroup <- substr(colnames(BAAyr5),1,1);      labels <- metaData$ageGroup
mm <- model.matrix(~0 + labels);        y <- voom(BAAyr5, mm, plot = T);          fit <- lmFit(y, mm); 
contr <- makeContrasts(labelsA - labelsY, levels = colnames(coef(fit)));        tmp <- contrasts.fit(fit, contr);       tmp <- eBayes(tmp)
top.table5 <- topTable(tmp, sort.by = "T", n = Inf);      # write.csv(top.table5, file = "DifferentialExpression/ComparePublishedGeneSets/WistarBAA/diffExpr_Yr5_YvE.csv")
head(top.table5, 20); tail(top.table5, 20)


# ******************************************************  Plot all of the individual GSEA results  ***************************************



# ********************************************************************
# *******************  WISTAR BAA YEAR 2  ****************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetBAAyr2_YvE.GseaPreranked.1557855277477/")
BAAyr2gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855277477.xls"), sep="\t")
BAAyr2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t")
NES <- BAAyr2gseaFDRd$NES[grep("WISTAR",BAAyr2gseaFDRd$NAME)]
FDR <- BAAyr2gseaFDRd$FDR.q.val[grep("YOUTH",BAAyr2gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
  {   rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetBAAyr2_YvE.GseaPreranked.1557855277477.rpt"),row.names=NULL, stringsAsFactors = F);  FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3]) }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=BAAyr2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY622 - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr2_YOUTH.pdf", device="pdf", height=3.5, width=5)


BAAyr2gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855277477.xls"), sep="\t")
BAAyr2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t")
NES <- BAAyr2gseaFDR$NES[grep("WISTAR",BAAyr2gseaFDR$NAME)]
FDR <- BAAyr2gseaFDR$FDR.q.val[grep("AGING",BAAyr2gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
  {  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetBAAyr2_YvE.GseaPreranked.1557855277477.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3]) }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=BAAyr2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY622 - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  # annotate("text", x = -Inf, y = Inf, label = annotationInfo, hjust = 0, vjust = 1, parse = TRUE) + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr2_AGING.pdf", device="pdf", height=3.5, width=5)





# ********************************************************************
# *******************  WISTAR BAA YEAR 3  ****************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetBAAyr3_YvE.GseaPreranked.1557855251717/")
BAAyr3gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855251717.xls"), sep="\t")
BAAyr3gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
BAAyr3gsea <- rbind(zeroRankRow,BAAyr3gsea)
NES <- BAAyr3gseaFDRd$NES[grep("WISTAR",BAAyr3gseaFDRd$NAME)]
FDR <- BAAyr3gseaFDRd$FDR.q.val[grep("YOUTH",BAAyr3gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
  {  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetBAAyr3_YvE.GseaPreranked.1557855251717.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=BAAyr3gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY648 - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr3_YOUTH.pdf", device="pdf", height=3.5, width=5)


BAAyr3gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855251717.xls"), sep="\t")
BAAyr3gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855251717.xls"), sep="\t")); BAAyr3gsea <- rbind(BAAyr3gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- BAAyr3gseaFDR$NES[grep("WISTAR",BAAyr3gseaFDR$NAME)]
FDR <- BAAyr3gseaFDR$FDR.q.val[grep("AGING",BAAyr3gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
  {  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetBAAyr3_YvE.GseaPreranked.1557855251717.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=BAAyr3gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY648 - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr3_AGING.pdf", device="pdf", height=3.5, width=5)






# ********************************************************************
# *******************  WISTAR BAA YEAR 5  ****************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetBAAyr5_YvE.GseaPreranked.1557855086230/")
BAAyr5gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855086230.xls"), sep="\t")
BAAyr5gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
BAAyr5gsea <- rbind(zeroRankRow,BAAyr5gsea)
NES <- BAAyr5gseaFDRd$NES[grep("WISTAR",BAAyr5gseaFDRd$NAME)]
FDR <- BAAyr5gseaFDRd$FDR.q.val[grep("YOUTH",BAAyr5gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetBAAyr5_YvE.GseaPreranked.1557855086230.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=BAAyr5gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY819 - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr5_YOUTH.pdf", device="pdf", height=3.5, width=5)


BAAyr5gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855086230.xls"), sep="\t")
BAAyr5gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855086230.xls"), sep="\t")); BAAyr5gsea <- rbind(BAAyr5gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- BAAyr5gseaFDR$NES[grep("WISTAR",BAAyr5gseaFDR$NAME)]
FDR <- BAAyr5gseaFDR$FDR.q.val[grep("AGING",BAAyr5gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetBAAyr5_YvE.GseaPreranked.1557855086230.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=BAAyr5gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Immport SDY819 - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetBAAyr5_AGING.pdf", device="pdf", height=3.5, width=5)








# ********************************************************************
# *******************  GSE 79396 visit 1  ****************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetGSE79396_v1.GseaPreranked.1557854598532/")
GSE79396gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557854598532.xls"), sep="\t")
GSE79396gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# GSE79396gsea <- rbind(zeroRankRow,GSE79396gsea)
NES <- GSE79396gseaFDRd$NES[grep("WISTAR",GSE79396gseaFDRd$NAME)]
FDR <- GSE79396gseaFDRd$FDR.q.val[grep("YOUTH",GSE79396gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE79396_v1.GseaPreranked.1557854598532.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=GSE79396gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE79396 ZostaVax day 0 - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE79396_d0_YOUTH.pdf", device="pdf", height=3.5, width=5)



GSE79396gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557854598532.xls"), sep="\t")
GSE79396gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557854598532.xls"), sep="\t")); GSE79396gsea <- rbind(GSE79396gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- GSE79396gseaFDR$NES[grep("WISTAR",GSE79396gseaFDR$NAME)]
FDR <- GSE79396gseaFDR$FDR.q.val[grep("AGING",GSE79396gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE79396_v1.GseaPreranked.1557854598532.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=GSE79396gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE79396 ZostaVax day 0 - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE79396_d0_AGING.pdf", device="pdf", height=3.5, width=5)






# ********************************************************************
# *******************  GSE 123696 2013  ******************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetGSE123696_2013.GseaPreranked.1557855311404/")
GSE123696gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855311404.xls"), sep="\t")
GSE123696gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# GSE123696gsea <- rbind(zeroRankRow,GSE123696gsea)
NES <- GSE123696gseaFDRd$NES[grep("WISTAR",GSE123696gseaFDRd$NAME)]
FDR <- GSE123696gseaFDRd$FDR.q.val[grep("YOUTH",GSE123696gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE123696_2013.GseaPreranked.1557855311404.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=GSE123696gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123696 (2013) - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123696_2013_YOUTH.pdf", device="pdf", height=3.5, width=5)


GSE123696gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855311404.xls"), sep="\t")
GSE123696gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855311404.xls"), sep="\t")); GSE123696gsea <- rbind(GSE123696gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- GSE123696gseaFDR$NES[grep("WISTAR",GSE123696gseaFDR$NAME)]
FDR <- GSE123696gseaFDR$FDR.q.val[grep("AGING",GSE123696gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE123696_2013.GseaPreranked.1557855311404.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=GSE123696gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123696 (2013) - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123696_2013_AGING.pdf", device="pdf", height=3.5, width=5)







# ********************************************************************
# *******************  GSE 123687 2014  ******************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetGSE123687_2014.GseaPreranked.1557855332614/")
GSE123687gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855332614.xls"), sep="\t")
GSE123687gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# GSE123687gsea <- rbind(zeroRankRow,GSE123687gsea)
NES <- GSE123687gseaFDRd$NES[grep("WISTAR",GSE123687gseaFDRd$NAME)]
FDR <- GSE123687gseaFDRd$FDR.q.val[grep("YOUTH",GSE123687gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE123687_2014.GseaPreranked.1557855332614.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=GSE123687gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123687 (2014) - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123687_2014_YOUTH.pdf", device="pdf", height=3.5, width=5)



GSE123687gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855332614.xls"), sep="\t")
GSE123687gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855332614.xls"), sep="\t")); GSE123687gsea <- rbind(GSE123687gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- GSE123687gseaFDR$NES[grep("WISTAR",GSE123687gseaFDR$NAME)]
FDR <- GSE123687gseaFDR$FDR.q.val[grep("AGING",GSE123687gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE123687_2014.GseaPreranked.1557855332614.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=GSE123687gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123687 (2014) - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123687_2014_AGING.pdf", device="pdf", height=3.5, width=5)




# ********************************************************************
# *******************  GSE 123698 2015  ******************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetGSE123698_2015.GseaPreranked.1557855351693/")
GSE123698gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855351693.xls"), sep="\t")
GSE123698gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# GSE123698gsea <- rbind(zeroRankRow,GSE123698gsea)
NES <- GSE123698gseaFDRd$NES[grep("WISTAR",GSE123698gseaFDRd$NAME)]
FDR <- GSE123698gseaFDRd$FDR.q.val[grep("YOUTH",GSE123698gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE123698_2015.GseaPreranked.1557855351693.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=GSE123698gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123698 (2015) - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123698_2015_YOUTH.pdf", device="pdf", height=3.5, width=5)


GSE123698gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855351693.xls"), sep="\t")
GSE123698gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855351693.xls"), sep="\t")); GSE123698gsea <- rbind(GSE123698gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- GSE123698gseaFDR$NES[grep("WISTAR",GSE123698gseaFDR$NAME)]
FDR <- GSE123698gseaFDR$FDR.q.val[grep("AGING",GSE123698gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetGSE123698_2015.GseaPreranked.1557855351693.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=GSE123698gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE123698 (2015) - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSE123698_2015_AGING.pdf", device="pdf", height=3.5, width=5)






# ********************************************************************
# *******************  GSE 123698 2015  ******************************
# ********************************************************************

path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetMilieauInterieur.GseaPreranked.1557855965199/")
GSEMilieauInterieurgseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855965199.xls"), sep="\t")
GSEMilieauInterieurgsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
GSEMilieauInterieurgsea <- rbind(zeroRankRow,GSEMilieauInterieurgsea)
NES <- GSEMilieauInterieurgseaFDRd$NES[grep("WISTAR",GSEMilieauInterieurgseaFDRd$NAME)]
FDR <- GSEMilieauInterieurgseaFDRd$FDR.q.val[grep("YOUTH",GSEMilieauInterieurgseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetMilieauInterieur.GseaPreranked.1557855965199.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=GSEMilieauInterieurgsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("MilieauInterieur - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSEMilieauInterieurgsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


GSEMilieauInterieurgseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855965199.xls"), sep="\t")
GSEMilieauInterieurgsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855965199.xls"), sep="\t")); GSEMilieauInterieurgsea <- rbind(GSEMilieauInterieurgsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- GSEMilieauInterieurgseaFDR$NES[grep("WISTAR",GSEMilieauInterieurgseaFDR$NAME)]
FDR <- GSEMilieauInterieurgseaFDR$FDR.q.val[grep("AGING",GSEMilieauInterieurgseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetMilieauInterieur.GseaPreranked.1557855965199.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=GSEMilieauInterieurgsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("MilieauInterieur - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetGSEMilieauInterieurgsea_AGING.pdf", device="pdf", height=3.5, width=5)



# ********************************************************************
# *********************  +/+ YvE day 7  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetHiHi_v2_YvE.GseaPreranked.1557855723962/")
HiHiv2gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855723962.xls"), sep="\t")
HiHiv2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
HiHiv2gsea <- rbind(zeroRankRow,HiHiv2gsea)
NES <- HiHiv2gseaFDRd$NES[grep("WISTAR",HiHiv2gseaFDRd$NAME)]
FDR <- HiHiv2gseaFDRd$FDR.q.val[grep("YOUTH",HiHiv2gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetHiHi_v2_YvE.GseaPreranked.1557855723962.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=HiHiv2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS+CD38+ cTfh day 7 - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetHiHiv2gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


HiHiv2gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855723962.xls"), sep="\t")
HiHiv2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855723962.xls"), sep="\t")); HiHiv2gsea <- rbind(HiHiv2gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- HiHiv2gseaFDR$NES[grep("WISTAR",HiHiv2gseaFDR$NAME)]
FDR <- HiHiv2gseaFDR$FDR.q.val[grep("AGING",HiHiv2gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetHiHi_v2_YvE.GseaPreranked.1557855723962.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=HiHiv2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS+CD38+ cTfh day 7 - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetHiHiv2gsea_AGING.pdf", device="pdf", height=3.5, width=5)





# ********************************************************************
# *********************  +/+ YvE day 0  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetHiHi_v1_YvE.GseaPreranked.1557855402299/")
HiHiv1gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855402299.xls"), sep="\t"); HiHiv1gseaFDRd <- HiHiv1gseaFDRd[-grep("AGING", HiHiv1gseaFDRd$NAME, value = F),]
HiHiv1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# HiHiv1gsea <- rbind(zeroRankRow,HiHiv1gsea)
NES <- HiHiv1gseaFDRd$NES[grep("WISTAR",HiHiv1gseaFDRd$NAME)]
FDR <- HiHiv1gseaFDRd$FDR.q.val[grep("YOUTH",HiHiv1gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetHiHi_v1_YvE.GseaPreranked.1557855402299.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=HiHiv1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS+CD38+ cTfh day 0 - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetHiHiv1gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


HiHiv1gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855402299.xls"), sep="\t"); HiHiv1gseaFDR <- HiHiv1gseaFDR[-grep("YOUTH", HiHiv1gseaFDR$NAME, value = F),]
HiHiv1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855402299.xls"), sep="\t")); HiHiv1gsea <- rbind(HiHiv1gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- HiHiv1gseaFDR$NES[grep("WISTAR",HiHiv1gseaFDR$NAME)]
FDR <- HiHiv1gseaFDR$FDR.q.val[grep("AGING",HiHiv1gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetHiHi_v1_YvE.GseaPreranked.1557855402299.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=HiHiv1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS+CD38+ cTfh day 0 - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetHiHiv1gsea_AGING.pdf", device="pdf", height=3.5, width=5)




# ********************************************************************
# *********************  -/- YvE day 7  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetLoLo_v2_YvE.GseaPreranked.1557855890936/")
LoLov2gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855890936.xls"), sep="\t")
LoLov2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# LoLov2gsea <- rbind(zeroRankRow,LoLov2gsea)
NES <- LoLov2gseaFDRd$NES[grep("WISTAR",LoLov2gseaFDRd$NAME)]
FDR <- LoLov2gseaFDRd$FDR.q.val[grep("YOUTH",LoLov2gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetLoLo_v2_YvE.GseaPreranked.1557855890936.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=LoLov2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS-CD38- cTfh day 7 - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetLoLov2gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


LoLov2gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855890936.xls"), sep="\t")
LoLov2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855890936.xls"), sep="\t")); LoLov2gsea <- rbind(LoLov2gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- LoLov2gseaFDR$NES[grep("WISTAR",LoLov2gseaFDR$NAME)]
FDR <- LoLov2gseaFDR$FDR.q.val[grep("AGING",LoLov2gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetLoLo_v2_YvE.GseaPreranked.1557855890936.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=LoLov2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS-CD38- cTfh day 7 - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetLoLov2gsea_AGING.pdf", device="pdf", height=3.5, width=5)

# ********************************************************************
# *********************  -/- YvE day 0  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetLoLo_v1_YvE.GseaPreranked.1557855843503/")
LoLov1gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855843503.xls"), sep="\t"); LoLov1gseaFDRd <- LoLov1gseaFDRd[-grep("AGING", LoLov1gseaFDRd$NAME, value = F),]
LoLov1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# LoLov1gsea <- rbind(zeroRankRow,LoLov1gsea)
NES <- LoLov1gseaFDRd$NES[grep("WISTAR",LoLov1gseaFDRd$NAME)]
FDR <- LoLov1gseaFDRd$FDR.q.val[grep("YOUTH",LoLov1gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetLoLo_v1_YvE.GseaPreranked.1557855843503.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=LoLov1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS-CD38- cTfh day 0 - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetLoLov1gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


LoLov1gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855843503.xls"), sep="\t"); LoLov1gseaFDR <- LoLov1gseaFDR[-grep("YOUTH", LoLov1gseaFDR$NAME, value = F),]
LoLov1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855843503.xls"), sep="\t")); LoLov1gsea <- rbind(LoLov1gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- LoLov1gseaFDR$NES[grep("WISTAR",LoLov1gseaFDR$NAME)]
FDR <- LoLov1gseaFDR$FDR.q.val[grep("AGING",LoLov1gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetLoLo_v1_YvE.GseaPreranked.1557855843503.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=LoLov1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("ICOS-CD38- cTfh day 0 - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetLoLov1gsea_AGING.pdf", device="pdf", height=3.5, width=5)



# ********************************************************************
# *********************  Nav YvE day 7  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetNaive_v2_YvE.GseaPreranked.1557855934146/")
Naivev2gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855934146.xls"), sep="\t")
Naivev2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# Naivev2gsea <- rbind(zeroRankRow,Naivev2gsea)
NES <- Naivev2gseaFDRd$NES[grep("WISTAR",Naivev2gseaFDRd$NAME)]
FDR <- Naivev2gseaFDRd$FDR.q.val[grep("YOUTH",Naivev2gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetNaive_v2_YvE.GseaPreranked.1557855934146.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=Naivev2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Naive CD4 day 7 - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetNaivev2gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


Naivev2gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855934146.xls"), sep="\t")
Naivev2gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855934146.xls"), sep="\t")); Naivev2gsea <- rbind(Naivev2gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- Naivev2gseaFDR$NES[grep("WISTAR",Naivev2gseaFDR$NAME)]
FDR <- Naivev2gseaFDR$FDR.q.val[grep("AGING",Naivev2gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetNaive_v2_YvE.GseaPreranked.1557855934146.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=Naivev2gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Naive CD4 day 7 - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetNaivev2gsea_AGING.pdf", device="pdf", height=3.5, width=5)


# ********************************************************************
# *********************  Nav YvE day 0  ******************************
# ********************************************************************
path <- c("DifferentialExpression/GSEA/AgingSignature/GSEAresults/probeBAAyr4_targetNaive_v1_YvE.GseaPreranked.1557855915997/")
Naivev1gseaFDRd <- read.csv( paste0(path, "gsea_report_for_na_neg_1557855915997.xls"), sep="\t")
Naivev1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_YOUTH.xls"), sep="\t");   
# zeroRankRow <- c(0,0,0,0,1,0,0,0,0)     # force start at zero for some datasets to make GSEA plot cleaner
# Naivev1gsea <- rbind(zeroRankRow,Naivev1gsea)
NES <- Naivev1gseaFDRd$NES[grep("WISTAR",Naivev1gseaFDRd$NAME)]
FDR <- Naivev1gseaFDRd$FDR.q.val[grep("YOUTH",Naivev1gseaFDRd$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetNaive_v1_YvE.GseaPreranked.1557855915997.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=Naivev1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Naive CD4 day 0 - Down") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetNaivev1gsea_YOUTH.pdf", device="pdf", height=3.5, width=5)


Naivev1gseaFDR <- read.csv( paste0(path, "gsea_report_for_na_pos_1557855915997.xls"), sep="\t")
Naivev1gsea <- read.csv( paste0(path, "WISTARBAA_YEAR4_AGING.xls"), sep="\t");   
numGenes <- nrow(read.csv( paste0(path, "ranked_gene_list_na_pos_versus_na_neg_1557855915997.xls"), sep="\t")); Naivev1gsea <- rbind(Naivev1gsea, c(0,0,0,0,numGenes,0,0,0,0))
NES <- Naivev1gseaFDR$NES[grep("WISTAR",Naivev1gseaFDR$NAME)]
FDR <- Naivev1gseaFDR$FDR.q.val[grep("AGING",Naivev1gseaFDR$NAME)]
if (FDR == 0)  # if FDR == 0, then will set the FDR to 1/permutations in the .rpt file since that is the upper bound on what it could be
{  rpt <- read.delim(file=paste0(path, "probeBAAyr4_targetNaive_v1_YvE.GseaPreranked.1557855915997.rpt"),row.names=NULL, stringsAsFactors = F);   FDR <- 1/as.numeric(rpt[ grep("nperm", rpt[,2]), 3])   }
annotationInfo <- paste0("NES: ", round(NES, 2), "\n", "FDR: ", formatC(FDR, format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.1, y=0.2, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=Naivev1gsea, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("Naive CD4 day 0 - Up") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/probeBAAyr4_targetNaivev1gsea_AGING.pdf", device="pdf", height=3.5, width=5)






# ****************************************************** Summary plots of the GSEA results  ***************************************





extMicroarrays <- rbind(BAAyr2gseaFDR, BAAyr3gseaFDR, BAAyr5gseaFDR, GSE79396gseaFDR, GSE123687gseaFDR, GSE123696gseaFDR, GSE123698gseaFDR, GSEMilieauInterieurgseaFDR,
                        BAAyr2gseaFDRd, BAAyr3gseaFDRd, BAAyr5gseaFDRd, GSE79396gseaFDRd, GSE123687gseaFDRd, GSE123696gseaFDRd, GSE123698gseaFDRd, GSEMilieauInterieurgseaFDRd)
extMicroarrays$NAME <- c("Immport SDY622 - Up", "Immport SDY648 - Up", "Immport SDY819 - Up", "GSE79396 day 0 - Up", "GSE123687 (2013) - Up", "GSE123696 (2014) - Up", "GSE123698 (2015) - Up", "Milieau Interieur - Up",
                         "Immport SDY622 - Down", "Immport SDY648 - Down", "Immport SDY819 - Down", "GSE79396 day 0 - Down", "GSE123687 (2013) - Down", "GSE123696 (2014) - Down", "GSE123698 (2015) - Down", "Milieau Interieur - Down")
extMicroarrays$NAME <- factor(extMicroarrays$NAME, levels = extMicroarrays$NAME[order(extMicroarrays$NES, decreasing = F)])

ggplot(extMicroarrays, aes(x=NAME, y=NES)) + geom_point(size=8) + geom_bar(stat="identity", width=0.1, fill="black") + theme_bw() + ylab("Normalized Enrichment Score") + 
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=16,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.y = element_blank()) + coord_flip() + scale_y_continuous(breaks = seq(-6,6,2))
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/externalStudies_wholeBloodArrays.pdf")



CD4subsetsAgingSig <- rbind(Naivev1gseaFDR, Naivev2gseaFDR, LoLov1gseaFDR, LoLov2gseaFDR, HiHiv1gseaFDR, HiHiv2gseaFDR,
                            Naivev1gseaFDRd, Naivev2gseaFDRd, LoLov1gseaFDRd, LoLov2gseaFDRd, HiHiv1gseaFDRd, HiHiv2gseaFDRd)
CD4subsetsAgingSig$NAME <- c("Naive CD4 d0 - Up", "Naive CD4 d7 - Up", "ICOS-CD38- cTfh d0 - Up", "ICOS-CD38- cTfh d7 - Up", "ICOS+CD38+ cTfh d0 - Up", "ICOS+CD38+ cTfh d7 - Up",
                             "Naive CD4 d0 - Down", "Naive CD4 d7 - Down", "ICOS-CD38- cTfh d0 - Down", "ICOS-CD38- cTfh d7 - Down", "ICOS+CD38+ cTfh d0 - Down", "ICOS+CD38+ cTfh d7 - Down"  )
CD4subsetsAgingSig$DAY <- c(0,7,0,7,0,7,0,7,0,7,0,7)
CD4subsetsAgingSig$NAME <- factor(CD4subsetsAgingSig$NAME, levels = CD4subsetsAgingSig$NAME)
#  "ICOS+CD38+ cTfh" ="#0D0887",  "ICOS-CD38- cTfh" = "#E16462", "Naive CD4" ="#F0F921"
# myColors <- c("#FFEB00","#9A8E00","#EA0600","#9A0400","#2019C9","#08056B");   names(myColors) <- levels(CD4subsetsAgingSig$NAME)

# CD4subsetsAgingSig$FDR.q.val <- CD4subsetsAgingSig$FDR.q.val + 0.00002    # because 1/permutations in order to not have FDR=0
CD4subsetsAgingSig$FDR.q.val <- -log10(CD4subsetsAgingSig$FDR.q.val)
ggplot(CD4subsetsAgingSig, aes(x=NES, y=FDR.q.val)) + geom_point(aes(size=FDR.q.val)) + theme_bw() + geom_text_repel(label=CD4subsetsAgingSig$NAME, point.padding = 0.75, size=5) +   
#   scale_color_manual(name = "NAME", values=myColors,) +   # and in the main ggplot call,          aes( , color=NAME)
  scale_size(range=c(3,10)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme(panel.border = element_blank()) + 
  theme(axis.text = element_text(size=16,hjust = 0.5), axis.title = element_text(size=16,hjust = 0.5), plot.title = element_text(size=24,hjust = 0.5), axis.ticks.y = element_blank()) + 
  xlab("Normalized Enrichment Score") + ylab("-log10 False Discovery Rate") + ggtitle("Aging signatures in CD4 subsets") + 
  theme(axis.title.x = element_text(margin=margin(t=5,r=0,b=0,l=0)), legend.position = "left", plot.title = element_text(margin=margin(l=0,r=0,t=0,b=20))) + 
  #geom_text(aes(label=NAME), hjust=0, nudge_x = 0.1, nudge_y=0.25, size=5 ) + 
  scale_y_continuous(limits=c(0,5)) + scale_x_continuous(limits=c(-3,3), breaks=seq(-3,3,1))
# ggsave(file="DifferentialExpression/GSEA/AgingSignature/Images/AgingSignature_CD4subsets.pdf", width=10, height=6)



# annotateHeatmap <- data.frame(row.names = colnames(deltaNES[,-4]), Subset = c("ICOS+CD38+ cTfh", "ICOS-CD38- cTfh", "Naive CD4"))
# ann_colors = list(  Subset = c("ICOS+CD38+ cTfh" ="#0D0887", "ICOS-CD38- cTfh" = "#E16462", "Naive CD4" ="#F0F921")  )
pheatmap(CD4subsetsAgingSig$NES, scale="none", cluster_col=F, # annotation_col = annotateHeatmap, annotation_colors = ann_colors, 
         show_rownames=T, cluster_row=F, main="deltaNES",
         fontsize_row = 15, color=inferno(100)
         #, filename = "DifferentialExpression/GSEA/Images/deltaNES_heatmap.pdf"
)



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
# write.csv(overlapResults, file = "DifferentialExpression/GSEA/AgingSignature/comparisonToMsigDB_C7.csv")

overlapHiHi_leadingEdge <- sapply(c7signatures, function (x) {  return (fisher.test(makeContingency( HiHiv2gsea$PROBE[which(HiHiv2gsea$CORE.ENRICHMENT == "Yes")] , x) )$p.value)})
overlapHiHi_leadingEdge <- as.data.frame(overlapHiHi_leadingEdge[order(overlapHiHi_leadingEdge, decreasing = F)])
# write.csv(overlapHiHi_leadingEdge, file = "DifferentialExpression/GSEA/AgingSignature/HiHiv2_leadingEdge_comparisonToMsigDB_C7.csv")


overlapResults <- sapply(hallmarksets, function (x) {    return( fisher.test( makeContingency(agingSignature[,2], x) )$p.value)      })
overlapResults <- as.data.frame(overlapResults[order(overlapResults, decreasing = F)] )
# write.csv(overlapResults, file = "DifferentialExpression/GSEA/AgingSignature/comparisonToMsigDB_C7.csv")









































hallmarkYvE_v1_Hi <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v1_HALLMARK.GseaPreranked.1531000213854/gsea_report_for_na_pos_1531000213854.xls", sep="\t", stringsAsFactors = F),
                           read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v1_HALLMARK.GseaPreranked.1531000213854/gsea_report_for_na_neg_1531000213854.xls", sep="\t", stringsAsFactors = F)
)
hallmarkYvE_v1_Hi$NAME <- toTitleCase(tolower(substr(hallmarkYvE_v1_Hi$NAME,start =10, stop=50)))
hallmarkYvE_v1_Hi$NAME <- factor(hallmarkYvE_v1_Hi$NAME, levels = hallmarkYvE_v1_Hi$NAME[order(hallmarkYvE_v1_Hi$NES, decreasing = T)])
rownames(hallmarkYvE_v1_Hi) <- hallmarkYvE_v1_Hi$NAME

hallmarkYvE_v2_Hi <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/gsea_report_for_na_pos_1531000206530.xls", sep="\t", stringsAsFactors = F),
                           read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/gsea_report_for_na_neg_1531000206530.xls", sep="\t", stringsAsFactors = F)
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

hallmarkYvE_v2_Nav <- rbind(read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v2_Hallmark.GseaPreranked.1540914108109/gsea_report_for_na_pos_1540914108109.xls", sep="\t", stringsAsFactors = F), 
                            read.csv("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v2_Hallmark.GseaPreranked.1540914108109/gsea_report_for_na_neg_1540914108109.xls", sep="\t", stringsAsFactors = F)
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
ggplot(deltaNESmelted, aes(x=variable, y=value)) + geom_violin() + theme_bw() + #facet_wrap(~variable) +
  # geom_sina(data=to_lower_ascii(deltaNESmelted), aes(x=variable, y=value), alpha=0.4, scale=F, method="density", maxwidth = .6) +   # fails in ggforce 0.2.1, awaiting update
  #  stat_summary(fun.y=median, geom="point", color="red", size=10) + geom_dotplot(stackdir='center', binaxis='y', dotsize=0.3, binwidth=0.1) +
  geom_jitter(width = 0.05) + 
  theme(axis.text = element_text(size=12,hjust = 0.5)) + theme(axis.title = element_text(size=14,hjust = 0.5)) + theme(plot.title = element_text(size=18,hjust = 0.5))+
  ggtitle("delta NES for hallmark genesets") + theme(axis.text.x = element_text(angle=45,hjust=1))
# ggsave(filename = "DifferentialExpression/GSEA/Images/deltaNES_violinplot.pdf")

# plot in waterfall-bar graph form
ggplot(deltaNESmelted, aes(x=geneset, y=value)) + geom_point() + facet_wrap(~variable) + theme_bw() + 
  theme(axis.text = element_text(size=12,hjust = 0.5)) + theme(axis.title = element_text(size=14,hjust = 0.5)) + theme(plot.title = element_text(size=18,hjust = 0.5))+
  ggtitle("delta NES for hallmark genesets") + theme(axis.text.x = element_text(angle=45,hjust=1))


# numGenes_deltaNES <- c(nrow(BAA_deltaNES_HiHi), nrow(BAA_deltaNES_LoLo), nrow(BAA_deltaNES_Naive))
# names(numGenes_deltaNES) <- c("ICOS+CD38+ cTfh", "ICOS-CD38- cTfh", "NaivWe CD4")



## *****************************       Construct Aging gene signature     ***********************************************

#  now need to construct the leading edge for these genes in the HiHi dataset from the day 7 leading edge genes


enrichedGenesHi16 <- returnGeneList(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), cutoff=1.6, "Hi_delta")
# write.csv(enrichedGenesHi16, file="DifferentialExpression/GSEA/deltaNES_leadingEdge_HiHi.csv")

enrichedGenesLo16 <- returnGeneList(c("DifferentialExpression/GSEA/GSEA_Results/Y_LoLo_vs_E_LoLo_v2_Hallmark.GseaPreranked.1540914059884/"), cutoff=1.6, "Lo_delta")
# write.csv(enrichedGenesLo16, file="DifferentialExpression/GSEA/deltaNES_leadingEdge_LoLo.csv")

enrichedGenesNav16 <- returnGeneList(c("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v2_Hallmark.GseaPreranked.1540914108109/"), cutoff=1.6, "Nav_delta")
# write.csv(enrichedGenesNav16, file="DifferentialExpression/GSEA/deltaNES_leadingEdge_Naive.csv")

enrichedGenesHi25 <- returnGeneList(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), cutoff=2.5, "Hi_delta")
write.csv(enrichedGenesHi25, file="DifferentialExpression/GSEA/deltaNES/deltaNES_leadingEdge_HiHi_cut2.5.csv")

enrichedGenesLo25 <- returnGeneList(c("DifferentialExpression/GSEA/GSEA_Results/Y_LoLo_vs_E_LoLo_v2_Hallmark.GseaPreranked.1540914059884/"), cutoff=2.5, "Lo_delta")
write.csv(enrichedGenesLo25, file="DifferentialExpression/GSEA/deltaNES/deltaNES_leadingEdge_LoLo_cut2.5.csv")

enrichedGenesNav25 <- returnGeneList(c("DifferentialExpression/GSEA/GSEA_Results/Y_naive_vs_E_naive_v2_Hallmark.GseaPreranked.1540914108109/"), cutoff=2.5, "Nav_delta")
write.csv(enrichedGenesNav25, file="DifferentialExpression/GSEA/deltaNES/deltaNES_leadingEdge_Naive_cut2.5.csv")



returnGeneList <- function(path, cutoff, cellSubset)  # assuming deltaNES is in the global scope
{
  enrichedGenes <- data.frame(); temp <- data.frame()
  a <- which(colnames(deltaNES) == cellSubset)  # +/+, -/-, or naive delta columns
  genesetsAging <- deltaNES[which(deltaNES[,a] > cutoff),"geneset"]; 
  genesetsAging <- paste0("HALLMARK_",toupper(genesetsAging),".xls")
  for (x in seq_along(genesetsAging)) {
    temp <- read.table(paste0(path, genesetsAging[x]),stringsAsFactors = F, sep="\t")    # read in the xls files for the gene sets of interest
    enrichedGenes <- rbind(enrichedGenes, temp[-1,])       }      # mash genesets all together into one 
  enrichedGenes <- enrichedGenes[which(enrichedGenes$V8 == "Yes"),]   # retain core enrichment genes that make up leading edge
  colnames(enrichedGenes) <- temp[1,] 
  enrichedGenes <- enrichedGenes[,-grep(paste(c("NAME","GENE SYMBOL", "GENE_TITLE","RUNNING ES","NA"), collapse="|"), colnames(enrichedGenes))]   # excess columns
  enrichedGenes <- enrichedGenes[which(!duplicated(enrichedGenes$PROBE)),]  # deduplicate the full list since some genes are in multiple genesets
  return(enrichedGenes)
}



#  ************************* test individual blocks in the hierarchical clustering of genesets *********************************   


returnGeneListBlocks <- function(path, pathways)
{
  enrichedGenes <- data.frame(); temp <- data.frame()
  genesetsAging <- paste0("HALLMARK_",toupper(pathways),".xls")
  for (x in seq_along(genesetsAging)) {
    temp <- read.table(paste0(path, genesetsAging[x]),stringsAsFactors = F, sep="\t")  
    enrichedGenes <- rbind(enrichedGenes, temp[-1,])       }  
  enrichedGenes <- enrichedGenes[which(enrichedGenes$V8 == "Yes"),]  # leading edge only
  colnames(enrichedGenes) <- temp[1,]
  enrichedGenes <- enrichedGenes[,-grep(paste(c("NAME","GENE SYMBOL", "GENE_TITLE","RUNNING ES","NA"), collapse="|"), colnames(enrichedGenes))]
  enrichedGenes <- enrichedGenes[which(!duplicated(enrichedGenes$PROBE)),]
  return(enrichedGenes)
}

annotateHeatmap <- data.frame(row.names = colnames(deltaNES[,-4]), Subset = c("ICOS+CD38+ cTfh", "ICOS-CD38- cTfh", "Naive CD4"))
ann_colors = list(  Subset = c("ICOS+CD38+ cTfh" ="#0D0887", "ICOS-CD38- cTfh" = "#E16462", "Naive CD4" ="#F0F921")  )
pheatmap(deltaNES[,-4], scale="none", cluster_col=F, annotation_col = annotateHeatmap, show_colnames=F, main="deltaNES",
         annotation_colors = ann_colors, cutree_rows = 7,
         fontsize_row = 9, color=inferno(100)
         #, filename = "DifferentialExpression/GSEA/Images/deltaNES_heatmap.pdf"
)


genesetNames <- rownames(hallmark_subsets)
blockA <- genesetNames[grep(paste(c("Reactive","gamma","P53","Peroxisome","Bile"), collapse="|"), genesetNames, value=F)]
blockB <- genesetNames[grep(paste(c("Tnfa","Oxidative","Unfolded","alpha","Apical_surface"), collapse="|"), genesetNames, value=F)]
blockC <- genesetNames[grep(paste(c("Uv_response_dn","Spermatogenesis"), collapse="|"), genesetNames, value=F)]
blockD <- genesetNames[grep(paste(c("Mitotic","Fatty","Kras_signaling_dn","Dna_repair","G2m","Wnt"), collapse="|"), genesetNames, value=F)]
blockE <- genesetNames[grep(paste(c("Androgen","Myc_targets_v1","Notch","Protein_secretion","Estrogen_response_early"), collapse="|"), genesetNames, value=F)]
blockF <- genesetNames[grep(paste(c("Pi3k","Heme","Adipogenesis","Tgf","Glycolysis","Myc_targets_v2","Apical_junction"), collapse="|"), genesetNames, value=F)]
blockG <- genesetNames[grep(paste(c("Coagulation","Complement","Kras_signaling_up","Angiogenesis","Inflammatory","Hedgehog","Allograft","Epithelial",
                                    "Mtorc1","Apoptosis","Myogenesis","Cholesterol","Il2","Xenobiotic","Estrogen_response_late","Il6","Hypoxia",
                                    "E2f","Uv_response_up"), collapse="|"), genesetNames, value=F)]
setdiff(c(blockA,blockB,blockC,blockD,blockE,blockF, blockG), genesetNames)
setdiff(genesetNames, c(blockA,blockB,blockC,blockD,blockE,blockF,blockG))

blockAenriched <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), blockA)
blockBenriched <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), blockB)
blockCenriched <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), blockC)
blockDenriched <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), blockD)
blockEenriched <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), blockE)
blockFenriched <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), blockF)
blockGenriched <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), blockG)
mergeBlocks <- cbind(blockAenriched$PROBE, blockBenriched$PROBE, blockCenriched$PROBE, blockDenriched$PROBE,blockEenriched$PROBE,blockFenriched$PROBE,blockGenriched$PROBE)
# write.csv(mergeBlocks, file="DifferentialExpression/GSEA/deltaNES/mergeBlocks.csv")

#  external program work:   GSEA performed using Broad GSEA java applet 

mergeBlocksResults <- read.csv(file="DifferentialExpression/GSEA/deltaNES/mergeBlocks_mergedStudies.csv", stringsAsFactors = F)
ggplot(mergeBlocksResults, aes(x=NAME, y=NES)) + theme_bw()   +  geom_boxplot(width=0.5) + 
  geom_point(size=6, pch=21, fill="black", colour="white") + # facet_wrap(~COHORT) + 
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5)) +
  ggtitle("NES by Geneset Grouping") + theme(plot.title = element_text(size=24,hjust = 0.5)) # + geom_line(aes( group=COHORT))


#  *************************     d7 enrichment aging signatures  *********************************   


d7signaturesHiHi <- rownames(hallmark_subsets[order(hallmark_subsets$Hi_v2, decreasing = T), ])[c(1:15,17:20)]  # only FDR <0.05 through 20th row not including TGFbeta (FDR 0.051)
d7signaturesLoLo <- rownames(hallmark_subsets[order(hallmark_subsets$Lo_v2, decreasing = T), ])[c(1:2)]  # only FDR <0.05 through 2nd row, #3 is myogenesis which is FDR 0.068 
d7signaturesNaive <- rownames(hallmark_subsets[order(hallmark_subsets$Nav_v2, decreasing = T), ])[c(1:9, 12)]  # only FDR <0.05 through 2nd row, #3 is myogenesis which is FDR 0.068

d7signaturesHiHigenes <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), d7signaturesHiHi)  # using HiHi but doesn't matter
d7signaturesLoLogenes <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), d7signaturesLoLo)  # which directory since 
d7signaturesNaivegenes <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), d7signaturesNaive)# gene lists are all the same
setdiff(d7signaturesHiHigenes$PROBE, d7signaturesLoLogenes$PROBE)
mergeDay7 <- cbind(d7signaturesHiHigenes$PROBE, d7signaturesLoLogenes$PROBE, d7signaturesNaivegenes$PROBE, c(blockEenriched$PROBE, blockFenriched$PROBE, blockGenriched$PROBE) )  # separate column for deltaNES genes from hihi
colnames(mergeDay7) <- c("d7_hihi", "d7_lolo", "d7_naive", "blockEFG")
# write.csv(mergeDay7, file="DifferentialExpression/GSEA/AgingSignature/mergeDay7signatures.csv")



d0signaturesHiHi <- rownames(hallmark_subsets[order(hallmark_subsets$Hi_v1, decreasing = T), ])[c(1:2,4)]  # only FDR <0.05 through 20th row not including TGFbeta (FDR 0.051)
d0signaturesLoLo <- rownames(hallmark_subsets[order(hallmark_subsets$Lo_v1, decreasing = T), ])[c(1:3)]  # only FDR <0.05 through 2nd row, #3 is myogenesis which is FDR 0.068 
d0signaturesNaive <- rownames(hallmark_subsets[order(hallmark_subsets$Nav_v1, decreasing = T), ])[c(1)]  # only FDR <0.05 through 2nd row, #3 is myogenesis which is FDR 0.068

d0signaturesHiHigenes <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), d0signaturesHiHi)  # using HiHi but doesn't matter
d0signaturesLoLogenes <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), d0signaturesLoLo)  # which directory since 
d0signaturesNaivegenes <- returnGeneListBlocks(c("DifferentialExpression/GSEA/GSEA_Results/Y_hihi_vs_E_hihi_v2_HALLMARK.GseaPreranked.1531000206530/"), d0signaturesNaive)# gene lists are all the same

mergeDay0 <- cbind(d0signaturesHiHigenes$PROBE, d0signaturesLoLogenes$PROBE, d0signaturesNaivegenes$PROBE )  
colnames(mergeDay0) <- c("d0_hihi", "d0_lolo", "d0_naive")
# write.csv(mergeDay0, file="DifferentialExpression/GSEA/AgingSignature/mergeDay0signatures.csv")





#  ************************* Differential expression for Ab hi vs low responders, controlling for subject *********************************   DID NOT YIELD USEFUL RESULT


# ****************************************************** Subdivide cohort into neutAb strong vs weak responders  ***************************************

# phenotypeGSVA <- read.csv(file = "phenotypeGSVA.csv", stringsAsFactors = F)
a <- cor.test(phenotypeGSVA$H3N2.nAb.FCd28, phenotypeGSVA$H1N1.nAb.FCd28, use="complete")
annotationInfo <- paste0("r = ", round(a$estimate,2), "\n", "p = ", formatC(a$p.value, format="e", digits=1))
my_grob1 = grobTree(textGrob(annotationInfo, x=0.65,  y=0.15, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(phenotypeGSVA, aes(x=H3N2.nAb.FCd28, y=H1N1.nAb.FCd28)) + geom_point(size=5) + theme_bw() + 
  theme(axis.text = element_text(size=16,hjust = 0.5))+theme(axis.title = element_text(size=28,hjust = 0.5)) +
  ggtitle("Fold-change in neutAb") + theme(plot.title = element_text(size=24,hjust = 0.5)) + # xlim(-1,1) + ylim(-1,1)
  geom_smooth(method=lm, se=F, fullrange=T, size=2, alpha=0.1)  + annotation_custom(my_grob1)   

# now split cohort into above-median and below-median neutAb responses  


neutAbHigh <- rownames(phenotypeGSVA[which(phenotypeGSVA$H3N2.nAb.FCd28 > median(phenotypeGSVA$H3N2.nAb.FCd28)),])
neutAbLow <- rownames(phenotypeGSVA[which(phenotypeGSVA$H3N2.nAb.FCd28 < median(phenotypeGSVA$H3N2.nAb.FCd28)),])


#  ************************* Differential expression for Ab hi vs low responders, controlling for subject *********************************   DID NOT YIELD USEFUL RESULT

metaData_Abresponders <- data.frame(row.names=colnames(bestDataCounts));  metaData_Abresponders$condition <- "empty"
#  1:36 are young, then 37:83 are elderly

metaData_Abresponders$ageGroup <- c(rep("Y",36),rep("E",47))
metaData_Abresponders$subject <- substr(rownames(metaData_Abresponders), start=2,stop=7)
metaData_Abresponders$condition <- substr(rownames(metaData_Abresponders),start=9, stop=20)
metaData_Abresponders$responders <- NA
metaData_Abresponders$responders[grep(x=metaData_Abresponders$subject, pattern=paste(neutAbHigh, collapse="|"))] <- "Hi"
metaData_Abresponders$responders[grep(x=metaData_Abresponders$subject, pattern=paste(neutAbLow, collapse="|"))] <- "Low"
metaData_Abresponders$subgroup <- paste0(metaData_Abresponders$condition, "_", metaData_Abresponders$responders)

fullDataset_Abresponders <- DESeqDataSetFromMatrix(countData=bestDataCounts, colData=metaData_Abresponders, design= ~subgroup)
dds <- estimateSizeFactors(fullDataset_Abresponders)
temp <- counts(dds,normalized=TRUE)
idx <- rowSums( counts(dds, normalized=TRUE) >= 20 ) >= 9  # filter for genes with at least 20 counts in 25% of the remaining samples
fullDataset_Abresponders <- fullDataset_Abresponders[idx,]

DESdata_Abresponders <- DESeq(fullDataset_Abresponders, parallel=TRUE)

diffExpr_Abresponders <- list(
  HiHi_v2_LowvHi = results(DESdata_Abresponders, contrast=c("subgroup","HiHi_v2_Low","HiHi_v2_Hi"), parallel = T), 
  LoLo_v2_LowvHi = results(DESdata_Abresponders, contrast=c("subgroup","LoLo_v2_Low","LoLo_v2_Hi"), parallel = T),
  Naive_v2_LowvHi = results(DESdata_Abresponders, contrast=c("subgroup","Naive_v2_Low","Naive_v2_Hi"), parallel = T)     )

# for (i in 1:length(diffExpr_Abresponders)) { write.csv(diffExpr_Abresponders[[i]], file=paste0("DifferentialExpression/Abresponders_",names(diffExpr_Abresponders[i]),".csv")) }




## *****************************       GSEA results: external Genesets - GSE79396 zoster vaccine and aging  VISIT 1  ***********************************************

# all together at baseline
ExternalGenesetsPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/GSE79396_ZosterVacc_YvE_visit1.GseaPreranked.1555447501664/gsea_report_for_na_pos_1555447501664.xls", sep="\t")
ExternalGenesetsNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/GSE79396_ZosterVacc_YvE_visit1.GseaPreranked.1555447501664/gsea_report_for_na_neg_1555447501664.xls", sep="\t")
ExternalGenesets <- rbind(ExternalGenesetsPos, ExternalGenesetsNeg)

BAA_deltaNES_HiHi <- read.csv("DifferentialExpression/GSEA/GSEA_Results/GSE79396_ZosterVacc_YvE_visit1.GseaPreranked.1555447501664/BAA_DELTANES_LEADEDGE_HIHI.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(ExternalGenesets$NES[grep("HIHI",ExternalGenesets$NAME)],2), "\n", "FDR: ", formatC(ExternalGenesets$FDR.q.val[grep("HIHI",ExternalGenesets$NAME)], format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.68,  y=0.87, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=BAA_deltaNES_HiHi, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE79396 Zoster Vac: ICOS+CD38+ cTfh") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
ggsave(file="DifferentialExpression/GSEA/Images/GSE79396_deltaNES_HiHi.pdf", device="pdf", height=3.5, width=5)

BAA_deltaNES_LoLo <- read.csv("DifferentialExpression/GSEA/GSEA_Results/GSE79396_ZosterVacc_YvE_visit1.GseaPreranked.1555447501664/BAA_DELTANES_LEADEDGE_LOLO.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(ExternalGenesets$NES[grep("LOLO",ExternalGenesets$NAME)],2), "\n", "FDR: ", formatC(ExternalGenesets$FDR.q.val[grep("LOLO",ExternalGenesets$NAME)], format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.68,  y=0.87, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=BAA_deltaNES_LoLo, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE79396 Zoster Vac: ICOS-CD38- cTfh") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
ggsave(file="DifferentialExpression/GSEA/Images/GSE79396_deltaNES_LoLo.pdf", device="pdf", height=3.5, width=5)

BAA_deltaNES_Naive <- read.csv("DifferentialExpression/GSEA/GSEA_Results/GSE79396_ZosterVacc_YvE_visit1.GseaPreranked.1555447501664/BAA_DELTANES_LEADEDGE_NAIVE.xls", sep="\t")
annotationInfo <- paste0("NES: ", round(ExternalGenesets$NES[grep("NAIVE",ExternalGenesets$NAME)],2), "\n", "FDR: ", formatC(ExternalGenesets$FDR.q.val[grep("NAIVE",ExternalGenesets$NAME)], format="e", digits=1))
my_grob = grobTree(textGrob(annotationInfo, x=0.68,  y=0.87, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(data=BAA_deltaNES_Naive, aes(x=`RANK.IN.GENE.LIST`, y=`RUNNING.ES`) ) + geom_line(color="black", size=1) + geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
  ggtitle("GSE79396 Zoster Vac: Naive CD4") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=12,hjust = 0.5))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  annotation_custom(my_grob) + geom_hline(yintercept = 0)
ggsave(file="DifferentialExpression/GSEA/Images/GSE79396_deltaNES_Naive.pdf", device="pdf", height=3.5, width=5)



ggplot(ExternalGenesets, aes(x=NAME, y=NES)) +  theme_bw() + #geom_point() #+ 
  geom_bar(stat="identity", width=0.3) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme(axis.text = element_text(size=12))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  ggtitle("NES of subsetes for leading edges")  
ggsave(filename = "DifferentialExpression/GSEA/Images/GSE79396_deltaNES_visit1_barGraph.pdf")




## *****************************       GSEA results: external Genesets - GSE79396 zoster vaccine and aging  VISIT 4  ***********************************************


# all together at baseline
ExternalGenesetsPos <- read.csv("DifferentialExpression/GSEA/GSEA_Results/GSE79396_ZosterVacc_YvE_visit4.GseaPreranked.1555447506183/gsea_report_for_na_pos_1555447506183.xls", sep="\t")
ExternalGenesetsNeg <- read.csv("DifferentialExpression/GSEA/GSEA_Results/GSE79396_ZosterVacc_YvE_visit4.GseaPreranked.1555447506183/gsea_report_for_na_neg_1555447506183.xls", sep="\t")
ExternalGenesets <- rbind(ExternalGenesetsPos, ExternalGenesetsNeg)

ggplot(ExternalGenesets, aes(x=NAME, y=NES)) +  theme_bw() + #geom_point() #+ 
  geom_bar(stat="identity", width=0.3) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme(axis.text = element_text(size=12))+theme(axis.title = element_text(size=14,hjust = 0.5))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  ggtitle("NES of subsetes for leading edges")  
ggsave(filename = "DifferentialExpression/GSEA/Images/GSE79396_deltaNES_visit4_barGraph.pdf")



