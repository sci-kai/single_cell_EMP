library(Seurat)

theme_jb <- function(base_size=11) {
  (theme_bw() +
     theme(
       axis.title = element_text(size=11),
       axis.text = element_text(size=11),
       legend.text = element_text(size=11),
       #optional rotation of x-axis labels
       axis.text.x = element_text(angle=90, hjust=1,vjust=.5)
     ))
}
theme_jb_nogrid <- function(base_size=11) {
  (theme_bw() +
     theme(
       axis.title = element_text(size=11),
       axis.text = element_text(size=11),
       legend.text = element_text(size=11),
       #optional rotation of x-axis labels
       axis.text.x = element_text(angle=90, hjust=1,vjust=.5),
       panel.grid = element_blank(),
       panel.border = element_blank(),
       axis.line = element_line()
     ))
}
 
############
#Prepare cellphonedb input
############

#load all cells
OUTPUT = "~/tscr/kai/projects/01_SCC_scRNA/results/samples/cohort/D19200/rp/"
load(paste0(OUTPUT,"RData/D19200-filt.RData"))

#get annotation from fibroblasts and tumor subsets
load(paste0(OUTPUT,"RData/D19200-fibroendo.RData"))
load(paste0(OUTPUT,"RData/D19200-tumornofibmono.RData"))

#exclude those cells misclustered in fibroblasts/tumorcells
removecells <- c(names(D19200.filt$CellType[D19200.filt$CellType %in% c("Fibroblasts","Endothelial Cells") & !names(D19200.filt$CellType) %in% names(D19200.ef$orig.ident)]),
                 names(D19200.filt$CellType[D19200.filt$CellType=="Tumor" & !names(D19200.filt$CellType) %in% names(D19200.tumor$orig.ident)]))
keepcells <- names(D19200.filt$CellType)[!names(D19200.filt$CellType) %in% removecells]
D19200.filt2 <- subset(D19200.filt,cells = keepcells)

#renormalize RNA assay
#RNA normalization without regressing, necessary for example for differential gene expression not based on SCT.
DefaultAssay(D19200.filt2) <- "RNA"
#log-normalization for slot "data"
D19200.filt2 <- NormalizeData(D19200.filt2,
                              normalization.method = "LogNormalize", #default
                              scale.factor = 10000, #default, scale factor 10000 is recommendation from seurat tutorial
                              assay = "RNA",
                              margin = 1 # default; normalizes across features
)
#scaling for slot "scale.data"
allgenes <- rownames(D19200.filt2)
D19200.filt2 <- ScaleData(D19200.filt2,
                          features = allgenes,
                          do.scale = T,
                          do.center = T,
                          scale.max = 10,
                          assay = "RNA")

#rewrite dendritic and endothelial cells
D19200.filt2$CellType <- as.character(D19200.filt2$CellType)
D19200.filt2$CellType[D19200.filt2$CellType=="Dendritic Cells"] <- "DendriticCells"
D19200.filt2$CellType[D19200.filt2$CellType=="Endothelial Cells"] <- "EndothelialCells"

#add cluster annotation from tumor
D19200.tumor$trajectorycluster <- as.character(D19200.tumor$seurat_clusters_tumor_SCT_PCA)
D19200.tumor$trajectorycluster[D19200.tumor$trajectorycluster=="0"] <- "6 temp"
D19200.tumor$trajectorycluster[D19200.tumor$trajectorycluster=="1"] <- "3 temp"
D19200.tumor$trajectorycluster[D19200.tumor$trajectorycluster=="2"] <- "4 temp"
D19200.tumor$trajectorycluster[D19200.tumor$trajectorycluster=="3"] <- "2 temp"
D19200.tumor$trajectorycluster[D19200.tumor$trajectorycluster=="4"] <- "1 temp"
D19200.tumor$trajectorycluster[D19200.tumor$trajectorycluster=="5"] <- "7 temp"
D19200.tumor$trajectorycluster[D19200.tumor$trajectorycluster=="6"] <- "5 temp"
D19200.tumor$trajectorycluster <- sub(" temp","",D19200.tumor$trajectorycluster)
D19200.tumor$trajectorycluster <- factor(D19200.tumor$trajectorycluster)

D19200.filt2$CellTypeclusters <- as.character(D19200.filt2$CellType)
D19200.filt2$CellTypeclusters[names(D19200.tumor$trajectorycluster)] <- paste0("Tumor",as.character(D19200.tumor$trajectorycluster))

#save TXT and CSV files to run with cellphonedb
#save log-normalized RNA assay for counts
RNAdata <- as.data.frame(D19200.filt2@assays$RNA@data)
write.table(RNAdata,file = paste0(OUTPUT,"cellphonedb/RNAdata.txt"),quote = F,sep = "\t")
#manually add "gene" as first column name to output!
#save annotation only of celltypes
annotation <- data.frame(Cell=names(D19200.filt2$CellType),cell_type=D19200.filt2$CellType)
write.table(annotation,file = paste0(OUTPUT,"cellphonedb/celltype_annotation.txt"),quote = F,sep = "\t",row.names = F)
#annotation with clusters separately
annotationclusters <- cbind(names(D19200.filt2$CellTypeclusters),D19200.filt2$CellTypeclusters)
write.table(annotationclusters,file = paste0(OUTPUT,"cellphonedb/celltypeclusters_annotation.txt"),quote = F,sep = "\t",row.names = F)

############
#run cellphonedb with python environment
############

############
#ANALYSIS with all celltype 
############

means <- read.delim("~/tscr/kai/projects/01_SCC_scRNA/results/samples/cohort/D19200/rp/cellphonedb/celltype/means.txt")
pvalues <- read.delim("~/tscr/kai/projects/01_SCC_scRNA/results/samples/cohort/D19200/rp/cellphonedb/celltype/pvalues.txt")

#change 0 p-values to close to zero (as in plot_dot script from Teichlab)
pvalues[pvalues==0] = 0.0009
#in this script they change mean interactions from 0 to 1 to be 0 in log2-conversion. I do not agree with this practice, since this changes the values significantly.
#instead, change value to be the minimum non-zero value (so the smallest value)
means[means==0] <- min(means[,-(1:11)][means[,-(1:11)]!=0])

#Heatmap with number of interactions
#calculate number of interactions: those interactions with p-value less than 0.05 are counted as interacting
countint <- data.frame(colSums(pvalues[,-(1:11)]<=0.05))

rownames(countint) <- gsub("Endothelial.Cells","EndothelialCells",rownames(countint))
rownames(countint) <- gsub("Dendritic.Cells","DendriticCells",rownames(countint))
names(countint) <-  "count"

countint.long <- cbind(data.frame(do.call(rbind,strsplit(rownames(countint),split = "[.]"))),countint)
names(countint.long) <- c("partnerA","partnerB","count")

for(x in unique(countint.long$partnerA)) {
  pairs1 <- countint.long[countint.long$partnerB==x & countint.long$partnerA!=x,]
  pairs2 <- countint.long[countint.long$partnerA==x & countint.long$partnerB!=x,]
  
  if(sum(pairs1$partnerA == pairs2$partnerB)==4) {
    countint.long[countint.long$partnerB==x & countint.long$partnerA!=x,]$count <- pairs1$count + pairs2$count
    countint.long[countint.long$partnerA==x & countint.long$partnerB!=x,]$count <- pairs1$count + pairs2$count
  } else {
    warning("Exit.")
  }
}

countint.wide <- reshape2::dcast(countint.long,partnerA ~ partnerB,value.var = "count")
rownames(countint.wide) <- countint.wide$partnerA
countint.wide <- countint.wide[,-1]
#rows are the first, columns the second partner
#first partner: expressed
library(ComplexHeatmap)
pdf(file = paste0(OUTPUT,"figs-D19200/D19200-cpdb-intcounts.pdf"),width=10,height=10)
Heatmap(countint.wide,
        col = circlize::colorRamp2(breaks = c(0,max(countint.wide)),colors = c("white","black")),
        cluster_rows = F,
        cluster_columns = F,
        row_order = c("DendriticCells","Tumor","EndothelialCells","Fibroblasts","Macrophages"),
        column_order = c("DendriticCells","Tumor","EndothelialCells","Fibroblasts","Macrophages"),
        #graphics parameter
        row_names_gp = gpar(fontsize=11),
        column_names_gp = gpar(fontsize=11),
        row_title_gp = gpar(fontsize=11),
        heatmap_legend_param = list(title_gp = gpar(fontsize=11), labels_gp = gpar(fontsize=11)),
        border = T,
        use_raster = F,
        width = unit(5,"in"),
        height = unit(5,"in")
)
dev.off()

#Interaction Dotplots
#filter for those interactions with at at least one p-value below 0.05
pvalues.sig <- pvalues[rowSums(pvalues[,-(1:11)]>=0.05) != ncol(pvalues[,-(1:11)]>=0.05),]
means.sig <- means[rowSums(pvalues[,-(1:11)]>=0.05) != ncol(pvalues[,-(1:11)]>=0.05),]

##reformat from wide to long format
means.m <- reshape2::melt(means.sig[,-c(1,3:11)])
pvalues.m <- reshape2::melt(pvalues.sig[,-c(1,3:11)])

intact <- as.data.frame(cbind(interacting_pair=means.m$interacting_pair,
                              partner=as.character(means.m$variable),
                              mean=means.m$value,
                              pval=pvalues.m$value))
intact$mean <- as.numeric(intact$mean)
intact$pval <- as.numeric(intact$pval)
intact$partner <- gsub("Endothelial.Cells","EndothelialCells",intact$partner)
intact$partner <- gsub("Dendritic.Cells","DendriticCells",intact$partner)

library(ggplot2)
#Interaction Dotplot
intact.map <- ggplot(data=intact,aes(x=interacting_pair, y=partner, col=log2(mean), size=-log10(pval))) + 
  geom_point() + 
  scale_y_discrete(position = "right") +
  theme_jb() + 
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))
ggsave(filename = paste0(OUTPUT,"figs-D19200/D19200-cpdb-intmap.svg"),intact.map,width=30,height=10)

#special look at EGFR interactions
intact.map.specific <- ggplot(data=intact[unique(c(grep("EGFR",intact$interacting_pair),grep("VEGFA",intact$interacting_pair),grep("TGFB",intact$interacting_pair))),],
                              aes(x=interacting_pair, y=partner, col=log2(mean), size=-log10(pval))) + 
  geom_point() + 
  scale_y_discrete(position = "right") +
  theme_jb() + 
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))
ggsave(filename = paste0(OUTPUT,"figs-D19200/D19200-cpdb-intmap-specific.svg"),intact.map.specific,width=10,height=10)

############
#ANALYSIS with tumor populations
############

#analysis
means.sub <- read.delim("~/tscr/kai/projects/01_SCC_scRNA/results/samples/cohort/D19200/rp/cellphonedb/allclusters/means.txt")
pvalues.sub <- read.delim("~/tscr/kai/projects/01_SCC_scRNA/results/samples/cohort/D19200/rp/cellphonedb/allclusters/pvalues.txt")

#change 0 p-values to close to zero (as in plot_dot script from Teichlab)
pvalues.sub[pvalues.sub==0] = 0.0009
#in this script they change mean interactions from 0 to 1 to be 0 in log2-conversion. I do not agree with this practice, since this changes the values significantly.
#instead, change value to be the minimum non-zero value (so the smallest value)
means.sub[means.sub==0] <- min(means.sub[,-(1:11)][means.sub[,-(1:11)]!=0])

#filter for those interactions with at least one p-value less than 1
pvalues.sub.sig <- pvalues.sub[rowSums(pvalues.sub[,-(1:11)]>=0.05) != ncol(pvalues.sub[,-(1:11)]>=0.05),]
means.sub.sig <- means.sub[rowSums(pvalues.sub[,-(1:11)]>=0.05) != ncol(pvalues.sub[,-(1:11)]>=0.05),]

means.sub.m <- reshape2::melt(means.sub.sig[,-c(1,3:11)])
pvalues.sub.m <- reshape2::melt(pvalues.sub.sig[,-c(1,3:11)])

intact.sub <- as.data.frame(cbind(interacting_pair=means.sub.m$interacting_pair,
                                  partner=as.character(means.sub.m$variable),
                                  mean=means.sub.m$value,
                                  pval=pvalues.sub.m$value))

intact.sub$mean <- as.numeric(intact.sub$mean)
intact.sub$pval <- as.numeric(intact.sub$pval)

intact.sub$partner <- gsub("Endothelial.Cells","EndothelialCells",intact.sub$partner)
intact.sub$partner <- gsub("Dendritic.Cells","DendriticCells",intact.sub$partner)

intact.sub.map <- ggplot(data=intact.sub,aes(x=interacting_pair, y=partner, col=log2(mean), size=-log10(pval))) + 
  geom_point() + 
  scale_y_discrete(position = "right") +
  theme_jb() + 
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))
ggsave(filename = paste0(OUTPUT,"figs-D19200/D19200-cpdb-allclusters-intmap.png"),intact.sub.map,width=30,height=30)

#special look at specific interactions
intact.sub.map.specific <- ggplot(data=intact.sub[unique(c(grep("EGFR",intact.sub$interacting_pair),grep("VEGFA",intact.sub$interacting_pair),grep("TGFB",intact.sub$interacting_pair))),],
                                  aes(x=interacting_pair,y=partner, col=log2(mean), size=-log10(pval))) + 
  geom_point() + 
  scale_y_discrete(position = "right") +
  theme_jb() + 
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))
ggsave(filename = paste0(OUTPUT,"figs-D19200/D19200-cpdb-allclusters-intmap-specific.svg"),intact.sub.map.specific,width=10,height=30)

#Count Interactions and show as lineplot
#calculate number of interactions: those interactions with p-value less than 0.05 are counted as interacting
countint.sub <- data.frame(colSums(pvalues.sub[,-(1:11)]<=0.05))
rownames(countint.sub) <- gsub("Endothelial.Cells","EndothelialCells",rownames(countint.sub))
rownames(countint.sub) <- gsub("Dendritic.Cells","DendriticCells",rownames(countint.sub))
rownames(countint.sub) <- gsub("Myoblast.like","Myoblast-like",rownames(countint.sub))
rownames(countint.sub) <- gsub("Transition.CAF1.2","Transition-CAF1-2",rownames(countint.sub))
names(countint.sub) <-  "count"

countint.sub.long <- cbind(data.frame(do.call(rbind,strsplit(rownames(countint.sub),split = "[.]"))),countint.sub)
names(countint.sub.long) <- c("partnerA","partnerB","count")

for(x in unique(countint.sub.long$partnerA)) {
  pairs1.sub <- countint.sub.long[countint.sub.long$partnerB==x & countint.sub.long$partnerA!=x,]
  pairs2.sub <- countint.sub.long[countint.sub.long$partnerA==x & countint.sub.long$partnerB!=x,]
  
  if(sum(pairs1.sub$partnerA == pairs2.sub$partnerB)==length(pairs1.sub$partnerA)) {
    countint.sub.long[countint.sub.long$partnerB==x & countint.sub.long$partnerA!=x,]$count <- pairs1.sub$count + pairs2.sub$count
    countint.sub.long[countint.sub.long$partnerA==x & countint.sub.long$partnerB!=x,]$count <- pairs1.sub$count + pairs2.sub$count
  } else {
    warning("Exit.")
  }
}

countint.sub.wide <- reshape2::dcast(countint.sub.long,partnerA ~ partnerB,value.var = "count")
rownames(countint.sub.wide) <- countint.sub.wide$partnerA
countint.sub.wide <- countint.sub.wide[,-1]

countint.sub.long.TME <- reshape2::melt(cbind(partnerA=rownames(countint.sub.wide[grep("Tumor",rownames(countint.sub.wide),invert = T),grep("Tumor",rownames(countint.sub.wide),invert = F)]),
                     countint.sub.wide[grep("Tumor",rownames(countint.sub.wide),invert = T),grep("Tumor",rownames(countint.sub.wide),invert = F)]))
names(countint.sub.long.TME) <- c("partnerA","partnerB","count")

countint.sub.long.TME$partnerB <- sub("Tumor","",as.character(countint.sub.long.TME$partnerB))

countint.sub.TME.plot <- ggplot(data=countint.sub.long.TME,
       aes(x=partnerB,y=count)) +
  geom_point(aes(col=partnerA)) +
  geom_line(aes(group=partnerA,col=partnerA)) +
  scale_y_continuous(breaks=seq(0,200,10)) +
  theme_jb()

ggsave(filename = paste0(OUTPUT,"figs-D19200/D19200-cpdb-celltypeclusters-intcounts-TumortoTME-lineplot.svg"),countint.sub.TME.plot,width=10,height=5)
ggsave(filename = paste0(OUTPUT,"figs-D19200/D19200-cpdb-celltypeclusters-intcounts-TumortoTME-lineplot-noleg.svg"),countint.sub.TME.plot + theme(legend.position = "none") + ggtitle(""),width=10,height=5)

#sessionInfo
writeLines(capture.output(sessionInfo()), paste0(OUTPUT,"RData/18-D19200-cellphonedb-sessionInfo.txt"))
