#define general functions for usage in different scripts
#ggplot theme
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
#standard QC plots for scRNA processing
PlotUMIvsGENE <- function(barcode_metrics,col,colname,xintercept=NA,yintercept=NA,title=NULL) {
  ggExtra::ggMarginal(ggplot(data=barcode_metrics)+ 
                        geom_point(aes(x=log10(nCount),
                                       y=log10(nFeature),
                                       col=col),
                                   alpha=.5) + 
                                   {if(!is.na(xintercept)) geom_vline(xintercept = xintercept)} +
                                   {if(!is.na(yintercept)) geom_hline(yintercept = yintercept)} +
                                   {if(!is.null(title)) ggtitle(title)} +
                        scale_x_continuous(breaks=seq(0,10,0.5)) +
                        scale_y_continuous(breaks=seq(0,10,0.5)) +
                        labs(col=colname) +
                        theme_jb() +
                        theme(legend.position = "bottom"),
                      type="densigram", bins=100)
}

PlotUMIvsMITO <- function(barcode_metrics,col,colname,xintercept=NA,yintercept=NA,title=NULL) {
  ggExtra::ggMarginal(ggplot(data=barcode_metrics) + 
                        geom_point(aes(x=log10(nCount),
                                       y=percent.mito,
                                       col=col),
                                   alpha=0.5) + 
                                   {if(!is.na(xintercept)) geom_vline(xintercept = xintercept)} +
                                   {if(!is.na(yintercept)) geom_hline(yintercept = yintercept)} +
                                   {if(!is.null(title)) ggtitle(title)} +
                        scale_x_continuous(breaks=seq(0,10,0.5)) +
                        scale_y_continuous(breaks=seq(0,100,5), limits = c(0,100)) +
                        labs(col=colname) +
                        theme_jb() +
                        theme(legend.position="bottom"),
                      type="densigram",
                      xparams = list(bins=100),
                      yparams = list(binwidth=1))
}

PlotUMIvsHK <- function(barcode_metrics,col,colname,xintercept=NA,yintercept=NA,title=NULL) {
  ggExtra::ggMarginal(ggplot(data=barcode_metrics) + 
                        geom_point(aes(x=log10(nCount),
                                       y=n.exp.hkgenes,
                                       col=col),
                                   alpha=0.5) + 
                                   {if(!is.na(xintercept)) geom_vline(xintercept = xintercept)} +
                                   {if(!is.na(yintercept)) geom_hline(yintercept = yintercept)} +
                                   {if(!is.null(title)) ggtitle(title)} +
                        scale_x_continuous(breaks=seq(0,10,0.5)) +
                        scale_y_continuous(breaks=seq(0,100,5),limits = c(0,100)) +
                        labs(col=colname) +
                        theme_jb() +
                        theme(legend.position="bottom"),
                      type="densigram",
                      xparams = list(bins=100),
                      yparams = list(binwidth=1))
}

PlotHKvsMITO <- function(barcode_metrics,col,colname,xintercept=NA,yintercept=NA,title=NULL) {
  ggExtra::ggMarginal(ggplot(data=barcode_metrics) + 
                        geom_point(aes(x=n.exp.hkgenes,
                                       y=percent.mito,
                                       col=col),
                                   alpha=0.5) + 
                                   {if(!is.na(xintercept)) geom_vline(xintercept = xintercept)} +
                                   {if(!is.na(yintercept)) geom_hline(yintercept = yintercept)} +
                                   {if(!is.null(title)) ggtitle(title)} +
                        scale_y_continuous(breaks=seq(0,100,5),limits = c(0,100)) +
                        scale_x_continuous(breaks=seq(0,100,5),limits = c(0,100)) +
                        labs(col=colname) +
                        theme_jb() +
                        theme(legend.position="bottom"),
                      type="densigram",
                      xparams = list(binwidth=1),
                      yparams = list(binwidth=1))
}

#easier FeaturePlots
FeaturePlot_markers_comb <- function(SeuratObject,genes,reduction,ncol=2) {
  DefaultAssay(SeuratObject) <- "RNA"
  markerplot <- FeaturePlot(SeuratObject, 
                            features = genes, 
                            cols = c("grey95", "blue"),
                            pt.size = 1,
                            order = T,
                            combine = T,
                            ncol = ncol,
                            reduction = reduction,
                            slot = "data") + 
    labs(x="UMAP_1",y="UMAP_2") + 
    theme_jb_nogrid()
  markerplot
}

#celltype markers
celltypemarkers <- list()
#lymphocytes
#B-Cells & Plasma cells
celltypemarkers[["BCells"]] <- c("CD79A","MS4A1","CD19")
celltypemarkers[["BCellsubtypes"]] <- c("IGHD","IGHM","CD19","CD27","CD38","CD24","CR2")
celltypemarkers[["plasmablasts"]] <- c("IGHG1","IGHG2","IGHG3","IGHG4","IGKC")

#T-Cells
celltypemarkers[["TCells"]] <- c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B")
celltypemarkers[["TCellsNaive"]] <- c("TCF7")
celltypemarkers[["TCellsRegulatory"]] <- c("CD4","FOXP3","IL32","CTLA4")
celltypemarkers[["TCellsCytotoxic"]] <- c("CD8A","CD8B","GZMB","GZMH","GNLY","TNF","IFNG","NKG7")
celltypemarkers[["TCellsHelper"]] <- c("CD4","KLRB1","CCR4","CCR8","CXCR5","CXCR3")
celltypemarkers[["TCellsGammaDelta"]] <- c("TRGV9")
celltypemarkers[["NKCells"]] <- c("GNLY","NKG7","CD16")

#myeloid cells
celltypemarkers[["mastcells"]] <- c("KIT","TPSAB1","CPA3")

celltypemarkers[["cDCs"]] <- c("HLA-DQA1","FCER1A")
celltypemarkers[["cDC1"]] <- c("CLEC9A","FLT3","IDO1","XCER1","CADM1")
celltypemarkers[["cDC2"]] <- c("CD1C","FCER1A","HLA-DQA1","CD1A","SIRPA")
celltypemarkers[["cDCmature"]] <- c("LAMP3","CCR7","FSCN1")

celltypemarkers[["pDCs"]] <- c("LILRA4","GZMB","IL3RA","IRF4","SOX4","JCHAIN")

celltypemarkers[["monocytes"]] <- c("CD14","FCGR3A","CD68","CD163")
celltypemarkers[["nonclassicalmonocytes"]] <- c("CD14","FCGR3A","LST1","LILRB2")
celltypemarkers[["macrophages"]] <- c("CD68","CD163")

celltypemarkers[["RBCs"]] <- c("HBA1","HBA2","HBB")
celltypemarkers[["platelets"]] <- c("PPBP","ITGA2B")

#Stroma
celltypemarkers[["ECs"]] <- c("VWF","IFITM3")
celltypemarkers[["fib"]] <- c("VIM","DCN","PDPN","ACTA2","FAP","TAGLN","PDGFRB","FN1")
celltypemarkers[["fibECM"]] <- c("COL1A1","MMP1","MMP2","MMP3","MMP9","MMP10","MMP11","MMP14","MMP19")
celltypemarkers[["fibImmune"]] <- c("PLA2G2A","CXCL14","CXCL12","C3","C7","CCL19","CCL21")
celltypemarkers[["myf"]] <- c("ACTA2","MCAM","MYH11","MYLK")
celltypemarkers[["peri"]] <- c("RGS5")
celltypemarkers[["epithelial"]] <- c(paste0("KRT",c("6A","6B","6C","14","16","DAP")),"S100A7","S100A8","S100A9","FABP5")

# We only evaluate the foldchanges and PCT differences for differential expression and do not apply any p-value threshold, we do not calculate p-values.
calculate_foldchanges <- function(SeuratObject,idents) {
  Idents(SeuratObject) <- idents
  SeuratObject.clus <- lapply(unique(Idents(SeuratObject)), function(clus) {
    temp <- FoldChange(SeuratObject, 
                       features = NULL,
                       pseudocount.use = 1,
                       ident.1 = clus, 
                       ident.2 = NULL,
                       assay = "RNA",
                       slot = "data")
    temp$gene <- rownames(temp)
    temp$cluster <- clus
    temp
  })
  names(SeuratObject.clus) <- unique(Idents(SeuratObject))
  SeuratObject.clus <- do.call(rbind,SeuratObject.clus)
  return(SeuratObject.clus)
}

#volcano plots
ggvolcano <- function(data,cutoff,onlypos=F,pctdiff=F) {
  if(pctdiff == F) {
    if (onlypos == T) {
      if(length(cutoff) != 1) {stop("Please give only one positive cutoff if only positive values are shown.")}
      data <- data[data$avg_log2FC>0,]
      ggplot(data = data,
             aes(x=avg_log2FC,
                 y=pct.1-pct.2,
                 col=avg_log2FC>cutoff)) + 
        geom_point() + 
        geom_text_repel(aes(label=ifelse(avg_log2FC>cutoff,as.character(gene),"")),max.overlaps = 100) +
        scale_color_manual(values = c("TRUE"="red", "FALSE"="grey")) +
        scale_x_continuous(breaks=seq(-10,10,0.5)) +
        scale_y_continuous(breaks=seq(0,1,0.1)) +
        coord_cartesian(ylim = c(0,1), xlim=c(0,max(data$avg_log2FC)+0.1)) + #symmetric x-axis and y-axis with short buffer
        theme_jb() +
        theme(legend.position = "none", panel.grid.minor = element_blank())
    } else {
      cutoff_low <- cutoff[1]
      cutoff_up <- cutoff[2]
      ggplot(data = data,
             aes(x=avg_log2FC,
                 y=pct.1-pct.2,
                 col=avg_log2FC>cutoff_up | avg_log2FC<cutoff_low)) + 
        geom_point() + 
        geom_text_repel(aes(label=ifelse(avg_log2FC>cutoff_up | avg_log2FC<cutoff_low,as.character(gene),"")),max.overlaps = 100) +
        scale_color_manual(values = c("TRUE"="red", "FALSE"="grey")) +
        scale_x_continuous(breaks=seq(-10,10,0.5)) +
        scale_y_continuous(breaks=seq(-1,1,0.1)) +
        coord_cartesian(ylim = c(-1,1), xlim=c(-max(abs(c(min(data$avg_log2FC),max(data$avg_log2FC))))-0.1,max(abs(c(min(data$avg_log2FC),max(data$avg_log2FC))))+0.1)) + #symmetric x-axis and y-axis with short buffer
        geom_hline(yintercept = 0, col="black", size=0.1) +
        geom_vline(xintercept = 0, col="black", size=0.1) +
        theme_jb() +
        theme(legend.position = "none", panel.grid.minor = element_blank())
    }
  } else {
    if (onlypos == T) {
      if(length(cutoff) != 1) {stop("Please give only one positive cutoff if only positive values are shown.")}
      data <- data[data$avg_log2FC>0,]
      ggplot(data = data,
             aes(x=avg_log2FC,
                 y=pct.1-pct.2,
                 col=pct.1-pct.2>cutoff)) + 
        geom_point() + 
        geom_text_repel(aes(label=ifelse(pct.1-pct.2>cutoff,as.character(data$gene),"")),max.overlaps = 100) +
        scale_color_manual(values = c("TRUE"="red", "FALSE"="grey")) +
        scale_x_continuous(breaks=seq(0,10,0.5)) +
        scale_y_continuous(breaks=seq(0,1,0.1)) +
        coord_cartesian(ylim = c(0,1), xlim=c(0,max(data$avg_log2FC)+0.1)) + #symmetric x-axis and y-axis with short buffer
        theme_jb() +
        theme(legend.position = "none", panel.grid.minor = element_blank())
    } else {
      cutoff_low <- cutoff[1]
      cutoff_up <- cutoff[2]
      ggplot(data = data,
             aes(x=avg_log2FC,
                 y=pct.1-pct.2,
                 col=pct.1-pct.2>cutoff_up | pct.1-pct.2<cutoff_low)) + 
        geom_point() + 
        geom_text_repel(aes(label=ifelse(pct.1-pct.2>cutoff_up | pct.1-pct.2<cutoff_low,as.character(data$gene),"")),
                        max.overlaps = 100) +
        scale_color_manual(values = c("TRUE"="red", "FALSE"="grey")) +
        scale_x_continuous(breaks=seq(-10,10,0.5)) +
        scale_y_continuous(breaks=seq(-1,1,0.1)) +
        coord_cartesian(ylim = c(-1,1), xlim=c(-max(abs(c(min(data$avg_log2FC),max(data$avg_log2FC))))-0.1,max(abs(c(min(data$avg_log2FC),max(data$avg_log2FC))))+0.1)) + #symmetric x-axis and y-axis with short buffer
        geom_hline(yintercept = 0, col="black", size=0.1) +
        geom_vline(xintercept = 0, col="black", size=0.1) +
        theme_jb() +
        theme(legend.position = "none", panel.grid.minor = element_blank())
    }
  }
}

GSEAvis <- function(GSEAdf,
                    npathways=10,
                    pvalthres=NULL,
                    NESthres=NULL,
                    sort="padj",
                    positive=T) {
  if (!is.null(pvalthres) | !is.null(NESthres) | !is.null(npathways)) {
    if (!is.null(pvalthres)) {
      if (min(GSEAdf$padj) < pvalthres) {
        GSEAdf <- GSEAdf[GSEAdf$padj < pvalthres,]
      } else {
        stop("no gene sets below p-value threshold.")
      }
    }
    if (!is.null(NESthres)) {
      if(NESthres>0) {
        if (max(GSEAdf$NES) > NESthres) {
          GSEAdf <- GSEAdf[GSEAdf$NES > NESthres,]
        } else {
          stop("no gene sets above specified enrichment score threshold.")
        }
      } else {
        if (min(GSEAdf$NES) < NESthres) {
          GSEAdf <- GSEAdf[GSEAdf$NES < NESthres,]
        } else {
          stop("no gene sets above specified enrichment score threshold.")
        }
      }
    }
    if (positive==T) {
      GSEAdf <- GSEAdf[GSEAdf$NES>0,]
    } else {
      GSEAdf <- GSEAdf[GSEAdf$NES<0,]
    }
    #sort
    if(sort=="padj") {
      GSEAdf <- GSEAdf[order(GSEAdf$padj),]
    } 
    if(sort=="NES") {
      if(positive==T) {
        GSEAdf <- GSEAdf[order(GSEAdf$NES, decreasing = T),]
      } else {
        GSEAdf <- GSEAdf[order(GSEAdf$NES, decreasing = F),]
      }
    }
    if (nrow(GSEAdf)>npathways) {
      GSEAdf <- GSEAdf[1:npathways,]
    } else {
      warning(paste0("Less than ",npathways," pathways were present."))
    }
    
    #plotting
    GSEAdf$pathway <- factor(GSEAdf$pathway, levels=rev(GSEAdf$pathway))
    ggplot(data=GSEAdf) +
      geom_bar(aes(x=pathway,y=NES,fill=-log10(padj)),
               stat="identity") +
      scale_x_discrete(position = "top") +
      coord_flip() +
      #scale_y_continuous(breaks = seq(0,1,0.1)) +
      theme_jb() +
      theme(axis.text.y = element_text(size = 9))
    
  } else {
    warning("You have to specify either a p-value threshold, enrichment Score threshold or npathways to plot.")
  }  
}
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  intersection/union
}

colorfunc <- function(types,brewer=F) {
  n <- length(types)
  if (brewer==T & n<=9) {
    ifelse(n>9,
           colors <- rainbow(n = n),
           colors <- brewer.pal(n,"Set1"))
    
  } else {
    colors <- scales::hue_pal()(n)
  }
  names(colors) <- types
  colors
}
histbarplot <- function(anno1,anno2, relative=F) {
  anno.count <- table(anno1=anno1, anno2=anno2)
  anno.count.melt <-reshape2::melt(anno.count)
  
  if(relative==F) {
    #absolute barplot
    ggplot(data=anno.count.melt) +
      geom_bar(aes(x=anno1,
                   y=value,
                   fill=anno2),
               stat="identity",
               position="stack") +
      #scale_y_continuous(breaks = seq(0,1,0.1)) +
      labs(x="",y="fraction of cells",fill="types") +
      coord_flip() +
      theme_jb()
  } else {
    #relative barplot
    ggplot(data=anno.count.melt) +
      geom_bar(aes(x=anno1,
                   y=value,
                   fill=anno2),
               stat="identity",
               position=position_fill(reverse = T)) +
      scale_y_continuous(breaks = seq(0,1,0.1)) +
      labs(x="",y="fraction of cells",fill="types") +
      coord_flip() +
      theme_jb()
  }
}
piechart <- function(anno1,anno2) {
  
  anno.count <- table(anno1=anno1, anno2=anno2)
  anno.count.melt <-reshape2::melt(anno.count)
  
  p1 <- list()
  anno.count.melt.sub <- list()
  for (cl in unique(anno1)) {
    anno.count.melt.sub[[cl]] <- anno.count.melt[anno.count.melt$anno1==cl,]
    
    p1[[cl]] <- ggplot(data=anno.count.melt.sub[[cl]],aes(x="",y=value,fill=anno2)) +
      geom_bar(stat="identity", color="white") +
      coord_polar("y",start=0) +
      ggtitle(cl) +
      theme_void()
  }
  p1
}
