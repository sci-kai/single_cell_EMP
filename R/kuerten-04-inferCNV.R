library(Seurat)
library(infercnv)
library(ggpubr)
library(ComplexHeatmap)

OUTPUT = "~/tscr/kai/projects/01_SCC_scRNA/results/pubdatasets/kuertendata/HPVneg_CD45n_analysis/"
OUTPUTINFERCNV <- paste0(OUTPUT,"inferCNV/")
load(paste0(OUTPUT,"RData/HNSCC-filt.RData"))

#For every patient
cell.names <- list()
counts_matrix <- list()
annotation_df <- list()
infer_cnvs_obj <- list()
infer_cnvs <- list()
for (p in unique(HNSCC.filt$patients)[-1:-2]) {
  set.seed(42)
  cell.names[[p]] <-  names(HNSCC.filt$CellType_SCTccmito_PCA[HNSCC.filt$patients==p])
  counts_matrix[[p]] <- GetAssayData(HNSCC.filt, assay = "RNA", slot = "counts")[,cell.names[[p]]]
  annotation_df[[p]] <- HNSCC.filt@meta.data[cell.names[[p]], "CellType_SCTccmito_PCA", drop = F]
  write.table(annotation_df[[p]], paste0(OUTPUTINFERCNV,p,"-annotation.txt"), sep = "\t", col.names = F, row.names = T, quote = F)
  
  infer_cnvs_obj[[p]] <- CreateInfercnvObject(raw_counts_matrix = as.matrix(counts_matrix[[p]]),
                                              annotations_file = paste0(OUTPUTINFERCNV,p,"-annotation.txt"),
                                              delim = "\t",
                                              gene_order_file = "~/tscr/Peiffer/2 Projects/10x/hg19_genes_noDuplicates.bed",
                                              ref_group_names = as.character(unique(HNSCC.filt$CellType_SCTccmito_PCA[HNSCC.filt$patients==p])[unique(HNSCC.filt$CellType_SCTccmito_PCA[HNSCC.filt$patients==p])!="Epithelial"]))
  dir.create(paste0(OUTPUTINFERCNV,"/",p,"-infercnv"))
  infer_cnvs[[p]] <- infercnv::run(infer_cnvs_obj[[p]],
                                   cutoff = 0.1,
                                   out_dir = paste0(OUTPUTINFERCNV,p,"-infercnv"),
                                   cluster_by_groups = T,
                                   denoise = T,
                                   HMM = F,
                                   analysis_mode='samples')
}

#infer_cnvs <- list()
celltypepat <- list()
plotcells <- list()
infercnv.data <- list()
infercnv.data.auto <- list() 
infercnv.all.hm <- list()
#for (p in unique(HNSCC.filt$patients)[-1]) {
  #infer_cnvs[[p]] <- readRDS(paste0(OUTPUTINFERCNV,p,"-infercnv/","run.final.infercnv_obj"))
#}  
for (p in unique(HNSCC.filt$patients)) {
  celltypepat[[p]] <- as.character(HNSCC.filt$CellType_SCTccmito_PCA[HNSCC.filt$patients==p])
  names(celltypepat[[p]]) <- names(HNSCC.filt$CellType_SCTccmito_PCA[HNSCC.filt$patients==p])
  set.seed(42)
  plotcells[[p]] <- c(celltypepat[[p]][celltypepat[[p]]=="Epithelial"],sample(celltypepat[[p]][celltypepat[[p]]=="T-Cells"],replace = F,size = min(c(sum(celltypepat[[p]]=="Epithelial"),sum(celltypepat[[p]]=="T-Cells")))))
  #all cells
  infercnv.data[[p]] <- t(infer_cnvs[[p]]@expr.data)
  # Only use autosomes
  infercnv.data.auto[[p]] <- infercnv.data[[p]][,rownames(infer_cnvs[[p]]@gene_order)[!as.character(infer_cnvs[[p]]@gene_order$chr) %in% c("MT")]]
  # Heatmap with standard hierarchical clustering
  
  #Heatmap with selected cells
  infercnv.all.hm[[p]] <- Heatmap(infercnv.data.auto[[p]][names(plotcells[[p]]),],
                                  name="log2",
                                  col = circlize::colorRamp2(breaks = c(min(infercnv.data.auto[[p]]),0.85,1,1.15,max(infercnv.data.auto[[p]])),
                                                             colors = c("darkblue","blue","white","red","darkred")),
                                  show_column_names = F,
                                  show_row_names = F,
                                  #column clustering and split
                                  column_split = factor(paste0("chr",as.character(infer_cnvs[[p]]@gene_order$chr)[!as.character(infer_cnvs[[p]]@gene_order$chr) %in% c("MT")]), levels=paste0("chr",c(as.character(1:22),"X","Y"))),
                                  column_order = rownames(infer_cnvs[[p]]@gene_order)[!as.character(infer_cnvs[[p]]@gene_order$chr) %in% c("MT")],
                                  column_title_rot = 90,
                                  column_title_side = "bottom",
                                  cluster_columns = F,
                                  cluster_column_slices = F,
                                  #row clustering and split
                                  row_split = plotcells[[p]],
                                  row_title_rot = 0,
                                  cluster_rows = T,
                                  cluster_row_slices = F,
                                  clustering_method_rows =  "ward.D2",
                                  clustering_distance_rows = "euclidean",
                                  #graphic options
                                  row_names_gp = gpar(fontsize=11),
                                  column_names_gp = gpar(fontsize=11),
                                  row_title_gp = gpar(fontsize=11),
                                  heatmap_legend_param = list(title_gp = gpar(fontsize=11), labels_gp = gpar(fontsize=11)),
                                  #row_dend_width = unit("1","in"),
                                  border = T,
                                  border_gp = gpar(lwd=0.5),
                                  use_raster = T,
                                  raster_device = "CairoPNG",
                                  raster_quality = 4,
                                  width = unit("4.375","in"),
                                  height = unit("5","in"))
  
  pdf(paste0(OUTPUT,"figs/HNSCC-all-inferCNVauto-",p,"-allhm.pdf"), width=20, height=20)
  print(infercnv.all.hm[[p]])
  dev.off()
  
}

load("~/tscr/kai/projects/01_SCC_scRNA/results/alldatasets/RData/alltumor-perpatients-perpatient.RData")

phentype <- list()
infercnv.all.hm.pt <- list()
for (p in grep("HN",names(SCC.p),value = T)) {
  phentype[[p]] <- SCC.p[[p]]$pattumortype[unlist(lapply(strsplit(names(SCC.p[[p]]$pattumortype),split = "_"),function(x) gsub(", ","_",toString(x[-length(x)])))) %in% names(plotcells[[p]][plotcells[[p]]=="Epithelial"])]
  names(phentype[[p]]) <- unlist(lapply(strsplit(names(phentype[[p]]),split = "_"),function(x) gsub(", ","_",toString(x[-length(x)]))))
  #add T-Cells and only EMP-related cells
  phentype[[p]] <- c(phentype[[p]][c(grep("epi",phentype[[p]]),grep("pEMT",phentype[[p]]))],plotcells[[p]][plotcells[[p]]=="T-Cells"])
  
  #Heatmap with selected cells
  infercnv.all.hm.pt[[p]] <- Heatmap(infercnv.data.auto[[p]][names(phentype[[p]]),],
                                     name="log2",
                                     col = circlize::colorRamp2(breaks = c(min(infercnv.data.auto[[p]]),0.85,1,1.15,max(infercnv.data.auto[[p]])),
                                                                colors = c("darkblue","blue","white","red","darkred")),
                                     show_column_names = F,
                                     show_row_names = F,
                                     #column clustering and split
                                     column_split = factor(paste0("chr",as.character(infer_cnvs[[p]]@gene_order$chr)[!as.character(infer_cnvs[[p]]@gene_order$chr) %in% c("MT")]), levels=paste0("chr",c(as.character(1:22),"X","Y"))),
                                     column_order = rownames(infer_cnvs[[p]]@gene_order)[!as.character(infer_cnvs[[p]]@gene_order$chr) %in% c("MT")],
                                     column_title_rot = 90,
                                     column_title_side = "bottom",
                                     cluster_columns = F,
                                     cluster_column_slices = F,
                                     #row clustering and split
                                     row_split = phentype[[p]],
                                     row_title_rot = 0,
                                     cluster_rows = T,
                                     cluster_row_slices = F,
                                     clustering_method_rows =  "ward.D2",
                                     clustering_distance_rows = "euclidean",
                                     #graphic options
                                     row_names_gp = gpar(fontsize=11),
                                     column_names_gp = gpar(fontsize=11),
                                     row_title_gp = gpar(fontsize=11),
                                     heatmap_legend_param = list(title_gp = gpar(fontsize=11), labels_gp = gpar(fontsize=11)),
                                     #row_dend_width = unit("1","in"),
                                     border = T,
                                     border_gp = gpar(lwd=0.5),
                                     use_raster = T,
                                     raster_device = "CairoPNG",
                                     raster_quality = 4,
                                     width = unit(4.375,"in"),
                                     height = unit(5,"in"))
  
  pdf(paste0(OUTPUT,"figs/HNSCC-all-inferCNVauto-",p,"-pt-small.pdf"), width=20, height=20)
  print(infercnv.all.hm.pt[[p]])
  dev.off()
}