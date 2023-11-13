suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(enrichR)
    library(rafalib)
    library(Matrix)
    library(edgeR)
    library(MAST)
})


# set data in/output
in_dir <- ''
out_dir <- ""

# Read in Seurat 
scrna.combined <- readRDS(paste0(in_dir,".rds"))
#scrna.combined$treatment <- gsub("^.*?Control/p35", "Control_p35", scrna.combined$treatment)
#scrna.combined$treatment <- gsub("^.*?TauR406W/p35", "TauR406W_p35", scrna.combined$treatment)
DefaultAssay(scrna.combined) <- "RNA"
# Split Seurat 
scrna.list <- SplitObject(scrna.combined, split.by = "treatment")
scrna.QTau_Control<- merge(scrna.list$QTau, y = scrna.list$Control, project = "QTau_Control")
# select genes that are expressed in at least 2 treatments or 2 ctrls with > 2 # reads.
# Test on one celltype 
# select all cells in cluster
scrna.QTau_Control$celltype <- Idents(scrna.QTau_Control)
celltype <- unique(scrna.QTau_Control$celltype)
selected_cells <- list()
cell_selection <- list()
counts         <- list()
log_counts     <- list()
fData          <- list()
m              <- list()
cData          <- list()
sca            <- list()
treatment      <- list()
zlmCond        <- list()
summaryCond    <- list()
fcHurdle       <- list()
summaryDt      <- list()
fcHurdle1      <- list()
top50          <- list()
mastN          <- list()
p1             <- list()
p2             <- list()
sig_genes      <- list()
for (cell in celltype){
     selected_cells[[cell]] <- WhichCells(scrna.QTau_Control, idents = cell)
     cell_selection[[cell]] <- subset(scrna.QTau_Control, cells = selected_cells[[cell]])
     counts[[cell]] <- cell_selection[[cell]]@assays$RNA@counts
     log_counts[[cell]] <- log(counts[[cell]] + 1) / log(2)
     fData[[cell]] <- data.frame(primerid = rownames(log_counts[[cell]]))
     m[[cell]] = cell_selection[[cell]]@meta.data
     m[[cell]]$wellKey = rownames(m[[cell]])
     m[[cell]]$sample_id = factor(m[[cell]]$orig.ident)
     m[[cell]]$treatment = factor(m[[cell]]$treatment)
     cData[[cell]] <- m[[cell]]
     sca[[cell]] <- FromMatrix(as.matrix(log_counts[[cell]]), cData[[cell]], fData[[cell]])
     colData(sca[[cell]])$cngeneson <- scale(colSums(assay(sca[[cell]]) > 0))
     treatment<- factor(colData(sca[[cell]])$treatment)
     zlmCond[[cell]] <- suppressMessages(MAST::zlm(~treatment + cngeneson, sca[[cell]], method = "bayesglm", ebayes = T))
     dir.create(file.path(paste0(out_dir,cell,'_QTau_vs_Control')), recursive = TRUE)
     saveRDS(zlmCond[[cell]],paste0(out_dir,cell,"_QTau_vs_Control","/",cell,"_treatmentQTau.vs.Control_zlmCond.rds"))
     summaryCond[[cell]] <- suppressMessages(MAST::summary(zlmCond[[cell]], doLRT = "treatmentQTau"))
     saveRDS(summaryCond[[cell]],paste0(out_dir,cell,"_QTau_vs_Control","/",cell,"_treatmentQTau.vs.Control_summaryCond.rds"))
     #summaryCond[[cell]] <- readRDS(paste0(in_dir,cell,"_treatmentTauR406W.vs.Control_p35_summaryCond.rds"))
     summaryDt[[cell]] <- summaryCond[[cell]]$datatable
     fcHurdle[[cell]] <- merge(summaryDt[[cell]][summaryDt[[cell]]$contrast == "treatmentQTau" & summaryDt[[cell]]$component ==
                               "logFC", c(1, 7, 5, 6, 8)], summaryDt[[cell]][summaryDt[[cell]]$contrast == "treatmentQTau" &
                               summaryDt[[cell]]$component == "H", c(1, 4)], by = "primerid")
     fcHurdle1[[cell]] <- stats::na.omit(as.data.frame(fcHurdle[[cell]]))
     top50[[cell]] = head(fcHurdle1[[cell]][order(fcHurdle1[[cell]]$`Pr(>Chisq)`), ], 50)
     fcHurdle1[[cell]]$pval = fcHurdle1[[cell]]$`Pr(>Chisq)`
     fcHurdle1[[cell]]$dir = ifelse(fcHurdle1[[cell]]$z > 0, "up", "down")
     fcHurdle1[[cell]] %>% group_by(dir) %>% top_n(-50, pval) %>% arrange(z) -> mastN[[cell]]
     write.csv(fcHurdle1[[cell]],paste0(out_dir,cell,"_QTau_vs_Control","/","gene_list_QTau_vs._Control.csv"))
     sig_genes[[cell]] <- fcHurdle1[[cell]] %>% filter(pval <= 0.05)
     sig_genes[[cell]] <- sig_genes[[cell]][order(sig_genes[[cell]]$pval, decreasing = FALSE), ] 
     write.csv(sig_genes[[cell]],paste0(out_dir,cell,"_QTau_vs_Control","/","sig_gene_list_QTau_vs._Control.csv"))
     mastN[[cell]] = mastN[[cell]]$primerid
}    

for (cell in cells){
pdf(paste0(out_dir,cell,"_QTau_vs_Control","/","DotPlot_QTau_vs.Control.pdf"),height=30,width=20)
print(DotPlot(cell_selection[[cell]], features = mastN[[cell]], group.by = "sample_id", assay = "RNA") +
    coord_flip() + RotatedAxis() + ggtitle(paste0("MAST ", cell, " QTau vs. Control")))
dev.off()
}

for (cell in cells){
pdf(paste0(out_dir,cell,"_QTau_vs_Control","/","QTau_vs.Control_VolcanoPlot.pdf"),height=15,width=15)
print(EnhancedVolcano(fcHurdle1[[cell]], lab = fcHurdle1[[cell]]$primerid, x = 'coef', y = 'pval', title = paste0("MAST ",cell, " QTau vs. Control"), FCcutoff = 0.1))
dev.off()
}
