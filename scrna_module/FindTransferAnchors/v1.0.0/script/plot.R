library(Seurat)
library(ggplot2)

sto <- readRDS("/work/share/acuhtwkcu9/taoxiao/04_Project/P2024050717K4YH7T/taoxiao/Analysis/sto/12_HnEC_reAnalysis/result/05_yingshe/FindTransferAnchors/HnEC/HnEC_spot_cell.rds")
#SpatialDimPlot(sto)

pancreas.query <- sto
Idents(pancreas.query) <- pancreas.query$predicted.id
#SpatialDimPlot(pancreas.query, label.size=1,, stroke = NA) + theme(legend.position = "right")

pdf("SpatialDimPlot.pdf")
SpatialDimPlot(pancreas.query, label.size=1,, stroke = NA) + theme(legend.position = "right")
dev.off()