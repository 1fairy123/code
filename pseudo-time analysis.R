rm(list=ls())
setwd("F:\\HLQ\\pseudo_time\\")
library(Seurat)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(dplyr)
library(ggplot2)
library(tidyr)
library(piano)
library(msigdbr)
library(monocle)

## load seurat object
MP <- readRDS("cov.B.rds")
data <- as(as.matrix(MP@assays$RNA@counts), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = MP@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

diff.wilcox = FindAllMarkers(MP)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
diff.genes <- subset(all.markers,p_val_adj<0.01)$gene

mycds <- setOrderingFilter(mycds, diff.genes)
p1 <- plot_ordering_genes(mycds)

var.genes <- VariableFeatures(MP)
mycds <- setOrderingFilter(mycds, var.genes)
p2 <- plot_ordering_genes(mycds)

disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
p1|p2|p3
plot_ordering_genes

library(reduceDimension)
library(Rcpp)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)
p1 = plot_cell_trajectory(mycds4, color_by = "State")
p2 = plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~State, nrow = 1)
p1|p2

p1 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters",cell_name_size=20,cell_size=1.5,cell_link_size=1.5)
p2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters",cell_name_size=10,cell_size=1,cell_link_size=1.5,cex.axis=5,cex.lab=5) + facet_wrap(~seurat_clusters, nrow = 1)
p1/p2
ggsave('1.1.pdf',p1/p2,width=7,height=6)

p1 <- plot_cell_trajectory(mycds, color_by = "celltype",cell_size=1.5,cell_link_size=1.5)
p2 <- plot_cell_trajectory(mycds, color_by = "celltype",cell_size=1,cell_link_size=1.5) + facet_wrap(~seurat_clusters, nrow = 1)
p1/p2
ggsave('1.2.pdf',p1/p2,width=7,height=6)

test = orderCells(mycds,root_state = 1) 
plot4 <- plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size=1,cell_link_size=1.5)
ggsave('1.3.pdf',plot4,width=7,height=6)
plot3 <- plot_cell_trajectory(mycds, color_by = "celltype")

disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

aa<-read.csv("15darkgenes.csv")
aaa<-aa[,1]
p3<-plot_pseudotime_heatmap(mycds[aaa,], num_clusters=4,
                        show_rownames=T, return_heatmap=T)
ggsave('1.4.pdf',p3,width=7,height=6)
