#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')
source('scripts/_include_data.R')

library(monocle3)
library(tidyverse)

dataset = arguments[1]
if (length(arguments) > 1) {
	marker.gene.dataset = arguments[2]
} else {
	marker.gene.dataset = marker_genes
}

cds = readRDS(paste0('checkpoints/',dataset,'_reclustered.rds'))

system(paste0('mkdir -p figures/',dataset,'_marker_genes'))

marker.gene.df = get(marker.gene.dataset)



# p = plot_cells(cds, color_cells_by=predictor, group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
# 	geom_point(aes(color=cell_color,alpha=assigned_cell_type),size=0.2,stroke=0.1) + coord_equal() +
# 	scale_alpha_manual(values=rep(0.25,length(unique(cds$assigned_cell_type)))) +
# 	viridis::scale_color_viridis(option='C') +
#  	guides(color = FALSE, alpha = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
#  	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
# ggsave(p,file=paste0('figures/',dataset,'_clustering_',predictor,'_legend.png'),height=5)

for (i in 1:nrow(marker.gene.df)) {
	this = marker.gene.df[i,]
	cat(paste0('Now processing ',this$gene,'.\n'))
	if (this$gene %in% rowData(cds)$gene_short_name) {
		p = plot_cells(cds,genes=this$gene,label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
			coord_equal() +
			viridis::scale_color_viridis(option='viridis', 
                name = expression(log[10]*'(Expression)'), na.value = 'grey80', 
                end = 0.8, alpha = 0.25) +
			guides() +
			theme_classic(base_size=16) +
			theme(axis.text=element_blank(),axis.ticks=element_blank())
		suppressMessages(ggsave(p,file=gsub(' ','_',paste0('figures/',dataset,'_marker_genes/',marker.gene.dataset,'_',this$type,'_',this$gene,'.png')),height=5))
	} else {
		warning('Gene ',this$gene,' not found in dataset.')
	}
}
