#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)

dataset = arguments[1]

cds = readRDS(paste0('checkpoints/',dataset,'_reclustered.rds'))

system(paste0('mkdir -p figures/',dataset,'_clusters'))

for (i in 1:nlevels(clusters(cds))) {
	colors.vector = rep('#CCCCCC',nlevels(clusters(cds)))
	colors.vector[i] = '#FF0000'

	p = plot_cells(cds, color_cells_by='cluster', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
		coord_equal() +
		guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
		scale_color_manual(values=colors.vector) +
		theme_classic(base_size=16) + theme(legend.position='none',axis.text=element_blank(),axis.ticks=element_blank())
	ggsave(p,file=paste0('figures/',dataset,'_clusters/cluster_',dataset,'_',formatC(i,width=3,flag=0),'.png'))
}
