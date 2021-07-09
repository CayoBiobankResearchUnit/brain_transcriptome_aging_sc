#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)

dataset = arguments[1]
if (length(arguments) > 1) {
	doublet.rate = as.numeric(arguments[2])
} else {
	doublet.rate = NULL
}

cds = readRDS(paste0('checkpoints/cell_data_sets/',dataset,'.rds'))

write(paste0(nrow(colData(cds)),' cells to start.'),file=paste0('reports/',dataset,'_initial_qc.txt'))

cds = cds %>%
	{ scrub_low_quality(.,umi.min.cutoff=umi.min,mt.max.cutoff=mt.max) } %>%
	{ find_doublets(.,collision.rate=doublet.rate) }

write(paste0(nrow(colData(cds)),' cells pass with minimum UMI ',umi.min,' and maximum mitochondrial percentage ',mt.max,'%.'),file=paste0('reports/',dataset,'_initial_qc.txt'),append=TRUE)

system(paste0('mkdir -p figures/',dataset,'_doublets'))

plot_doublets(cds,score.column='Scr.kNN',outdir=paste0('figures/',dataset,'_doublets'))

plot_doublets(cds,score.column='scrublet_score',outdir=paste0('figures/',dataset,'_doublets'))

p = ggplot(as.data.frame(colData(cds)),aes(scrublet_score,Scr.kNN)) +
	geom_point(size=0.2,alpha=0.2,shape=21,fill=NA,stroke=0.1) +
	facet_wrap(~sample,nrow=sample.rows) +
	theme_classic() +
	coord_equal() +
	xlab('BBI') +
	ylab('Manual')
ggsave(p,file=paste0('figures/',dataset,'_doublets/manual_bbi_scrublet_comparison.pdf'),useDingbats=FALSE)

saveRDS(cds,file=paste0('checkpoints/',dataset,'_doublets_included.rds'))
