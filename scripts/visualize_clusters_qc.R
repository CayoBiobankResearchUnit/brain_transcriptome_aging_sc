#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)
library(RColorBrewer)

dataset = arguments[1]

# Plot initial thresholding

cds = readRDS(paste0('checkpoints/cell_data_sets/',dataset,'.rds'))

p = ggplot(as.data.frame(colData(cds)),aes(log10(perc_mitochondrial_umis + 0.01))) +
	geom_histogram(bins=50) +
	geom_path(aes(x,y),data=data.frame(x=c(log10(mt.max + 0.01),log10(mt.max + 0.01)),y=c(0,5000)),color='blue') +
	theme_classic(base_size=24) +
	xlab(expression(log[10](percentage~mt~genes + 0.01))) + ylab('Count')
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_mitochondrial_cutoff.pdf'),useDingbats=FALSE)

p = ggplot(within(unique(as.data.frame(colData(cds))[c('id','total.umi')]),{
		id = factor(id,levels=id[order(total.umi)])
	}),aes(id,total.umi)) +
	geom_bar(stat='identity') +
	coord_flip() +
	scale_y_continuous(
		breaks=seq(0,(round(max(unique(as.data.frame(colData(cds))[c('id','total.umi')])$total.umi)/1e7) * 1e7),5e6),
		labels=seq(0,(round(max(unique(as.data.frame(colData(cds))[c('id','total.umi')])$total.umi)/1e7) * 1e7),5e6) / 1e6
	) +
	theme_classic(base_size=24) +
#	theme(axis.text.y=element_blank(),axis.ticks.y = element_blank()) +
	xlab('Sample') + ylab('Total UMI (millions)')
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_umi_distribution.pdf'),useDingbats=FALSE)

p = ggplot(within(reshape2::melt(table(colData(cds)$id)),{id=factor(Var1,levels=Var1[order(value)]); n.indexes=value}),aes(id,n.indexes)) +
	geom_bar(stat='identity') +
	coord_flip() +
	scale_y_continuous(
		breaks=seq(0,round(max(table(colData(cds)$id))/1e3),1) * 1e3,
		labels=seq(0,round(max(table(colData(cds)$id))/1e3),1)
	) +
	theme_classic(base_size=24) +
#	theme(axis.text.y=element_blank(),axis.ticks.y = element_blank()) +
	xlab('Sample') + ylab('Indexes (thousands)')
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_index_distribution.pdf'),useDingbats=FALSE)

# Plot doublet thresholding

cds = readRDS(paste0('checkpoints/',dataset,'_doublets_included.rds'))

# Read in manually set thresholds
doublet.thresholds = readRDS(paste0('checkpoints/doublet_thresholds_',dataset,'.rds'))


p = ggplot(as.data.frame(colData(cds)),aes(Scr.kNN)) +
	geom_density(size=0.2) +
#	geom_vline(aes(xintercept=threshold),data=doublet.thresholds) +
	geom_path(aes(threshold,y),data=do.call(rbind,lapply(1:nrow(doublet.thresholds),function(i) { 
		rbind(within(doublet.thresholds[i,],{y=0}),within(doublet.thresholds[i,],{y=3}))
	})),color='blue',size=1) +
	facet_wrap(~sample,ncol=sample.rows) +
	theme_classic(base_size=16) +
	theme(strip.background=element_blank()) +
	xlab('Scrublet score') + ylab('Density')
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_doublet_score_cutoffs.pdf'),useDingbats=FALSE)


p = ggplot(as.data.frame(colData(cds)),aes(Scr.kNN)) +
	geom_density(size=0.5) +
#	geom_vline(aes(xintercept=threshold),data=doublet.thresholds) +
	geom_path(aes(threshold,y),data=do.call(rbind,lapply(1:nrow(doublet.thresholds),function(i) { 
		rbind(within(doublet.thresholds[i,],{y=0}),within(doublet.thresholds[i,],{y=3}))
	})),color='blue',size=1) +
	facet_wrap(~sample,ncol=sample.rows*3/4) +
	theme_classic(base_size=24) +
	theme(axis.text=element_blank(),strip.text=element_blank()) +
	xlab('Scrublet score') + ylab('Density')
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_doublet_score_cutoffs_no_labels.pdf'),useDingbats=FALSE,height=7)

# Plot after removing doublets

cds = readRDS(paste0('checkpoints/',dataset,'_doublets_removed.rds'))

# Plot clusters by partition
p = plot_cells(cds, color_cells_by='partition', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title='Partition',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_partitions.pdf'),useDingbats=FALSE)

# By cluster
p = plot_cells(cds, color_cells_by='cluster', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_clusters.pdf'),useDingbats=FALSE)

library(ggrastr)

# By cluster
p = plot_cells(cds, color_cells_by='cluster', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
 	guides(color = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(legend.position='none',axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_clusters_article.pdf'),useDingbats=FALSE,height=5,width=10)


# By doublet score
p = plot_cells(cds, color_cells_by='Scr.kNN', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title='Scrublet score',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_scrublet_custom.pdf'),useDingbats=FALSE)

library(viridis)
# By doublet score
p = plot_cells(cds, color_cells_by='Scr.kNN', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
	scale_color_viridis(option='D',name='Scrublet score') +
# 	guides(color = guide_legend(title='Scrublet score')) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_scrublet_custom_article.pdf'),useDingbats=FALSE,height=5,width=10)

# # Age
# # By doublet score
# p = plot_cells(cds, color_cells_by='age', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
# 	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
# 	coord_equal() +
# 	scale_color_viridis(option='C',name='Age') +
# # 	guides(color = guide_legend(title='Scrublet score')) +
#  	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
# ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_age_custom_article.pdf'),useDingbats=FALSE)


# By doublet score 2
p = plot_cells(cds, color_cells_by='scrublet_score', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title='Scrublet score',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_scrublet_bbi.pdf'),useDingbats=FALSE)

# Plot doublet score distribution by cluster
p = ggplot(within(as.data.frame(colData(cds)),{
		cluster = factor(clusters(cds),levels=names(sort(tapply(colData(cds)$Scr.kNN,clusters(cds),mean),decreasing=TRUE)))
	}),aes(cluster,Scr.kNN)) +
	geom_violin(draw_quantiles=0.5) +
#	scale_color_manual(values=cell.colors) +
#	facet_wrap(~assigned_cell_type,ncol=1,scales='free_x') +
	theme_classic() +
#	theme(legend.position='bottom') +
	theme(axis.text.x=element_text(angle=-45,hjust=0)) +
	xlab('Cluster') +
	ylab('Scrublet score')
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_scrublet_by_cluster_custom.pdf'),useDingbats=FALSE,width=20)

# Plot cluster sizes
p = ggplot(within(as.data.frame(colData(cds)),{
		cluster = factor(clusters(cds),levels=names(sort(tapply(colData(cds)$Scr.kNN,clusters(cds),mean),decreasing=TRUE)))
	}),aes(cluster)) +
	geom_bar(stat='count') +
#	scale_color_manual(values=cell.colors) +
#	facet_wrap(~assigned_cell_type,ncol=1,scales='free_x') +
	theme_classic() +
#	theme(legend.position='bottom') +
	theme(axis.text.x=element_text(angle=-45,hjust=0)) +
	xlab('Cluster') +
	ylab('Count')
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_scrublet_by_cluster_custom_counts.pdf'),useDingbats=FALSE,width=20)

# Plot doublet score distribution by cluster
p = ggplot(within(as.data.frame(colData(cds)),{
		cluster = factor(clusters(cds),levels=names(sort(tapply(colData(cds)$scrublet_score,clusters(cds),mean),decreasing=TRUE)))
	}),aes(cluster,scrublet_score)) +
	geom_violin(draw_quantiles=0.5) +
#	scale_color_manual(values=cell.colors) +
#	facet_wrap(~assigned_cell_type,ncol=1,scales='free_x') +
	theme_classic() +
#	theme(legend.position='bottom') +
	theme(axis.text.x=element_text(angle=-45,hjust=0)) +
	xlab('Cluster') +
	ylab('Scrublet score')
ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_scrublet_by_cluster_bbi.pdf'),useDingbats=FALSE,width=20)

# # Plot manual cell type assignments
# p = plot_cells(cds, color_cells_by='assigned_cell_type', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
# 	coord_equal() +
# 	scale_color_manual(name='Cell type',values=cell.prediction.colors[cell.prediction.levels %in% cds$assigned_cell_type],drop=TRUE) +
# 	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
# 	theme_classic(base_size=16) +
# 	theme(axis.text=element_blank(),axis.ticks=element_blank())
# ggsave(p,file=paste0('figures/',dataset,'_clustering_qc_manual_cell_type.pdf'),useDingbats=FALSE,height=5)
# 
