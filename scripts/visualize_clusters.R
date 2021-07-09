#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)
library(RColorBrewer)
library(ggrastr)

dataset = arguments[1]

cds = readRDS(paste0('checkpoints/',dataset,'_classified_manual.rds'))

# Plot clusters by partition
p = plot_cells(cds, color_cells_by='partition', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title='Partition',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_partitions.pdf'),useDingbats=FALSE)

# Plot clusters by partition
p = plot_cells(cds, color_cells_by='partition', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
 	guides(color = guide_legend(title='Partition',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_partitions_article.pdf'),useDingbats=FALSE)

# By cluster
p = plot_cells(cds, color_cells_by='cluster', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_clusters.pdf'),useDingbats=FALSE)

# By cluster
p = plot_cells(cds, color_cells_by='cluster', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	coord_equal() +
 	guides(color = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_clusters_article.pdf'),useDingbats=FALSE,height=5)

# By original cluster
p = plot_cells(cds, color_cells_by='monocle_cluster', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_original_clusters.pdf'),useDingbats=FALSE)

# By partitions 2
p = plot_cells(cds, color_cells_by='partition', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color,alpha=assigned_cell_type),size=0.2,stroke=0.1) + coord_equal() +
	scale_alpha_manual(values=rep(0.1,length(unique(cds$assigned_cell_type)))) +
	scale_color_manual(values=scales::hue_pal(c(0, 360), 100,65)(nlevels(partitions(cds)))[c(seq(2,nlevels(partitions(cds)),2),seq(1,nlevels(partitions(cds)),2))]) +
 	guides(color = FALSE, alpha = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_partitions_legend.pdf'),height=5)

# By cluster 2
p = plot_cells(cds, color_cells_by='cluster', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color,alpha=assigned_cell_type),size=0.2,stroke=0.1) + coord_equal() +
	scale_alpha_manual(values=rep(0.1,length(unique(cds$assigned_cell_type)))) +
	scale_color_manual(values=scales::hue_pal(c(0,360)+15,100,65,0,1)(nlevels(clusters(cds)))[c(which(as.logical((1:nlevels(clusters(cds))) %% 2)),which(!as.logical((1:nlevels(clusters(cds))) %% 2)))]) +
 	guides(color = FALSE, alpha = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_clusters_legend.pdf'),height=5)


	
# p = plot_cells(foo, color_cells_by='Scr.kNN', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
# 	geom_point(aes(color=Scr.kNN),size=0.2,stroke=0.1) + coord_equal() + scale_alpha_manual(values=rep(0.25,8)) + viridis::scale_color_viridis(breaks=c(seq(0,0.16,0.02),0.15),labels=c('0.00','0.02','0.04','0.06','0.08','0.10','0.12','0.14','0.16','Oligodendrocyte precursor cells'),option='C') +
#  	guides(alpha=FALSE,color=guide_legend(title='Scrublet score',override.aes = list(size=2,alpha=1))) +
#  	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
 
# By predictor
p = plot_cells(cds, color_cells_by=predictor, group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),size=0.2,alpha=0.02) +
	scale_color_viridis(option='C',name=predictor.label) +
	coord_equal() +
# 	guides(color = guide_legend(title=predictor.label,override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_',predictor,'_article.pdf'),useDingbats=FALSE,height=5)

p = plot_cells(cds, color_cells_by=predictor, group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title=predictor.label,override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_',predictor,'.pdf'),useDingbats=FALSE)

p = plot_cells(cds, color_cells_by=predictor, group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color,alpha=assigned_cell_type),size=0.2,stroke=0.1) + coord_equal() +
	scale_alpha_manual(values=rep(0.1,length(unique(cds$assigned_cell_type)))) +
	viridis::scale_color_viridis(option='C') +
 	guides(color = FALSE, alpha = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_',predictor,'_legend.pdf'),height=5)

p = plot_cells(cds, color_cells_by=predictor, group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.1,size=0.2,stroke=0.1) + coord_equal() +
	viridis::scale_color_viridis(name=predictor.label,option='C') +
 	guides(color = guide_legend(title=predictor.label,override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_',predictor,'_nolegend.pdf'),height=5)

p = plot_cells(cds, color_cells_by=predictor, group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color,alpha=assigned_cell_type),size=0.2,stroke=0.1) + coord_equal() +
	scale_alpha_manual(values=rep(0.1,length(unique(cds$assigned_cell_type)))) +
	viridis::scale_color_viridis(option='C') +
 	guides(color = guide_legend(title='Age',override.aes = list(size=2,alpha=1)), alpha = guide_legend(title='Cluster',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_',predictor,'_legend2.pdf'),height=5)

# By batch
p = plot_cells(cds, color_cells_by=batch.variable, group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title='Batch',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_batch.pdf'),useDingbats=FALSE)

# By doublet score
p = plot_cells(cds, color_cells_by='Scr.kNN', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title='Scrublet score',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_scrublet_custom.pdf'),useDingbats=FALSE)

p = plot_cells(cds, color_cells_by='Scr.kNN', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color,alpha=assigned_cell_type),size=0.2,stroke=0.1) +
	coord_equal() +
	viridis::scale_color_viridis(option='C') +
	scale_alpha_manual(values=rep(0.1,nlevels(cds$assigned_cell_type))) +
 	guides(color=FALSE,alpha = guide_legend(title='Scrublet score',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_scrublet_custom_simple.pdf'),height=5)

### Scrublet
p = plot_cells(cds, color_cells_by='Scr.kNN', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.02,size=0.2) +
	coord_equal() +
	viridis::scale_color_viridis(option='D',name='Scrublet score') +
# 	guides(color=FALSE,alpha = guide_legend(title='Scrublet score',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_scrublet_custom_article.pdf'),height=5)


cds$logUMI = log2(cds$n.umi)
### Scrublet
p = plot_cells(cds, color_cells_by='logUMI', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.02,size=0.2) +
	coord_equal() +
	viridis::scale_color_viridis(option='D',name=expression(log[2]('UMIs'))) +
# 	guides(color=FALSE,alpha = guide_legend(title='Scrublet score',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_umi_custom_article.pdf'),height=5)

p = plot_cells(cds, color_cells_by='perc_mitochondrial_umis', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.02,size=0.2) +
	coord_equal() +
	viridis::scale_color_viridis(option='D',name='Perc. mt reads') +
# 	guides(color=FALSE,alpha = guide_legend(title='Scrublet score',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_mt_custom_article.pdf'),height=5)

p = plot_cells(cds, color_cells_by='extraction_batch', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.02,size=0.2) +
	coord_equal() +
	scale_color_brewer(palette='Set1') +
#	viridis::scale_color_viridis(option='D',name='Perc. mt reads') +
 	guides(color = guide_legend(title='Batch',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_batch_custom_article.pdf'),height=5)

p = plot_cells(cds, color_cells_by='group', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.02,size=0.2) +
	coord_equal() +
	scale_color_brewer(palette='Set1') +
#	viridis::scale_color_viridis(option='D',name='Perc. mt reads') +
 	guides(color = guide_legend(title='Group',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_group_custom_article.pdf'),height=5)


# By doublet score 2
p = plot_cells(cds, color_cells_by='scrublet_score', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
 	guides(color = guide_legend(title='Scrublet score',override.aes = list(size=2,alpha=1))) +
 	theme_classic(base_size=16) + theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_scrublet_bbi.pdf'),useDingbats=FALSE)

# Plot doublet score distribution by cluster
p = ggplot(within(as.data.frame(colData(cds)),{
		cluster = factor(garnett_cluster,levels=sort(unique(garnett_cluster)))
	}),aes(factor(garnett_cluster),Scr.kNN,color=assigned_cell_type)) +
	geom_violin(draw_quantiles=0.5) +
	scale_color_manual(values=cell.colors) +
#	facet_wrap(~assigned_cell_type,ncol=1,scales='free_x') +
	theme_classic() +
	theme(legend.position='bottom') +
	xlab('Cluster') +
	ylab('Scrublet score')
ggsave(p,file=paste0('figures/',dataset,'_clustering_scrublet_by_cluster_custom.pdf'),useDingbats=FALSE,width=20)

# Plot doublet score distribution by cluster
p = ggplot(within(as.data.frame(colData(cds)),{
		cluster = factor(garnett_cluster,levels=sort(unique(garnett_cluster)))
	}),aes(factor(garnett_cluster),scrublet_score,color=assigned_cell_type)) +
	geom_violin(draw_quantiles=0.5) +
	scale_color_manual(values=cell.colors) +
#	facet_wrap(~assigned_cell_type,ncol=1,scales='free_x') +
	theme_classic() +
	theme(legend.position='bottom') +
	xlab('Cluster') +
	ylab('Scrublet score')
ggsave(p,file=paste0('figures/',dataset,'_clustering_scrublet_by_cluster_bbi.pdf'),useDingbats=FALSE,width=20)

# Plot garnett cell type predictions
p = plot_cells(cds, color_cells_by='cell_type', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
	scale_color_manual(name='Cell type',values=cell.prediction.colors[cell.prediction.levels %in% cds$cell_type],drop=TRUE) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	theme_classic(base_size=16) +
	theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_garnett_cell_type.pdf'),useDingbats=FALSE,height=5)

# Plot garnett cluster extended predictions
p = plot_cells(cds, color_cells_by='cluster_ext_type', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
	coord_equal() +
	scale_color_manual(name='Cell type',values=cell.prediction.colors[cell.prediction.levels %in% cds$cluster_ext_type],drop=TRUE) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	theme_classic(base_size=16) +
	theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_garnett_cluster_ext_type.pdf'),useDingbats=FALSE,height=5)

# Plot manual cell type assignments
p = plot_cells(cds, color_cells_by='assigned_cell_type', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0,group_label_size=3) +
	geom_point_rast(aes(color=cell_color),alpha=0.1,size=0.2,stroke=0.1) + coord_equal() +
	coord_equal() +
	scale_color_manual(name='Cell type',values=cell.prediction.colors[cell.prediction.levels %in% cds$assigned_cell_type],drop=TRUE) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	theme_classic(base_size=16) +
	theme(axis.text=element_blank(),axis.ticks=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_clustering_manual_cell_type.pdf'),useDingbats=FALSE,height=5)

# Plot again in PNG to control file size
ggsave(p,file=paste0('figures/',dataset,'_clustering_manual_cell_type.png'),height=5)

for (i in cell.levels) {
	# Plot manual cell type assignments
	this.colors = rep('#eeeeee',length(cell.prediction.levels))
	this.colors[match(i,cell.prediction.levels)] = cell.prediction.colors[match(i,cell.prediction.levels)]
	p = plot_cells(cds, color_cells_by='assigned_cell_type', group_cells_by='partition',label_cell_groups=F,cell_size=0.2,alpha=0.25,group_label_size=3) +
		coord_equal() +
		scale_color_manual(name='Cell type',values=this.colors[cell.prediction.levels %in% cds$assigned_cell_type],drop=TRUE) +
		guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
		theme_classic(base_size=16) +
		theme(axis.text=element_blank(),axis.ticks=element_blank())
	ggsave(p,file=paste0('figures/',dataset,'_clustering_manual_cell_type_',gsub(' ','_',tolower(i)),'.pdf'),useDingbats=FALSE,height=5)

	# Plot again in PNG to control file size
	ggsave(p,file=paste0('figures/',dataset,'_clustering_manual_cell_type_',gsub(' ','_',tolower(i)),'.png'),height=5)
}



# Plot cell type proportions by individual against predictor

cds.split = split(subset(as.data.frame(colData(cds)),assigned_cell_type %in% cell.levels),cds$animal_id)
sample.cell.proportions = do.call(rbind,lapply(cds.split,function(x) {
	nums = table(x$assigned_cell_type) / nrow(x)
	out = data.frame(id=unique(x$animal_id),group=unique(x$group),predictor=unique(x[[predictor]]),type=as.character(names(nums)), prop = as.numeric(nums) * 100,stringsAsFactors=FALSE)
	out
}))

sample.cell.proportions$type = factor(sample.cell.proportions$type,levels=cell.levels)
sample.cell.proportions$type.label = sample.cell.proportions$type
levels(sample.cell.proportions$type.label) = paste0(levels(sample.cell.proportions$type),unlist(lapply(split(sample.cell.proportions,sample.cell.proportions$type),function(x) ifelse (coef(summary(lm(prop~predictor+group,data=x)))['predictor','Pr(>|t|)'] < 0.05/nlevels(sample.cell.proportions$type),' (*)',''))))
p.adjust(unlist(lapply(split(sample.cell.proportions,sample.cell.proportions$type),function(x) coef(summary(lm(prop~predictor+group,data=x)))['predictor','Pr(>|t|)'])),'fdr')
lapply(split(sample.cell.proportions,sample.cell.proportions$type),function(x) coef(summary(lm(prop~predictor+group,data=x)))['predictor',])

sample.cluster.proportions = do.call(rbind,lapply(cds.split,function(x) {
	nums = table(x$assigned_cell_cluster) / nrow(x)
	out = data.frame(id=unique(x$animal_id),group=unique(x$group),predictor=unique(x[[predictor]]),cluster_type=as.character(names(nums)),type=gsub(' [0-9]+$','',as.character(names(nums))), prop = as.numeric(nums) * 100,stringsAsFactors=FALSE)
	out
}))

sample.cluster.proportions$type = factor(sample.cluster.proportions$type,levels=cell.levels)
sample.cluster.proportions$cluster_type = factor(sample.cluster.proportions$cluster_type,levels=levels(cds$assigned_cell_cluster))
sample.cluster.proportions$type.label = sample.cluster.proportions$cluster_type
levels(sample.cluster.proportions$type.label) = paste0(levels(sample.cluster.proportions$cluster_type),unlist(lapply(split(sample.cluster.proportions,sample.cluster.proportions$cluster_type),function(x) ifelse (coef(summary(lm(prop~predictor+group,data=x)))['predictor','Pr(>|t|)'] < 0.05/nlevels(sample.cluster.proportions$cluster_type) ,' (*)',''))))
p.adjust(unlist(lapply(split(sample.cluster.proportions,sample.cluster.proportions$cluster_type),function(x) coef(summary(lm(prop~predictor+group,data=x)))['predictor','Pr(>|t|)'])),'fdr')

lapply(split(sample.cluster.proportions,sample.cluster.proportions$cluster_type),function(x) coef(summary(lm(prop~predictor+group,data=x)))['predictor',])

saveRDS(sample.cell.proportions,file=paste0('checkpoints/',dataset,'_cell_proportions.rds'))

saveRDS(sample.cluster.proportions,file=paste0('checkpoints/',dataset,'_cluster_proportions.rds'))

p = ggplot(sample.cell.proportions,aes(predictor,prop,color=type)) +
	geom_point() + geom_smooth(method=lm,se=FALSE) +
	scale_color_manual(name='Cell type',values=cell.prediction.colors[cell.prediction.levels %in% sample.cell.proportions$type]) +
	scale_x_continuous(limits=c(0,ceiling(max(sample.cell.proportions$predictor)))) +
#	scale_y_continuous(limits=c(0,100)) +
	facet_wrap(~type.label,ncol=1,scales='free_y',strip.position='right') +
#	coord_fixed(ratio=0.2) +
	xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) + ylab('Percentage of cells') +
	theme_classic(base_size=14) +
	theme(legend.position='none',strip.background.y=element_blank(),strip.text.y=element_text(hjust=0,angle = 0),axis.text.y=element_text(size=6))
# ggsave(p,file=paste0('figures/',dataset,'_clustering_cell_proportions_',predictor,'.pdf'),useDingbats=FALSE,height=5)

gt = ggplot_gtable(ggplot_build(p))
gt$layout$clip[grepl('panel',gt$layout$name)] = 'off'

pdf(file=paste0('figures/',dataset,'_clustering_cell_proportions_',predictor,'.pdf'),height=5,useDingbats=FALSE)
	grid::grid.draw(gt)
dev.off()

p = ggplot(sample.cell.proportions,aes(predictor,prop,color=type)) +
	geom_point() + geom_smooth(method=lm,se=FALSE) +
	scale_color_manual(name='Cell type',values=cell.prediction.colors[cell.prediction.levels %in% sample.cell.proportions$type]) +
	scale_x_continuous(limits=c(0,ceiling(max(sample.cell.proportions$predictor)))) +
#	scale_y_continuous(limits=c(0,100)) +
	facet_wrap(~type.label,ncol=2,scales='free_y') +
#	coord_fixed(ratio=0.2) +
	xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) + ylab('Percentage of cells') +
	theme_article(base_size=14) +
	theme(legend.position='none',strip.background=element_blank(),axis.text.y=element_text())

gt = ggplot_gtable(ggplot_build(p))
gt$layout$clip[grepl('panel',gt$layout$name)] = 'off'

pdf(file=paste0('figures/',dataset,'_clustering_cell_proportions_',predictor,'_no_se.pdf'),height=5,useDingbats=FALSE)
	grid::grid.draw(gt)
dev.off()

p = ggplot(sample.cell.proportions,aes(predictor,prop,color=type,fill=type)) +
	geom_point() + geom_smooth(method=lm,se=TRUE) +
	scale_color_manual(name='Cell type',values=cell.prediction.colors[cell.prediction.levels %in% sample.cell.proportions$type]) +
	scale_fill_manual(name='Cell type',values=cell.prediction.colors[cell.prediction.levels %in% sample.cell.proportions$type]) +
	scale_x_continuous(limits=c(0,ceiling(max(sample.cell.proportions$predictor)))) +
#	scale_y_continuous(limits=c(0,100)) +
	facet_wrap(~type.label,ncol=2,scales='free_y') +
#	coord_fixed(ratio=0.2) +
	xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) + ylab('Percentage of cells') +
	theme_article(base_size=14) +
	theme(legend.position='none',strip.background=element_blank(),axis.text.y=element_text())

gt = ggplot_gtable(ggplot_build(p))
gt$layout$clip[grepl('panel',gt$layout$name)] = 'off'

pdf(file=paste0('figures/',dataset,'_clustering_cell_proportions_',predictor,'_article.pdf'),height=5,useDingbats=FALSE)
	grid::grid.draw(gt)
dev.off()

library(egg)

p = ggplot(sample.cluster.proportions,aes(predictor,prop,color=type)) +
	geom_point() + geom_smooth(method=lm,se=FALSE) +
	scale_color_manual(name='Cell type',values=cell.prediction.colors[cell.prediction.levels %in% sample.cell.proportions$type]) +
	scale_x_continuous(limits=c(0,ceiling(max(sample.cell.proportions$predictor)/5) * 5)) +
#	scale_y_continuous(limits=c(0,100)) +
	facet_wrap(~type.label,scales='free_y',ncol=4,dir='h') +
#	coord_fixed(ratio=0.2) +
	xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) + ylab('Percentage of cells') +
	theme_article() +
	theme(legend.position='none',axis.text.y=element_text(size=6),strip.text=element_text(size=6),panel.spacing=unit(0,'lines'))
ggsave(p,file=paste0('figures/',dataset,'_clustering_cluster_proportions_',predictor,'_article.pdf'),useDingbats=FALSE,height=5)

if (FALSE) {
# Generate top markers

top.markers = top_markers(cds,group_cells_by='assigned_cell_cluster',genes_to_test_per_group=50,cores=4)
top.markers$cell_group = factor(top.markers$cell_group,levels=levels(colData(cds)$assigned_cell_cluster))
top.markers = top.markers[order(top.markers$cell_group,-top.markers$pseudo_R2),]

saveRDS(top.markers,file=paste0('checkpoints/',dataset,'_clusters_top_markers.rds'))

p = plot_genes_by_group(cds,
                    top.markers %>% filter(fraction_expressing >= 0.10) %>% group_by(cell_group) %>% top_n(3, pseudo_R2) %>% pull(gene_short_name) %>% unique,
                    group_cells_by='assigned_cell_cluster',
                    ordering_type='none',
                    max.size=3) +
	viridis::scale_color_viridis(name = 'Mean expression') +
	scale_size(name = 'Percent expressed',range=c(0,3)) +
	theme(axis.title=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_marker_genes_clusters.pdf'),useDingbats=FALSE,height=7)




cds.en = cds[,cds$assigned_cell_type == 'Excitatory neurons']
top.markers = top_markers(cds.en,group_cells_by='assigned_cell_cluster',genes_to_test_per_group=50,cores=4)
top.markers$cell_group = droplevels(factor(top.markers$cell_group,levels=levels(colData(cds)$assigned_cell_cluster)))
top.markers = top.markers[order(top.markers$cell_group,-top.markers$pseudo_R2),]
cds.en$assigned_cell_cluster = droplevels(cds.en$assigned_cell_cluster)

p = plot_genes_by_group(cds.en,
                    top.markers %>% filter(fraction_expressing >= 0.10) %>% group_by(cell_group) %>% top_n(5, pseudo_R2) %>% pull(gene_short_name) %>% unique,
                    group_cells_by='assigned_cell_cluster',
                    ordering_type='none',
                    max.size=3) +
	viridis::scale_color_viridis(name = 'Mean expression') +
	scale_size(name = 'Percent expressed',range=c(0,3)) +
	theme(axis.title=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_marker_genes_exc_clusters.pdf'),useDingbats=FALSE,height=7)




cds.in = cds[,cds$assigned_cell_type == 'Inhibitory neurons']
top.markers = top_markers(cds.in,group_cells_by='assigned_cell_cluster',genes_to_test_per_group=50,cores=4)
top.markers$cell_group = droplevels(factor(top.markers$cell_group,levels=levels(colData(cds)$assigned_cell_cluster)))
top.markers = top.markers[order(top.markers$cell_group,-top.markers$pseudo_R2),]
cds.in$assigned_cell_cluster = droplevels(cds.in$assigned_cell_cluster)

p = plot_genes_by_group(cds.in,
                    top.markers %>% filter(fraction_expressing >= 0.10) %>% group_by(cell_group) %>% top_n(5, pseudo_R2) %>% pull(gene_short_name) %>% unique,
                    group_cells_by='assigned_cell_cluster',
                    ordering_type='none',
                    max.size=3) +
	viridis::scale_color_viridis(name = 'Mean expression') +
	scale_size(name = 'Percent expressed',range=c(0,3)) +
	theme(axis.title=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_marker_genes_inh_clusters.pdf'),useDingbats=FALSE,height=7)



top.markers = top_markers(cds,group_cells_by='assigned_cell_type',genes_to_test_per_group=50,cores=4)
top.markers$cell_group = factor(top.markers$cell_group,levels=levels(colData(cds)$assigned_cell_type))
top.markers = top.markers[order(top.markers$cell_group,-top.markers$pseudo_R2),]

saveRDS(top.markers,file=paste0('checkpoints/',dataset,'_celltypes_top_markers.rds'))

p = plot_genes_by_group(cds,
                    top.markers %>% filter(fraction_expressing >= 0.10) %>% group_by(cell_group) %>% top_n(5, pseudo_R2) %>% pull(gene_short_name) %>% unique,
                    group_cells_by='assigned_cell_type',
                    ordering_type='none',
                    max.size=3) +
	viridis::scale_color_viridis(name = 'Mean expression') +
	scale_size(name = 'Percent expressed',range=c(0,3)) +
	theme(axis.title=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_marker_genes_cells.pdf'),useDingbats=FALSE,height=7)
}