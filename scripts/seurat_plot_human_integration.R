#!/usr/bin/env Rscript

options(future.globals.maxSize = 100000 * 1024^2)

arguments = commandArgs(trailing=TRUE)

dataset = arguments[1]

library(Seurat)
library(Matrix)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

cds = readRDS(paste0('checkpoints/',dataset,'_macmul_homsap_integrated.rds'))

cds = RunPCA(object = cds, verbose = FALSE)

cds = RunUMAP(object = cds, dims = 1:30)
saveRDS(cds,file=paste0('checkpoints/',dataset,'_macmul_homsap_integrated_umap.rds'))

# Grab UMAP
a.umap = do.call(rbind,lapply(rownames(cds@meta.data),function(i) cds@reductions$umap[[i]]))

a.umap = data.frame(as.data.frame(a.umap),cds@meta.data)

a.umap$dataset = 'human'
a.umap$dataset[!is.na(a.umap$assigned_cell_type)] = 'rhesus'

saveRDS(a.umap,file=paste0('checkpoints/',dataset,'_macmul_homsap_integrated_umap_only.rds'))

library(ggplot2)
library(ggrastr)
library(egg)
library(RColorBrewer)
library(viridis)

a.umap$assigned_cell_type = factor(a.umap$assigned_cell_type,levels=c(cell.levels,'Unknown'))
a.umap$class_label[!is.na(a.umap$class_label) & !nchar(a.umap$class_label)] = NA

a.umap.test = a.umap[sample(1:nrow(a.umap),5000),]

a.umap$class2 = a.umap$class_label
a.umap$class2[a.umap$class_label %in% 'Non-Neuronal'] = a.umap$subclass_label[a.umap$class_label %in% 'Non-Neuronal']

a.umap$class2 = factor(a.umap$class2,levels=c(with(subset(a.umap,!is.na(class2)),c(sort(intersect(unique(class_label),unique(class2))),sort(setdiff(unique(class2),unique(class_label)))))))

p = ggplot() +
	geom_point_rast(aes(UMAP_1,UMAP_2),data=subset(a.umap,dataset=='human'),color='#cccccc',alpha=0.05,size=0.01,show.legend=FALSE) +
	geom_point_rast(aes(UMAP_1,UMAP_2,color=assigned_cell_type),data=subset(a.umap,dataset=='rhesus'),alpha=0.05,size=0.01) +
	coord_equal() +
	scale_color_manual(name='Cell type',values=cell.colors,drop=TRUE) +
	theme_classic(base_size=16) +
	theme(axis.ticks=element_blank(),axis.text=element_blank()) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	xlab('UMAP 1') + ylab('UMAP 2')
ggsave(p,file=paste0('figures/human_macaque_clustering_assigned_cell_type.pdf'),useDingbats=FALSE,height=5,width=10)


p = ggplot() +
	geom_point_rast(aes(UMAP_1,UMAP_2,color=class_label),data=subset(a.umap,dataset=='human' & !is.na(class_label)),alpha=0.05,size=0.01,show.legend=TRUE) +
	geom_point_rast(aes(UMAP_1,UMAP_2),data=subset(a.umap,dataset=='rhesus'),color='#cccccc',alpha=0.05,size=0.01) +
	coord_equal() +
#	scale_color_manual(name='Cell type',values=cell.colors,drop=TRUE) +
	scale_color_brewer(name='Class',palette='Set1') +
	theme_classic(base_size=16) +
	theme(axis.ticks=element_blank(),axis.text=element_blank()) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	xlab('UMAP 1') + ylab('UMAP 2')
ggsave(p,file=paste0('figures/human_macaque_clustering_cell_class.pdf'),useDingbats=FALSE,height=5,width=10)


p = ggplot() +
	geom_point_rast(aes(UMAP_1,UMAP_2,color=class2),data=subset(a.umap,dataset=='human' & !is.na(class2)),alpha=0.05,size=0.01,show.legend=TRUE) +
	geom_point_rast(aes(UMAP_1,UMAP_2),data=subset(a.umap,dataset=='rhesus'),color='#cccccc',alpha=0.05,size=0.01) +
	coord_equal() +
#	scale_color_manual(name='Cell type',values=cell.colors,drop=TRUE) +
	scale_color_brewer(name='Class',palette='Accent') +
	theme_classic(base_size=16) +
	theme(axis.ticks=element_blank(),axis.text=element_blank()) +
	guides(color = guide_legend(override.aes = list(size=2,alpha=1))) +
	xlab('UMAP 1') + ylab('UMAP 2')
ggsave(p,file=paste0('figures/human_macaque_clustering_cell_subclass.pdf'),useDingbats=FALSE,height=5,width=10)

