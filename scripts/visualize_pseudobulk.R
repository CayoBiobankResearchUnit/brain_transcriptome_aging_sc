#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(tidyverse)
library(RColorBrewer)
library(umap)

dataset = arguments[1]

pseudobulk = readRDS(paste0('checkpoints/',dataset,'_pseudobulk.rds'))
keep.genes = readRDS(paste0('checkpoints/',dataset,'_keep_genes.rds'))
meta = readRDS(paste0('checkpoints/',dataset,'_metadata.rds'))

meta = subset(meta,sample %in% Reduce(union,lapply(pseudobulk,colnames)))

pseudobulk.libraries = do.call(cbind,lapply(cell.levels,function(i) {
	x = pseudobulk[[i]]
	colnames(x) = paste(colnames(x),i,sep='.')
	x
}))

meta.libraries = do.call(rbind,lapply(cell.levels,function(i) {
	m = meta
	m$library_id = paste(m$sample,i,sep='.')
	m$cell = factor(i,levels=cell.levels)
	m
}))
rownames(meta.libraries) = meta.libraries$library_id

# Calculate total number of reads (UMI)
meta.libraries$umi = colSums(pseudobulk.libraries)

p = ggplot(meta.libraries,aes_string(predictor,'umi',color='cell')) +
	geom_point(size=1.25) +
	geom_smooth(method=lm,se=FALSE) +
	scale_color_manual(name='Cell type',values=cell.colors) +
	facet_wrap(~cell,scales='free_y',ncol=1) +
	scale_y_continuous() +
	theme_classic() +
	xlab(predictor.label) +
	ylab('Number of UMI') +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_',predictor,'_vs_umi.pdf'),useDingbats=FALSE)


set.seed(1)
a = umap(t(pseudobulk.libraries), n_neighbors = 50, min_dist = 0.5)
a.umap = data.frame(as.data.frame(a$layout),meta.libraries)

p = ggplot(a.umap,aes_string('V1','V2',color='cell')) +
	geom_point(size=1.25) +
	scale_color_manual(name='Cell type',values=cell.colors) +
	theme_classic(base_size=18) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_umap.pdf'),useDingbats=FALSE,height=3)

p = ggplot(a.umap,aes_string('V1','V2',color='age')) +
	geom_point(size=1.25) +
#	scale_color_manual(name='Cell type',values=cell.colors) +
	theme_classic(base_size=18) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_umap_age.pdf'),useDingbats=FALSE,height=3)

b = prcomp(cor(pseudobulk.libraries[,!colnames(pseudobulk.libraries) %in% names(which(colSums(pseudobulk.libraries)==0))]))
b.pca = data.frame(as.data.frame(b$rotation),meta.libraries[!colnames(pseudobulk.libraries) %in% names(which(colSums(pseudobulk.libraries)==0)),])
b.pca.long = data.frame(reshape2::melt(b$rotation),meta.libraries[reshape2::melt(b$rotation)$Var1,])
b.pca.long = do.call(rbind,lapply(split(b.pca.long,b.pca.long$Var2),function(x) {
	pc = unique(x$Var2)
	x$pc.label = paste0(pc,' (',gsub(' ','',format(round(summary(b)$importance['Proportion of Variance',pc]*100,1),nsmall=1)),'%)')
	x
}))

p = ggplot(b.pca,aes_string('PC1','PC2',color='cell')) +
	geom_point(size=1.25) +
	scale_color_manual(name='Cell type',values=cell.colors) +
	theme_classic(base_size=18) +
	xlab(paste0('PC1',' (',round(summary(b)$importance['Proportion of Variance',paste0('PC1')]*100,1),'%)')) +
	ylab(paste0('PC2',' (',round(summary(b)$importance['Proportion of Variance',paste0('PC2')]*100,1),'%)')) +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_pca.pdf'),useDingbats=FALSE,height=3)

p = ggplot(subset(b.pca.long,Var2 %in% paste0('PC',1:length(cell.levels))),aes_string('umi','value',color='cell')) +
	geom_point(size=1.25) +
	scale_color_manual(name='Cell type',values=cell.colors) +
	facet_wrap(~pc.label,nrow=cell.rows) +
	scale_x_continuous(trans='log10',breaks=10^(0:ceiling(log10(max(b.pca.long$umi)))),labels=paste0('1E',log10(10^(0:ceiling(log10(max(b.pca.long$umi))))))) +
	theme_classic(base_size=18) +
	xlab('Number of UMI') +
	ylab('PC') +
	theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_pca_pc_vs_umi.pdf'),useDingbats=FALSE,width=11)

# PCA for each region
pca = do.call(rbind,lapply(cell.levels,function(i) {
	genes = keep.genes[[i]]
	ids = setdiff(rownames(subset(meta.libraries,cell == i)),names(which(colSums(pseudobulk.libraries) == 0)))

	b = prcomp(cor(pseudobulk.libraries[genes,ids]))
	b.pca = do.call(rbind,lapply(1:length(keep.genes),function(x) {
		data.frame(
			meta.libraries[ids,],
			pc = x,
			pc.label = paste0('PC',formatC(x,width=2,flag=0)),
			loading = b$rotation[,paste0('PC',x)],
			region = i,
			prop.variance = summary(b)$importance['Proportion of Variance',paste0('PC',x)],
			cum.variance = summary(b)$importance['Cumulative Proportion',paste0('PC',x)],
			label.cell = paste0(
				i, ' (',
				format(round(summary(b)$importance['Proportion of Variance',paste0('PC',x)] * 100,2),nsmall=2),
				'%)'
			),
			label.pc = paste0(
				paste0('PC',formatC(x,width=2,flag=0)), ' (',
				format(round(summary(b)$importance['Proportion of Variance',paste0('PC',x)] * 100,2),nsmall=2),
				'%)'
			),
			stringsAsFactors=FALSE
		)
	}))
	b.pca
}))

for (i in 1:length(cell.levels)) {
	p = ggplot(subset(pca,pc == i),aes_string(predictor,'loading')) +
		geom_point() +
		geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~label.cell,scales='free_y',nrow=cell.rows) +
		theme_classic() +
		theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) + 
		xlab(predictor.label) +
		ylab(paste0('PC',i))
	ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_',tolower(predictor.label),'_vs_pc',formatC(i,width=2,flag=0),'.pdf'),useDingbats=FALSE,height=7,width=11)
}
for (r in cell.levels) {
	p = ggplot(subset(pca,cell == r & pc <= 16),aes_string(predictor,'loading')) +
		geom_point() +
		geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~label.pc,scales='free_y',nrow=cell.rows) +
		theme_classic() +
		theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) + 
		xlab(predictor.label) +
		ylab(paste0('PC'))
	ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_',tolower(predictor.label),'_vs_pcs_cell_',gsub(' ','_',tolower(r)),'.pdf'),useDingbats=FALSE,height=7,width=11)
}

library(GMPR)

keep.genes = readRDS(paste0('checkpoints/',dataset,'_keep_genes.rds'))


pseudobulk.gmpr = lapply(cell.levels,function(x) {
	m = pseudobulk[[x]]
	size.factor = suppressWarnings(GMPR(t(m),min_ct=1,intersect_no=1))
	g = t(t(m) / size.factor)
	g
#	v = voom(pseudobulk[[x]])
#	v$E[keep.genes[[x]],]
})

e.keep = do.call(cbind,pseudobulk.gmpr)[Reduce(union,keep.genes),]

colnames(e.keep) = with(expand.grid(meta$sample,cell.levels),gsub(' ','.',paste(Var1,Var2,sep='_')))

remove.terms = function(model,covariates) {
	remove.term = function(model,term) {
		# Find term in model matrix/output
		search.string = gsub('\\.','\\\\.',paste0('^',term))
	
		# Number of classes of that term
		n.classes = sum(grepl(search.string,names(coef(model))))

		# Extract betas
		betas = coef(model)[grep(search.string,names(coef(model)))]
	
		# Extract model terms
		mat = model.matrix(model)[,grep(search.string,colnames(model.matrix(model)))]
		if (n.classes > 1) {
			names(betas) = colnames(mat) = gsub(search.string,'',names(betas))
			regress.term = rowSums(t(betas * t(mat)))
		} else {
			regress.term = betas * mat
		}
		as.numeric(replace(regress.term,which(is.na(regress.term)),0))
	}
	Reduce(`+`,lapply(covariates,function(x) {
		remove.term(model,x)
	}))
}

e.meta = data.frame(id=colnames(e.keep))
e.meta$sample_id = factor(substr(e.meta$id,1,5),levels=meta$sample)
e.meta$cell_type = factor(gsub('\\.',' ',gsub('NSM[0-9]{2}_','',e.meta$id)),levels=cell.levels)

e.meta = merge(e.meta,meta,by.x='sample_id',by.y='sample')
rownames(e.meta) = e.meta$id

clus = makeCluster(n.cores)
registerDoParallel(cores=n.cores)
clusterExport(clus,varlist=c('e.keep','e.meta','remove.terms','model.covariates','batch.variable'),envir=environment())

e.regressed = t(parApply(clus,e.keep,1,function(y) {
	
	# "this" is a data frame with the metadata and the expression for a single gene
	this = e.meta
	this$e = y


	model = lm(as.formula(paste('e',paste(model.covariates,collapse=' + '),sep=' ~ ')),data=this)

	partial.resid = y - remove.terms(model,batch.variable)

	return(partial.resid)
	# return(resid(model)) # Return residuals of model NOT including biological terms of interest
}))

stopCluster(clus)

e.regressed = e.regressed[,apply(e.regressed,2,function(x) !any(is.na(x)))]

e.regressed = e.regressed[,rownames(subset(e.meta,cell_type != 'Pericytes'))]
#e.meta = e.meta[colnames(e.regressed),]

library(tidyverse)
library(ggplot2)
library(umap)
library(RColorBrewer)

set.seed(3)

a = umap(t(e.regressed), n_neighbors = 50, min_dist = 0.5)
a.umap = data.frame(as.data.frame(a$layout),e.meta[rownames(a$layout),c(model.covariates,'animal_id','cell_type')])

p = ggplot(a.umap,aes_string('V1','V2',color='cell_type')) +
	geom_point(size=1.5) +
	scale_color_manual(values=cell.colors,name='Cell type') +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_gmpr_umap.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(a.umap,aes_string('V1','V2',color='group')) +
	geom_point(size=1.5) +
	scale_color_brewer(palette='Set1',name='Group') +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_gmpr_umap_group.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(a.umap,aes_string('V1','V2',color='age')) +
	geom_point(size=1.5) +
	scale_color_viridis(option='C',name='Age') +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_gmpr_umap_age.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(a.umap,aes_string('V1','V2',color='animal_id')) +
	geom_point(size=1.5) +
	scale_color_discrete(name='ID') +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_gmpr_umap_id.pdf'),useDingbats=FALSE,width=12,height=7)


b = prcomp(cor(e.regressed))
b.pca = data.frame(as.data.frame(b$rotation),e.meta[rownames(a$layout),c(model.covariates,'animal_id','cell_type')])
# b.pca.long = data.frame(reshape2::melt(b$rotation),e.meta[reshape2::melt(b$rotation)$Var1,])
# b.pca.long = do.call(rbind,lapply(split(b.pca.long,b.pca.long$Var2),function(x) {
# 	pc = unique(x$Var2)
# 	x$pc.label = paste0(pc,' (',gsub(' ','',format(round(summary(b)$importance['Proportion of Variance',pc]*100,1),nsmall=1)),'%)')
# 	x
# }))


p = ggplot(b.pca,aes_string('PC1','PC2',color='cell_type')) +
	geom_point(size=1.5) +
	scale_color_manual(values=cell.colors,name='Cell type') +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',round(summary(b)$importance['Proportion of Variance',1] * 100,2),'%)')) +
	ylab(paste0('PC2 (',round(summary(b)$importance['Proportion of Variance',2] * 100,2),'%)')) +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_gmpr_pca.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(b.pca,aes_string('PC1','PC2',color='group')) +
	geom_point(size=1.5) +
	scale_color_brewer(palette='Set1',name='Group') +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',round(summary(b)$importance['Proportion of Variance',1] * 100,2),'%)')) +
	ylab(paste0('PC2 (',round(summary(b)$importance['Proportion of Variance',2] * 100,2),'%)')) +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_gmpr_pca_group.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(b.pca,aes_string('PC1','PC2',color='age')) +
	geom_point(size=1.5) +
	scale_color_viridis(option='C',name='Age') +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',round(summary(b)$importance['Proportion of Variance',1] * 100,2),'%)')) +
	ylab(paste0('PC2 (',round(summary(b)$importance['Proportion of Variance',2] * 100,2),'%)')) +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_gmpr_pca_age.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(b.pca,aes_string('PC1','PC2',color='animal_id')) +
	geom_point(size=1.5) +
	scale_color_discrete(name='ID') +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',round(summary(b)$importance['Proportion of Variance',1] * 100,2),'%)')) +
	ylab(paste0('PC2 (',round(summary(b)$importance['Proportion of Variance',2] * 100,2),'%)')) +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_pseudobulk_gmpr_pca_id.pdf'),useDingbats=FALSE,width=12,height=7)



keep.genes = readRDS(paste0('checkpoints/',dataset,'_keep_genes.rds'))

for (j in c('gmpr','voom')) {
	message(j)
	pseudobulk.keep = readRDS(paste0('checkpoints/',dataset,'_pseudobulk_normalized_',j,'.rds'))

	pseudobulk.ordered = lapply(names(pseudobulk.keep),function(i) {
		x = pseudobulk.keep[[i]]
		x[keep.genes[[i]],meta$sample[order(meta$age)]]
	})
	names(pseudobulk.ordered) = names(pseudobulk)

	for (i in names(pseudobulk.keep)) {
		message(i)
		pdf(file=paste0('figures/',dataset,'_pseudobulk_heatmap_',j,'_',tolower(gsub(' ','_',i)),'.pdf'))
			heatmap(pseudobulk.ordered[[i]],Colv=NA,labRow=FALSE)
		dev.off()
	}
}