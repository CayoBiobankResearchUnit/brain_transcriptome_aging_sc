#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(tidyverse)
library(mashr)
library(reshape2)
library(egg)
library(ggrastr)

dataset = arguments[1]

if (length(arguments) > 1) {
	model.method = arguments[2]
} else {
	model.method = 'emma'
}

model.output = readRDS(paste0('checkpoints/',dataset,'_model_',model.method,'.rds'))

model.stats = bind_model_outputs(model.output,predictor)

model.stats[[paste('pval',predictor,sep='.')]][model.stats$converged %in% FALSE] = NA

p = ggplot(model.stats,aes_string(paste('pval',predictor,sep='.'),fill='cell')) +
	geom_histogram(binwidth=0.005) +
	facet_wrap(~cell,nrow=cell.rows) +
	scale_fill_manual(values=cell.colors) +
	scale_x_continuous(
		limits = c(0,1),
		breaks = seq(0,1,0.1),
		labels = c('0.0',rep('',4),'0.5',rep('',4),'1.0')
	) +
	xlab('p value') + ylab('Count') +
	theme_classic() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_histogram.pdf'),useDingbats=FALSE,height=5)

# QQ-plot the stat_qq way
# nlog_trans = function (base = exp(1)) {
#     force(base)
#     trans = function(x) -log(x, base)
#     inv = function(x) base^(-x)
#     scales::trans_new(paste0("log-", format(base)), trans, inv, scales::log_breaks(base = base), 
#         domain = c(1e-100, Inf))
# }
# nlog10_trans = nlog_trans(10)
# 
# p = ggplot(
# 		model.stats,
# 		aes_string(sample=paste('pval',predictor,sep='.'),color='cell')
# 	) +
# 	stat_qq(distribution=stats::qunif) +
# 	geom_abline(slope=1,col='black',size=0.2) +
# 	facet_wrap(~cell) +
# 	scale_x_continuous(trans=nlog10_trans) +
# 	scale_y_continuous(trans=nlog10_trans) +
# 	scale_color_manual(values=cell.colors) +
# # 	scale_x_continuous(breaks=seq(0,1,0.05),labels=c(0,rep('',19),1)) +
# # 	scale_y_continuous(breaks=seq(0,1,0.05),labels=c(0,rep('',19),1)) +
# 	xlab('Expected') + ylab('Observed') + theme_classic() +
# 	theme(legend.position='none')
# ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_qqplot.pdf'),useDingbats=FALSE)

# Custom code for QQ plot
p = ggplot(
		subset(do.call(rbind,lapply(split(model.stats,model.stats$cell),function(x) {
			data.frame(
				cell=factor(unique(x$cell),levels=cell.levels,labels=cell.short.levels),
				expected=-log10(seq(1/nrow(x),1,1/nrow(x))),
				observed=-log10(quantile(x[[paste('pval',predictor,sep='.')]],seq(1/nrow(x),1,1/nrow(x)),na.rm=TRUE))
			)
		})),!is.na(observed)),
		aes(expected,observed,color=cell)
	) +
#	geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
	geom_point_rast(size=0.5) +
	geom_abline(slope=1,col='black',size=0.2) +
#	facet_wrap(~cell,nrow=cell.rows,strip.position='left') +
	facet_wrap(~cell,nrow=1,strip.position='top',drop=TRUE) +
	coord_fixed() +
#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
	scale_color_manual(values=cell.colors) +
 	scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(model.stats$cell))-1)))),1)) +
 	scale_y_continuous(breaks=seq(0,ceiling(max(-log10(model.stats[[paste('pval',predictor,sep='.')]]),na.rm=TRUE)),1)) +
# 	scale_y_continuous(breaks=seq(0,1,0.05),labels=c(0,rep('',19),1)) +
	xlab(expression(-log[10] * ('Expected'~italic(p)))) +
	ylab(expression(-log[10] * ('Observed'~italic(p)))) +
	theme_article(base_size=16) +
	theme(legend.position='none')
#	theme(legend.position='none',strip.text=element_text(angle = 0))
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_qqplot.pdf'),useDingbats=FALSE,height=5)

# If reference bulk dataset found
if (file.exists(paste0('checkpoints/',dataset,'_bulk_reference_',model.method,'.rds'))) {
	# Compare super-pseudobulk data to reference bulk dataset

	superbulk = readRDS(paste0('checkpoints/',dataset,'_model_',model.method,'_superbulk.rds'))
	referencebulk = readRDS(paste0('checkpoints/',dataset,'_bulk_reference_',model.method,'.rds'))

	superbulk = superbulk[,colnames(referencebulk)]
	superbulk = cbind(superbulk,p.adjust(superbulk[,paste('pval',predictor,sep='.')]))
	referencebulk = cbind(referencebulk,p.adjust(referencebulk[,paste('pval',predictor,sep='.')]))
	colnames(superbulk)[ncol(superbulk)] = colnames(referencebulk)[ncol(referencebulk)] = paste('qval',predictor,sep='.')


	colnames(superbulk) = paste('super',colnames(superbulk),sep='.')
	colnames(referencebulk) = paste('reference',colnames(referencebulk),sep='.')
	compare.bulk = as.data.frame(cbind(
		superbulk[intersect(rownames(superbulk),rownames(referencebulk)),],
		referencebulk[intersect(rownames(superbulk),rownames(referencebulk)),]
	))
	compare.bulk$significant = factor((compare.bulk[[paste('super.pval',predictor,sep='.')]] < 0.05) %in% TRUE & (compare.bulk[[paste('reference.pval',predictor,sep='.')]] < 0.05) %in% TRUE,levels=c('TRUE','FALSE'))
	levels(compare.bulk$significant) = c('p < 0.05','p ≥ 0.05')

	compare.bulk = compare.bulk[complete.cases(compare.bulk),]

	p = ggplot() +
		geom_hline(yintercept=0) +
		geom_vline(xintercept=0) +
		geom_smooth(data=compare.bulk,aes_string(paste('reference','beta',predictor,sep='.'),paste('super','beta',predictor,sep='.'),color='significant',alpha='significant'),method=lm,show.legend=FALSE) +
		geom_point(data=compare.bulk,aes_string(paste('reference','beta',predictor,sep='.'),paste('super','beta',predictor,sep='.'),color='significant',alpha='significant'),size=0.25) +
#		geom_text(
#			data = data.frame(label = lm_eqn(data.frame(x=compare.bulk[[paste('reference','beta',predictor,sep='.')]], y=compare.bulk[[paste('super','beta',predictor,sep='.')]])[(compare.bulk[[paste('reference.pval',predictor,sep='.')]] < 0.05) %in% TRUE & (compare.bulk[[paste('super.pval',predictor,sep='.')]] < 0.05) %in% TRUE,])),
#			aes(
#				x=-1.5*max(abs(c(min(compare.bulk[[paste('reference','beta',predictor,sep='.')]],na.rm=TRUE),max(compare.bulk[[paste('reference','beta',predictor,sep='.')]],na.rm=TRUE)))) + max(abs(c(min(compare.bulk[[paste('reference','beta',predictor,sep='.')]],na.rm=TRUE),max(compare.bulk[[paste('reference','beta',predictor,sep='.')]],na.rm=TRUE)))) * 0.5,
#				y=min(compare.bulk[[paste('super','beta',predictor,sep='.')]],na.rm=TRUE) + (max(compare.bulk[[paste('super','beta',predictor,sep='.')]],na.rm=TRUE) - min(compare.bulk[[paste('super','beta',predictor,sep='.')]],na.rm=TRUE)) * 0.95,
#				label=label
#			),
#			size=4,
#			parse=TRUE
#		) +
		scale_color_manual(name='Significance',values=c('#ff0000','#cccccc'),labels=c(expression(italic(p) < 0.05),expression(italic(p) >= 0.05)),na.value=NA) +
		scale_alpha_manual(values=c(1,0.1),na.value=NA) +
#		scale_x_continuous(limits=1.5*c(-max(abs(c(min(compare.bulk[[paste('reference','beta',predictor,sep='.')]],na.rm=TRUE),max(compare.bulk[[paste('reference','beta',predictor,sep='.')]],na.rm=TRUE)))),max(abs(c(min(compare.bulk[[paste('reference','beta',predictor,sep='.')]],na.rm=TRUE),max(compare.bulk[[paste('reference','beta',predictor,sep='.')]],na.rm=TRUE)))))) +
#		scale_y_continuous(limits=c(-max(abs(c(min(compare.bulk[[paste('super','beta',predictor,sep='.')]],na.rm=TRUE),max(compare.bulk[[paste('super','beta',predictor,sep='.')]],na.rm=TRUE)))),max(abs(c(min(compare.bulk[[paste('super','beta',predictor,sep='.')]],na.rm=TRUE),max(compare.bulk[[paste('super','beta',predictor,sep='.')]],na.rm=TRUE)))))) +
		coord_fixed() +
		xlab(expression(beta ~ (true~bulk))) +
		ylab(expression(beta ~ (pseudobulk))) +
		ggtitle(parse(text=lm_eqn(data.frame(x=compare.bulk[[paste('reference','beta',predictor,sep='.')]], y=compare.bulk[[paste('super','beta',predictor,sep='.')]])[(compare.bulk[[paste('reference.pval',predictor,sep='.')]] < 0.05) %in% TRUE & (compare.bulk[[paste('super.pval',predictor,sep='.')]] < 0.05) %in% TRUE,]))) +
		theme_classic(base_size=18) +
		theme(axis.ticks=element_blank(),axis.text=element_blank(),plot.title=element_text(size=14)) +
		guides(color = guide_legend(override.aes=list(size=3,alpha=1)),alpha=FALSE)
	ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_pseudobulk_vs_referencebulk.pdf'),width=7,height=5,useDingbats=FALSE)

	p = ggplot(compare.bulk,aes_string(paste('reference','beta',predictor,sep='.'),paste('super','beta',predictor,sep='.'),color='significant',alpha='significant')) +
		geom_hline(yintercept=0) +
		geom_vline(xintercept=0) +
		geom_smooth(method=lm,show.legend=FALSE,alpha=0,fill=NA,color=NA) +
		geom_point(size=0.05,alpha=0) +
		scale_color_manual(name='Significance',values=c('#ff0000','#cccccc'),labels=c(expression(italic(p) < 0.05),expression(italic(p) >= 0.05)),na.value=NA) +
		scale_alpha_manual(values=c(1,0.1),na.value=NA) +
		coord_fixed() +
		xlab(expression(beta ~ (true~bulk))) +
		ylab(expression(beta ~ (pseudobulk))) +
		ggtitle(parse(text=lm_eqn(data.frame(x=compare.bulk[[paste('reference','beta',predictor,sep='.')]], y=compare.bulk[[paste('super','beta',predictor,sep='.')]])[(compare.bulk[[paste('reference.pval',predictor,sep='.')]] < 0.05) %in% TRUE & (compare.bulk[[paste('super.pval',predictor,sep='.')]] < 0.05) %in% TRUE,]))) +
		theme_classic(base_size=18) +
		theme(axis.ticks=element_blank(),axis.text=element_blank(),plot.title=element_text(size=14,color=NA)) +
		guides(color = guide_legend(override.aes=list(size=3,alpha=1)),alpha=FALSE)
	ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_pseudobulk_vs_referencebulk_blank.pdf'),width=7,height=5,useDingbats=FALSE)

}

# Use bulk age predictors on superbulk expression matrix to predict age
bulk.best.predictors = readRDS('checkpoints/bulk_best_predictors.rds')
pseudobulk = readRDS(paste0('checkpoints/',dataset,'_pseudobulk.rds'))

meta = read.delim(paste0('data/metadata_',dataset,'.txt'),stringsAsFactors=FALSE)
meta = subset(meta,sample %in% Reduce(union,lapply(pseudobulk,colnames)))

e.superbulk = Reduce(`+`,pseudobulk)

e.superbulk.keep = e.superbulk[rowMeans(e.superbulk/colSums(e.superbulk) * 1e6) > tpm.cutoff,]

e.superbulk.scaled = t(apply(e.superbulk.keep,1,scale))
colnames(e.superbulk.scaled) = colnames(e.superbulk.keep)

superbulk.gene.predictions = unlist(lapply(colnames(e.superbulk.scaled),function(x) sum(e.superbulk.scaled[intersect(names(bulk.best.predictors),rownames(e.superbulk.scaled)),x] * bulk.best.predictors[intersect(names(bulk.best.predictors),rownames(e.superbulk.scaled))])))
names(superbulk.gene.predictions) = colnames(e.superbulk.scaled)
superbulk.standardized.predictions = mean(meta[[predictor]]) + (superbulk.gene.predictions - mean(superbulk.gene.predictions)) * sd(meta[[predictor]]) / sd(superbulk.gene.predictions)

superbulk.predictions = data.frame(meta,known=meta[[predictor]],predicted=superbulk.standardized.predictions[meta$sample])

p = ggplot(superbulk.predictions,aes(known,predicted)) +
	geom_abline(slope=1,intercept=0,col='black',size=0.5) +
	geom_point() +
	geom_text(data = data.frame(label = lm_eqn(with(superbulk.predictions,data.frame(x=known, y=predicted)))), aes(x=(ceiling(with(superbulk.predictions,max(c(known,predicted))) / 5) * 5) / 4, y=(ceiling(with(superbulk.predictions,max(c(known,predicted))) / 5) * 5) * 0.95, label=label), size=6, parse=TRUE) +
	geom_smooth(method=lm,se=TRUE) +
	coord_fixed(xlim=c(0,ceiling(with(superbulk.predictions,max(c(known,predicted))) / 5) * 5),ylim=c(0,ceiling(with(superbulk.predictions,max(c(known,predicted))) / 5) * 5)) +
	scale_x_continuous(limits=c(0,ceiling(with(superbulk.predictions,max(c(known,predicted))) / 5) * 5)) +
	scale_y_continuous(limits=c(floor(with(superbulk.predictions,min(c(known,predicted))) / 5) * 5,ceiling(with(superbulk.predictions,max(c(known,predicted))) / 5) * 5)) +
	theme_classic(base_size=24) + xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_bulk_pseudobulk_predictions.pdf'),useDingbats=FALSE)

# Visualize mashr

mash.results = readRDS(paste0('checkpoints/',dataset,'_model_',model.method,'_mashr.rds'))

mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)

model.stats[[paste('qval',predictor,sep='.')]] = p.adjust(model.stats[[paste('pval',predictor,sep='.')]],'fdr')

p = ggplot(model.stats[model.stats[[paste('qval',predictor,sep='.')]] < fsr.cutoff,],aes_string(paste('beta',predictor,sep='.'),color='cell')) +
	scale_color_manual(name='Cell type',values=cell.colors) +
	geom_density() +
	theme_classic(base_size=12) +
	guides(color = guide_legend(ncol = 1)) +
	xlab(expression(italic(beta))) +
	ylab('Density')
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_beta_density.pdf'),width=7,height=3,useDingbats=FALSE)

p = ggplot(melt(tapply(model.stats[[paste('qval',predictor,sep='.')]],model.stats$cell,function(x) sum(x < fsr.cutoff,na.rm=TRUE))),aes(Var1,value,fill=Var1)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Cell type',values=cell.colors) +
	theme_classic(base_size=12) +
	theme(legend.position='none',axis.text.x=element_text(angle=-45,hjust=0,vjust=1)) +
	xlab('Cell types') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_beta_counts.pdf'),width=7,height=3,useDingbats=FALSE)

b.mash = data.frame(
	expand.grid(rownames(mash.beta),colnames(mash.beta)),
	qval=as.numeric(mash.lfsr),
	beta=as.numeric(mash.beta),
	stringsAsFactors=FALSE
)

b.mash$Var2 = factor(b.mash$Var2,levels=cell.levels)

p = ggplot(subset(b.mash,qval < fsr.cutoff),aes(beta,color=Var2)) +
	scale_color_manual(name='Cell type',values=cell.colors[cell.levels %in% b.mash$Var2]) +
	geom_density() +
	theme_classic(base_size=12) +
	guides(color = guide_legend(ncol = 1)) +
	xlab(expression(italic(beta))) +
	ylab('Density')
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_mashr_beta_density.pdf'),width=7,height=3,useDingbats=FALSE)

p = ggplot(melt(tapply(b.mash$qval,b.mash$Var2,function(x) sum(x < fsr.cutoff,na.rm=TRUE))),aes(Var1,value,fill=Var1)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Cell type',values=cell.colors[cell.levels %in% b.mash$Var2]) +
	theme_classic(base_size=12) +
	theme(legend.position='none',axis.text.x=element_text(angle=-45,hjust=0,vjust=1)) +
	xlab('Cell types') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_mashr_beta_counts.pdf'),width=7,height=3,useDingbats=FALSE)

b.glmm = model.stats[c('gene','cell',paste('qval',predictor,sep='.'),paste('beta',predictor,sep='.'))]
names(b.glmm) = c('Var1','Var2','qval','beta')


b.glmm$qval.signed = with(b.glmm,ifelse(beta>0,1,-1) * qval)
b.mash$qval.signed = with(b.mash,ifelse(beta>0,1,-1) * qval)

# Put together longform data frame with gene counts (up vs. downregulated, split by method)
model.counts.combined = rbind(
	within(melt(tapply(b.glmm$qval.signed,b.glmm$Var2,function(x) -sum(abs(x) < fsr.cutoff & x < 0,na.rm=TRUE))), {method=as.character(model.method.labels[model.method]); direction='down'} ),
	within(melt(tapply(b.glmm$qval.signed,b.glmm$Var2,function(x) sum(abs(x) < fsr.cutoff & x >= 0,na.rm=TRUE))), {method=as.character(model.method.labels[model.method]); direction='up'} ),
	within(melt(tapply(b.mash$qval.signed,b.mash$Var2,function(x) -sum(abs(x) < fsr.cutoff & x < 0,na.rm=TRUE))), {method='MASH'; direction='down'} ),
	within(melt(tapply(b.mash$qval.signed,b.mash$Var2,function(x) sum(abs(x) < fsr.cutoff & x >= 0,na.rm=TRUE))), {method='MASH'; direction='up'} )
)

ylimit = ceiling(with(model.counts.combined,max(abs(value),na.rm=TRUE))/100) * 100

if (packageVersion('ggplot2') < 3.3) warning('Some plots may not display correctly with ggplot2 version < 3.3.0')

p = ggplot(model.counts.combined,aes(Var1,value,fill=Var1,alpha=direction)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Cell type',values=cell.colors) +
	scale_alpha_manual(values=c(0.75,1)) +
	scale_y_continuous(
		limits = c(-ylimit,ylimit),
		breaks = c(-ylimit,-ylimit*0.5,0,ylimit*0.5,ylimit),
		labels = c(formatC(ylimit,width=5,flag=' '),'Decrease',formatC(0,width=5,flag=' '),'Increase',formatC(ylimit,width=5,flag=' '))
	) +
	theme_classic(base_size=12) +
	theme(
		axis.ticks.y = element_line(linetype=c(1,0,1,0,1)),
		axis.title.x = element_blank(),
		axis.title.y = element_text(),
		axis.text.x = element_blank(),
		axis.text.y=element_text(
			face = c('plain','bold','plain','bold','plain'),
#			size = axis.text.size * c(1,2,1,2,1),
			angle = c(0,90,0,90,0), hjust=0.5
		)
	) +
	facet_wrap(~method,nrow=2) +
	guides(fill = guide_legend(nrow = 3), alpha=FALSE) +
	theme(legend.position='bottom') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_mashr_comparison_beta_counts.pdf'),width=7,height=6,useDingbats=FALSE)


p = ggplot(mutate(subset(model.counts.combined,!is.na(value) & value != 0),Var1=factor(Var1,levels=cell.levels,labels=cell.short.levels)),aes(Var1,value,fill=Var1,alpha=direction)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Cell type',values=cell.colors) +
	scale_alpha_manual(values=c(0.75,1)) +
	scale_y_continuous(
		limits = c(-ylimit,ylimit),
		breaks = c(-ylimit,-ylimit*0.5,0,ylimit*0.5,ylimit),
		labels = c(formatC(ylimit,width=5,flag=' '),'Decrease',formatC(0,width=5,flag=' '),'Increase',formatC(ylimit,width=5,flag=' '))
	) +
	theme_article(base_size=12) +
	theme(
		axis.ticks.y = element_line(linetype=c(1,0,1,0,1)),
		axis.title.x = element_blank(),
		axis.title.y = element_text(),
		axis.text.x = element_text(
			angle = -45, hjust = 0, vjust = 1
		),
		axis.text.y=element_text(
			face = c('plain','bold','plain','bold','plain'),
#			size = axis.text.size * c(1,2,1,2,1),
			angle = c(0,90,0,90,0), hjust=0.5
		)
	) +
	facet_wrap(~method,nrow=2) +
	guides(fill = guide_legend(nrow = 3), alpha=FALSE) +
	theme(legend.position='none') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_mashr_comparison_beta_counts_presentation.pdf'),width=7,height=5,useDingbats=FALSE)


ylimit = ceiling(with(subset(model.counts.combined,method=='MASH'),max(abs(value),na.rm=TRUE))/100) * 100

p = ggplot(droplevels(subset(model.counts.combined,method=='MASH')),aes(Var1,value,fill=Var1,alpha=direction)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Cell type',values=cell.colors) +
	scale_alpha_manual(values=c(0.75,1)) +
	scale_y_continuous(
		limits = c(-ylimit,ylimit),
		breaks = c(-ylimit,-ylimit*0.5,0,ylimit*0.5,ylimit),
		labels = c(formatC(ylimit,width=5,flag=' '),'Decrease',formatC(0,width=5,flag=' '),'Increase',formatC(ylimit,width=5,flag=' '))
	) +
	theme_classic(base_size=12) +
	theme(
		axis.ticks.y = element_line(linetype=c(1,0,1,0,1)),
		axis.title.x = element_blank(),
		axis.title.y = element_text(),
		axis.text.x = element_text(angle=-45,hjust=0,vjust=1),
		axis.text.y=element_text(
			face = c('plain','bold','plain','bold','plain'),
#			size = axis.text.size * c(1,2,1,2,1),
			angle = c(0,90,0,90,0), hjust=0.5
		),
		legend.position='none'
	) +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_mashr_mash_beta_counts.pdf'),width=7,height=6,useDingbats=FALSE)

library(tidyverse)

p = ggplot(droplevels(subset(within(model.counts.combined,{levels(Var1)=cell.short.levels}),method=='MASH')),aes(Var1,value,fill=Var1,alpha=direction)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Cell type',values=cell.colors) +
	scale_alpha_manual(values=c(0.75,1)) +
	scale_y_continuous(
		limits = c(-ylimit,ylimit),
		breaks = c(-ylimit,-ylimit*0.5,0,ylimit*0.5,ylimit),
		labels = c(formatC(ylimit,width=5,flag=' '),'Decrease',formatC(0,width=5,flag=' '),'Increase',formatC(ylimit,width=5,flag=' '))
	) +
	theme_classic(base_size=24) +
	theme(
		axis.ticks.y = element_line(linetype=c(1,0,1,0,1)),
		axis.title.x = element_blank(),
		axis.title.y = element_text(),
		axis.text.x = element_text(angle=-45,hjust=0,vjust=1),
		axis.text.y=element_text(
			face = c('plain','bold','plain','bold','plain'),
#			size = axis.text.size * c(1,2,1,2,1),
			angle = c(0,90,0,90,0), hjust=0.5
		),
		legend.position='none'
	) +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_mashr_mash_beta_counts_short.pdf'),width=7,height=6,useDingbats=FALSE)


b.mash.split.genes = split(b.mash,b.mash$Var1)

# Tally up sig cell-combos, as well as up-regulated and down-regulated
cell.combinations = table(unlist(lapply(b.mash.split.genes,function(x) {
	paste(subset(x,qval < fsr.cutoff & beta != 0)$Var2,collapse='-')
})))

# Get rid of blanks
cell.combinations = cell.combinations[as.logical(nchar(names(cell.combinations)))]

# Sort in same order (by total)
cell.combinations = cell.combinations[order(cell.combinations,decreasing=TRUE)]

cell.combinations.inc = cell.combinations.dec = integer(length(cell.combinations))
names(cell.combinations.inc) = names(cell.combinations.dec) = names(cell.combinations)

# Ensure all have the same order
cell.combinations.inc[names(cell.combinations)] = table(unlist(lapply(b.mash.split.genes,function(x) {
	paste(subset(x,qval < fsr.cutoff & beta > 0)$Var2,collapse='-')
})))[names(cell.combinations)]

cell.combinations.dec[names(cell.combinations)] = table(unlist(lapply(b.mash.split.genes,function(x) {
	paste(subset(x,qval < fsr.cutoff & beta < 0)$Var2,collapse='-')
})))[names(cell.combinations)]

# Set NAs to 0
cell.combinations.inc[is.na(cell.combinations.inc)] = 0
cell.combinations.dec[is.na(cell.combinations.dec)] = 0

# The sum of up-regulated genes and down-regulated genes do not add up to the sum of significant genes
# This is because genes can have different directions, causing their cells to be tallied differently.
# Thus, resort, based on the sum of up- and down-regulated cell combos
cell.combinations.sum = cell.combinations.inc + cell.combinations.dec

cell.combinations.sum = sort(cell.combinations.sum,decreasing=TRUE)
cell.combinations = cell.combinations[names(cell.combinations.sum)]
cell.combinations.inc = cell.combinations.inc[names(cell.combinations.sum)]
cell.combinations.dec = cell.combinations.dec[names(cell.combinations.sum)]

# Set NAs to 0
cell.combinations.inc[is.na(cell.combinations.inc)] = 0
cell.combinations.dec[is.na(cell.combinations.dec)] = 0

width.of.bars = 0.8
cell.combinations.results = do.call(rbind,lapply(1:length(cell.combinations),function(i) {
	x = names(cell.combinations)[i]
	n = unlist(lapply(strsplit(names(cell.combinations),'-'),length))[i]
	count.all = as.integer(cell.combinations[x])
	count.inc = as.integer(cell.combinations.inc[x])
	count.dec = as.integer(cell.combinations.dec[x])
	out = integer(length(cell.levels))
	names(out) = cell.levels
	out[unlist(strsplit(x,split='-'))] = 1
	out = rbind(
		data.frame(
			combination=i,
			n_cells = n,
			share_cell = n > fraction.shared.cutoff * length(cell.levels),
			cell=factor(cell.levels,levels=cell.levels),
			cell_sig = factor(ifelse(as.logical(out),names(out),NA),levels=cell.levels),
			value = 1,
			xmin = seq(1,length(cell.levels)) - (width.of.bars)/2,
			xmax = seq(1,length(cell.levels)) + (width.of.bars)/2,
			ymin = i - (width.of.bars)/2,
			ymax = i + (width.of.bars)/2,
			chart = 'meta'
		),
		data.frame(
			combination=i,
			n_cells = n,
			share_cell = n > fraction.shared.cutoff * length(cell.levels),
			cell=NA,
			cell_sig=NA,
			value=count.all,
			xmin = NA,
			xmax = NA,
			ymin = NA,
			ymax = NA,
			chart='count_all'
		),
		data.frame(
			combination=i,
			n_cells = n,
			share_cell = n > fraction.shared.cutoff * length(cell.levels),
			cell=NA,
			cell_sig=NA,
			value=count.inc,
			xmin = NA,
			xmax = NA,
			ymin = NA,
			ymax = NA,
			chart='count_increase'
		),
		data.frame(
			combination=i,
			n_cells = n,
			share_cell = n > fraction.shared.cutoff * length(cell.levels),
			cell=NA,
			cell_sig=NA,
			value=count.dec,
			xmin = NA,
			xmax = NA,
			ymin = NA,
			ymax = NA,
			chart='count_decrease'
		)
	)
	out
}))

max.plot = length(cell.levels) * 2
base.color = '#000000'
fade.color = '#000000'
highlight.color = '#ff0000'
base.size = 14

meta.combinations.results.plot = subset(cell.combinations.results,combination < max.plot)

meta.combinations.results.plot$cell = factor(meta.combinations.results.plot$cell,levels=cell.levels,labels=cell.short.levels)
meta.combinations.results.plot$cell_sig = factor(meta.combinations.results.plot$cell_sig,levels=cell.levels,labels=cell.short.levels)

p1 = ggplot(subset(meta.combinations.results.plot,chart=='meta')) +
	geom_blank(aes(x = cell,y=combination)) +
	geom_rect(
		aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=cell,alpha=cell_sig,color=cell_sig),
		linetype=1
	) +
	scale_fill_manual(values=cell.colors) +
	scale_alpha_manual(values=rep(1,length(cell.colors)),na.value=0) +
	scale_color_manual(values=rep('black',length(cell.colors)),na.value=0) +
	scale_y_continuous(trans='reverse',expand=c(0,0)) +
	scale_x_discrete(position='top',expand=c(0,0)) +
	coord_equal() +
	theme_classic(base_size=base.size) + 
	theme(
		legend.position='none',
		axis.line=element_blank(),
		axis.title=element_blank(),
		axis.text.x=element_text(angle=45,vjust=0,hjust=0,margin=margin(t=-1)),
		axis.text.y=element_blank(),
		axis.ticks=element_blank()
	)
p2 = ggplot(subset(meta.combinations.results.plot,chart=='count_increase')) +
	geom_bar(aes(combination,value),stat='identity',fill=base.color,width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0),limits=c(0,ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100)) +
	coord_flip() +
	ylab('Upregulated') +
	theme_classic(base_size=base.size) +
	theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
p3 = ggplot(subset(meta.combinations.results.plot,chart=='count_decrease')) +
	geom_bar(aes(combination,value),stat='identity',fill=base.color,width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(trans='reverse',expand=c(0,0),limits=c(ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100,0)) +
	coord_flip() +
	ylab('Downregulated') +
	theme_classic(base_size=base.size) +
	theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

pdf(file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_mashr_comparison_beta_counts_meta.pdf'),useDingbats=FALSE,height=7,width=11)
	ggarrange(p3,p1,p2,ncol=3,nrow=1,widths=c(1.5,1,1.5),newpage=FALSE)
dev.off()

p4 = ggplot(subset(meta.combinations.results.plot,chart=='count_increase')) +
	geom_bar(aes(combination,value,fill=share_cell),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0),limits=c(0,ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Upregulated') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
p5 = ggplot(subset(meta.combinations.results.plot,chart=='count_decrease')) +
	geom_bar(aes(combination,value,fill=share_cell),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(trans='reverse',expand=c(0,0),limits=c(ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100,0)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Downregulated') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

pdf(file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_mashr_comparison_beta_counts_meta_shared.pdf'),useDingbats=FALSE,height=7,width=11)
	ggarrange(p5,p1,p4,ncol=3,nrow=1,widths=c(1.5,1,1.5),newpage=FALSE)
dev.off()

p6 = ggplot(subset(meta.combinations.results.plot,chart=='count_increase')) +
	geom_bar(aes(combination,value,fill=n_cells==1),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0),limits=c(0,ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Upregulated') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
p7 = ggplot(subset(meta.combinations.results.plot,chart=='count_decrease')) +
	geom_bar(aes(combination,value,fill=n_cells==1),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(trans='reverse',expand=c(0,0),limits=c(ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100,0)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Downregulated') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

pdf(file=paste0('figures/',dataset,'_model_results_',predictor,'_',model.method,'_mashr_comparison_beta_counts_meta_single.pdf'),useDingbats=FALSE,height=7,width=11)
	ggarrange(p7,p1,p6,ncol=3,nrow=1,widths=c(1.5,1,1.5),newpage=FALSE)
dev.off()






# Inputs
# mash.lfsr: matrix with regions as columns, genes as rows, and mashr LFSRs as values
# mash.beta: matrix with regions as columns, genes as rows, and mashr betas as values
# region.levels: vector with regions in preferred order
# region.colors: vector with color choices in hex RGB.
threshold.range = seq(0.01,0.2,0.01)
# Loop through a series of cutoffs
which.cell.tally = do.call(rbind,lapply(threshold.range,function(threshold) {
	# raw counts
	# unique counts
	# Create matrices on only those genes that are significant in that direction in only one cell
	this.inc.lfsr = mash.lfsr[rowSums(mash.lfsr < threshold & mash.beta > 0) == 1,]
	this.inc.beta = mash.beta[rowSums(mash.lfsr < threshold & mash.beta > 0) == 1,]
	this.dec.lfsr = mash.lfsr[rowSums(mash.lfsr < threshold & mash.beta < 0) == 1,]
	this.dec.beta = mash.beta[rowSums(mash.lfsr < threshold & mash.beta < 0) == 1,]
	out = rbind(
		# Counts/proportions of total significant genes with positive betas in each cell
		data.frame(
			threshold,
			cell = colnames(mash.lfsr),
			count = as.integer(colSums(mash.lfsr < threshold)),
			proportion = as.numeric(
				colSums(mash.lfsr < threshold) /
					sum(mash.lfsr < threshold)
			),
			direction = 'Increase or decrease',
			type = 'total'
		),
		data.frame(
			threshold,
			cell = colnames(mash.lfsr),
			count = as.integer(colSums(mash.lfsr < threshold & mash.beta > 0)),
			proportion = as.numeric(
				colSums(mash.lfsr < threshold & mash.beta > 0) /
					sum(mash.lfsr < threshold & mash.beta > 0)
			),
			direction = 'Increase',
			type = 'total'
		),
		# Counts/proportions of total significant genes with negative betas in each cell
		data.frame(
			threshold,
			cell = colnames(mash.lfsr),
			count = as.integer(colSums(mash.lfsr < threshold & mash.beta < 0)),
			proportion = as.numeric(
				colSums(mash.lfsr < threshold & mash.beta < 0) /
					sum(mash.lfsr < threshold & mash.beta < 0)
			),
			direction = 'Decrease',
			type = 'total'
		),
		# Counts/proportions of unique significant genes with positive betas in each cell
		data.frame(
			threshold,
			cell = colnames(mash.lfsr),
			count = as.integer(colSums(this.inc.lfsr < threshold)),
			proportion = as.numeric(
				colSums(this.inc.lfsr < threshold) /
					sum(this.inc.lfsr < threshold)
			),
			direction = 'Increase or decrease',
			type = 'unique'
		),
		# Counts/proportions of unique significant genes with positive betas in each cell
		data.frame(
			threshold,
			cell = colnames(mash.lfsr),
			count = as.integer(colSums(this.inc.lfsr < threshold & this.inc.beta > 0)),
			proportion = as.numeric(
				colSums(this.inc.lfsr < threshold & this.inc.beta > 0) /
					sum(this.inc.lfsr < threshold & this.inc.beta > 0)
			),
			direction = 'Increase',
			type = 'unique'
		),
		# Counts/proportions of unique significant genes with negative betas in each cell
		data.frame(
			threshold,
			cell = colnames(mash.lfsr),
			count = as.integer(colSums(this.inc.lfsr < threshold & this.inc.beta < 0)),
			proportion = as.numeric(
				colSums(this.inc.lfsr < threshold & this.inc.beta < 0) /
					sum(this.inc.lfsr < threshold & this.inc.beta < 0)
			),
			direction = 'Decrease',
			type = 'unique'
		)
	)
	out$direction = factor(out$direction,levels=c('Increase or decrease','Increase','Decrease'))
	out$cell = factor(out$cell,levels=cell.levels)
	out
}))

# The y-axis below can be switched between proportion and count below. The rank order will stay the same
# Plot for total number of significant genes per cell
p = ggplot(droplevels(subset(which.cell.tally,direction %in% c('Increase','Decrease') & type=='total')),aes(threshold,count,color=cell)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=cell.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_counts_total_',tolower(predictor.label),'_cell_tally.pdf'),useDingbats=FALSE,height=5)

# Plot of unique significant genes per cell
p = ggplot(droplevels(subset(which.cell.tally,direction %in% c('Increase','Decrease') & type=='unique')),aes(threshold,count,color=cell)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=cell.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_counts_unique_',tolower(predictor.label),'_cell_tally.pdf'),useDingbats=FALSE,height=5)

# Plot for total number of significant genes per cell
p = ggplot(droplevels(subset(which.cell.tally,direction %in% c('Increase','Decrease') & type=='total')),aes(threshold,proportion,color=cell)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=cell.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_proportions_total_',tolower(predictor.label),'_cell_tally.pdf'),useDingbats=FALSE,height=5)

# Plot of unique significant genes per cell
p = ggplot(droplevels(subset(which.cell.tally,direction %in% c('Increase','Decrease') & type=='unique')),aes(threshold,proportion,color=cell)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=cell.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_proportions_unique_',tolower(predictor.label),'_cell_tally.pdf'),useDingbats=FALSE,height=5)




# The y-axis below can be switched between proportion and count below. The rank order will stay the same
# Plot for total number of significant genes per cell
p = ggplot(droplevels(subset(which.cell.tally,direction %in% c('Increase or decrease','Increase','Decrease') & type=='total')),aes(threshold,count,color=cell)) +
	geom_line() +
	facet_wrap(~direction,scales='free_y',ncol=1) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=cell.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_counts_total_',tolower(predictor.label),'_cell_tally_withtotal.pdf'),useDingbats=FALSE,height=5)

# Plot of unique significant genes per cell
p = ggplot(droplevels(subset(which.cell.tally,direction %in% c('Increase or decrease','Increase','Decrease') & type=='unique')),aes(threshold,count,color=cell)) +
	geom_line() +
	facet_wrap(~direction,scales='free_y',ncol=1) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=cell.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_counts_unique_',tolower(predictor.label),'_cell_tally_withtotal.pdf'),useDingbats=FALSE,height=5)

# Plot for total number of significant genes per cell
p = ggplot(droplevels(subset(which.cell.tally,direction %in% c('Increase or decrease','Increase','Decrease') & type=='total')),aes(threshold,proportion,color=cell)) +
	geom_line() +
	facet_wrap(~direction,scales='free_y',ncol=1) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=cell.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_proportions_total_',tolower(predictor.label),'_cell_tally_withtotal.pdf'),useDingbats=FALSE,height=5)

# Plot of unique significant genes per cell
p = ggplot(droplevels(subset(which.cell.tally,direction %in% c('Increase or decrease','Increase','Decrease') & type=='unique')),aes(threshold,proportion,color=cell)) +
	geom_line() +
	facet_wrap(~direction,scales='free_y',ncol=1) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=cell.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_proportions_unique_',tolower(predictor.label),'_cell_tally_withtotal.pdf'),useDingbats=FALSE,height=5)


p = ggplot(droplevels(subset(which.cell.tally,direction %in% c('Increase or decrease','Increase','Decrease') & type=='total')),aes(threshold,count,color=cell)) +
	geom_line() +
	facet_wrap(~direction,scales='free_y',ncol=1) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=cell.colors) +
	theme_article(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_counts_total_',tolower(predictor.label),'_cell_tally_withtotal_presentation.pdf'),useDingbats=FALSE,height=5)




# Number of genes in n cells
n.cell.tally = do.call(rbind,lapply(threshold.range,function(threshold) {
	# raw counts
	table(rowSums(mash.lfsr < threshold & mash.beta > 0))
	table(rowSums(mash.lfsr < threshold & mash.beta < 0))
	out = rbind(
		data.frame(
			threshold,
			n.cells = as.integer(names(table(rowSums(mash.lfsr < threshold & mash.beta > 0)))),
			count = as.integer(table(rowSums(mash.lfsr < threshold & mash.beta > 0))),
			proportion = as.numeric(
				table(rowSums(mash.lfsr < threshold & mash.beta > 0)) / nrow(mash.lfsr)
			),
			direction = 'Increase'
		),
		data.frame(
			threshold,
			n.cells = as.integer(names(table(rowSums(mash.lfsr < threshold & mash.beta < 0)))),
			count = as.integer(table(rowSums(mash.lfsr < threshold & mash.beta < 0))),
			proportion = as.numeric(
				table(rowSums(mash.lfsr < threshold & mash.beta < 0)) / nrow(mash.lfsr)
			),
			direction = 'Decrease'
		)
	)
	out$direction = factor(out$direction,levels=c('Increase','Decrease'))
	out
}))
p = ggplot(droplevels(subset(n.cell.tally,n.cells>0)),aes(threshold,count,color=factor(n.cells))) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
#	scale_color_manual(values=c('#000000',cell.colors)) +
	scale_color_manual(values=cell.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_counts_',tolower(predictor.label),'_n_cells_tally.pdf'),useDingbats=FALSE,height=5)

p = ggplot(droplevels(subset(n.cell.tally,n.cells>0)),aes(threshold,proportion,color=factor(n.cells))) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
#	scale_color_manual(values=c('#000000',cell.colors)) +
	scale_color_manual(values=cell.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/',dataset,'_model_results_sig_proportions_',tolower(predictor.label),'_n_cells_tally.pdf'),useDingbats=FALSE,height=5)








go.enrichment.results.all = readRDS(paste0('checkpoints/',dataset,'_topgo_results_',model.method,'_mashr.rds'))
# disease.enrichment.results.all = readRDS('checkpoints/disease_enrichment_results.rds')

i = 'kst'

go.enrichment.results = readRDS(paste0('checkpoints/',dataset,'_topgo_results_',model.method,'_mashr.rds'))

# Visualize GO results
p = ggplot(
	do.call(rbind,lapply(split(go.enrichment.results,list(go.enrichment.results$direction, go.enrichment.results$cell)),function(x) {
		within(data.frame(unique(x[c('direction','cell')]),expected=-log10(seq(1/nrow(x),1,1/nrow(x))),observed=-log10(quantile(x$pval,seq(1/nrow(x),1,1/nrow(x)),na.rm=TRUE))),{
			cell=factor(cell,levels=cell.levels,labels=cell.short.levels)
			direction=factor(direction,levels=c('increase','decrease'),labels=c('Increase','Decrease'))
		})
	})),
	aes(expected,observed,color=direction)) +
#		geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
	geom_point_rast(size=0.1) +
	geom_abline(slope=1,col='black',size=0.2) +
	facet_wrap(~cell,nrow=cell.rows) +
	scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(go.enrichment.results$direction,go.enrichment.results$cell))-1)))),1)) +
	scale_y_continuous(breaks=seq(0,ceiling(max(-log10(go.enrichment.results$pval))),1)) +
	scale_color_manual(name='Direction',values=c('#386cb0','#f0027f')) +
#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
	xlab(expression(-log[10] * ('Expected'~italic(p)))) +
	ylab(expression(-log[10] * ('Observed'~italic(p)))) +
	theme_classic() +
	theme() +
	guides(color=guide_legend(override.aes=list(size=3)))
ggsave(p,file=paste0('figures/',dataset,'_model_results_pval_enrichment_topgo_',tolower(predictor.label),'_qqplot.pdf'),useDingbats=FALSE)
