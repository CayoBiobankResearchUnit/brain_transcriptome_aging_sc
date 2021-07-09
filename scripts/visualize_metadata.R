#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)
library(RColorBrewer)

dataset = arguments[1]

cds = readRDS(paste0('checkpoints/cell_data_sets/',dataset,'.rds'))

meta = read.delim(paste0('data/metadata_',dataset,'.txt'),stringsAsFactors=FALSE)

meta = meta[meta$sample %in% cds$sample,]

p = ggplot(meta,aes_string(age.variable,fill='sex')) +
	geom_histogram(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,0.5)) +
	facet_wrap(~group,ncol=1) +
	scale_fill_manual(values=sex.colors) +
	scale_x_continuous(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,5)) +
	scale_y_continuous(limits=c(0,2),breaks=seq(0,2,1)) +
	xlab('Age (years)') + ylab('Count') +
	theme_classic(base_size=24) +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/',dataset,'_metadata_age_histogram.pdf'),useDingbats=FALSE,height=5)

p = ggplot(mutate(meta,group=paste('group',group)),aes_string(age.variable,fill='sex')) +
	geom_histogram(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,0.5)) +
	facet_wrap(~group,ncol=1) +
	scale_fill_manual(values=sex.colors) +
	scale_x_continuous(breaks=seq(0,30,10)) +
	scale_y_continuous(limits=c(0,3),breaks=seq(0,3,1)) +
	xlab('Age (years)') + ylab('Count') +
	theme_classic(base_size=24) +
	theme(legend.position='none',strip.background=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_metadata_age_histogram_article.pdf'),useDingbats=FALSE,height=5)

p = ggplot(meta,aes_string(age.variable,fill='sex')) +
	geom_histogram(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,0.5)) +
	facet_wrap(~group,ncol=1) +
	scale_fill_manual(values=sex.colors) +
	scale_x_continuous(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,5)) +
	scale_y_continuous(limits=c(0,2),breaks=seq(0,2,1)) +
	xlab('Age (years)') + ylab('Count') +
	theme_classic(base_size=16) +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/',dataset,'_metadata_age_histogram_small.pdf'),useDingbats=FALSE,height=4)

p = ggplot(meta,aes_string(age.variable,'1',color='sex')) +
	geom_point() +
	facet_wrap(~group,ncol=1) +
	scale_color_manual(values=sex.colors) +
	scale_x_continuous(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,5)) +
	scale_y_continuous(limits=c(0,2),breaks=seq(0,2,1)) +
	xlab('Age (years)') +
	theme_classic(base_size=24) +
	theme(legend.position='none',axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_metadata_age_dotplot.pdf'),useDingbats=FALSE,height=5)

p = ggplot(meta,aes_string(age.variable,'1',color='sex')) +
	geom_point() +
	facet_wrap(~group,ncol=1) +
	scale_color_manual(values=sex.colors) +
	scale_x_continuous(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,5)) +
	scale_y_continuous(limits=c(0,2),breaks=seq(0,2,1)) +
	xlab('Age (years)') +
	theme_classic(base_size=16) +
	theme(legend.position='none',axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_metadata_age_dotplot_small.pdf'),useDingbats=FALSE,height=4)

p = ggplot(meta,aes_string('group',age.variable,color='sex')) +
	geom_point() +
	scale_color_manual(values=sex.colors) +
	scale_y_continuous(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,5)) +
	ylab('Age (years)') +
	theme_classic(base_size=24) +
	theme(legend.position='none',axis.title.x=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_metadata_age_dotplot2.pdf'),useDingbats=FALSE,width=5)

p = ggplot(meta,aes_string(1,age.variable,color='sex')) +
	geom_point() +
	scale_color_manual(values=sex.colors) +
	scale_y_continuous(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,5)) +
	ylab('Age (years)') +
	theme_classic(base_size=24) +
	theme(legend.position='none',axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
ggsave(p,file=paste0('figures/',dataset,'_metadata_age_dotplot3.pdf'),useDingbats=FALSE,width=5)

p = ggplot(meta,aes_string(age.variable,fill='sex')) +
	geom_histogram(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,0.5)) +
	scale_fill_manual(values=sex.colors) +
	scale_x_continuous(breaks=seq(0,ceiling(max(meta[[age.variable]])/5)*5,5)) +
	coord_flip() +
	xlab('Count') +
	ylab('Age (years)') +
	theme_classic(base_size=24) +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/',dataset,'_metadata_age_histogram2.pdf'),useDingbats=FALSE,width=5)
