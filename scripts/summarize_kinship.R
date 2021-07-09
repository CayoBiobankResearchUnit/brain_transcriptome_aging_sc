#!/usr/bin/env RScript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(monocle3)
library(tidyverse)

dataset = arguments[1]
if (length(arguments) > 1) {
	kinship.file = arguments[2]
} else {
	kinship.file = paste0('data/',dataset,'_ngsrelate_results.txt')
}

cds = readRDS(paste0('checkpoints/',dataset,'_classified_manual.rds'))

ngsrelate = read.delim(kinship.file,stringsAsFactors=FALSE)

rmatrix = matrix(
	NA,
	nrow=length(unique(colData(cds)$sample)),
	ncol=length(unique(colData(cds)$sample)),
	dimnames=list(
		sort(unique(colData(cds)$sample)),
		sort(unique(colData(cds)$sample))
	)
)

diag(rmatrix) = 1
rmatrix[as.matrix(ngsrelate[,c('ida','idb')])] = ngsrelate$rab
rmatrix[as.matrix(ngsrelate[,c('idb','ida')])] = ngsrelate$rab

k = rmatrix

z = matrix(0,nrow=length(unique(colData(cds)$sample)),ncol=nrow(k))
rownames(z) = sort(unique(colData(cds)$sample))
colnames(z) = rownames(k)

diag(z) = 1

saveRDS(list('k'=k,'z'=z),file=paste0('checkpoints/',dataset,'_relatedness.rds'))
