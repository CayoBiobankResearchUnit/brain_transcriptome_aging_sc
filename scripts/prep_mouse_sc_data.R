#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(Matrix)

dataset = arguments[1]
n.subsets = as.integer(arguments[2])

# If a checkpoint does not exist, make the mouse ortholog matrix
if (!file.exists('checkpoints/musmus_mmul_matrix.rds')) {
	# If the raw matrix does not exist in binary form, make it
	if (!file.exists('checkpoints/musmus_matrix.rds')) {
		# If the CSV  does not exist, download it
		if (!file.exists('data/musmus_matrix.csv')) download.file('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/matrix.csv',destfile='data/musmus_matrix.csv')

		musmus.matrix = read.csv('data/musmus_matrix.csv',row.names=1)
		musmus.matrix = as.matrix(musmus.matrix)
		musmus.matrix = t(musmus.matrix)

		saveRDS(musmus.matrix,file='checkpoints/musmus_matrix.rds')
	} else {
		musmus.matrix = readRDS('checkpoints/musmus_matrix.rds')
	}

	# Subset only 1:1 orthologous genes
	library(biomaRt)
	mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='mmulatta_gene_ensembl')

	mmul.mmus = getBM(
		attributes=c('ensembl_gene_id','external_gene_name','mmusculus_homolog_ensembl_gene','mmusculus_homolog_associated_gene_name','mmusculus_homolog_orthology_type'),
		mart = mmul)
	mmul.mmus = subset(mmul.mmus,mmusculus_homolog_orthology_type == 'ortholog_one2one')

	# Read in list of genes that are in the macaque dataset
	mmul.genes = readRDS(paste0('checkpoints/',dataset,'_mmul_genes.rds'))

	mmul.mmus = subset(mmul.mmus,ensembl_gene_id %in% mmul.genes & mmusculus_homolog_associated_gene_name != '' & (! mmusculus_homolog_associated_gene_name %in% names(table(mmul.mmus$mmusculus_homolog_associated_gene_name)[table(mmul.mmus$mmusculus_homolog_associated_gene_name) > 1])))

	musmus.ortho = musmus.matrix[rownames(musmus.matrix)[rownames(musmus.matrix) %in% mmul.mmus$mmusculus_homolog_associated_gene_name],]

	# Convert gene names
	rownames(mmul.mmus) = mmul.mmus$mmusculus_homolog_associated_gene_name
	rownames(musmus.ortho) = mmul.mmus[rownames(musmus.ortho),'ensembl_gene_id']

	# Sample 200,000 cells
	set.seed(seed)
	musmus.ortho = musmus.ortho[,sample(1:ncol(musmus.ortho),200000)]

	saveRDS(musmus.ortho,file='checkpoints/musmus_mmul_matrix.rds')
} else {
	musmus.ortho = readRDS('checkpoints/musmus_mmul_matrix.rds')
}

# If the mouse metadata checkpoint does not exist, create it
if (!file.exists('checkpoints/musmus_metadata.rds')) {
	if (!file.exists('data/musmus_metadata.csv')) download.file('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/metadata.csv',destfile='data/musmus_metadata.csv')
	musmus.meta = read.csv('data/musmus_metadata.csv')
	rownames(musmus.meta) = musmus.meta[[1]]
	musmus.meta = musmus.meta[colnames(musmus.ortho),]
	saveRDS(musmus.meta,file='checkpoints/musmus_metadata.rds')
} else {
	musmus.meta = readRDS('checkpoints/musmus_metadata.rds')
}

# Ignore cells that have no reference labels
musmus.meta = subset(musmus.meta,class_label != '')

# Equalize cell orders in the gene x cell matrix and the metadata
musmus.ortho = musmus.ortho[,rownames(musmus.meta)]

# Save macaque dataset as Seurat

if (!file.exists(paste0('checkpoints/',dataset,'_seurat_sct.rds'))) {
	library(monocle3)
	macmul.cds = readRDS(paste0('checkpoints/',dataset,'_classified_manual.rds'))

	library(Seurat)
	macmul.seurat = CreateSeuratObject(counts=assays(macmul.cds)$counts,meta.data=as.data.frame(colData(macmul.cds)),project=dataset)
	macmul.sct = SCTransform(macmul.seurat,verbose=FALSE)

	saveRDS(macmul.sct,file=paste0('checkpoints/',dataset,'_seurat_sct.rds'))
}

library(Seurat)
musmus.seurat = CreateSeuratObject(counts=musmus.ortho,meta.data=musmus.meta,project='mouse')
musmus.sct = SCTransform(musmus.seurat,verbose=FALSE)

saveRDS(musmus.sct,file=paste0('checkpoints/musmus_seurat_sct.rds'))
