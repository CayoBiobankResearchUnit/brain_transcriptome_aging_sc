#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(Matrix)

dataset = arguments[1]

# If a checkpoint does not exist, make the human ortholog matrix
if (!file.exists('checkpoints/homsap_mmul_matrix.rds')) {
	# If the raw matrix does not exist in binary form, make it
	if (!file.exists('checkpoints/homsap_matrix.rds')) {
		# If the CSV  does not exist, download it
		if (!file.exists('data/homsap_matrix.csv')) download.file('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv',destfile='data/homsap_matrix.csv')

		homsap.matrix = read.csv('data/homsap_matrix.csv',row.names=1)
		homsap.matrix = as.matrix(homsap.matrix)
		homsap.matrix = t(homsap.matrix)

		saveRDS(homsap.matrix,file='checkpoints/homsap_matrix.rds')
	} else {
		homsap.matrix = readRDS('checkpoints/homsap_matrix.rds')
	}

	# Subset only 1:1 orthologous genes
	library(biomaRt)
	mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='mmulatta_gene_ensembl')

	mmul.hsap = getBM(
		attributes=c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_associated_gene_name','hsapiens_homolog_orthology_type'),
		mart = mmul)
	mmul.hsap = subset(mmul.hsap,hsapiens_homolog_orthology_type == 'ortholog_one2one')

	# Read in list of genes that are in the macaque dataset
	mmul.genes = readRDS(paste0('checkpoints/',dataset,'_mmul_genes.rds'))

	mmul.hsap = subset(mmul.hsap,ensembl_gene_id %in% mmul.genes & hsapiens_homolog_associated_gene_name != '' & (! hsapiens_homolog_associated_gene_name %in% names(table(mmul.hsap$hsapiens_homolog_associated_gene_name)[table(mmul.hsap$hsapiens_homolog_associated_gene_name) > 1])))

	homsap.ortho = homsap.matrix[rownames(homsap.matrix)[rownames(homsap.matrix) %in% mmul.hsap$hsapiens_homolog_associated_gene_name],]

	# Convert gene names
	rownames(mmul.hsap) = mmul.hsap$hsapiens_homolog_associated_gene_name
	rownames(homsap.ortho) = mmul.hsap[rownames(homsap.ortho),'ensembl_gene_id']

	saveRDS(homsap.ortho,file='checkpoints/homsap_mmul_matrix.rds')
} else {
	homsap.ortho = readRDS('checkpoints/homsap_mmul_matrix.rds')
}

# If the human metadata checkpoint does not exist, create it
if (!file.exists('checkpoints/homsap_metadata.rds')) {
	if (!file.exists('data/homsap_metadata.csv')) download.file('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv',destfile='data/homsap_metadata.csv')
	homsap.meta = read.csv('data/homsap_metadata.csv')
	rownames(homsap.meta) = homsap.meta[[1]]
	homsap.meta = homsap.meta[colnames(homsap.ortho),]
	saveRDS(homsap.meta,file='checkpoints/homsap_metadata.rds')
} else {
	homsap.meta = readRDS('checkpoints/homsap_metadata.rds')
}

# Ignore cells that have no reference labels
homsap.meta = subset(homsap.meta,class_label != '')

# Equalize cell orders in the gene x cell matrix and the metadata
homsap.ortho = homsap.ortho[,rownames(homsap.meta)]

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
homsap.seurat = CreateSeuratObject(counts=homsap.ortho,meta.data=homsap.meta,project='human')
homsap.sct = SCTransform(homsap.seurat,verbose=FALSE)

saveRDS(homsap.sct,file=paste0('checkpoints/homsap_seurat_sct.rds'))
