#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(tidyverse)
library(reshape2)
library(mashr)

dataset = arguments[1]

if (length(arguments) > 1) {
	model.method = arguments[2]
} else {
	model.method = 'emma'
}

model.output = readRDS(paste0('checkpoints/',dataset,'_model_',model.method,'.rds'))

model.stats = bind_model_outputs(model.output,predictor)

model.stats[[paste('pval',predictor,sep='.')]][model.stats$converged %in% FALSE] = NA

model.stats = do.call(rbind,lapply(split(model.stats,model.stats$cell),function(x) {
	x[[paste('qval',predictor,sep='.')]] = p.adjust(x[[paste('pval',predictor,sep='.')]],'fdr')
	x
}))
rownames(model.stats) = NULL

ensembl.gene.names = sort(unique(model.stats$gene))

model.beta = do.call(cbind,lapply(cell.levels,function(x) {
	bind_model_statistics(x,model.stats,paste('beta',predictor,sep='.'))
}))
model.bvar = do.call(cbind,lapply(cell.levels,function(x) {
	bind_model_statistics(x,model.stats,paste('bvar',predictor,sep='.'))
}))
model.qval = do.call(cbind,lapply(cell.levels,function(x) {
	bind_model_statistics(x,model.stats,paste('qval',predictor,sep='.'))
}))
model.sbeta = model.beta / sqrt(model.bvar)

if (ignore.checkpoints || !file.exists(paste0('checkpoints/',dataset,'_rnaseq_genes_go.rds'))) {
	library(biomaRt)

	mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='mmulatta_gene_ensembl')

	mmul.go = getBM(
		attributes=c('ensembl_gene_id','go_id'),
		filters = 'ensembl_gene_id',
		values = ensembl.gene.names,
		mart = mmul)
	saveRDS(mmul.go,file=paste0('checkpoints/',dataset,'_rnaseq_genes_go.rds'))
} else {
	message('Checkpoint found!\nLoading GO annotations from file.')

	mmul.go = readRDS(paste0('checkpoints/',dataset,'_rnaseq_genes_go.rds'))
}

gene2go = lapply(
	unique(mmul.go$ensembl_gene_id),
	function(x) {
		out = sort(mmul.go[mmul.go$ensembl_gene_id == x,'go_id'])
		out[out != '']
	}
)

names(gene2go) = unique(mmul.go$ensembl_gene_id)

library(topGO)

# Fetch names for GO terms identified as relations by topGO
if (ignore.checkpoints || !file.exists(paste0('checkpoints/',dataset,'_rnaseq_go_names.rds'))) {
	# Create dummy vector for topGO
	all.genes = numeric(length(gene2go))
	names(all.genes) = names(gene2go)

	# Initialize topGO object
	go.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=all.genes,geneSelectionFun=function(x) x>0,nodeSize = 10,annotationFun = annFUN.gene2GO,gene2GO = gene2go)

	# Put together all relevant GO terms
	go.ids = sort(union(mmul.go$go_id,go.topgo@graph@nodes))

	# Exclude blank GO terms
	go.ids = go.ids[nchar(go.ids) == 10]

	# Grab names and namespaces from the GO database
	library(GO.db)

	go.terms = data.frame(
		go_id = go.ids,
		go_namespace = Ontology(go.ids),
		go_name = Term(go.ids),
		stringsAsFactors=FALSE
	)

	# Sometimes, there are names that are not in the GO database. Fetch these names from AmiGO

	library(XML)

	to.do = go.terms$go_id[!complete.cases(go.terms)]

	if (length(to.do)) {
		message('Fetching names for ',length(to.do),' GO terms')
		# Fetch missing GO metadata from AmiGO
		for (i in to.do) {
			message(i)
			search.url = paste0('http://amigo.geneontology.org/amigo/term/',i)

			amigo.data = xpathApply(htmlParse(search.url), '//dd', xmlValue)
			go.terms$go_namespace[match(i,go.terms$go_id)] = unlist(lapply(strsplit(amigo.data[[3]],'_'),function(x) paste(toupper(substr(x,1,1)),collapse='')))
			go.terms$go_name[match(i,go.terms$go_id)] = amigo.data[[2]]
		}
	}
	rownames(go.terms) = go.terms$go_id
	saveRDS(go.terms,file=paste0('checkpoints/',dataset,'_rnaseq_go_names.rds'))
} else {
	# If checkpoint is found, use it
	message('Checkpoint found!\nLoading GO metadata from file.')

	go.terms = readRDS(paste0('checkpoints/',dataset,'_rnaseq_go_names.rds'))
}

keep.cells = names(which(unlist(lapply(split(model.stats,model.stats$cell),function(x) mean(sort(x[[paste('pval',predictor,sep='.')]],decreasing=TRUE)[1:ceiling(nrow(x)/100)]))) > qq.cutoff))

# 
single.cell.results = do.call(rbind,lapply(intersect(cell.levels,keep.cells),function(i) {

#out = vector(mode='list',length=length(intersect(cell.levels,keep.cells)))
#names(out) = intersect(cell.levels,keep.cells)
#
#for (i in intersect(cell.levels,keep.cells)) {
#message(i)
	this.cell = model.sbeta[,i]
	this.cell = this.cell[!is.na(this.cell)]

	this.cell.inc.topgo = new('topGOdata',
		description='Simple session',
		ontology='BP',
		allGenes=-this.cell,
		geneSelectionFun=function(x) x > 0, # ignored
		nodeSize = 10,
		annotationFun = annFUN.gene2GO,
		gene2GO = gene2go
	)

	this.cell.dec.topgo = new('topGOdata',
		description='Simple session',
		ontology='BP',
		allGenes=this.cell,
		geneSelectionFun=function(x) x < 0, # ignored
		nodeSize = 10,
		annotationFun = annFUN.gene2GO,
		gene2GO = gene2go
	)

	this.cell.inc.test = runTest(this.cell.inc.topgo,algorithm='weight01',statistic='KS')
	this.cell.dec.test = runTest(this.cell.dec.topgo,algorithm='weight01',statistic='KS')

	this.cell.inc.results = reshape2::melt(this.cell.inc.test@score)
	this.cell.dec.results = reshape2::melt(this.cell.dec.test@score)

	go.i.inc = gsub('\\..*$','',rownames(this.cell.inc.results))

	this.cell.inc.results = data.frame(
		go_id = go.i.inc,
		go.terms[go.i.inc,setdiff(colnames(go.terms),'go_id')],
		pval = this.cell.inc.results$value,
		qval = p.adjust(this.cell.inc.results$value,'fdr'),
		direction = 'increase',
		set = 'single',
		cell = i,
		model = model.method,
		stringsAsFactors=FALSE
	)

	go.i.dec = gsub('\\..*$','',rownames(this.cell.dec.results))

	this.cell.dec.results = data.frame(
		go_id = go.i.dec,
		go.terms[go.i.dec,setdiff(colnames(go.terms),'go_id')],
		pval = this.cell.dec.results$value,
		qval = p.adjust(this.cell.dec.results$value,'fdr'),
		direction = 'decrease',
		set = 'single',
		cell = i,
		model = model.method,
		stringsAsFactors=FALSE
	)

#	out[[i]] = rbind(this.cell.inc.results,this.cell.dec.results)
#}

	rbind(this.cell.inc.results,this.cell.dec.results)
}))

single.cell.results = single.cell.results[order(
	match(single.cell.results$cell,cell.levels),
	match(single.cell.results$direction,c('increase','decrease')),
	single.cell.results$pval
),]

# Placeholder (will bind other data later)
go.enrichment.results = single.cell.results

rownames(go.enrichment.results) = NULL

saveRDS(go.enrichment.results,file=paste0('checkpoints/',dataset,'_topgo_results_',model.method,'.rds'))

if (file.exists(paste0('checkpoints/',dataset,'_model_',model.method,'_mashr.rds'))) {
	mash.results = readRDS(paste0('checkpoints/',dataset,'_model_',model.method,'_mashr.rds'))

	mash.beta = get_pm(mash.results)
	mash.berr = get_psd(mash.results)
	mash.lfsr = get_lfsr(mash.results)
	mash.sbet = mash.beta / mash.berr

	mashr.genes = rownames(mash.beta)
	names(mashr.genes) = mashr.genes

	single.cell.results = do.call(rbind,lapply(intersect(cell.levels,keep.cells),function(i) {
		message(i)
		this.cell = mash.sbet[,i]
		this.cell = this.cell[!is.na(this.cell)]

		this.cell.inc.topgo = new('topGOdata',
			description='Simple session',
			ontology='BP',
			allGenes=-this.cell,
			geneSelectionFun=function(x) x > 0, # ignored
			nodeSize = 10,
			annotationFun = annFUN.gene2GO,
			gene2GO = gene2go
		)

		this.cell.dec.topgo = new('topGOdata',
			description='Simple session',
			ontology='BP',
			allGenes=this.cell,
			geneSelectionFun=function(x) x < 0, # ignored
			nodeSize = 10,
			annotationFun = annFUN.gene2GO,
			gene2GO = gene2go
		)

		setMethod('GOKSTest',signature(c(object='classicScore')),get.ks.pval)

		this.cell.inc.test = runTest(this.cell.inc.topgo,algorithm='weight01',statistic='KS')
		this.cell.dec.test = runTest(this.cell.dec.topgo,algorithm='weight01',statistic='KS')

		this.cell.inc.results = reshape2::melt(this.cell.inc.test@score)
		this.cell.dec.results = reshape2::melt(this.cell.dec.test@score)

		go.i.inc = gsub('\\..*$','',rownames(this.cell.inc.results))

		setMethod('GOKSTest',signature(c(object='classicScore')),get.ks.score)

		this.cell.inc.results = data.frame(
			go_id = go.i.inc,
			go.terms[go.i.inc,setdiff(colnames(go.terms),'go_id')],
			score = runTest(this.cell.inc.topgo,algorithm='weight01',statistic='KS')@score,
			pval = this.cell.inc.results$value,
			qval = p.adjust(this.cell.inc.results$value,'fdr'),
			direction = 'increase',
			set = 'single',
			cell = i,
			model = 'mash',
			stringsAsFactors=FALSE
		)

		go.i.dec = gsub('\\..*$','',rownames(this.cell.dec.results))

		this.cell.dec.results = data.frame(
			go_id = go.i.dec,
			go.terms[go.i.dec,setdiff(colnames(go.terms),'go_id')],
			score = runTest(this.cell.dec.topgo,algorithm='weight01',statistic='KS')@score,
			pval = this.cell.dec.results$value,
			qval = p.adjust(this.cell.dec.results$value,'fdr'),
			direction = 'decrease',
			set = 'single',
			cell = i,
			model = 'mash',
			stringsAsFactors=FALSE
		)

		rbind(this.cell.inc.results,this.cell.dec.results)
	}))

	single.cell.results = single.cell.results[order(
		match(single.cell.results$cell,cell.levels),
		match(single.cell.results$direction,c('increase','decrease')),
		single.cell.results$pval
	),]

#	go.enrichment.results = rbind(go.enrichment.results,single.cell.results)

	saveRDS(single.cell.results,file=paste0('checkpoints/',dataset,'_topgo_results_',model.method,'_mashr.rds'))

	mash.enrichment.results = single.cell.results

	n.permutations = 1000

	get.cell.specificity = function(enrichment.results,effect.direction) {
		this.df = subset(enrichment.results,direction == effect.direction)
		this.df = subset(this.df,complete.cases(this.df))
		this.wide = tidyr::pivot_wider(droplevels(this.df),id_cols='go_id',names_from='cell',values_from='score')
		this.mat = matrix(as.matrix(this.wide[,2:ncol(this.wide)]),nrow=nrow(this.wide),ncol=ncol(this.wide)-1,dimnames=list(this.wide$go_id,colnames(this.wide)[2:ncol(this.wide)]))
	# 	apply(this.mat,1,function(x) max(x) - quantile(x,0.75))
		specificity = apply(this.mat,1,function(x) sort(x)[length(x)] - sort(x)[length(x)-1])
		specificity.null = do.call(cbind,parallel::mclapply(1:n.permutations,function(i) {
			apply(matrix(apply(this.mat,2,sample),nrow=nrow(this.mat),dimnames=dimnames(this.mat)),1,function(x) sort(x)[length(x)] - sort(x)[length(x)-1])
		},mc.cores=n.cores))
		data.frame(
			go_id = rownames(this.mat),
			go_name = as.character(unique(subset(this.df,select=c('go_id','go_name')))$go_name),
			direction = effect.direction,
			top.cell = colnames(this.mat)[apply(this.mat,1,which.max)],
			specificity,
			specificity.pval = rowMeans(specificity < specificity.null),
			specificity.qval = p.adjust(rowMeans(specificity < specificity.null),'fdr')
		)
	}

	go.kst.inc.spec = get.cell.specificity(mash.enrichment.results,'increase')
	go.kst.dec.spec = get.cell.specificity(mash.enrichment.results,'decrease')

	go.specificity.results = rbind(
		go.kst.inc.spec,
		go.kst.dec.spec
	)

	go.specificity.results$top.cell.score = apply(go.specificity.results,1,function(x) {
		subset(mash.enrichment.results,go_name == x['go_name'] & cell == x['top.cell'] & direction == x['direction'])$score
	})
	go.specificity.results$top.cell.pval = apply(go.specificity.results,1,function(x) {
		subset(mash.enrichment.results,go_name == x['go_name'] & cell == x['top.cell'] & direction == x['direction'])$pval
	})
	go.specificity.results$top.cell.qval = apply(go.specificity.results,1,function(x) {
		subset(mash.enrichment.results,go_name == x['go_name'] & cell == x['top.cell'] & direction == x['direction'])$qval
	})

	saveRDS(go.specificity.results,file='checkpoints/topgo_specificity_results.rds')

}