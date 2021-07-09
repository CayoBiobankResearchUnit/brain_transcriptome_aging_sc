#!/usr/bin/env Rscript

# Function for converting marker file in one species to another

convert.marker.file = function(in.file,species.old,species.new,out.file,gene_id_type_in='SYMBOL',gene_id_type_out='ENSEMBL') {

	marker.file = suppressMessages(scan(file=in.file,what='',sep='\n',quiet=TRUE))

	# Clean up stray tabs and trailing spaces
	marker.file = marker.file %>% { gsub('\t',' ',.) } %>% { gsub(' *$','',.) }

	# Find cell type classes via pattern matching
	marker.cells = marker.file[grep('^>',marker.file)]

	# Infer start:end line numbers associated with each cell type definition
	marker.lines = apply(cbind(grep('^>',marker.file),c(grep('^>',marker.file)[2:sum(grepl('^>',marker.file))],length(marker.file)+1)),1,function(x) x[1]:(x[2]-1))

	# Find subtype definitions via pattern matching (return one line per cell type
	#      definition to preserve geometry)
	marker.subtypes = unlist(lapply(1:length(marker.cells),function(i) {
		out = marker.file[marker.lines[[i]]]
		out = out[grep('^subtype',out)]
		if (length(out)) out else NA
	}))
	# Find expressed marker gene definitions and return a list with one vector of
	#      parsed marker genes per cell type definition
	marker.genes = marker.file[grep('^ *expressed:',marker.file)] %>%
						{ gsub('#.*','',.) } %>%
						{ gsub('\\t','',.) } %>% { gsub(' *','',.) } %>%
						{ gsub('^expressed:','',.) } %>%
						strsplit(',')

	# Load mouse biomaRt database
	mart = biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset=paste0(species.old,'_gene_ensembl'))

	lookup.column = if (gene_id_type_in == 'SYMBOL') 'external_gene_name' else if (gene_id_type_in == 'ENSEMBL') 'ensembl_gene_id'
	output.column = if (gene_id_type_out == 'SYMBOL') paste0(species.new,'_homolog_associated_gene_name') else if (gene_id_type_out == 'ENSEMBL') paste0(species.new,'_homolog_ensembl_gene')
	species.lookup = biomaRt::getBM(attributes=c('ensembl_gene_id','external_gene_name',
										paste0(species.new,'_homolog_ensembl_gene'),
										paste0(species.new,'_homolog_associated_gene_name'),
										paste0(species.new,'_homolog_orthology_type')),
							filters = lookup.column,
							values = unique(unlist(marker.genes)),
							mart = mart)

	# Return only one2one orthologs
	species.match = species.lookup[species.lookup[[paste0(species.new,'_homolog_orthology_type')]] == 'ortholog_one2one',]
	species.match[[paste0(species.new,'_gene')]] = species.match[[output.column]]

	# If gene name is blank, substitute the ensembl gene ID
	species.match[[paste0(species.new,'_gene')]][species.match[[paste0(species.new,'_gene')]] == ''] = species.match[[paste0(species.new,'_homolog_ensembl_gene')]][species.match[[paste0(species.new,'_gene')]] == '']

	# Calculate named vector to convert mouse genes to macaque genes
	old2new = species.match[[paste0(species.new,'_gene')]]
	names(old2new) = species.match[[lookup.column]]

	# Reformat converted macaque genes based on the expressed genes line of
	#       the marker gene file
	marker.genes.new = lapply(marker.genes,function(x) {
		out = as.character(old2new[x])
		out = out[!is.na(out)]
		if (length(out)) {
			paste('expressed:',paste(out,collapse=', '))
		} else {
			NA
		}
	})

	# Create broken links vector as those cell types that lack genes due to missing
	#     macaque orthologs
	broken.link = unlist(lapply(marker.genes.new,is.na))

	marker.file[setdiff(grep('^ *expressed:',marker.file),unlist(marker.lines[broken.link]))] = paste0(unlist(lapply(marker.genes.new[!broken.link],function(x) { x[is.na(x)] = ''; x})),'\t# ',marker.file[setdiff(grep('^ *expressed:',marker.file),unlist(marker.lines[broken.link]))])

	# For cell types that have been excluded due to missing orthologs, need to
	#     fix broken subtype references. Otherwise, these will trigger errors
	#     in garnett

	# Calculate the lines that indicate subtype relationships
	marker.subtype.lines = unlist(lapply(1:length(marker.lines),function(i) {
		these.lines = marker.lines[[i]]
		if (sum(grepl('^subtype of: *',marker.file[these.lines]))) {
			these.lines[grep('^subtype of: *',marker.file[these.lines])]
		} else {
			NA
		}
	}))

	# Save the original so that we can track changes later
	marker.subtypes.original = marker.file[marker.subtype.lines]

	if (any(!is.na(marker.subtype.lines))) {
		# Iteratively correct broken subtypes references by assigning the parent class of
		#      the broken subtype (sometimes the parent is also missing, hence the loop)
		while(!all(broken.link | is.na(marker.file[marker.subtype.lines]) | gsub('^subtype of: *','',marker.file[marker.subtype.lines]) %in% gsub('^>','',marker.cells[!broken.link]))) {
			temp.subtypes = unique(marker.file[marker.subtype.lines][!is.na(marker.file[marker.subtype.lines]) & !broken.link])
			bad.subtypes = temp.subtypes[!gsub('^subtype of: *','',temp.subtypes) %in% gsub('^>','',marker.cells[!broken.link])]
			temp.replacements = do.call(rbind,lapply(bad.subtypes,function(i) {
				cbind(i,marker.file[marker.subtype.lines][grep(gsub('^subtype of: *','',i),gsub('^>','',marker.cells))])
			}))
			for (i in 1:nrow(temp.replacements)) {
				marker.file[setdiff(which(marker.file %in% temp.replacements[i,1]),unlist(marker.lines[broken.link]))] = temp.replacements[i,2]
			}
		}

		marker.file[is.na(marker.file)] = ''

		# Wherever the lines were overwritten, replace with the updated subtype line and
		#     append the old line as a comment
		marker.file[marker.subtype.lines[!is.na(marker.subtype.lines) & marker.file[marker.subtype.lines] != marker.subtypes.original]] = paste0(marker.file[marker.subtype.lines][!is.na(marker.subtype.lines) & marker.file[marker.subtype.lines] != marker.subtypes.original],'\t# ',marker.subtypes.original[!is.na(marker.subtype.lines) & marker.file[marker.subtype.lines] != marker.subtypes.original])
	}
	# The lines of the marker file that match unbroken cell type definitions are
	#       finally returned and written to a new marker gene file
	marker.file[unlist(marker.lines[broken.link])] = paste0('#\t',marker.file[unlist(marker.lines[broken.link])])

	# # Iteratively check to see if any cell types in the broken.links vector are
	# #       referenced by other cell subtypes. These links are then added to the
	# #       broken links vector until no missing references remain.
	# while (any(!(gsub('^subtype of: *','',marker.subtypes[!is.na(marker.subtypes) & !broken.link]) %in% gsub('^>','',marker.cells[!broken.link])))) {
	# 	broken.link = broken.link | (!is.na(marker.subtypes) & !(gsub('^subtype of: *','',marker.subtypes) %in% gsub('^>','',marker.cells[!broken.link])))
	# }

	# Make final file prettier
	marker.file = marker.file[!grepl('^#* *\\t*$',marker.file)]

	marker.file[grep('^#*\\t*>',marker.file)[2:sum(grepl('^#*\\t*>',marker.file))]] = paste0('\n',marker.file[grep('^#*\\t*>',marker.file)[2:sum(grepl('^#*\\t*>',marker.file))]])

	# There is a bug that if subtype links were broken and no parent cell types were found,
	#      the cell type definition was interrupted by a commented-out line. To prevent
	#      this, move the comment to the preceding line to preserve a record of changes.
	marker.file.out = paste(marker.file,collapse='\n')
	marker.file.out = gsub('\\n(\\t*# subtype of:)','\\1',marker.file.out)

	# Write finished file
	write(marker.file.out,file=out.file,sep='\n')

}

# Functions
scrub_low_quality = function(cds,umi.min.cutoff=100,mt.max.cutoff=10) {
	cds = cds[,cds$n.umi >= umi.min.cutoff & cds$perc_mitochondrial_umis <= mt.max.cutoff]
	cds
}

reduce_dim_and_cluster = function(cds,k=10) {
	cds = preprocess_cds(cds,method='PCA',num_dim=100)
	cds = reduce_dimension(cds,reduction_method='UMAP',preprocess_method='PCA')
	cds = cluster_cells(cds, cluster_method='leiden',k=k)
	cds
}

run.DoubletFinder = function(cell_data_set, id, collision.rate=NULL) {
	require(R.devices)

	require(Seurat)
	cds.monocle = cell_data_set[,cell_data_set$sample %in% id]
	cds = CreateSeuratObject(exprs(cds.monocle))

	cds = NormalizeData(cds)
	cds = ScaleData(cds)
	cds = FindVariableFeatures(cds, selection.method = 'vst', nfeatures = 2000)

	cds = RunPCA(cds)

	cds = JackStraw(cds, num.replicate=200, verbose=FALSE)
	cds = ScoreJackStraw(cds,dims=1:20)
	pcs_sig = 1:(which(cds@reductions$pca@jackstraw$overall.p.values[,'Score'] > 0.05)[1] - 1)

	require(DoubletFinder)
	sweep.res.list = paramSweep_v3(cds, PCs = pcs_sig, sct = FALSE)
	sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
	bcmvn = suppressGraphics(find.pK(sweep.stats))
	pK_best = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

	if (is.null(collision.rate)) collision.rate = 0.075

	homotypic.prop = modelHomotypic(partitions(cds.monocle))
	nExp_poi = round( collision.rate * ncol(cds) )
	nExp_poi.adj = round( nExp_poi * (1 - homotypic.prop) )

	cds = doubletFinder_v3(cds, PCs = pcs_sig, pN = 0.25, pK = pK_best, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	cds = doubletFinder_v3(cds, PCs = pcs_sig, pN = 0.25, pK = pK_best, nExp = nExp_poi.adj, reuse.pANN = colnames(cds@meta.data)[grep('^pANN_',colnames(cds@meta.data))], sct = FALSE)
	out = cds@meta.data[,4:6]
	colnames(out) = c('DF.pANN','DF.low.conf','DF.high.conf')
	out
}

run.Scrublet = function(cell_data_set, id, collision.rate=NULL) {

	cds.monocle = cell_data_set[,cell_data_set$sample %in% id]

	require(reticulate)
	scr = import('scrublet')

	count.matrix = assay(cds.monocle)

	if (is.null(collision.rate)) {
		scrub = scr$Scrublet(t(count.matrix))
	} else {
		scrub = scr$Scrublet(t(count.matrix),expected_doublet_rate=collision.rate)
	}

#	out = scrub$scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=as.integer(30))

	out = scrub$scrub_doublets()
	out = do.call(data.frame,out)
	names(out) = c('Scr.kNN','Scr.call')
	out[[2]] = c('Singlet','Doublet')[as.integer(out[[2]]) + 1]

	rownames(out) = colnames(cds.monocle)
	
	out
}

find_doublets = function(cds,detection_method='Scrublet',sample.column='sample',collision.rate=NULL) {
	sample.list = unique(colData(cds)[[sample.column]])
	cell.order = rownames(colData(cds))
	if (detection_method == 'Scrublet') {
		doublet_data = do.call(rbind,lapply(sample.list,function(x) run.Scrublet(cds,x,collision.rate=collision.rate)))
		doublet_calls = doublet_data[cell.order,'Scr.call']
	} else if (detection_method == 'DoubletFinder') {
		doublet_data = do.call(rbind,lapply(sample.list,function(x) run.DoubletFinder(cds,x,collision.rate=collision.rate)))
		doublet_calls = doublet_data[cell.order,'DF.high.conf']
	}
	colData(cds) = cbind(colData(cds),doublet_data[cell.order,])

	message('Detected ',sum(doublet_calls != 'Singlet'),' doublets\n')
#	cds = cds[,doublet_calls == 'Singlet']

	cds
}

plot_doublets = function(cds,detection_method='Scrublet',score.column='Scr.kNN',sample.column='sample',outdir=NULL) {
	for (i in sort(unique(colData(cds)[[sample.column]]))) {
		cds.sample = cds[,cds[[sample.column]] == i]
		p = ggplot(as.data.frame(colData(cds.sample)),aes_string(score.column)) + geom_density() + scale_x_continuous(breaks=seq(floor(min(colData(cds.sample)[[score.column]]) * 100)/100,ceiling(max(colData(cds.sample)[[score.column]]) * 100) / 100,0.01),minor_breaks=seq(floor(min(colData(cds.sample)[[score.column]]) * 1000)/1000,ceiling(max(colData(cds.sample)[[score.column]]) * 1000) / 1000,0.001)) + theme_minimal() + theme(axis.text.x=element_text(size=8,angle=-45,hjust=0)) + ggtitle(i)
		if (is.null(outdir)) outdir = '.'
		score.suffix = if (score.column == 'Scr.kNN') 'manual' else if (score.column == 'scrublet_score') 'bbi' else score.column
		ggsave(p,file=paste0(outdir,'/doublet_scores_',i,'_',score.suffix,'.pdf'),useDingbats=FALSE,width=7,height=3)
	}
}

mark_doublets = function(cds,thresholds,detection_method='Scrublet',sample.column='sample',collision.rate=NULL) {
	sample.list = unique(colData(cds)[[sample.column]])
	cell.order = rownames(colData(cds))

	if (detection_method == 'Scrublet') {
		score.column = 'Scr.kNN'
	} else if (detection_method == 'DoubletFinder') {
		score.column = 'DF.pANN'
	}

	if (is.null(collision.rate)) collision.rate = 0.1

	cds.out = cds
	cds.out$doublet.call = 'Singlet'

	cdata = colData(cds.out)

	doublet.cells = unlist(lapply(sample.list,function(i) {
		threshold = thresholds$threshold[thresholds[[sample.column]] == i]
		rownames(subset(cdata,sample ==i))[subset(cdata,sample ==i)$Scr.kNN > threshold]
	}))

	cds.out$doublet.call[colnames(cds.out) %in% doublet.cells] = 'Doublet'

	cds.out
}

remove_doublets = function(cds) {
	cds[,cds$doublet.call == 'Singlet']
}

merge_metadata = function(cds,metadata.file,sample.column='sample') {
	meta = read.delim(metadata.file,stringsAsFactors=FALSE)
	rownames(meta) = meta[[sample.column]]
	meta[[sample.column]] = NULL
	colData(cds) = cbind(colData(cds),meta[colData(cds)[[sample.column]],])
	cds
}

correct_batch = function(cds,batch_column='extraction_batch') {
	cds = align_cds(cds,alignment_group=batch_column)
	cds = reduce_dimension(cds,reduction_method='UMAP',preprocess_method='Aligned')
	cds = cluster_cells(cds, resolution=c(10^seq(-6,-1)),k=k)
	cds
}

# # DEFUNCT (see monocle3's combine_cds)
# merge_cds = function(cds1,cds2) {
# 	require(Seurat)
# 
# 	cds = Reduce(merge,list(CreateAssayObject(exprs(cds1)),CreateAssayObject(exprs(cds2))))
# 	headings.cols = intersect(colnames(colData(cds1)),colnames(colData(cds2)))
# 	headings.rows = intersect(colnames(rowData(cds1)),colnames(rowData(cds2)))
# 	new_cell_data_set(
# 		cds@counts,
# 		rbind(colData(cds1)[,headings.cols],colData(cds2)[,headings.cols]),
# 		unique(rbind(rowData(cds1)[,headings.rows],rowData(cds2)[,headings.rows]))
# 	)
# }

# Function for taking model outputs as a list and binding them into a long-form data frame
bind_model_outputs = function(output,predictor) {
	out = data.frame(
		gene=as.character(unlist(lapply(output,rownames))),
		cell=factor(unlist(lapply(names(output),function(x) { foo = output[[x]]; rep(x,nrow(foo)) })),levels=names(output)),
		pval=as.numeric(unlist(lapply(output,function(y) y[,paste('pval',predictor,sep='.')]))),
		beta=as.numeric(unlist(lapply(output,function(y) y[,paste('beta',predictor,sep='.')]))),
		bvar=as.numeric(unlist(lapply(output,function(y) y[,paste('bvar',predictor,sep='.')]))),
		converged=as.logical(Reduce(union,lapply(output,function(y) if ('converged' %in% names(y)) y[,'converged'] else NA))),
		stringsAsFactors=FALSE
	)
	names(out) = c('gene','cell',paste('pval',predictor,sep='.'),paste('beta',predictor,sep='.'),paste('bvar',predictor,sep='.'),'converged')
	out
}

# Function for taking model statistics and binding them into cell x gene matrices
bind_model_statistics = function(cell.type,stats.df,stat='pval') {
	genes = sort(unique(stats.df$gene))
	stats.df[[stat]][stats.df$converged %in% FALSE] = NA
	val = stats.df[stats.df$cell %in% cell.type,][[stat]]
	names(val) = stats.df[stats.df$cell %in% cell.type,]$gene
	r = rep(NA,length(genes))
	names(r) = genes
	r[names(val)] = val
	matrix(r,ncol=1,dimnames=list(names(r),cell.type))
}

# Credit: https://groups.google.com/d/msg/ggplot2/bVQw8KfXZO8/StxSHQWHBxAJ
qnlog10 = function(p) -log10(qunif(rev(p)))

# Credit: http://goo.gl/K4yh
lm_eqn = function(df){
    m = lm(y ~ x, df);
    eq = substitute(italic(y) == b * italic(x) + a*","~~italic(r)^2~"="~r2, 
            list(   a = as.character(format(coef(m)[1], digits = 2)), 
                    b = as.character(format(coef(m)[2], digits = 2)), 
                    r2 = as.character(format(summary(m)$r.squared, digits = 3))))
    as.character(as.expression(eq));                 
}

# Tweaked Matrix::readMM function to allow for numbers greater than 2147483648 (.Machine$integer.max)
readmm = function (file) 
{
	require(Matrix)
    if (is.character(file)) 
        file <- if (file == "") 
            stdin()
        else file(file)
    if (!inherits(file, "connection")) 
        stop("'file' must be a character string or connection")
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    scan1 <- function(what, ...) scan(file, nmax = 1, what = what, 
        quiet = TRUE, ...)
    if (scan1(character()) != "%%MatrixMarket") 
        stop("file is not a MatrixMarket file")
    if (!(typ <- tolower(scan1(character()))) %in% "matrix") 
        stop(gettextf("type '%s' not recognized", typ), domain = NA)
    if (!(repr <- tolower(scan1(character()))) %in% c("coordinate", 
        "array")) 
        stop(gettextf("representation '%s' not recognized", repr), 
            domain = NA)
    elt <- tolower(scan1(character()))
    if (!elt %in% c("real", "complex", "integer", "pattern")) 
        stop(gettextf("element type '%s' not recognized", elt), 
            domain = NA)
    sym <- tolower(scan1(character()))
    if (!sym %in% c("general", "symmetric", "skew-symmetric", 
        "hermitian")) 
        stop(gettextf("symmetry form '%s' not recognized", sym), 
            domain = NA)
    nr <- scan1(integer(), comment.char = "%")
    nc <- scan1(integer())
    nz <- scan1(numeric()) # initialize as numeric to avoid R's integer limit
    checkIJ <- function(els) {
        if (any(els$i < 1 | els$i > nr)) 
            stop("readMM(): row\t values 'i' are not in 1:nr", 
                call. = FALSE)
        if (any(els$j < 1 | els$j > nc)) 
            stop("readMM(): column values 'j' are not in 1:nc", 
                call. = FALSE)
    }
    if (repr == "coordinate") {
        switch(elt, real = , integer = {
            els <- scan(file, nmax = nz, quiet = TRUE, what = list(i = integer(), 
                j = integer(), x = numeric()))
            checkIJ(els)
            switch(sym, general = {
                new("dgTMatrix", Dim = c(nr, nc), i = els$i - 
                  1L, j = els$j - 1L, x = els$x)
            }, symmetric = {
                new("dsTMatrix", uplo = "L", Dim = c(nr, nc), 
                  i = els$i - 1L, j = els$j - 1L, x = els$x)
            }, `skew-symmetric` = {
                stop("symmetry form 'skew-symmetric' not yet implemented for reading")
                new("dgTMatrix", uplo = "L", Dim = c(nr, nc), 
                  i = els$i - 1L, j = els$j - 1L, x = els$x)
            }, hermitian = {
                stop("symmetry form 'hermitian' not yet implemented for reading")
            }, stop(gettextf("symmetry form '%s' is not yet implemented", 
                sym), domain = NA))
        }, pattern = {
            els <- scan(file, nmax = nz, quiet = TRUE, what = list(i = integer(), 
                j = integer()))
            checkIJ(els)
            switch(sym, general = {
                new("ngTMatrix", Dim = c(nr, nc), i = els$i - 
                  1L, j = els$j - 1L)
            }, symmetric = {
                new("nsTMatrix", uplo = "L", Dim = c(nr, nc), 
                  i = els$i - 1L, j = els$j - 1L)
            }, `skew-symmetric` = {
                stop("symmetry form 'skew-symmetric' not yet implemented for reading")
                new("ngTMatrix", uplo = "L", Dim = c(nr, nc), 
                  i = els$i - 1L, j = els$j - 1L)
            }, hermitian = {
                stop("symmetry form 'hermitian' not yet implemented for reading")
            }, stop(gettextf("symmetry form '%s' is not yet implemented", 
                sym), domain = NA))
        }, complex = {
            stop("element type 'complex' not yet implemented")
        }, stop(gettextf("'%s()' is not yet implemented for element type '%s'", 
            "readMM", elt), domain = NA))
    }
    else stop(gettextf("'%s()' is not yet implemented for  representation '%s'", 
        "readMM", repr), domain = NA)
}