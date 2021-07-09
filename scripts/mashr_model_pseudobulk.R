#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

source('scripts/_include_functions.R')
source('scripts/_include_options.R')

library(mashr)

dataset = arguments[1]

if (length(arguments) > 1) {
	model.method = arguments[2]
} else {
	model.method = 'emma'
}

model.output = readRDS(paste0('checkpoints/',dataset,'_model_',model.method,'.rds'))

model.stats = bind_model_outputs(model.output,predictor)

model.stats = cbind(model.stats,unlist(tapply(model.stats[,paste('pval',predictor,sep='.')],model.stats[,'cell'],p.adjust,method='fdr')))
colnames(model.stats)[ncol(model.stats)] = paste('qval',predictor,sep='.')
rownames(model.stats) = NULL

saveRDS(model.stats,file=paste0('checkpoints/',dataset,'_model_',model.method,'_results.rds'))

# model.array = do.call(cbind,lapply(cell.levels,function(x) {
# 	bind_model_statistics(x,model.stats,paste('pval',predictor,sep='.'))
# }))

Bhat = do.call(cbind,lapply(cell.levels,function(x) {
	bind_model_statistics(x,model.stats,paste('beta',predictor,sep='.'))
}))

Shat =  do.call(cbind,lapply(cell.levels,function(x) {
	sqrt(bind_model_statistics(x,model.stats,paste('bvar',predictor,sep='.')))
}))

# Filter out cell types that show a global enrichment of low p-values

keep.cells = names(which(unlist(lapply(split(model.stats,model.stats$cell),function(x) mean(sort(x[[paste('pval',predictor,sep='.')]],decreasing=TRUE)[1:ceiling(nrow(x)/100)]))) > qq.cutoff))

Bhat = Bhat[,keep.cells]
Shat = Shat[,keep.cells]

# For missing data, set beta to 0 and standard error to 1000
# See https://github.com/stephenslab/mashr/issues/17#issuecomment-330628338
# See also pmid:31320509
Bhat[is.na(Bhat)] = 0
Shat[is.na(Shat)] = 1000

# Create the mashr data object
mash.data = mash_set_data(Bhat,Shat)

# Compute canonical covariance matrices
U.c = cov_canonical(mash.data)  

# # Apply multivariate adaptive shrinkage method in initial model ("naive" run)
# now = Sys.time()
# m.c = mash(mash.data, U.c)
# print(Sys.time() - now) # print time elapsed cause this sometimes takes awhile
# 
# message('Log likelihood (canonical covariance matrix): ',format(get_loglik(m.c),digits=10))
# 
# # Save pairwise sharing matrix
# mash.shared = get_pairwise_sharing(m.c, factor = 0.5)

m.1by1 = mash_1by1(mash.data)
strong.subset = get_significant_results(m.1by1, thresh = 0.05)

# # Get significant results
# strong.subset = get_significant_results(m.c, thresh = 0.05, sig_fn = ashr::get_lfdr)

# # The code below is the preferred method (computes covariance matrix on a null dataset to learn correlation structure among null tests)

# Get random subset (random choose half of all genes)
set.seed(seed)
random.subset = sample(1:nrow(Bhat),ceiling(nrow(Bhat)/2))

# Set temporary objects in order to estimate null correlation structure
temp = mash_set_data(Bhat[random.subset,],Shat[random.subset,])
temp.U.c = cov_canonical(temp)
Vhat = estimate_null_correlation(temp,temp.U.c)
rm(list=c('temp','temp.U.c'))

mash.random = mash_set_data(Bhat[random.subset,],Shat[random.subset,],V=Vhat)
mash.strong = mash_set_data(Bhat[strong.subset,],Shat[strong.subset,], V=Vhat)

# Perform PCA and extreme deconvolution to obtain data-driven covariances
U.pca = cov_pca(mash.strong,5)
U.ed = cov_ed(mash.strong, U.pca)

# Fit mash model
U.c = cov_canonical(mash.random)

now = Sys.time()
m.r = mash(mash.random, Ulist = c(U.ed,U.c), outputlevel = 1)
print(Sys.time() - now)

now = Sys.time()
m = mash(mash.data, g=get_fitted_g(m.r), fixg=TRUE)
print(Sys.time() - now)

message('Log likelihood (final model): ',format(get_loglik(m),digits=10))


# # The code below is an alternative method (does not compute covariance matrix on a null dataset)
# # Perform PCA and extreme deconvolution
# U.pca = cov_pca(mash.data,5,subset=strong.subset)
# U.ed = cov_ed(mash.data, U.pca,subset=strong.subset)
# 
# # Run with data-driven covariance matrix
# now = Sys.time()
# m.ed = mash(mash.data, U.ed)
# print(Sys.time() - now)
# 
# message('Log likelihood (data-driven covariance matrix): ',format(get_loglik(m.ed),digits=10))
# 
# # Finally, run with both data-driven and canonical covariance matrices
# now = Sys.time()
# m = mash(mash.data, c(U.c,U.ed))
# print(Sys.time() - now)
# 
# message('Log likelihood (canonical + data-driven covariance matrices): ',format(get_loglik(m),digits=10))
# 
# # Get number of significant genes in each region
# apply(get_lfsr(m),2,function(x) sum(x < fsr.cutoff))

saveRDS(m,file=paste0('checkpoints/',dataset,'_model_',model.method,'_mashr.rds'))