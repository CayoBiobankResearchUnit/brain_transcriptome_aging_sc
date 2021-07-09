#!/usr/bin/env Rscript

# Options

# Number of cores to use for parallel computing (set to 0 to use all detected cores)
n.cores = 4

# Option to redo steps even when checkpoints exist
ignore.checkpoints = FALSE

# Transcripts-per-million cutoff
tpm.cutoff = 10

# False sign rate cutoff (mashr)
fsr.cutoff = 0.2

# Fraction of cells to be considered universal (0.8 ensures 6 out of 7)
fraction.shared.cutoff = 0.8

# Fraction of cells to be considered unique (0.2 ensures 1 out of 7)
fraction.unique.cutoff = 0.2

# Minimum age for modeling
min.age = 0

# Drop cells if the highest 1% of p values are below a threshold.
# Theoretical expected mean p-value for the first 1% is 0.995
qq.cutoff = 0.975

# Minimum UMI to call a nucleus
umi.min = 100

# Maximum percentage mitochondrial reads
mt.max = 10

# Clustering k
cluster.k.initial = 5
cluster.k.final = 20
cluster.k.recluster = 5

# Set column name for predictor variable
predictor = 'age'

# Set label for predictor variable
predictor.label = 'Age'

# Set units for predictor
predictor.units = 'years'

# Set column name for age variable
age.variable = 'age'

# Set column name for primary batch variable
batch.variable = 'extraction_batch'

# Covariates to include in model(s)
model.covariates = c('age','group','extraction_batch')

# Define factor levels
cell.levels = c('Excitatory neurons','Inhibitory neurons','Astrocytes','Oligodendrocytes','Oligodendrocyte precursor cells','Microglia','Endothelial cells','Pericytes')

# Define short factor levels
cell.short.levels = c('Exc.','Inh.','Astr.','Olig.','OPC','Micr.','End.','Per.')

# Define factor levels for predicted cells (including unknown)
cell.prediction.levels = c('Neurons','Excitatory neurons','Inhibitory neurons','Astrocytes','Oligodendrocytes','Oligodendrocyte precursor cells','Microglia','Endothelial cells','Pericytes','Unknown')

# Define cell colors
cell.colors = c('#1b9e77','#66a61e','#d95f02','#7570b3','#e7298a','#e6ab02','#a6761d','#666666')

# Define cell colors for predicted cells (including unknown)
cell.prediction.colors = c('#666666','#1b9e77','#66a61e','#d95f02','#7570b3','#e7298a','#e6ab02','#a6761d','#666666','#eeeeee')

# Rows for facets (cells)
cell.rows = 2

# Rows for facets (samples)
sample.rows = 4

# Sex colors (all females in this dataset, but this will keep visuals in line with other projects)
sex.colors = c('#4daf4a','#984ea3')

# Model methods
model.method.labels = c('emma' = 'EMMA', 'pqlseq' = 'PQLseq')

# Set random seed
seed = 42

# Set number of chunks (for assembling a matrix in pieces
n.chunks = 100

# Set number of subsets (with selection of mouse cells per set)
n.subsets = 5

get.ks.pval = function (object) 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(ks.test(x.a, seq_len(N)[-x.a], alternative = "greater")$p.value)
}

get.ks.score = function (object) 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(ks.test(x.a, seq_len(N)[-x.a], alternative = "greater")$statistic)
}

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
#                             End configurations
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

# If not set, sets n.cores to detected cores
if (!n.cores) n.cores = ifelse('future' %in% rownames(installed.packages()),future::availableCores(methods='mc.cores') + 1,parallel::detectCores(logical=FALSE))
