# Brain single-nucleus RNA-seq analysis pipeline

**Repository for the aging brain transcriptome project from Cayo Santiago macaques (single-nucleus RNA-seq)**

This repository contains scripts used in the analysis of age effects in single-nucleus RNA-seq data for the Cayo Santiago rhesus macaque population. Corresponding code use in the analysis of age effects in bulk-tissue RNA-seq data can be found in the respository [CayoBiobankResearchUnit/brain\_transriptome\_aging\_bulk](https://github.com/CayoBiobankResearchUnit/brain_transcriptome_aging_bulk).

Note that we ran most steps on the University of Washington ([Mox](https://wiki.cac.washington.edu/display/hyakusers/Hyak+mox+Overview)) and Arizona State University ([Agave](https://cores.research.asu.edu/research-computing/user-guide)) high-performance computing clusters. We have aimed to generalize the code here by removing system-specific references to installed software and modules. Instead, we document required software and version numbers below (excluding standard Unix programs and R). For HPC systems, the required scripts and binaries must be in the PATH. The easiest way to do this is to use an existing module or to install your own. In these cases, the modules should be loaded prior to running the appropriate code below.

As Mox and Agave use the [slurm](https://slurm.schedmd.com/documentation.html) scheduler, most code below should run on slurm systems with little or no modification. For non-slurm HPC systems, slurm scripts and environmental variables will need to be adjusted, though hopefully without too much hassle.

We ran most analysis steps using [R](https://cran.r-project.org) (v4.1). We recommend the following utility or visualization packages to extend base R's functionality for single-cell RNA-seq analysis.

| Package                                                              | Description                                        |
| -----------                                                          | -----------                                        |
| [tidyverse](https://www.tidyverse.org/)                              | utilities for data manipulation and visualization  |
| [reshape2](https://cran.r-project.org/web/packages/reshape2)         | data manipulation                                  |
| [XML](https://cran.r-project.org/web/packages/XML)                   | parsing XML files                                  |
| [ggrastr](https://cran.r-project.org/web/packages/ggrastr)           | rasterizing big data visualizations                |
| [ggrepel](https://cran.r-project.org/web/packages/ggrepel)           | avoiding overplotting of labels                    |
| [egg](https://cran.r-project.org/web/packages/egg)                   | advanced layouts for visualization                 |
| [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer) | color advice for plots                             |
| [viridis](https://cran.r-project.org/web/packages/viridis)           | perceptually uniform color gradients               |
| [doParallel](https://cran.r-project.org/web/packages/doParallel)     | support for parallel computing                     |
| [future](https://cran.r-project.org/web/packages/future)             | support for parallel computing                     |
| [Matrix](https://cran.r-project.org/web/packages/Matrix)             | support for sparse matrices                        |
| [monocle3](https://cole-trapnell-lab.github.io/monocle3)             | tools for single-cell genomic analysis             |
| [Seurat](https://satijalab.org/seurat)                               | tools for single-cell genomic analysis             |

More specialized R packages are listed with their specific scripts below.

# Inputs

This pipeline picks up from the outputs of [bbi-lab/bbi-dmux](https://github.com/bbi-lab/bbi-dmux) and [bbi-lab/bbi-sci](https://github.com/bbi-lab/bbi-sci). The following files are expected:

* A monocle3 cell data set (cds) object saved as an RDS file in `checkpoints/cell_data_sets/dlpfc.rds`

* Demultiplexed, dedepulicated BAM files should be located in the `bam/` folder with the naming convention `<library ID>.bam`. Each bam file contains all mapped reads for an animal, with duplicate sci-RNA-seq barcodes and UMIs removed. BAM headers contain barcode information, which is not needed once the cds is generated.

* A library metadata file should be placed in `data/metadata_dlpfc.txt`

# Pipeline

## Compute relatedness with ANGSD and ngsRelate

* ***Required software***: [R](https://cran.r-project.org) (v3.6.3), [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) (v0.9.21), [ngsRelate](https://github.com/ANGSD/NgsRelate) (v1.0)

```
# Prepare regions and bam lists for ANGSD
sbatch scripts/angsd_setup.sh

# Calculate genotype likelihoods with ANGSD
sbatch --array=1-20 scripts/angsd_genotype.sh dlpfc

# Calculate kinship with ngsRelate
sbatch scripts/angsd_ngsrelate.sh dlpfc
```

## Visualize metadata

```
scripts/visualize_metadata.R dlpfc
```

## Find multiplets

A collision rate (`0.045`) is specified below, estimated from a mouse-human barnyard experiment included in the same sci-RNA-seq3 run.

```
scripts/monocle_find_doublets.R dlpfc 0.045
```

* After running Scrublet, visualize Scrublet kNN score distributions and manually assign thresholds. Save thresholds as a data frame with columns `sample` (sample ID) and `threshold` to the RDS file: `checkpoints/doublet_thresholds_dlpfc.rds`

## Preprocess cell data set

* ***Required input***: `checkpoints/doublet_thresholds_dlpfc.rds` (file with data frame manually saved in the previous step).

```
scripts/monocle_preprocess.R dlpfc
```

## Visualize initial UMAP and clustering results
```
scripts/visualize_clusters_qc.R dlpfc
```

* After clustering, view Scublet kNN scores per cluster and save a vector of failed clusters (clusters with noticeably high score distributions) to the RDS file: `checkpoints/doublet_clusters_dlpfc.rds`

## Remove clusters with noticeably high multiplet scores

* ***Required input***: `checkpoints/doublet_clusters_dlpfc.rds` (file with vector manually saved in the previous step).

```
scripts/manual_qc_remove_doublets.R dlpfc
```

## Visualize marker gene expression

```
scripts/visualize_marker_genes.R dlpfc darmanis
scripts/visualize_marker_genes.R dlpfc zhu
```

## Recluster data

* ***Required input***: `checkpoints/lumped_clusters_dlpfc.rds` (file with vector of cluster IDs to recluster with a lower *k* setting, resulting in higher granularity).

```
scripts/monocle_recluster.R dlpfc
```

## Visualize clusters individually to facilitate reclassification

```
scripts/visualize_cluster_membership.R dlpfc
```

* After visualizing, assign clusters to cell types in a data frame with columns `clusters` (comma-separated cluster IDs) and `cell_type` to the RDS file: `checkpoints/cell_assignments_dlpfc.rds`

## Assign cell types

* ***Required input***: `checkpoints/cell_assignments_dlpfc.rds` (file with data frame manually saved in the previous step).

```
scripts/manual_reclassify.R dlpfc
```

## Visualize final cell data set
```
scripts/visualize_clusters.R dlpfc
```

## Integrate with Allen Brain Map mouse data

* ***Key libraries***: [biomaRt](https://doi.org/doi:10.18129/B9.bioc.biomaRt), [SingleR](https://doi.org/doi:10.18129/B9.bioc.SingleR)

```
# Download and prep CDS from mouse data
scripts/prep_mouse_sc_data.R dlpfc

# Integrate data
scripts/seurat_integrate_mouse_macaque.R dlpfc

# Plot integrated data
scripts/seurat_plot_mouse_integration.R dlpfc

# Transfer labels
scripts/singleR_transfer_mouse_labels.R
```

## Integrate with Allen Brain Map human data

* ***Key libraries***: [biomaRt](https://doi.org/doi:10.18129/B9.bioc.biomaRt), [SingleR](https://doi.org/doi:10.18129/B9.bioc.SingleR)

```
# Download and prep CDS from human data
scripts/prep_human_sc_data.R dlpfc

# Integrate data
scripts/seurat_integrate_human_macaque.R dlpfc

# Plot integrated data
scripts/seurat_plot_human_integration.R dlpfc

# Transfer labels
scripts/singleR_transfer_human_labels.R
```

## Make pseudobulk libraries

* ***Key libraries***: [umap](https://cran.r-project.org/web/packages/umap), [GMPR](https://github.com/jchen1981/GMPR)

```
# Combine single-cell transcriptomes into cell-type-specific pseudobulk transcriptomes
scripts/make_pseudobulk.R dlpfc

# Visualize pseudobulk libraries
scripts/visualize_pseudobulk.R dlpfc

```

## Model pseudobulk data

* ***Key libraries***: [limma](https://doi.org/doi:10.18129/B9.bioc.limma), [GMPR](https://github.com/jchen1981/GMPR), [EMMREML](https://cran.r-project.org/web/packages/EMMREML), [mashr](https://cran.r-project.org/web/packages/mashr)

```
# Write final kinship matrix
scripts/summarize_kinship.R dlpfc

# Run EMMA model
scripts/emma_model_pseudobulk.R dlpfc

# Run MASH model
scripts/mashr_model_pseudobulk.R dlpfc emma
```

## Run enrichment analysis

* ***Key libraries***: [biomaRt](https://doi.org/doi:10.18129/B9.bioc.biomaRt), [topGO](https://doi.org/doi:10.18129/B9.bioc.topGO), [GO.db](https://doi.org/doi:10.18129/B9.bioc.GO.db)

```
scripts/topgo_enrichment.R dlpfc
```

## Compare pseudobulk model results and bulk-tissue RNA-seq results

* ***Required inputs***: `checkpoints/dlpfc_bulk_reference_emma.rds` and `checkpoints/bulk_best_predictors.rds`. These files can be generated from the [CayoBiobankResearchUnit/brain\_transriptome\_aging\_bulk](https://github.com/CayoBiobankResearchUnit/brain_transcriptome_aging_bulk) pipeline and transferred with the following code (assumes that the two repositories share a parent directory):

```
dataset = 'dlpfc'
model.method='emma'
source('../brain_transcriptome_aging_bulk/scripts/_include_options.R')
emma.results = readRDS('../brain_transcriptome_aging_bulk/checkpoints/emma_results.rds')
this.results = emma.results[,paste(c('beta','bvar','pval'),predictor,sep='.'),region.levels[match(dataset,tolower(region.levels))]]
source('scripts/_include_options.R')
colnames(this.results) = paste(c('beta','bvar','pval'),predictor,sep='.')
saveRDS(this.results,paste0('checkpoints/',dataset,'_bulk_reference_',model.method,'.rds'))

this.predictors = readRDS('../brain_transcriptome_aging_bulk/checkpoints/predictions_best_predictors.rds')
saveRDS(this.predictors,'checkpoints/bulk_best_predictors.rds')
```

Once the above files are generated, the following script summarizes and visualizes all model results:

```
scripts/visualize_model_results.R dlpfc
```
