# scPropensity


In single-cell sequencing datasets, cells are often assigned meta data information (_e.g._, batch, sample ID, condition...). 
scPropensity provides a score to evaluate whether there is any relationship between meta data assigment and molecular similarity across cells. It only requires a cell-cell similarity matrix and a label assignment, thus it is cluster-free.

A notable application of scPropensity concerns the analysis of cancer clones: cells are assigned a clone identity, which may be inferred by copy-number profile, SNV, or the readout of lineage-tracing assays, and we may ask whether two clones have a more or less similar transcriptional profile than expected by chance. scPropensity computes a clone-clone transcriptional similarity, which can then be used to classify clones into distinct transcriptional groups. This has implications on the evaluation of the transcriptional heterogeneity of cancer clones, which is linked to their differentiation potential.

## Description

scPropensity is inspired from a concept in structural bioinformatics called [statistical potential](https://en.wikipedia.org/wiki/Statistical_potential), which is useful to evaluate the likelihood a protein complex model via a pseudo-energy function computed from a database of experimental protein structures. 

In particular, pair propensity scores are derived from amino acid pairings at the protein-protein interface and are defined as _p_(_x_,_y_) = _F_(_x_,_y_)/_G_(_x_,_y_), where _F_(_x_,_y_) is the observed frequency of pair (_x_,_y_) and _G_(_x_,_y_) is the expected frequency of pair (_x_,_y_). Depending on the value of _p_(_x_,_y_), _x_ and _y_ are more (> 1), less (< 1) or equally (= 1) likely to be in contact with each other than expected by chance.

Here, _x_ and _y_ are cell labels.
A cell-cell similarity measure is derived from the assay (gene expression, chromatin accessiblity state...) and is used to build a _k_-nn graph, where nodes are cells and a directed edge connects cell _i_ with cell _j_ if and only if _j_ is one of the closest _k_ cells to _i_ according to this measure. 

Hence, _F_(_x_,_y_) is the number of edges (_i_,_j_) in the _k_-nn graph such that _i_ is labelled with _x_ and _j_ is labelled with _j_; 
_G_(_x_,_y_) is the expected number of edges labelled with (_x_,_y_) given the neighbourhood size _k_ and the number of cells labelled with _x_ and _y_, respectively  (see [[Nadalin _et al._]](https://doi.org/10.1101/2023.06.28.546923) for details).
Therefore, _p_(_x_,_y_) tells whether cells labelled with _x_ tend to be more (> 1), less (< 1) or equally (= 1) similar to the cells labelled with _y_ than expected by chance. 


## Requirements

+ R v4.0.3
+ Seurat v4.0.5

## Instructions

scPropensity is implemented in R, it takes as input a Seurat object and a meta data field ID.
To compute the pair propensity score on `object.Rds` with respect to `sample.name`, run:

```
scPropensity(object.file = "object.Rds", slot = "sample.name", outdir = "dir")
```

The above function builds a _k_-nn graph, computes the pair propensities of the labels in ``object@meta.data$sample.name` and creates the output folder `dir`. It contains a _n_ x _n_ matrix _M_, where _n_ is the number of distict values in sample.name`, and _M_\[_x_,_y_\] is the log pair propensity of (_x_,_y_). It also contains the L2-normalised version of _M_.

## Citing

If you find this software useful, please cite:

Nadalin _et al._ Multi-omic lineage tracing predicts the transcriptional, epigenetic and genetic determinants of cancer evolution. biorxiv. https://doi.org/10.1101/2023.06.28.546923