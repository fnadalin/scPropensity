To compute the pair propensity on a Seurat object, run:

```
scPropensity(object.file = "object.Rds", slot = "sample.name", outdir = "dir")
```
The function builds a k-nn graph, computes the pair propensities of the unique labels in object@meta.data$sample.name, and creates a folder "dir" containing output files. 
