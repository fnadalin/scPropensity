
#' PrepareSeuratObject
#'
#' Compute highly variable genes, PCA and k-nn graph 
#' @param object The input Seurat object
#' @param cells The cell ID to subset the object on (default: use all cells)
#' @param selection.method The method to use for computing HVGs (default: use the already computed hvgs)
#' @param reduction The pre-computed dimensional reduction to use (default: run PCA)
#' @param dims The number of top dimensions to use (default: 10)
#' @param k The number of neighbours in the k-nn graph (if not already computed)
#' @return The object with updated hvg, PC, and k-nn
#' @export
PrepareSeuratObject <- function(object, cells = NULL, reduction = "", dims = 10, k = 30, ...) {

    if (length(cells) > 0) {
        object <- subset(object, subset = CELL_SUBSET) 
    }
    
    assay <- object@active.assay
    nn_graph <- paste0(assay, "_nn")
    
    if (reduction != "") {
        if (!(reduction %in% Reductions(object))) {
            stop("The object does not contain the dimensional reduction", call. = FALSE)
        }
    } else {
        reduction = "pca"
        var_features <- object@assays[[assay]]@var.features
        if (length(var_features) == 0) {
            object <- FindVariableFeatures(object, ...)
        }
        object <- RunPCA(object, ...)
    }
    object <- FindNeighbors(object, k.param = k, reduction = reduction, dims = 1:dims, ...) 
    
    return(object)
}

#' ObservedLabelPairs
#'
#' Compute the matrix containing the observed count of paired labels
#' @param object The input Seurat object
#' @param slot The slot in object@meta.data containing the labels assigned to the cells
#' @param label.list external list of labels to consider (default: all labels in the metadata slot)
#' @return The observed matrix
#' @export
ObservedLabelPairs <- function(object, slot, label.list = NULL) {

    assay <- object@active.assay
    nn_graph <- paste0(assay, "_nn")
    if(!(nn_graph %in% Graphs(object))) {
        stop("The object does not contain the graph", call. = FALSE)
    }
    
    M <- object@graphs[[nn_graph]]
    class <- object@meta.data[[slot]]
    
    if (length(label.list) == 0) {
        label.list = unique(class)
    }
    idx <- which(label.list %in% class)
    label.list <- label.list[idx]
    label.idx <- match(class, label.list, nomatch = 0)
    
    # compute the number of direct edges between all possible label pairs (non-symmetric matrix)
    C_obs <- matrix(0, nrow = length(label.list), ncol = length(label.list))
    colnames(C_obs) <- rownames(C_obs) <- label.list
    n <- 0
    y <- 0
    for (x in M@i) {
        n <- n+1
        i <- label.idx[x+1]
        if (i > 0) {
            while (y+1 <= length(M@p) & M@p[y+1] < n) {
                y <- y+1
            }
            if (y != x+1) { # exclude loops
                j <- label.idx[y]
                if (j > 0) {
                    C_obs[i,j] <- C_obs[i,j]+1
                }
            }
        }
    }
        
    return(C_obs)
}

#' ExpectedLabelPairs
#'
#' Compute the matrix containing the expected count of paired labels, given the number of label occurrences and the number of the neighbours
#' @param object The input Seurat object
#' @param slot The slot in object@meta.data containing the labels assigned to the cells
#' @param label.list external list of labels to consider (default: all labels in the metadata slot)
# '@param reduction dimensional reduction (default: "pca")
#' @return The expected matrix
#' @export
ExpectedLabelPairs <- function(object, slot, label.list = NULL, reduction = "pca") {

    assay <- object@active.assay
    nn_graph <- paste0(assay, "_nn")
    if(!(nn_graph %in% Graphs(object))) {
        stop("The object does not contain the graph", call. = FALSE)
    }
    
    M <- object@graphs[[nn_graph]]
    class <- object@meta.data[[slot]]
    
    field <- paste("FindNeighbors", assay, reduction, sep = ".")
    if (!(field %in% names(object@commands))) {
        stop(paste0("The object does not contain the graph information (", field, ")"), call. = FALSE)
    }
    K <- object@commands[[field]]$k.param
    
    if (length(label.list) == 0) {
        label.list = unique(class)
    }
    idx <- which(label.list %in% class)
    label.list <- label.list[idx]
    label.idx <- match(class, label.list, nomatch = 0)
    
    # compute the expected number of direct edges 
    C_exp <- matrix(0, nrow = length(label.list), ncol = length(label.list))
    colnames(C_exp) <- rownames(C_exp) <- label.list
    N <- length(class)
    for (i in 1:nrow(C_exp)) {
        n_i <- sum(class == label.list[i])
        for (j in 1:ncol(C_exp)) {
            n_j <- sum(class == label.list[j]) - (i == j) # exclude node i from the total count
            C_exp[i,j] <- n_i*K*n_j / (N-1) 
        }
    }
        
    return(C_exp)
}

#' Propensity
#'
#' Compute the pair propensity matrix
#' @param object The input Seurat object
#' @param slot The slot in object@meta.data containing the labels assigned to the cells
#' @param label.list external list of labels to consider (default: all labels in the metadata slot)
#' @param reduction dimensional reduction (default: "pca")
#' @return The propensity matrix
#' @export
Propensity <- function(object, slot, label.list = NULL, reduction = "pca") {

    C_obs <- ObservedLabelPairs(object, slot, label.list)
    C_exp <- ExpectedLabelPairs(object, slot, label.list, reduction)
    
    C_prop <- matrix(0, nrow = nrow(C_obs), ncol = ncol(C_obs))
    colnames(C_prop) <- rownames(C_prop) <- rownames(C_obs)

    for (i in 1:nrow(C_exp)) {
        for (j in 1:ncol(C_exp)) {
            C_prop[i,j] <- C_obs[i,j] / C_exp[i,j]
        }
    }
    
    return(C_prop)
}

#' Symmetrise
#'
#' Compute the symmetric version of the pair propensity matrix (L2 norm)
#' @param C The propensity matrix
#' @return The simmetric propensity matrix
#' @export
Symmetrise <- function(C) {

    C_symm <- t(C)*C
    C_symm <- apply(C_symm, 2, sqrt)
    
    return(C_symm)
}

#' WriteMatrix
#'
#' Write output propensity matrices (before and after L2 normalisation
#' @param C The propensity matrix
#' @param filename tsv filename where to print the matrix
#' @return The simmetric propensity matrix
#' @export
WriteMatrix <- function(C, filename) {

    outdir <- dirname(filename)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    write.table(C, file = filename, quote = FALSE, sep = "\t")
}

#' scPropensity
#'
#' Wrapper function for pair propensity matrix calculation
#' @param object The input Seurat object
#' @param slot The slot in object@meta.data containing the labels assigned to the cells
#' @param label.list external list of labels to consider (default: all labels in the metadata slot)
#' @param outdir Folder where to print the propensity matrices
#' @return NA
#' @export
scPropensity <- function(object.file, slot, outdir, label.file = "", ...) {

    cat("Read object...")
    object <- readRDS(object.file)
    label.list = NULL
    if (label.file != "") {
        label.list <- read.table(label.file, sep = "\t")[,1]
    }
    cat("done.\n")
    
    cat("Prepare Seurat Object...")
    object <- PrepareSeuratObject(object, ...)
    cat("done.\n")
    
    cat("Compute propensity and symmetrise...")
    C_prop <- Propensity(object, slot, label.list = label.list)
    C_prop_symm <- Symmetrise(C_prop)
    cat("done.\n")
    
    cat("Write to file...")
    WriteMatrix(C_prop, file.path(outdir, "pair_propensity.tsv"))
    WriteMatrix(C_prop_symm, file.path(outdir, "pair_propensity_symm.tsv"))
    cat("done.\n")
}



