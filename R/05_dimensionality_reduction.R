#' Wrapper to perform and export UMAP dimensionality reduction on the backbone
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param extra_args_UMAP Named list of arguments to pass to uwot:umap. Defaults to list(n_neighbors=15L,min_dist=0.2,metric="euclidean",verbose=verbose,n_epochs=1000L)
#' @param chans vector of backbone channels' names
#' @param preds matrix of imputed data
#' @param verbose Verbosity
#' @noRd
#' @importFrom RcppHNSW hnsw_knn
perform_UMAP_dimensionality_reduction <- function(
                                               paths,
                                               extra_args_UMAP,
                                               chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                                               verbose=TRUE,
                                               annot=read.table(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE)
                                               )
{
    if(verbose){
        message("Performing dimensionality reduction")
    }
    xp_backbone <- vector("list", nrow(annot))
    for(i in seq_len(nrow(annot))){
        index <- h5read(file = paths["h5"], name = paste0("sampling/predictions/", i)) == 1L
            
        xp_backbone[[i]] <- h5read(file = paths["h5"],name = paste0("/input/expression_transformed_scaled/", i))[index, ]
        colnames(xp_backbone[[i]]) <- h5readAttributes(file = paths["h5"], name = paste0("/input/expression/", i))$colnames
    }
    ns <- vapply(xp_backbone, nrow, integer(1))
    xp_backbone <- do.call(rbind, xp_backbone)
    
    umap <- do.call(uwot::umap,c(list(X=xp_backbone[,chans]),extra_args_UMAP))
    colnames(umap) <- c("UMAP1","UMAP2")
    saveRDS(umap,file=file.path(paths["rds"],"umap.Rds"))

    umap <- split_matrix(umap, rep(seq_along(ns), ns))
    for(i in seq_along(umap)){
        h5write(umap[[i]], file = paths["h5"], name = paste0("/umap/backbone/", i))
    }
    
    invisible()
}
