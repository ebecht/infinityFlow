#' Wrapper to perform and export UMAP dimensionality reduction on the backbone
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param extra_args_UMAP Named list of arguments to pass to uwot:umap. Defaults to list(n_neighbors=15L,min_dist=0.2,metric="euclidean",verbose=verbose,n_epochs=1000L)
perform_UMAP_dimensionality_reduction=function(
                                               paths,
                                               extra_args_UMAP,
                                               chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                                               preds=readRDS(file.path(paths["rds"],"predictions.Rds")),
                                               verbose=TRUE
                                               )
{
    if(verbose){
        message("Performing dimensionality reduction")
    }

    umap=do.call(uwot::umap,c(list(X=preds[[1]][,chans]),extra_args_UMAP))
    colnames(umap)=c("UMAP1","UMAP2")

    saveRDS(umap,file=file.path(paths["rds"],"umap.Rds"))
    invisible()
}
