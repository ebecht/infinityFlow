#' Plot the UMAP embedding of the backbone with color-coded intensities for backbone and predicted measurements
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param chop_quantiles removes the top and bottom \emph{chop_quantiles} of the intensity scale for each marker when mapping intensities to colors.

plot_results=function(
                      paths,
                      chop_quantiles=0.005,
                      chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                      events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
                      preds=readRDS(file.path(paths["rds"],"predictions_cbound.Rds")),
                      sampling=readRDS(file.path(paths["rds"],"sampling_preds.Rds")),
                      prediction_colnames=readRDS(file.path(paths["rds"],"prediction_colnames.Rds")),
                      a=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE),
                      verbose=TRUE,
                      file_name=file.path(paths["output"],"umap_plot_annotated.pdf"),
                      global_palette=FALSE
                      ){
    if(verbose){
        message("Plotting")
    }
    
    ## chop_quantiles=0.005;
    ## chans=readRDS(file.path(paths["rds"],"chans.Rds"));
    ## umap=readRDS(file.path(paths["rds"],"umap.Rds"));
    ## events.code=readRDS(file.path(paths["rds"],"pe.Rds"));
    ## preds=readRDS(file.path(paths["rds"],"predictions_cbound.Rds"));
    ## sampling=readRDS(file.path(paths["rds"],"sampling_preds.Rds"));
    ## a=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE)
    ## verbose=TRUE;
    ## prediction_colnames=readRDS(file.path(paths["rds"],"prediction_colnames.Rds"))
    
    a=setNames(as.character(a[,"target",]),a[,"file"])
    a[is.na(a)]=paste0("Autofluorescence",1:sum(is.na(a)))
    
    if(verbose){
        message("\tChopping off the top and bottom ",chop_quantiles," quantiles")
    }
    for(col in prediction_colnames){
        q=quantile(preds[,col],c(chop_quantiles,1-chop_quantiles))
        preds[,col][preds[,col]<=q[1]]=q[1]
        preds[,col][preds[,col]>=q[2]]=q[2]
        preds[,col]=preds[,col]
    }

    if(verbose){
        message("\tShuffling the order of cells (rows)")
    }
    scrbl=sample(1:nrow(preds))
    preds=preds[scrbl,]
    
    colnames(preds)=gsub("/","-",colnames(preds))
    channels.code=setNames(colnames(preds),colnames(preds))
    ## preds=cbind(umap,preds[,!w],preds[,sort(colnames(preds)[w])])

    if(verbose){
        message("\tProducing plot")
    }
    color_biplot_by_channels(
        preds,
        x_axis="UMAP1",
        y_axis="UMAP2",
        global_across_channels=global_palette,
        file_name=file_name,
        palette=jet.colors(100),
        pch=16,
        cex=min(1,1.1-0.15*log10(nrow(preds))),
        res=72,
        raster.height=360*4,
        raster.width=360*4
    )
}
