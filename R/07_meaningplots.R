#' Plot the UMAP embedding of the backbone with color-coded intensities for backbone and predicted measurements
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param chop_quantiles removes the top and bottom \emph{chop_quantiles} of the intensity scale for each marker when mapping intensities to colors.

plot_results=function(
                      paths,
                      chop_quantiles=0.005,
                      chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                      umap=readRDS(file.path(paths["rds"],"umap.Rds")),
                      events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
                      preds=readRDS(file.path(paths["rds"],"predictions.Rds")),
                      sampling=readRDS(file.path(paths["rds"],"sampling_preds.Rds")),
                      a=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE)
                      ){
    a=setNames(as.character(a[,"target",]),a[,"file"])
    a[is.na(a)]=paste0("Autofluorescence",1:sum(is.na(a)))

    scrbl=sample(1:nrow(preds))
    umap=umap[scrbl,]

    w=colnames(preds)%in%names(a)
    colnames(preds)[w]=paste0(a[colnames(preds)[w]],".predicted")
    colnames(preds)=gsub("/","_",colnames(preds))
    colnames(preds)=make.unique(colnames(preds))

    for(col in colnames(preds)){
        q=quantile(preds[,col],c(chop_quantiles,1-chop_quantiles))
        preds[,col][preds[,col]<=q[1]]=q[1]
        preds[,col][preds[,col]>=q[2]]=q[2]
        preds[,col]=preds[,col][scrbl]
    }

    channels.code=setNames(colnames(preds),colnames(preds))
    preds=cbind(umap,preds[,!w],preds[,sort(colnames(preds)[w])])

    color_biplot_by_channels(
        preds,
        x_axis="UMAP1",
        y_axis="UMAP2",
        global_across_channels=FALSE,
        file_name=file.path(paths["output"],"umap_plot_annotated.pdf"),
        palette=jet.colors(100),
        pch=16,
        cex=min(1,1.1-0.15*log10(nrow(preds))),
        res=72,
        raster.height=360*4,
        raster.width=360*4
    )
}
