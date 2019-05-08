#' Scaling of the backbone measurements
#' @param yvar name of the exploratory measurement
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param scale_function Scaling function to apply. Should apply to a matrix (events x backbone measurements) and return a matrix of similar size. Defaults to scale(,center=TRUE,scale=FALSE)
standardize_backbone_data_across_wells=function(
                                                yvar,
                                                paths,
                                                scale_function=function(x){scale(x,center=TRUE,scale=TRUE)}
                                                ){
    env=environment()
    sapply(
        c("chans"),
        function(object){
            assign(object,value=readRDS(file.path(paths["rds"],paste0(object,".Rds"))),envir=env)
            invisible()
        }
    )
    events.code=readRDS(file.path(paths["rds"],"pe.Rds"))

    xp=readRDS(file.path(paths["rds"],"xp_transformed.Rds"))

    xp=split(as.data.frame(xp),events.code)
    xp=lapply(xp,as.matrix)
    xp=lapply(xp,function(x){
        x[,chans]=scale_function(x[,chans])
        x
    })
    xp=do.call(rbind,xp)
    saveRDS(xp,file=file.path(paths["rds"],"xp_transformed_scaled.Rds"))
    invisible()
}
