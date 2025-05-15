#' Scaling of the backbone measurements
#' @param yvar name of the exploratory measurement
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param scale_function Scaling function to apply. Should apply to a matrix (events x backbone measurements) and return a matrix of similar size. Defaults to scale(,center=TRUE,scale=FALSE)
#' @param xp Logicle-transformed backbone expression matrix
#' @param chans vector of backbone channels' names
#' @param events.code vector of length nrow(xp) specifying from which well each event originates
#' @param verbose Verbosity
#' @importFrom rhdf5 h5readAttributes
#' @noRd
standardize_backbone_data_across_wells <- function(
                                                   yvar,
                                                   paths,
                                                   scale_function=function(x){scale(x,center=TRUE,scale=TRUE)},
                                                   xp=readRDS(file.path(paths["rds"],"xp_transformed.Rds")),
                                                   chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                                                   annot=read.table(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE),
                                                   verbose=TRUE
                                                   ){
    if(verbose){
        message("Harmonizing backbone data")
        message("\tScaling expression matrices")
    }
    
    for(i in seq_len(nrow(annot))){
        xp = h5read(file = paths["h5"], name = paste0("/input/expression_transformed/", i))
        colnames(xp) = h5readAttributes(file = paths["h5"], name = paste0("/input/expression/", i))$colnames
        for(chan in c(chans, yvar)){
            xp[, chan] = scale_function(xp[, chan])
        }
        h5write(obj = xp, file = paths["h5"], name = paste0("/input/expression_transformed_scaled/", i))
    }

    invisible()
}
