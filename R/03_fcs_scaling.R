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
        ## scaling_parameters <- matrix(ncol = ncol(xp), nrow = 2, dimnames = list(c("center", "scale"), colnames(xp)))
        ## scaling_parameters["center", ] = 0
        ## scaling_parameters["scale", ] = 1
        for(chan in chans){
            xp[, chan] <- scale_function(xp[, chan])
            ## scaling_parameters["center", chan] = attributes(scaled_data)$'scaled:center'
            ## scaling_parameters["scale", chan] = attributes(scaled_data)$'scaled:scale'
        }
        h5write(obj = xp, file = paths["h5"], name = paste0("/input/expression_transformed_scaled/", i))
        ## h5writeAttribute(attr = scaling_parameters, name = "scaling_parameters", h5obj = paths["h5"], h5loc = paste0("/input/expression_transformed_scaled/", i))
    }

    invisible()
}
## h5writeAttribute(attr = colnames(xp), h5obj = paths["h5"], name = "colnames", h5loc = paste0("input/expression/", i))
