#' Set up h5 dataset
#' @param input_events_downsampling See \code{\link{infinity_flow}}
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param extra_args_read_FCS passed to flowCore:read.FCS
#' @param verbose Verbosity
#' @noRd
#' @importFrom flowCore read.FCS
#' @importFrom flowCore read.FCSheader
#' @importFrom Biobase pData exprs
#' @importFrom rhdf5 h5createFile
#' @importFrom rhdf5 h5createGroup
#' @importFrom rhdf5 h5createDataset
#' @importFrom rhdf5 h5write

## Creating HDF5 dataset with fixed-sized matrices
create_h5_file <- function(
                           paths,
                           input_events_downsampling,
                           prediction_events_downsampling,
                           extra_args_read_FCS,
                           verbose=TRUE,
                           annot=read.table(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE)
                           ){
    message("Creating HDF5 dataset")
    message("\t", paths["h5"])

    if(file.exists(paths["h5"])){
        file.remove(paths["h5"])
    }
    h5createFile(paths["h5"])

    h5createGroup(paths["h5"], "input")
    h5createGroup(paths["h5"], "input/expression") ## Raw downsampled expression
    h5createGroup(paths["h5"], "input/expression_transformed") ## Transformed downsampled expression
    h5createGroup(paths["h5"], "input/expression_transformed_scaled") ## Transformed and scaled downsampled expression

    files <- file.path(paths["input"], annot$file)##list.files(paths["input"],full.names=TRUE,recursive=FALSE,include.dirs=FALSE,pattern=".fcs")
    ns <- as.integer(sapply(read.FCSheader(files), "[[", "$TOT"))
    ns[ns>input_events_downsampling]=input_events_downsampling
    nchans <- as.integer(read.FCSheader(files[1])[[1]]["$PAR"])
    ns_out <- ns - floor(ns/2)
    ns_out[ns_out > prediction_events_downsampling] <- prediction_events_downsampling
    for(i in seq_along(files)){
        suppressMessages(
            h5createDataset(
                paths["h5"],
                paste0("input/expression/", i),
                dims = c(ns[i], nchans)
            )
        )
    }
    for(i in seq_along(files)){
        suppressMessages(
            h5createDataset(
                paths["h5"],
                paste0("input/expression_transformed/", i),
                dims = c(ns[i], nchans)
            )
        )
    }
    for(i in seq_along(files)){
        suppressMessages(
            h5createDataset(
                paths["h5"],
                paste0("input/expression_transformed_scaled/", i),
                dims = c(ns[i], nchans)
            )
        )
    }

    h5createDataset(file = paths["h5"], dataset = "/dimensions", dims = c(length(files), 2))
    h5write(file = paths["h5"], name = "dimensions", obj = cbind(ns, ns_out))
    
    ## Set-up vector of booleans to select events for predictions
    h5createGroup(paths["h5"], "sampling")
    h5createGroup(paths["h5"], "sampling/fitting/")
    h5createGroup(paths["h5"], "sampling/predictions/")
    for(i in seq_along(files)){
        suppressMessages(
            h5createDataset(
                paths["h5"],
                paste0("sampling/fitting/", i),
                dims = ns[i],
                storage.mode = "integer",
                H5type = "H5T_STD_I8LE"
            )
        )
        suppressMessages(
            h5createDataset(
                paths["h5"],
                paste0("sampling/predictions/", i),
                dims = ns[i],
                storage.mode = "integer",
                H5type = "H5T_STD_I8LE"
            )
        )
    }
    
    h5createGroup(paths["h5"], "predictions")
    h5createGroup(paths["h5"], "predictions/raw") ## Raw predictions expression
    h5createGroup(paths["h5"], "predictions/background_corrected") ## Background-corrected expression
    h5createGroup(paths["h5"], "predictions/intra_well") ## Predictions for all events from an input file with the model trained on that file

    for(i in seq_along(files)){
        h5createDataset(
            paths["h5"],
            paste0("predictions/raw/", i),
            dims = c(ns_out[i], length(files)),
            chunk = c(ns_out[i], 1)
        )
        h5createDataset(
            paths["h5"],
            paste0("predictions/background_corrected/", i),
            dims = c(ns_out[i], length(files)),
            chunk = c(ns_out[i], 1)
        )
        suppressMessages(
            h5createDataset(
                paths["h5"],
                paste0("predictions/intra_well/", i),
                dims = ns[i]
            )
        )
    }

    h5createGroup(paths["h5"], "/umap/")
    h5createGroup(paths["h5"], "/umap/backbone/")
    for(i in seq_along(files)){
        suppressMessages(
            h5createDataset(
                paths["h5"],
                paste0("umap/backbone/", i),
                dims = c(ns_out[i], 2),
                chunk = c(ns_out[i], 2)
            )
        )
    }
}
