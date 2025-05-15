#' Parsing FCS files
#' @param input_events_downsampling See \code{\link{infinity_flow}}
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param extra_args_read_FCS passed to flowCore:read.FCS
#' @param name_of_PE_parameter Name of the exploratory measurement
#' @description This function reads the input FCS files, downsample to a user-defined number of events if necessary, harmonize the input data names, save a concatenated expression matrix and the corresponding vector mapping events to their file of origin.
#' @param verbose Verbosity
#' @noRd
#' @importFrom flowCore read.FCS
#' @importFrom flowCore read.FCSheader
#' @importFrom Biobase pData exprs
#' @importFrom rhdf5 h5writeAttribute
subsample_data <- function(
                           input_events_downsampling,
                           paths,
                           extra_args_read_FCS,
                           name_of_PE_parameter,
                           verbose=TRUE
                           ){
    ## Subsampling
    if(verbose){
        message("Parsing and subsampling input data")
        message("\tDownsampling to ",input_events_downsampling," events per input file")
    }
    files <- list.files(paths["input"],full.names=TRUE,recursive=TRUE,pattern=".fcs")
    invisible(
        lapply(
            files,
            function(file){
                res <- do.call(read.FCS,c(list(filename=file),extra_args_read_FCS))
                w <- sort(sample(seq_len(nrow(res)),min(input_events_downsampling,nrow(res))))
                res <- res[w,]
                write.FCS(res,file.path(paths["subset"],basename(file)))
            }
        )
    )
    
    ## convert to .h5
    
    if(verbose){
        message("\tSaving expression matrices to disk")
    }
    files <- list.files(paths["subset"],full.names=TRUE,recursive=FALSE,include.dirs=FALSE,pattern=".fcs")
    ns = integer(length(files))
    
    for(i in seq_along(files)){
        file = files[i]
        xp <- do.call(read.FCS,c(list(filename=file),extra_args_read_FCS))
        annot <- pData(xp@parameters)
        xp <- exprs(xp)
        targets <- annot$desc
        targets[is.na(targets)] <- annot$name[is.na(targets)]
        colnames(xp)[colnames(xp)!=name_of_PE_parameter] <- targets[colnames(xp)!=name_of_PE_parameter]
        colnames(xp)[colnames(xp)==name_of_PE_parameter] <- name_of_PE_parameter
        
        h5write(obj = xp, file = paths["h5"], name = paste0("input/expression/", i))
        h5writeAttribute(attr = colnames(xp), h5obj = paths["h5"], name = "colnames", h5loc = paste0("input/expression/", i))

        ns[i] = nrow(xp)
    }
    
    ## Map which events originate from which file.
    if(verbose){
        message("\tWriting to disk")
    }
    events.code <- unlist(
        lapply(
            names(ns),
            function(x){
                rep(tail(strsplit(x,"/")[[1]],1), ns[x])
            }
        )
    )

    saveRDS(events.code,file=file.path(paths["rds"],"pe.Rds"))
    invisible()
}
