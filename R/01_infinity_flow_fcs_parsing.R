#' Parsing FCS files
#' @param input_events_downsampling See \code{\link{infinity_flow}}
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param extra_args_read_FCS passed to flowCore:read.FCS
#' @param name_of_PE_parameter Name of the exploratory measurement
#' @description This function reads the input FCS files, downsample to a user-defined number of events if necessary, harmonize the input data names, save a concatenated expression matrix and the corresponding vector mapping events to their file of origin.
subsample_data=function(
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
    files=list.files(paths["input"],full.names=TRUE,recursive=TRUE,pattern=".fcs")
    invisible(sapply(files,function(file){
        res=do.call(read.FCS,c(list(filename=file),extra_args_read_FCS))
        w=sort(sample(1:nrow(res),min(input_events_downsampling,nrow(res))))
        res=res[w,]
        write.FCS(res,sub(paths["input"],paths["subset"],file,fixed=TRUE))
    }))

    ## convert to .Rds
    if(verbose){
        message("\tConcatenating expression matrices")
    }
    files=list.files(paths["subset"],full.names=TRUE,recursive=FALSE,include.dirs=FALSE,pattern=".fcs")
    ns=setNames(integer(length(files)),files)
    xp=sapply(files,function(file){
        xp=do.call(read.FCS,c(list(filename=file),extra_args_read_FCS))
        annot=pData(xp@parameters)
        ns[file]<<-nrow(xp)
        xp=exprs(xp)
        targets = annot$desc
        targets[is.na(targets)] = annot$name[is.na(targets)]
        colnames(xp)[colnames(xp)!=name_of_PE_parameter]=targets[colnames(xp)!=name_of_PE_parameter]
        colnames(xp)[colnames(xp)==name_of_PE_parameter]=name_of_PE_parameter
        xp
    },simplify=FALSE)
    xp=do.call(rbind,xp)

    ## Map which events originate from which file.
    if(verbose){
        message("\tWriting to disk")
    }
    events.code=unlist(unname(sapply(names(ns),function(x){rep(tail(strsplit(x,"/")[[1]],1),ns[x])},simplify=FALSE)))
    saveRDS(xp,file=file.path(paths["rds"],"xp.Rds"))
    saveRDS(events.code,file=file.path(paths["rds"],"pe.Rds"))
    invisible()
}
