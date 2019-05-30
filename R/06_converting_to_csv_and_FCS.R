#' Exporting results
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param FCS_export if FALSE, no FCS export. if "concatenated", export one FCS file with all the data. if "split", export the data in the result folder under the subfolder FCS, with each file corresponding to a (subsampled) input file.
export_data=function(
                     paths,
                     FCS_export,
                     CSV_export,
                     chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                     transforms_chan=readRDS(file.path(paths["rds"],"transforms_chan.Rds")),
                     transforms_pe=readRDS(file.path(paths["rds"],"transforms_pe.Rds")),
                     xp=readRDS(file.path(paths["rds"],"xp.Rds")),
                     xp_scaled=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds")),
                     umap=readRDS(file.path(paths["rds"],"umap.Rds")),
                     events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
                     preds=readRDS(file.path(paths["rds"],"predictions.Rds")),
                     sampling=readRDS(file.path(paths["rds"],"sampling_preds.Rds")),
                     a=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE),
                     verbose=TRUE
                     ){
    chans=readRDS(file.path(paths["rds"],"chans.Rds"));
    transforms_chan=readRDS(file.path(paths["rds"],"transforms_chan.Rds"));
    transforms_pe=readRDS(file.path(paths["rds"],"transforms_pe.Rds"));
    xp=readRDS(file.path(paths["rds"],"xp.Rds"));
    xp_scaled=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds"));
    umap=readRDS(file.path(paths["rds"],"umap.Rds"));
    events.code=readRDS(file.path(paths["rds"],"pe.Rds"));
    preds=readRDS(file.path(paths["rds"],"predictions.Rds"));
    sampling=readRDS(file.path(paths["rds"],"sampling_preds.Rds"));
    a=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE);
    verbose=TRUE

    if(verbose){
        message("Exporting results")
    }
    
    a[,"target"]=make.unique(a[,"target"])
    a=setNames(as.character(a[,"target",]),a[,"file"])
    a[is.na(a)]="Autofluorescence"
    if(any(a=="Autofluorescence")){
        a[a=="Autofluorescence"]=paste0("Autofluorescence",1:sum(a=="Autofluorescence"))
    }

    if(verbose){
        message("\tTransforming predictions back to a linear scale")
    }

    preds_raw=preds

    preds=lapply(preds,function(x){
        for(pe in names(transforms_pe)){
            ilgcl=inverseLogicleTransform(trans=transforms_pe[[pe]])
            x[,pe]=ilgcl(x[,pe])
        }
        x
    })

    preds=lapply(
        names(preds),
        function(x){
            w=colnames(preds[[x]])[colnames(preds[[x]])%in%names(a)]
            preds[[x]]=preds[[x]][,w]
            
            colnames(preds[[x]])=paste0(a[colnames(preds[[x]])],".",x)
            preds[[x]]
        }
    )
    preds_raw=lapply(
        names(preds_raw),
        function(x){
            w=colnames(preds_raw[[x]])[colnames(preds_raw[[x]])%in%names(a)]
            preds_raw[[x]]=preds_raw[[x]][,w]
            
            colnames(preds_raw[[x]])=paste0(a[colnames(preds_raw[[x]])],".",x)
            preds_raw[[x]]
        }
    )
    
    prediction_colnames=sort(do.call(c,lapply(preds,colnames)))
    preds=do.call(cbind,preds)[,prediction_colnames]
    preds_raw=do.call(cbind,preds_raw)[,prediction_colnames]
    
    unique_pes=unique(events.code)
    PE_id=sapply(events.code[sampling],match,table=unique_pes)

    write.csv(file=file.path(paths["output"],"Exploratory_Ab_ID_table.csv"),data.frame(file=unique_pes,target=a[unique_pes],PE_id=unique(PE_id)),row.names=FALSE)
    
    ## To make UMAP easier to plot for FlowJo users
    umap=apply(umap,2,function(x){
        (x-min(x))/(max(x)-min(x))*10000
    })

    preds=cbind(xp[sampling,],preds[,!colnames(preds)%in%colnames(xp)],Exploratory_Ab_ID=PE_id,umap)
    preds_raw=cbind(xp_scaled[sampling,],preds_raw[,!colnames(preds_raw)%in%colnames(xp_scaled)],Exploratory_Ab_ID=PE_id,umap)
    
    colnames(preds)=make.unique(colnames(preds))
    colnames(preds_raw)=make.unique(colnames(preds_raw))

    if(CSV_export){
        if(verbose){
            message("\t","Exporting as CSV")
        }
        write.csv(preds,file=file.path(paths["output"],"results.csv"),row.names=FALSE)
    }
    
    if(any(FCS_export=="split")){
        if(verbose){
            message("\t","Exporting FCS files (1 per well)")
        }
        preds_tmp=lapply(split(as.data.frame(preds),events.code[sampling]),as.matrix)
        dir.create(file.path(paths["output"],"FCS","split"),showWarnings=FALSE,recursive=TRUE)
        lapply(
            names(preds_tmp),
            function(file){
                FCS=flowFrame(preds_tmp[[file]])
                FCS@parameters$desc=as.character(FCS@parameters$desc)
                FCS@parameters$name=as.character(FCS@parameters$name)
                FCS=generate_description(FCS)
                invisible(write.FCS(FCS,file=file.path(paths["output"],"FCS","split",file)))
            }
        )
    }
    if(any(FCS_export=="concatenated")){
        if(verbose){
            message("\t","Exporting concatenated FCS file")
        }
        dir.create(file.path(paths["output"],"FCS","concatenated"),showWarnings=FALSE,recursive=TRUE)
        FCS=flowFrame(preds)
        FCS@parameters$desc=as.character(FCS@parameters$desc)
        FCS@parameters$name=as.character(FCS@parameters$name)
        FCS=generate_description(FCS)
        invisible(write.FCS(FCS,file=file.path(paths["output"],"FCS","concatenated","concatenated_results.fcs")))
    }

    saveRDS(preds_raw,file.path(paths["rds"],"predictions_cbound.Rds"))
    saveRDS(prediction_colnames,file.path(paths["rds"],"prediction_colnames.Rds"))
    preds
}
