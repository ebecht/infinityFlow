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
                     umap=readRDS(file.path(paths["rds"],"umap.Rds")),
                     events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
                     preds=readRDS(file.path(paths["rds"],"svms_predictions.Rds")),
                     sampling=readRDS(file.path(paths["rds"],"sampling_preds.Rds")),
                     a=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE)
                     ){
    ## env=environment()
    ## sapply(
    ##     c("transforms_chan","transforms_pe","xp","chans","umap"),
    ##     function(object){
    ##         assign(object,value=readRDS(file.path(paths["rds"],paste0(object,".Rds"))),envir=env)
    ##         invisible()
    ##     }
    ## )
    ## events.code=readRDS(file.path(paths["rds"],"pe.Rds"))
    ## preds=readRDS(file.path(paths["rds"],"svms_predictions.Rds"))
    ## sampling=readRDS(file.path(paths["rds"],"sampling_preds.Rds"))

    a[,"target"]=make.unique(a[,"target"])
    a=setNames(as.character(a[,"target",]),a[,"file"])
    a[is.na(a)]="Autofluorescence"
    if(any(a=="Autofluorescence")){
        a[a=="Autofluorescence"]=paste0("Autofluorescence",1:sum(a=="Autofluorescence"))
    }

    for(pe in names(transforms_pe)){
        ilgcl=inverseLogicleTransform(trans=transforms_pe[[pe]])
        preds[,pe]=ilgcl(preds[,pe])
    }
    colnames(preds)[colnames(preds)%in%names(a)]=paste0(a[colnames(preds)[colnames(preds)%in%names(a)]],".predicted")

    unique_pes=unique(events.code)
    PE_id=sapply(events.code[sampling],match,table=unique_pes)

    write.csv(file=file.path(paths["output"],"Exploratory_Ab_ID_table.csv"),data.frame(file=unique_pes,target=a[unique_pes],PE_id=unique(PE_id)),row.names=FALSE)
    
    ## To make UMAP easier to plot for FlowJo users
    umap=apply(umap,2,function(x){
        (x-min(x))/(max(x)-min(x))*10000
    })

    preds=cbind(xp[sampling,],preds[,!colnames(preds)%in%colnames(xp)],Exploratory_Ab_ID=PE_id,umap)
    colnames(preds)=make.unique(colnames(preds))

    if(CSV_export){
        write.csv(preds,file=file.path(paths["output"],"results.csv"),row.names=FALSE)
    }
    
    if(any(FCS_export=="split")){
        preds_tmp=lapply(split(as.data.frame(preds),events.code[sampling]),as.matrix)
        dir.create(file.path(paths["output"],"FCS","split"),showWarnings=FALSE)
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
        dir.create(file.path(paths["output"],"FCS","concatenated"),showWarnings=FALSE)
        FCS=flowFrame(preds)
        FCS@parameters$desc=as.character(FCS@parameters$desc)
        FCS@parameters$name=as.character(FCS@parameters$name)
        FCS=generate_description(FCS)
        invisible(write.FCS(FCS,file=file.path(paths["output"],"FCS","concatenated","concatenated_results.fcs")))
    }

    preds
}
