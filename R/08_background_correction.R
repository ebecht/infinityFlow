correct_background=function(
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
                            sampling=readRDS(file.path(paths["rds"],"sampling_preds.Rds")),
                            preds=readRDS(file.path(paths["rds"],"predictions.Rds")),
                            prediction_colnames=readRDS(file.path(paths["rds"],"prediction_colnames.Rds")),
                            a=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE),
                            verbose=TRUE
                            ){

    if(verbose){
        message("Background correcting")
    }

    a=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE)
    rownames(a)=a[,"file"]
    a=a[intersect(rownames(a),colnames(preds[[1]])),]
    
    preds_raw=preds
    
    preds_rawbgc=lapply(
        preds_raw,
        function(x){
            x[,a$file]
        }
    )
    
    for(i in seq_along(preds_rawbgc)){
        for(file in colnames(preds_rawbgc[[i]])){
            iso=a[file,"isotype"][1]
            x=preds_raw[[i]][,iso]
            y=preds_raw[[i]][,file]
            lm=lm(y~x)$coefficients
            intercept=lm[1]
            slope=lm[2]
            orthogonal_residuals=(-slope*x+y-intercept)/sqrt(slope^2+1)
            preds_rawbgc[[i]][,file]=orthogonal_residuals
        }
    }

    saveRDS(preds_rawbgc,file=file.path(paths["rds"],"predictions_backgroundcorrected.Rds"))
    
    ## ##################
    ## Exporting Phenograph / BGC / UMAPs
    ## ##################

    if(verbose){
        message("\tTransforming background-corrected predictions. (Use logarithm to visualize)")
    }
    ## Preparing data for display using log transform
    preds_rawbgc_linear=lapply(
        preds_rawbgc,function(x){
            apply(x,2,function(x){
                10^(x-min(x))
            })
        }
    )
    
    ## Adding uncorrected data (isotypes and autofluorescence)
    ## I don't think this does anything anymore
    for(i in seq_along(preds_rawbgc_linear)){
        preds_rawbgc_linear[[i]]=cbind(
            preds_rawbgc_linear[[i]],
            preds_raw[[i]][,setdiff(colnames(preds[[i]]),c(colnames(preds_rawbgc_linear[[i]]),chans)),drop=FALSE]
        )
        preds_rawbgc[[i]]=cbind(
            preds_rawbgc[[i]],
            preds_raw[[i]][,setdiff(colnames(preds[[i]]),c(colnames(preds_rawbgc[[i]]),chans)),drop=FALSE]
        )
    }
    
    preds_rawbgc_linear=lapply(
        preds_rawbgc_linear,
        function(x){
            x[,rownames(a)]
        }
    )
    preds_rawbgc=lapply(
        preds_rawbgc,
        function(x){
            x[,rownames(a)]
        }
    )
    for(x in names(preds_rawbgc_linear)){
        colnames(preds_rawbgc_linear[[x]])=paste0(a[colnames(preds_rawbgc_linear[[x]]),"target"],".",x,"_bgc")
        colnames(preds_rawbgc[[x]])=paste0(a[colnames(preds_rawbgc[[x]]),"target"],".",x,"_bgc")
        preds_rawbgc_linear[[x]]=preds_rawbgc_linear[[x]][,order(colnames(preds_rawbgc_linear[[x]]))]
        preds_rawbgc[[x]]=preds_rawbgc[[x]][,order(colnames(preds_rawbgc[[x]]))]
    }

    preds_rawbgc_linear=do.call(cbind,preds_rawbgc_linear)
    preds_rawbgc=do.call(cbind,preds_rawbgc)

    preds_rawbgc=preds_rawbgc[,sort(colnames(preds_rawbgc))]
    preds_rawbgc_linear=preds_rawbgc_linear[,sort(colnames(preds_rawbgc_linear))]
    
    preds_rawbgc_linear=cbind(xp[sampling,],preds_rawbgc_linear,minmax_scale(umap))
    preds_rawbgc=cbind(xp_scaled[sampling,],preds_rawbgc,minmax_scale(umap))
    
    unique_pes=unique(events.code)
    PE_id=sapply(events.code[sampling],match,table=unique_pes)

    preds_rawbgc_linear=cbind(preds_rawbgc_linear,PE_id=PE_id)
    preds_rawbgc=cbind(preds_rawbgc,PE_id=PE_id)

    if(CSV_export){
        if(verbose){
            message("\t","Exporting as CSV")
        }
        dir.create(file.path(paths["output"],"FCS_background_corrected/"),showWarnings=FALSE,recursive=TRUE)
        write.csv(as.data.frame(preds_rawbgc_linear),row.names=FALSE,file=file.path(paths["output"],"predicted_data_background_corrected.csv"))
    }

    if(any(FCS_export=="split")){
        if(verbose){
            message("\t","Exporting FCS files (1 per well)")
        }
        dir.create(file.path(paths["output"],"FCS_background_corrected/","split"),recursive=TRUE,showWarnings=FALSE)
        FCS_list=split_matrix(preds_rawbgc_linear,events.code[sampling])
        FCS_list=lapply(FCS_list,function(x){
            x=flowFrame(x)
            x@parameters$desc=as.character(x@parameters$desc)
            x@parameters$name=as.character(x@parameters$name)
            generate_description(x)
        })
        lapply(
            names(FCS_list),
            function(x){
                invisible(write.FCS(FCS_list[[x]],file=file.path(paths["output"],"FCS_background_corrected/","split",paste0(sub(".fcs","",x),"_target_",gsub("/","-",a[x,"target"]),".fcs"))))
            }
        )   
    }

    if(any(FCS_export=="concatenated")){
        if(verbose){
            message("\t","Exporting concatenated FCS file")
        }
        FCS=flowFrame(preds_rawbgc_linear)
        FCS@parameters$desc=as.character(FCS@parameters$desc)
        FCS@parameters$name=as.character(FCS@parameters$name)
        FCS=generate_description(FCS)
        dir.create(file.path(paths["output"],"FCS_background_corrected/","concatenated"),recursive=TRUE,showWarnings=FALSE)
        invisible(write.FCS(FCS,file=file.path(paths["output"],"FCS_background_corrected","concatenated","concatenated_results_background_corrected.fcs")))
    }

    saveRDS(preds_rawbgc,file.path(paths["rds"],"predictions_bgc_cbound.Rds"))
    preds_rawbgc
}
