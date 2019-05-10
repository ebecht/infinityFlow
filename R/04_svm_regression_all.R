#' Wrapper to SVM training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_svms
svm_fitter=function(x){
            try({
                w=sample(rep(c(TRUE,FALSE),times=c(floor(nrow(x)/2),nrow(x)-floor(nrow(x)/2))))
                svm=do.call(function(...){svm(...,x=x[w,chans],y=x[w,yvar])},params)
                pred=predict(svm,x[,chans])
                x=cbind(x,SVM=pred,train_set=ifelse(w,1,0))
                return(list(data=x,svm=svm))
            })
}

#' Wrapper to SVM predict. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_svms
svm_predictor=function(x){
    res=predict(x,xp)
}
            
#' Train SVM regressions
#' @param yvar name of the exploratory measurement
#' @param params Named list of arguments passed to e1071:svm()
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param cores Number of cores to use for parallel computing 
fit_svms=function(
                  yvar,
                  params=list(
                      type="eps-regression",
                      cost=1,
                      epsilon=0.5
                  ),
                  paths,
                  cores,
                  xp=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds")),
                  chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                  events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
                  transforms_chan=readRDS(file.path(paths["rds"],"transforms_chan.Rds")),
                  transforms_pe=readRDS(file.path(paths["rds"],"transforms_pe.Rds"))
                  )
{
    require(parallel)
    cl=makeCluster(cores)
    
    RNGkind("L'Ecuyer-CMRG")
    mc.reset.stream()

    env=environment()
    
    ## sapply(
    ##     c("chans","transforms_chan","transforms_pe"),
    ##     function(object){
    ##         assign(object,value=readRDS(file.path(paths["rds"],paste0(object,".Rds"))),envir=env)
    ##         invisible()
    ##     }
    ## )
    ## events.code=readRDS(file.path(paths["rds"],"pe.Rds"))
    ## xp=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds"))

    d.e=split_matrix(xp,events.code)
    rm(xp)

    clusterEvalQ(
        cl,
        library(e1071)
    )

    clusterExport(
        cl,
        c("yvar","chans","params"),
        envir=env
    )

    svms=parLapplyLB(
        X=d.e,
        fun=svm_fitter,
        cl=cl
    )

    stopCluster(cl)
    
    ## svms=pbmclapply(
    ##     d.e,
    ##     function(x){
    ##         try({
    ##             w=sample(rep(c(TRUE,FALSE),times=c(floor(nrow(x)/2),nrow(x)-floor(nrow(x)/2))))
    ##             svm=do.call(function(...){svm(...,x=x[w,chans],y=x[w,yvar])},params)
    ##             pred=predict(svm,x[,chans])
    ##             x=cbind(x,SVM=pred,train_set=ifelse(w,1,0))
    ##             return(list(data=x,svm=svm))
    ##         })
    ##     },
    ##     mc.cores=cores,
    ##     mc.preschedule=FALSE,
    ##     mc.set.seed=TRUE
    ## )
    saveRDS(svms,file=file.path(paths["rds"],"svms_models.Rds"))
    invisible()
}

#' Predict missing measurements using fitted SVM regressions
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param prediction_events_downsampling Number of events to predict data for per file
#' @param cores Number of cores to use for parallel computing

predict_svms=function(
                      paths,
                      prediction_events_downsampling,
                      cores,
                      chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                      events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
                      svms=readRDS(file.path(paths["rds"],"svms_models.Rds"))
                      )
{
    require(parallel)
    cl=makeCluster(cores)

    mc.reset.stream()
    
    env=environment()
    ## sapply(
    ##     c("chans"),
    ##     function(object){
    ##         assign(object,value=readRDS(file.path(paths["rds"],paste0(object,".Rds"))),envir=env)
    ##         invisible()
    ##     }
    ## )
    ## events.code=readRDS(file.path(paths["rds"],"pe.Rds"))
    ## svms=readRDS(file.path(paths["rds"],"svms_models.Rds"))
    
    xp=lapply(svms,"[[",1)
    xp=lapply(xp,function(x){x[x[,"train_set"]==0,]})
    xp=lapply(xp,function(x){x[sort(sample(1:nrow(x),min(prediction_events_downsampling,nrow(x)))),]})
    xp=do.call(rbind,xp)[,chans]
    svms=lapply(svms,"[[",2)

    clusterEvalQ(
        cl,
        library(e1071)
    )

    clusterExport(
        cl,
        c("xp"),
        envir=env
    )
    
    preds=do.call(
        cbind,
        setNames(
            parLapplyLB(
                cl=cl,
                X=svms,
                svm_predictor
            ),
            names(svms)
        )
    )
    
    stopCluster(cl)
    
    preds=cbind(xp,preds)
    sampling=as.numeric(rownames(preds))
    saveRDS(sampling,file=file.path(paths["rds"],"sampling_preds.Rds"))
    saveRDS(preds,file=file.path(paths["rds"],"svms_predictions.Rds"))
    invisible()
}

