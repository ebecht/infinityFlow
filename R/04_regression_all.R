#' Wrapper to SVM training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @export
fitter_svm=function(x){
    w=x[,"train_set"]==1
    model=do.call(function(...){svm(...,x=x[w,chans],y=x[w,yvar])},params)
    pred=predict(model,x[,chans])
    x=cbind(x,SVM=pred)
    return(list(data=x,svm=model))
}

#' Wrapper to XGBoost training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @export
fitter_xgboost=function(x){
    w=x[,"train_set"]==1
    model=do.call(function(...){xgboost(...,data=x[w,chans],label=x[w,yvar],nthread=1L)},params)
    pred=predict(model,x[,chans])
    x=cbind(x,XGBoost=pred)
    return(list(data=x,xgboost=model))
}

polynomial_formula=function(variables,degree){
    require(gtools)
    n=length(variables)
    polys=lapply(
        1:degree,
        function(deg){
            res=apply(gtools::combinations(n,deg,variables,repeats.allowed=TRUE),1,paste,collapse="*")
            res=paste0("I(",res,")")
        }
    )
    paste(do.call(c,lapply(polys,paste,sep="+")),collapse="+")
}

#' Wrapper to linear model training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @export
fitter_linear=function(x){
    w=x[,"train_set"]==1
    fmla=paste0(make.names(yvar),"~",polynomial_formula(variables=chans,degree=params$degree))
    model=lm(formula=fmla,data=as.data.frame(x[w,c(chans,yvar)]))
    pred=predict(model,as.data.frame(x[,chans]))
    x=cbind(x,LM=pred)
    return(list(data=x,lm=model))
}

#' Wrapper to SVM predict. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_svms
predict_wrapper=function(x){
    if(class(x)=="lm"){
        xp=as.data.frame(xp)
    }
    res=predict(x,xp)
}
            
#' Train SVM regressions
#' @param yvar name of the exploratory measurement
#' @param params Named list of arguments passed to regression fitting functions
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param cores Number of cores to use for parallel computing 
fit_regressions=function(
                  yvar,
                  params,
                  paths,
                  cores,
                  xp=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds")),
                  chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                  events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
                  transforms_chan=readRDS(file.path(paths["rds"],"transforms_chan.Rds")),
                  transforms_pe=readRDS(file.path(paths["rds"],"transforms_pe.Rds")),
                  regression_function
                  )
{
    require(parallel)
    cl=makeCluster(cores)
    
    RNGkind("L'Ecuyer-CMRG")
    mc.reset.stream()

    env=environment()

    colnames(xp)=make.names(colnames(xp))
    
    d.e=split_matrix(xp,events.code)
    rm(xp)
    d.e=lapply(
        d.e,
        function(x){
            w=sample(rep(c(TRUE,FALSE),times=c(floor(nrow(x)/2),nrow(x)-floor(nrow(x)/2))))
            x=cbind(x,train_set=ifelse(w,1,0))
            x
        }
    )
    
    clusterEvalQ(
        cl,
        {
            library(e1071)
            library(xgboost)
            library(gtools)
        }
    )

    clusterExport(
        cl,
        c("yvar","chans","params"),
        envir=env
    )

    clusterEvalQ(
        cl,
        {
            chans=make.names(chans)
            yvar=make.names(yvar)
        }
    )

    models=parLapplyLB(
        X=d.e,
        fun=regression_function,
        cl=cl
    )

    stopCluster(cl)
    saveRDS(models,file=file.path(paths["rds"],"regression_models.Rds"))
    invisible()
}

#' Predict missing measurements using fitted SVM regressions
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param prediction_events_downsampling Number of events to predict data for per file
#' @param cores Number of cores to use for parallel computing

predict_from_models=function(
                      paths,
                      prediction_events_downsampling,
                      cores,
                      chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                      events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
                      models=readRDS(file.path(paths["rds"],"regression_models.Rds"))
                      )
{
    require(parallel)
    cl=makeCluster(cores)

    mc.reset.stream()
    
    env=environment()
    
    xp=lapply(models,"[[",1)
    ## xp=lapply(xp,function(x){x[x[,"train_set"]==0,]})
    ## xp=lapply(xp,function(x){x[sort(sample(1:nrow(x),min(prediction_events_downsampling,nrow(x)))),]})
    spl=lapply(
        xp,
        function(x){
            spl=rep(FALSE,nrow(x))
            w=x[,"train_set"]==0
            spl[w][sample(1:sum(w),min(prediction_events_downsampling,sum(w)))]=TRUE
            spl
        }
    )
    sampling=which(do.call(c,spl))

    for(i in seq_along(xp)){
        xp[[i]]=xp[[i]][spl[[i]],]
    }
    
    xp=do.call(rbind,xp)[,make.names(chans)]

    models=lapply(models,"[[",2)

    clusterEvalQ(
        cl,
        {
            library(e1071)
            library(xgboost)
            library(gtools)
        }
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
                X=models,
                predict_wrapper
            ),
            names(models)
        )
    )
    stopCluster(cl)

    colnames(xp)=chans
    preds=cbind(xp,preds)
    preds=as.matrix(preds)
    
    saveRDS(sampling,file=file.path(paths["rds"],"sampling_preds.Rds"))   
    saveRDS(preds,file=file.path(paths["rds"],"predictions.Rds"))
    invisible()
}

