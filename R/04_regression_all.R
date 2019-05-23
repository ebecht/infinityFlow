#' Wrapper to SVM training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @export
fitter_svm=function(x,params){
    require(e1071)
    w=x[,"train_set"]==1
    model=do.call(function(...){svm(...,x=x[w,chans],y=x[w,yvar])},params)
    pred=predict(model,x[,chans])
    return(list(pred=pred,model=model))
}

#' Wrapper to XGBoost training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @export
fitter_xgboost=function(x,params){
    require(xgboost)
    w=x[,"train_set"]==1
    model=do.call(function(...){xgboost(...,data=x[w,chans],label=x[w,yvar],nthread=1L)},params)
    pred=predict(model,x[,chans])
    return(list(pred=pred,model=model))
}

#' Wrapper to linear model training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @export
fitter_linear=function(x,params){
    w=x[,"train_set"]==1
    fmla=paste0(make.names(yvar),"~",polynomial_formula(variables=chans,degree=params$degree))
    model=lm(formula=fmla,data=as.data.frame(x[w,c(chans,yvar)]))
    pred=predict(model,as.data.frame(x[,chans]))
    return(list(pred=pred,model=model))
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
                  regression_functions
                  )
{
    ## xp=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds"));
    ## chans=readRDS(file.path(paths["rds"],"chans.Rds"));
    ## events.code=readRDS(file.path(paths["rds"],"pe.Rds"));
    ## transforms_chan=readRDS(file.path(paths["rds"],"transforms_chan.Rds"));
    ## transforms_pe=readRDS(file.path(paths["rds"],"transforms_pe.Rds"));

                  
    require(parallel)
    cl=makeCluster(min(cores,length(unique(events.code))))
    
    RNGkind("L'Ecuyer-CMRG")
    mc.reset.stream()

    env=environment()
    
    colnames(xp)=make.names(colnames(xp))
    
    d.e=split_matrix(xp,events.code)
    d.e=lapply(
        d.e,
        function(x){
            w=sample(rep(c(TRUE,FALSE),times=c(floor(nrow(x)/2),nrow(x)-floor(nrow(x)/2))))
            x=cbind(x,train_set=ifelse(w,1,0))
            x
        }
    )
    train_set=matrix(ncol=1,dimnames=list(NULL,"train_set"),do.call(c,sapply(d.e,function(x){x[,"train_set"]},simplify=FALSE)))
    
    clusterExport(
        cl,
        c("yvar","chans"),
        envir=env
    )

    clusterEvalQ(
        cl,
        {
            chans=make.names(chans)
            yvar=make.names(yvar)
        }
    )
    
    models=list()
    timings=numeric()
    for(i in seq_along(regression_functions)){
        cat("\t",names(regression_functions)[i])
        t0=Sys.time()
        models[[i]]=parLapplyLB(
            X=d.e,
            fun=regression_functions[[i]],
            params=params[[i]],
            cl=cl
        )
        t1=Sys.time()
        dt=difftime(t1,t0,units="secs")
        cat("\t",dt," seconds","\n")
        timings=c(timings,dt)
    }
    names(models)=names(regression_functions)
    
    stopCluster(cl)
    saveRDS(models,file=file.path(paths["rds"],"regression_models.Rds"))
    saveRDS(train_set,file=file.path(paths["rds"],"train_set.Rds"))
    list(timings=timings)
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
                      models=readRDS(file.path(paths["rds"],"regression_models.Rds")),
                      xp=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds")),
                      train_set=readRDS(file.path(paths["rds"],"train_set.Rds"))
                      )
{
    ## chans=readRDS(file.path(paths["rds"],"chans.Rds"));
    ## events.code=readRDS(file.path(paths["rds"],"pe.Rds"));
    ## models=readRDS(file.path(paths["rds"],"regression_models.Rds"));
    ## xp=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds"));
    ## train_set=readRDS(file.path(paths["rds"],"train_set.Rds"))
    
    require(parallel)
    cl=makeCluster(cores)
    
    clusterEvalQ(
        cl,
        {
            library(xgboost)
            library(e1071)
        }
    )

    mc.reset.stream()
    
    env=environment()
    
    xp=cbind(xp,train_set)
    xp=split_matrix(xp,events.code)
    spl=lapply(
        xp,
        function(x){
            spl=rep(FALSE,nrow(x))
            w=x[,"train_set"]==0
            spl[w][sample(1:sum(w),min(prediction_events_downsampling,sum(w)))]=TRUE
            spl
        }
    )
    pred_set=which(do.call(c,spl))
    
    xp=do.call(rbind,xp)
    xp=xp[pred_set,]
    
    clusterExport(
        cl,
        c("xp","chans"),
        envir=env
    )
    invisible(clusterEvalQ(
        cl,
        {
            colnames(xp)=make.names(colnames(xp))
            xp=xp[,make.names(chans)]
        }
    ))
    
    for(i in seq_along(models)){
        models[[i]]=lapply(models[[i]],"[[",2)
    }
    
    preds=list()
    timings=numeric()
    t0=Sys.time()
    for(i in seq_along(models)){
        cat("\t",names(models)[i])
        preds[[i]]=do.call(
            cbind,
            parLapplyLB(
                cl=cl,
                X=models[[i]],
                predict_wrapper
            )
        )
        t1=Sys.time()
        dt=difftime(t1,t0,units="secs")
        cat("\t",dt," seconds","\n")
        timings=c(timings,dt)
    }
    
    stopCluster(cl)

    preds=lapply(preds,function(x){cbind(xp,x)})
    preds=lapply(preds,as.matrix)
    names(preds)=names(models)
    
    saveRDS(preds,file=file.path(paths["rds"],"predictions.Rds"))
    saveRDS(pred_set,file=file.path(paths["rds"],"sampling_preds.Rds"))   
    list(timings=timings)
}

