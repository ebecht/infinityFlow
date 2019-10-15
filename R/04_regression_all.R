#' Wrapper to SVM training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @export
fitter_svm=function(x,params){
    require(e1071)
    w=x[,"train_set"]==1
    model=do.call(function(...){svm(...,x=x[w,chans],y=x[w,yvar])},params)
    ## model=butcher(model)
    pred=predict(model,x[,chans])
    rm(list=setdiff(ls(),c("pred","model")))
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
    rm(list=setdiff(ls(),c("pred","model")))
    return(list(pred=pred,model=model))
}

#' Wrapper to linear model training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @export
fitter_linear=function(x,params){
    w=x[,"train_set"]==1
    fmla=paste0(make.names(yvar),"~",polynomial_formula(variables=chans,degree=params$degree))
    model=lm(formula=fmla,data=as.data.frame(x[w,c(chans,yvar)]))
    ## model=butcher(model)
    pred=predict(model,as.data.frame(x[,chans]))
    rm(list=setdiff(ls(),c("pred","model")))
    model$model=NULL ## Trim down for slimmer objects
    model$qr$qr=NULL ## Trim down for slimmer objects
    return(list(pred=pred,model=model))
}

#' Wrapper to glmnet. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @export
fitter_glmnet=function(x,params){
    w=x[,"train_set"]==1
    fmla=paste0(make.names(yvar),"~",polynomial_formula(variables=chans,degree=params$degree))
    params=params[setdiff(names(params),"degree")]
    params=c(
        params,
        list(
            formula=fmla,
            data=as.data.frame(x[w,c(chans,yvar)]),
            use.model.frame=TRUE
        )
    )
    ## MF=model.frame(fmla,data=as.data.frame(x[w,c(chans,yvar)]))
    ## model=glmnet(x=as.matrix(MF[,-match(yvar,colnames(MF))]),y=MF[,yvar],alpha=1)
    ## cvfit=cv.glmnet(as.matrix(MF[,-match(yvar,colnames(MF))]),y=MF[,yvar],alpha=1,type.measure="mse",nfolds=20)
    ## lambda.min = cvfit$lambda.min

    model=do.call(glmnetUtils:::cv.glmnet.formula,params)
    ## model=butcher:::axe_call.glmnet(model)
    model$call = NULL ## Slimming down object
    model$glmnet.fit$call = NULL ## Slimming down object
    attributes(model$terms)[[".Environment"]] = NULL ## Slimming down object
    pred=predict(model,as.data.frame(x[,chans]),s=model$lambda.min)
    
    ## attributes(model$terms)[".Environment"]=NULL ## Trim down for slimmer objects

    
    ## pred=predict(model,newx=as.matrix(MF[,-match(yvar,colnames(MF))]),s=lambda.min)[, 1]
    rm(list=setdiff(ls(),c("pred","model")))
    ## model$model=NULL ## Trim down for slimmer objects
    ## model$qr$qr=NULL ## Trim down for slimmer objects
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
    if("lm"%in%class(x)){
        xp=as.data.frame(xp)
    }
    if("cv.glmnet"%in%class(x)){
        xp=as.data.frame(xp)
        return(predict(x,xp,s=x$lambda.min)[,1])
    }
    if(class(x)=="raw"){
        require(keras)
        x=unserialize_model(x)
    }
    if(class(x)=="xgb.Booster"){
        x=xgb.Booster.complete(x)
        xgb.parameters(x) <- list(nthread = 1)
    }
    res=predict(x,xp)
}

#' Wrapper to Neural Network training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions. Defines model architecture
#' @export
fitter_nn=function(x,params){
    require(keras)
    model=unserialize_model(params$object)
    params=params[setdiff(names(params),"object")]

    early_stop=callback_early_stopping(monitor = "val_loss", patience = 20)
    callbacks=list(early_stop)
    
    w=x[,"train_set"]==1

    fit_history=do.call(
        function(...){
            fit(
                ...,
                object=model,
                x=x[w,chans],
                y=x[w,yvar],
                callbacks=callbacks
            )
        },
        params
    )
    
    pred = predict(model, x[, chans])
    rm(list=setdiff(ls(),c("pred","model","fit_history")))
    return(list(pred = pred, model = serialize_model(model), fit_history = fit_history))
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
                  regression_functions,
                  verbose=TRUE,
                  neural_networks_seed
                  )
{
    ## xp=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds"));
    ## chans=readRDS(file.path(paths["rds"],"chans.Rds"));
    ## events.code=readRDS(file.path(paths["rds"],"pe.Rds"));
    ## transforms_chan=readRDS(file.path(paths["rds"],"transforms_chan.Rds"));
    ## transforms_pe=readRDS(file.path(paths["rds"],"transforms_pe.Rds"));
    ## yvar = name_of_PE_parameter
    ## chans=make.names(chans)
    ## yvar=make.names(yvar)
    
    if(verbose){
        message("Fitting regression models")
    }
    
    require(parallel)
    cl=makeCluster(min(cores,length(unique(events.code))))
    
    RNGkind("L'Ecuyer-CMRG")
    mc.reset.stream()

    env=environment()
    
    colnames(xp)=make.names(colnames(xp))

    
    if(verbose){
        message("\tRandomly selecting 50% of the subsetted input files to fit models")
    }
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
        c("yvar","chans","neural_networks_seed"),
        envir=env
    )

    clusterEvalQ(
        cl,
        {
            chans=make.names(chans)
            yvar=make.names(yvar)
            library(tensorflow)
            library(keras)
            library(glmnetUtils)
            library(glmnet)
            library(butcher)
            if(!is.null(neural_networks_seed)){
                use_session_with_seed(neural_networks_seed) ## This will make results reproducible, disable GPU and CPU parallelism (which is good actually). Source: https://keras.rstudio.com/articles/faq.html#how-can-i-obtain-reproducible-results-using-keras-during-development
            }
        }
    )
    
    if(verbose){
        message("\tFitting...")
    }
    models=list()
    timings=numeric()
    for(i in seq_along(regression_functions)){
        cat("\t\t",names(regression_functions)[i],sep="")
        t0=Sys.time()
        models[[i]]=parLapplyLB(
            X=d.e,
            fun=regression_functions[[i]],
            params=params[[i]],
            cl=cl
        )
        t1=Sys.time()
        dt=difftime(t1,t0,units="secs")
        cat("\t",dt," seconds","\n",sep="")
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
                      train_set=readRDS(file.path(paths["rds"],"train_set.Rds")),
                      verbose=TRUE,
                      neural_networks_seed
                      )
{
    ## chans=readRDS(file.path(paths["rds"],"chans.Rds"));
    ## events.code=readRDS(file.path(paths["rds"],"pe.Rds"));
    ## models=readRDS(file.path(paths["rds"],"regression_models.Rds"));
    ## xp=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds"));
    ## train_set=readRDS(file.path(paths["rds"],"train_set.Rds"))
    
    if(verbose){
        message("Imputing missing measurements")
    }
    
    require(parallel)
    cl=makeCluster(cores)
    
    clusterEvalQ(
        cl,
        {
            library(xgboost)
            library(e1071)
            library(keras)
            library(tensorflow)
            library(glmnetUtils)
            library(glmnet)
        }
    )

    mc.reset.stream()
    
    env=environment()

    if(verbose){
        message("\tRandomly drawing events to predict from the test set")
    }

    spl=split(train_set[,1],events.code)
    spl=lapply(
        spl,
        function(x){
            res=rep(FALSE,length(x))
            w=x==0
            res[w][sample(1:sum(w),min(prediction_events_downsampling,sum(w)))]=TRUE
            res
        }
    )

    pred_set=which(do.call(c,spl))
    
    xp=xp[pred_set,]
    
    clusterExport(
        cl,
        c("xp","chans","neural_networks_seed"),
        envir=env
    )
    invisible(clusterEvalQ(
        cl,
        {
            colnames(xp)=make.names(colnames(xp))
            xp=xp[,make.names(chans)]
            if(!is.null(neural_networks_seed)){
                use_session_with_seed(neural_networks_seed) ## This will make results reproducible, disable GPU and CPU parallelism (which is good actually). Source: https://keras.rstudio.com/articles/faq.html#how-can-i-obtain-reproducible-results-using-keras-during-development
            }
        }
    ))
    
    for(i in seq_along(models)){
        models[[i]]=lapply(models[[i]],"[[",2)
    }

    if(verbose){
        message("\tImputing...")
    }
    preds=list()
    timings=numeric()
    for(i in seq_along(models)){
        cat("\t\t",names(models)[i],sep="")
        t0=Sys.time()
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
        cat("\t",dt," seconds","\n",sep="")
        timings=c(timings,dt)
        colnames(preds[[i]])=names(models[[i]])
    }
    
    stopCluster(cl)

    if(verbose){
        message("\tConcatenating predictions")
    }
    preds=lapply(preds,function(x){cbind(xp,x)})
    preds=lapply(preds,as.matrix)
    names(preds)=names(models)

    if(verbose){
        message("\tWriting to disk")
    }
    saveRDS(preds,file=file.path(paths["rds"],"predictions.Rds"))
    saveRDS(pred_set,file=file.path(paths["rds"],"sampling_preds.Rds"))   
    list(timings=timings)
}

