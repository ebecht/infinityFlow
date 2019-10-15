#' Wrapper to the Infinity Flow pipeline
#' @param path_to_fcs Path to the input directory where input FCS files are stored (one file per well). Will look for FCS files recursively in that directory.
#' @param path_to_output Path to the output directory where final results will be stored
#' @param path_to_intermediary_results Path to results to store temporary data. If left blank, will default to a temporary directory. It may be useful to store the intermediary results to further explore the data, tweak the pipeline or to resume computations.
#' @param backbone_selection_file If that argument is missing and R is run interactively, the user will be prompted to state whether each channel in the FCS file should be considered backbone measurement, exploratory measurement or ignored. Otherwise, the user should run \code{\link{select_backbone_and_exploratory_markers}} in an interactive R session, save its output using \emph{write.csv(row.names=FALSE)} and set this \emph{backbone_selection_file} parameter to the path of the saved output.
#' @param annotation Named character vector. Elements should be the targets of the exploratory antibodies, names should be the name of the FCS file where that exploratory antibody was measured.
#' @param annotation Named character vector. Elements should be the isotype used in each of the well and that (e.g. IgG2). The corresponding isotype should be present in \emph{annotation} (e.g. Isotype_IgG2, with this capitalization exactly). Autofluorescence measurements should be listed here as "Blank"
#' @param input_events_downsampling How many event should be kept per input FCS file. Default to no downsampling. In any case, half of the events will be used to train regression models and half to test the performance. Predictions will be made only on events from the test set.
#' @param cores Number of cores to use for parallel computing. Defaults to 1 (no parallel computing)
#' @param your_random_seed Set a seed for reproducible results. Defaults to 123
#' @param verbose Whether to print information about progress
#' @param extra_args_plotting list of named arguments to pass to plot_results. Defaults to list(chop_quantiles=0.005) which removes the top 0.05% and bottom 0.05% of the scale for each marker when mapping color palettes to intensities.
#' @param extra_args_read_FCS list of named arguments to pass to flowCore:read.FCS. Defaults to list(emptyValue=FALSE,truncate_max_range=FALSE,ignore.text.offset=TRUE) which in our experience avoided issues with data loading.
#' @param extra_args_UMAP list of named arguments to pass to uwot:umap. Defaults to list(n_neighbors=15L,min_dist=0.2,metric="euclidean",verbose=verbose,n_epochs=1000L)
#' @export

infinity_flow=function(
                       ## Input FCS files
                       path_to_fcs, ## Where the source FCS files are
                       path_to_output, ## Where the results will be stored
                       path_to_intermediary_results=tempdir(), ## Storing intermediary results. Default to a temporary directory. Can be a user-specified directory to store intermediary results (to resume interrupted computation)
                       backbone_selection_file=NULL, ## Define backbone and exploratory channels. If missing will be defined interactively and the selection will be saved in the output folder under the name backbone_selection_file.csv. To define it the first time you should call select_backbone_and_exploratory_markers(read.files(path_to_fcs,recursive=TRUE)) in an interactive R session 
                       

                       ## Annotation
                       annotation=NULL, ## Named vector with names = files. Use name of input files if missing
                       isotype=NULL, ## Named vector with names = files and values = which isotype this target maps to

                       ## Downsampling
                       input_events_downsampling=Inf, ## Number of cells to downsample to per input FCS file
                       prediction_events_downsampling=1000, ## Number of events to get predictions for per file
                       
                       ## Set-up multicore computing. When in doubt, set it to 1
                       cores=1L,
                       
                       ## Setting a random seed (for reproducibility despite stochasticity). We enforce reproducibility for the predictions even when doing multicore comparisons. 
                       your_random_seed=123,
                       verbose=TRUE,

                       extra_args_read_FCS=list(emptyValue=FALSE,truncate_max_range=FALSE,ignore.text.offset=TRUE),
                       regression_functions=list(SVM=fitter_svm,XGBoost=fitter_xgboost,LM2=fitter_linear,NN=fitter_nn,LASSO=fitter_glmnet)[1],
                       extra_args_regression_params=list(
                           list(type="eps-regression",cost=1,epsilon=0.5),
                           list(nrounds=10),
                           list(degree=2),
                           list(
                               object={
                                   model=keras_model_sequential()
                                   model %>%
                                       layer_dense(units = 8, activation = "relu", input_shape=14) %>%
                                       layer_dense(units = 8, activation = "relu") %>%
                                       layer_dense(units = 8, activation = "relu") %>%
                                       layer_dense(units = 8, activation = "relu") %>%
                                       layer_dense(units = 8, activation = "relu") %>%
                                       layer_dense(units = 8, activation = "relu") %>%
                                       layer_dense(units = 8, activation = "relu") %>%
                                       layer_dense(units = 8, activation = "relu") %>%
                                       layer_dense(units = 1, activation = "linear")
                                   model %>% compile(
                                                 loss = "mean_squared_error",
                                                 optimizer = optimizer_sgd(lr = 0.01)
                                             )
                                   serialize_model(model)
                               },
                               ##callbacks=list(callback_early_stopping(monitor = "val_loss", patience = 20)), ## This was moved inside fitter_nn because it causes crashes
                               epochs=1000,
                               validation_split=0.2,
                               verbose=0,
                               batch_size=128 ## If you want reproducibility with neural networks you should also use use_session_with_seed(your_random_seed + 3) before calling infinity flow. Set to NULL if you want to remove reproducibility of the NN (which can be good practice to ensure robustness)
                           ),
                           list(degree=2,alpha=1,nfolds=10)
                       )[1],
                       extra_args_UMAP=list(n_neighbors=15L,min_dist=0.2,metric="euclidean",verbose=verbose,n_epochs=1000L,n_threads=cores,n_sgd_threads=cores),
                       extra_args_export=list(FCS_export=c("split","concatenated","none")[1],CSV_export=FALSE),
                       extra_args_correct_background=list(FCS_export=c("split","concatenated","none")[1],CSV_export=FALSE),
                       extra_args_plotting=list(chop_quantiles=0.005),
                       neural_networks_seed=your_random_seed + 3 ## Set to NULL to disable reproducibility
                       ){
    ## Loading packages
    require(flowCore) ##Bioconductor
    require(e1071) ## For SVMs
    require(pbmcapply)
    require(matlab)
    require(png)
    require(raster)
    require(grid)
    require(uwot)
    require(xgboost)
    require(keras)
    require(glmnet)
    require(glmnetUtils)
    require(butcher)
    
    if(!is.null(neural_networks_seed)){
        use_session_with_seed(neural_networks_seed)
    }

    if(length(extra_args_regression_params)!=length(regression_functions)){
        stop("extra_args_regression_params and regression_functions should be lists of the same lengths")
    }
        
    ##/!\ Potentially add a check here to make sure parameters are consistent with FCS files

    settings=initialize(
        path_to_fcs=path_to_fcs,
        path_to_output=path_to_output,
        path_to_intermediary_results=path_to_intermediary_results,
        backbone_selection_file=backbone_selection_file,
        annotation=annotation,
        isotype=isotype,
        verbose=verbose,
        regression_functions=regression_functions
    )
    name_of_PE_parameter=settings$name_of_PE_parameter
    paths=settings$paths
    regression_functions=settings$regression_functions
    
    ## Subsample FCS files
    M=input_events_downsampling ## Number of cells to downsample to for each file
    set.seed(your_random_seed)
    subsample_data(
        input_events_downsampling=input_events_downsampling,
        paths=paths,
        extra_args_read_FCS=extra_args_read_FCS,
        name_of_PE_parameter=name_of_PE_parameter,
        verbose=verbose
    )
    
    ## Automatically scale backbone markers and each PE. ## /!\ This only accomodates a single exploratory measurement
    logicle_transform_input(
        yvar=name_of_PE_parameter,
        paths=paths,
        verbose=verbose
    )
    
    ## Z-score transform each backbone marker for each sample
    standardize_backbone_data_across_wells(
        yvar=name_of_PE_parameter,
        paths=paths,
        verbose=verbose
    )

    ## Regression models training and predictions
    set.seed(your_random_seed+1)
    timings_fit=fit_regressions(
        regression_functions=regression_functions,
        yvar=name_of_PE_parameter,
        paths=paths,
        cores=cores,
        params=extra_args_regression_params,
        verbose=verbose,
        neural_networks_seed=neural_networks_seed
        )

    set.seed(your_random_seed+2)
    timings_pred=predict_from_models(
        paths=paths,
        prediction_events_downsampling=prediction_events_downsampling,
        cores=cores,
        verbose=verbose,
        neural_networks_seed=neural_networks_seed
    )
    
    ## UMAP dimensionality reduction
    perform_UMAP_dimensionality_reduction(
        paths=paths,
        extra_args_UMAP=extra_args_UMAP,
        verbose=verbose
    )
    
    ## Export of the data and predicted data
    res=do.call(export_data,c(list(paths=paths,verbose=verbose),extra_args_export))
    
    ## Plotting
    do.call(plot_results,c(list(paths=paths,verbose=verbose),extra_args_plotting))

    ## Background correcting
    res_bgc=do.call(correct_background,c(list(paths=paths,verbose=verbose),extra_args_correct_background))

    ## Plotting background corrected

    ## Plotting
    do.call(
        plot_results,
        c(
            list(
                paths=paths,
                verbose=verbose,
                preds=readRDS(file.path(paths["rds"],"predictions_bgc_cbound.Rds")),
                file_name=file.path(paths["output"],"umap_plot_annotated_backgroundcorrected.pdf"),
                prediction_colnames = paste0(readRDS(file.path(paths["rds"],"prediction_colnames.Rds")),"_bgc")
            ),
            extra_args_plotting
        )
    )
        
    timings=cbind(fit=timings_fit$timings,pred=timings_pred$timings)
    rownames(timings)=names(regression_functions)
    saveRDS(timings,file.path(paths["rds"],"timings.Rds"))
    
    list(raw=res,bgc=res_bgc,timings=timings)
}

initialize=function(
                    path_to_fcs=path_to_fcs,
                    path_to_output=path_to_output,
                    path_to_intermediary_results=path_to_intermediary_results,
                    backbone_selection_file=backbone_selection_file,
                    annotation=annotation,
                    isotype=isotype,
                    verbose=TRUE,
                    regression_functions=regression_functions
                    ){
    
    ## The paths below have to point to directories. If they do not exist the script will create them to store outputs
    path_to_subsetted_fcs=file.path(path_to_intermediary_results,"subsetted_fcs") ## Subsetted here means for instance "gated on live singlets CD45+"
    path_to_rds=file.path(path_to_intermediary_results,"rds")
    ## A CSV file to map files to PE targets
    path_to_annotation_file=file.path(path_to_intermediary_results,"annotation.csv") ## Has to be a comma-separated csv file, with two columns. The first column has to be the name of the FCS files and the second the marker bound to the PE reporter. The first line (column names) have to be "file" and "target". You can leave the unlabelled PEs empty

    paths=c(
        input=path_to_fcs,
        intermediary=path_to_intermediary_results,
        subset=path_to_subsetted_fcs,
        rds=path_to_rds,
        annotation=path_to_annotation_file,
        output=path_to_output
    )
    
    if(path_to_intermediary_results==tempdir()){
        if(verbose){
            message("Using ",tempdir()," temporary directory to store intermediary results as no non-temporary directory has been specified")
        }
        sapply(paths[c("subset","rds")],dir.create,showWarnings=FALSE,recursive=TRUE)
    }
    if(!all(dir.exists(paths[-match("annotation",names(paths))]))){
        message(paste(paths[-match("annotation",names(paths))][!dir.exists(paths[-match("annotation",names(paths))])],collapse=" and "),": directories not found, creating directory(ies)")
        sapply(paths[-match("annotation",names(paths))][!dir.exists(paths[-match("annotation",names(paths))])],dir.create,recursive=TRUE,showWarnings=FALSE)
    }

    if(verbose){
        message("Using directories...")
        sapply(names(paths),function(x){message("\t",x,": ",paths[x])})
    }
    
    files=list.files(path_to_fcs,pattern="^.*.(fcs)$",ignore.case=TRUE,recursive=TRUE)
    if(missing(annotation)){
        annotation=setNames(files,files)
    }
    annotation=data.frame(file=names(annotation),target=annotation)
    if(!missing(isotype)){
        annotation=cbind(annotation,isotype=annotation$file[match(isotype,annotation$target)])
    }
    annotation$target=make.unique(as.character(annotation$target))
    write.csv(annotation,row.names=FALSE,file=paths["annotation"])
    
    files=list.files(path_to_fcs,pattern="^.*.(fcs)$",ignore.case=TRUE,recursive=TRUE,full.names=TRUE)
    if(missing(backbone_selection_file)){
        backbone_definition=select_backbone_and_exploratory_markers(files)
        write.csv(backbone_definition,file=file.path(path_to_output,"backbone_selection_file.csv"),row.names=FALSE)
    } else {
        backbone_definition=read.csv(backbone_selection_file,stringsAsFactors=FALSE)
    }
    backbone_definition$desc[is.na(backbone_definition$desc)]=backbone_definition$name[is.na(backbone_definition$desc)]
    chans=subset(backbone_definition,type=="backbone")
    chans=setNames(chans$desc,chans$name)
    name_of_PE_parameter=subset(backbone_definition,type=="exploratory")$name

    saveRDS(chans,file=file.path(path_to_rds,"chans.Rds"))

    if(is.null(names(regression_functions))){
        names(regression_functions)=paste0("Alg",seq_along(regression_functions))
    }
    w=is.na(names(regression_functions))|names(regression_functions)==""
    if(any(w)){
        names(regression_functions)[w]=paste0("Alg",seq_along(regression_functions))[w]
    }
    names(regression_functions)=make.unique(names(regression_functions))

    return(
        list(
            paths=paths,
            chans=chans,
            name_of_PE_parameter=name_of_PE_parameter,
            regression_functions=regression_functions
        )
    )
}
