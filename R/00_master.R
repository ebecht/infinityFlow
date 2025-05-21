#' Wrapper to the Infinity Flow pipeline
#' @param path_to_fcs Path to the input directory where input FCS files are stored (one file per well). Will look for FCS files recursively in that directory.
#' @param path_to_output Path to the output directory where final results will be stored
#' @param path_to_intermediary_results Path to results to store temporary data. If left blank, will default to a temporary directory. It may be useful to store the intermediary results to further explore the data, tweak the pipeline or to resume computations.
#' @param backbone_selection_file If that argument is missing and R is run interactively, the user will be prompted to state whether each channel in the FCS file should be considered backbone measurement, exploratory measurement or ignored. Otherwise, the user should run \code{\link{select_backbone_and_exploratory_markers}} in an interactive R session, save its output using \emph{write.csv(row.names=FALSE)} and set this \emph{backbone_selection_file} parameter to the path of the saved output.
#' @param annotation Named character vector. Elements should be the targets of the exploratory antibodies, names should be the name of the FCS file where that exploratory antibody was measured.
#' @param isotype Named character vector. Elements should be the isotype used in each of the well and that (e.g. IgG2). The corresponding isotype should be present in \emph{annotation} (e.g. Isotype_IgG2, with this capitalization exactly). Autofluorescence measurements should be listed here as "Blank"
#' @param regression_functions named list of fitter_* functions (see ls("package:infinityFlow") for the complete list). The names should be desired names for the different models. Each object of the list will correspond to a machine learning model to train. Defaults to list(XGBoost = fitter_xgboost).
#' @param extra_args_regression_params list of lists the same length as the regression_functions argument. Each element should be a named list, that will be passed as named arguments to the corresponding fitter_ function. Defaults to list(list(nrounds = 500, eta = 0.05)).
#' @param input_events_downsampling How many event should be kept per input FCS file. Default to no downsampling. In any case, half of the events will be used to train regression models and half to test the performance. Predictions will be made only on events from the test set, and downsampled according to prediction_events_downsampling.
#' @param prediction_events_downsampling How many event should be kept per input FCS file to output prediction for. Default to 1000.
#' @param cores Number of cores to use for parallel computing. Defaults to 1 (no parallel computing)
#' @param verbose Whether to print information about progress
#' @param extra_args_plotting list of named arguments to pass to plot_results. Defaults to list(chop_quantiles=0.005) which removes the top 0.05\% and bottom 0.05\% of the scale for each marker when mapping color palettes to intensities.
#' @param extra_args_read_FCS list of named arguments to pass to flowCore:read.FCS. Defaults to list(emptyValue=FALSE,truncate_max_range=FALSE,ignore.text.offset=TRUE) which in our experience avoided issues with data loading.
#' @param extra_args_UMAP list of named arguments to pass to uwot:umap. Defaults to list(n_neighbors=15L,min_dist=0.2,metric="euclidean",verbose=verbose,n_epochs=1000L)
#' @param extra_args_export Whether raw imputed data should be exported. Possible values are list(FCS_export = "split") to export one FCS file per input well, list(FCS_export = "concatenated") to export a single concatenated FCS file containing all the dataset, list(FCS_export = "csv") for a single CSV file containing all the dataset. You can export multiple modalities by using for instance extra_args_export = list(FCS_export = c("split", "concatenated", "csv"))
#' @param extra_args_correct_background Whether background-corrected imputed data should be exported. Possible values are list(FCS_export = "split") to export one FCS file per input well, list(FCS_export = "concatenated") to export a single concatenated FCS file containing all the dataset, list(FCS_export = "csv") for a single CSV file containing all the dataset. You can export multiple modalities by using for instance extra_args_export = list(FCS_export = c("split", "concatenated", "csv"))
#' @export

infinity_flow <- function(
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
                          verbose=TRUE,

                          extra_args_read_FCS=list(emptyValue=FALSE,truncate_max_range=FALSE,ignore.text.offset=TRUE),
                          regression_functions=list(
                              XGBoost=fitter_xgboost
                          ),
                          extra_args_regression_params=list(nrounds=100, eta = 0.1),
                          extra_args_UMAP=list(n_neighbors=15L,min_dist=0.2,metric="euclidean",verbose=verbose,n_epochs=1000L,n_threads=cores,n_sgd_threads=cores),
                          extra_args_export=list(FCS_export=c("split","concatenated","csv","none")[1]),
                          extra_args_correct_background=list(FCS_export=c("split","concatenated","csv","none")[1]),
                          extra_args_plotting=list(chop_quantiles=0.005)
                          ){

    ## Making sure that optional dependencies are installed if used.
    lapply(
        regression_functions,
        function(fun){
            fun(x = NULL, params = NULL)
        }
    )
    
    if(any(!isotype %in% annotation)){
        stop("The following values from the isotype argument are not matching any of the values from the annotation argument: ", paste0(isotype[!isotype %in% annotation], collapse = ", "))
    }
    ##/!\ Potentially add a check here to make sure parameters are consistent with FCS files

    settings <- initialize(
        path_to_fcs=path_to_fcs,
        path_to_output=path_to_output,
        path_to_intermediary_results=path_to_intermediary_results,
        backbone_selection_file=backbone_selection_file,
        annotation=annotation,
        isotype=isotype,
        verbose=verbose,
        regression_functions=regression_functions
    )
    name_of_PE_parameter <- settings$name_of_PE_parameter
    paths <- settings$paths
    regression_functions <- settings$regression_functions

    ## Create H5 dataset
    create_h5_file(
        paths=paths,
        input_events_downsampling=input_events_downsampling,
        prediction_events_downsampling=prediction_events_downsampling,
        extra_args_read_FCS=extra_args_read_FCS,
        verbose=verbose
    )

    ## Subsample FCS files
    M  <-  input_events_downsampling ## Number of cells to downsample to for each file
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
    timings_fit <- fit_regressions(
        regression_functions=regression_functions,
        yvar=name_of_PE_parameter,
        paths=paths,
        cores=cores,
        params=extra_args_regression_params,
        verbose=verbose
    )

    timings_pred <- predict_from_models(
        paths=paths,
        prediction_events_downsampling=prediction_events_downsampling,
        cores=cores,
        verbose=verbose
    )
    
    ## UMAP dimensionality reduction
    perform_UMAP_dimensionality_reduction(
        paths=paths,
        extra_args_UMAP=extra_args_UMAP,
        verbose=verbose
    )
    
    ## Export of the data and predicted data
    res <- do.call(export_data,c(list(paths=paths,verbose=verbose,yvar=name_of_PE_parameter),extra_args_export))
        
    ## Plotting
    do.call(plot_results,c(list(paths=paths,verbose=verbose,yvar=name_of_PE_parameter,cores=cores),extra_args_plotting))

    ## Background correcting
    res_bgc <- do.call(correct_background,c(list(paths=paths,verbose=verbose),extra_args_correct_background))

    ## Plotting background corrected
    transforms_tmp <- readRDS(file.path(paths["rds"], "transforms.Rds")) ## Disabling reverse transformation for background-corrected predictions
    transforms_tmp[[name_of_PE_parameter]]$backward = identity
    transforms_tmp[[name_of_PE_parameter]]$forward = identity
    do.call(
        plot_results,
        c(
            list(
                paths=paths,
                verbose=verbose,
                file_name=file.path(paths["output"],"umap_plot_annotated_backgroundcorrected.pdf"),
                predictions_group = "/predictions/background_corrected/",
                yvar=name_of_PE_parameter,
                transforms=transforms_tmp,
                cores=cores
            ),
            extra_args_plotting
        )
    )

    ## Generate quality controls for model fits
    quality_controls(
        paths = paths,
        yvar = name_of_PE_parameter,
        annot=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE),
        verbose=TRUE
    )
        
    timings <- cbind(fit=timings_fit$timings,pred=timings_pred$timings)
    rownames(timings) <- names(regression_functions)
    saveRDS(timings,file.path(paths["rds"],"timings.Rds"))
    
    list(raw=res,bgc=res_bgc,timings=timings)
}

initialize <- function(
                       path_to_fcs=path_to_fcs,
                       path_to_output=path_to_output,
                       path_to_intermediary_results=path_to_intermediary_results,
                       backbone_selection_file=backbone_selection_file,
                       annotation=annotation,
                       isotype=isotype,
                       verbose=TRUE,
                       regression_functions=regression_functions
                       ){
    
    
    if(path_to_intermediary_results==tempdir()){
        tmpdir = path_to_intermediary_results
        if(dir.exists(tmpdir)){
            i = 1
            while(dir.exists(file.path(tmpdir, i))){
                i = i + 1
            }
            path_to_intermediary_results = file.path(tmpdir, i)
        }
        if(verbose){
            message("Using ", tmpdir, " temporary directory to store intermediary results as no non-temporary directory has been specified")
        }
    }

    ## The paths below have to point to directories. If they do not exist the script will create them to store outputs
    path_to_subsetted_fcs <- file.path(path_to_intermediary_results,"subsetted_fcs") ## Subsetted here means for instance "gated on live singlets CD45+"
    path_to_rds <- file.path(path_to_intermediary_results,"rds")
    ## A CSV file to map files to PE targets
    path_to_annotation_file <- file.path(path_to_intermediary_results,"annotation.csv") ## Has to be a comma-separated csv file, with two columns. The first column has to be the name of the FCS files and the second the marker bound to the PE reporter. The first line (column names) have to be "file" and "target". You can leave the unlabelled PEs empty

    path_to_h5 = file.path(path_to_rds, "if.h5")
    paths <- c(
        input=path_to_fcs,
        intermediary=path_to_intermediary_results,
        subset=path_to_subsetted_fcs,
        rds=path_to_rds,
        h5=path_to_h5,
        annotation=path_to_annotation_file,
        output=path_to_output
    )
    paths <- vapply(paths, path.expand, "path")
    
    if(!all(dir.exists(paths[-match(c("annotation","h5"), names(paths))]))){
        message(
            paste(
                paths[-match("annotation",names(paths))][!dir.exists(paths[-match(c("annotation", "h5"),names(paths))])], collapse = " and "),
            ": directories not found, creating directory(ies)"
        )
        vapply(paths[-match(c("annotation", "h5"), names(paths))][!dir.exists(paths[-match(c("annotation", "h5"),names(paths))])],dir.create,recursive=TRUE,showWarnings=FALSE, FUN.VALUE = TRUE)
    }
    
    if(verbose){
        message("Using directories...")
        lapply(names(paths),function(x){message("\t",x,": ",paths[x])})
    }
    
    files <- list.files(path_to_fcs,pattern="^.*.(fcs)$",ignore.case=TRUE,recursive=TRUE)
    files_mismatch = names(annotation)[!names(annotation) %in% basename(list.files(paths["input"]))]
    if(length(files_mismatch) > 0){
        stop("The following files were not found in ", path_to_fcs, " :\n", paste0(files_mismatch, collapse = "\n"))
    }
    if(missing(annotation)){
        annotation <- setNames(files,files)
    }
    annotation <- data.frame(file=names(annotation),target=annotation)
    if(!missing(isotype)){
        annotation <- cbind(annotation,isotype=annotation$file[match(isotype,annotation$target)])
    }
    annotation$target <- make.unique(as.character(annotation$target))
    write.csv(annotation,row.names=FALSE,file=paths["annotation"])
    
    files <- list.files(path_to_fcs,pattern="^.*.(fcs)$",ignore.case=TRUE,recursive=TRUE,full.names=TRUE)
    if(missing(backbone_selection_file)){
        backbone_definition <- select_backbone_and_exploratory_markers(files)
        write.csv(backbone_definition,file=file.path(path_to_output,"backbone_selection_file.csv"),row.names=FALSE)
    } else {
        backbone_definition <- read.csv(backbone_selection_file,stringsAsFactors=FALSE)
    }
    backbone_definition$desc[is.na(backbone_definition$desc)] <- backbone_definition$name[is.na(backbone_definition$desc)]
    chans <- subset(backbone_definition,type=="backbone")
    chans <- setNames(chans$desc,chans$name)
    name_of_PE_parameter <- subset(backbone_definition,type=="exploratory")$name

    transforms = list()
    for(i in seq(1, nrow(backbone_definition), by = 1)){
        if(backbone_definition[i, "type" ]!="discard"){
            if(backbone_definition[i, "type"] == "backbone"){
                name <- backbone_definition[i, "desc"]
                if(is.na(name)){
                    name <- backbone_definition[i, "name"]
                }
            }
            if(backbone_definition[i, "type"] == "exploratory"){
                name = backbone_definition[i, "name"]
            }
            if(backbone_definition[i, "transformation"] == "identity"){
                transforms[[name]]$forward <- identity
                transforms[[name]]$backward <- identity
            } else if(backbone_definition[i, "transformation"] == "asinh"){
                transforms[[name]]$forward <- eval(bquote(function(x){asinh(x/.(backbone_definition[i, "cofactor"]))}))
                transforms[[name]]$backward <- eval(bquote(function(x){.(backbone_definition[i, "cofactor"])*sinh(x)}))
            } else {
                stop("'transformation' column in the backbone_selection_file should only be one of 'identity' or 'asinh'")
            }
        }
    }
    saveRDS(transforms, file=file.path(paths["rds"], "transforms.Rds"))
    saveRDS(chans,file=file.path(paths["rds"],"chans.Rds"))

    if(is.null(names(regression_functions))){
        names(regression_functions) <- paste0("Alg",seq_along(regression_functions))
    }
    w <- is.na(names(regression_functions))|names(regression_functions)==""
    if(any(w)){
        names(regression_functions)[w] <- paste0("Alg",seq_along(regression_functions))[w]
    }
    names(regression_functions) <- make.unique(names(regression_functions))

    return(
        list(
            paths=paths,
            chans=chans,
            name_of_PE_parameter=name_of_PE_parameter,
            regression_functions=regression_functions
        )
    )
}

utils::globalVariables(c("type"))
