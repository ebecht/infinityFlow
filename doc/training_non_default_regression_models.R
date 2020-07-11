## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- packages_installation, eval = FALSE-------------------------------------
#  if(!require(devtools)){
#  	install.packages("devtools")
#  }
#  if(!require(infinityFlow)){
#  	library(devtools)
#  	install_github("ebecht/infinityFlow")
#  }

## ---- preparation-------------------------------------------------------------
library(infinityFlow)

data(steady_state_lung)
data(steady_state_lung_annotation)
data(steady_state_lung_backbone_specification)

dir <- file.path(tempdir(), "infinity_flow_example")
input_dir <- file.path(dir, "fcs")
write.flowSet(steady_state_lung, outdir = input_dir)

write.csv(steady_state_lung_backbone_specification, file = file.path(dir, "backbone_selection_file.csv"), row.names = FALSE)

path_to_fcs <- file.path(dir, "fcs")
path_to_output <- file.path(dir, "output")
path_to_intermediary_results <- file.path(dir, "tmp")
backbone_selection_file <- file.path(dir, "backbone_selection_file.csv")

targets <- steady_state_lung_annotation$Infinity_target
names(targets) <- rownames(steady_state_lung_annotation)
isotypes <- steady_state_lung_annotation$Infinity_isotype
names(isotypes) <- rownames(steady_state_lung_annotation)

input_events_downsampling <- 1000
prediction_events_downsampling <- 500
cores = 1L

## ---- fitter_functions--------------------------------------------------------
print(grep("fitter_", ls("package:infinityFlow"), value = TRUE))

## ---- dependencies------------------------------------------------------------
       optional_dependencies <- c("glmnetUtils", "e1071")
       unmet_dependencies <- setdiff(optional_dependencies, rownames(installed.packages()))
       if(length(unmet_dependencies) > 0){
           install.packages(unmet_dependencies)
       }
       for(pkg in optional_dependencies){
       	   library(pkg, character.only = TRUE)
       }

## ---- regression_functions----------------------------------------------------
regression_functions <- list(
    XGBoost = fitter_xgboost, # XGBoost
    SVM = fitter_svm, # SVM
    LASSO2 = fitter_glmnet, # L1-penalized 2nd degree polynomial model
    LM = fitter_linear # Linear model
)

## ---- extra_args_regression_params--------------------------------------------
backbone_size <- table(read.csv(backbone_selection_file)[,"type"])["backbone"]
extra_args_regression_params <- list(
     ## Passed to the first element of `regression_functions`, e.g. XGBoost. See ?xgboost for which parameters can be passed through this list
    list(nrounds = 500, eta = 0.05),

    # ## Passed to the second element of `regression_functions`, e.g. neural networks through keras::fit. See https://keras.rstudio.com/articles/tutorial_basic_regression.html
    # list(
    #         object = { ## Specifies the network's architecture, loss function and optimization method
    #             model = keras_model_sequential()
    #             model %>%
    #                 layer_dense(units = backbone_size, activation = "relu", input_shape = backbone_size) %>%
    #                 layer_dense(units = backbone_size, activation = "relu", input_shape = backbone_size) %>%
    #                 layer_dense(units = 1, activation = "linear")
    #             model %>%
    #                 compile(loss = "mean_squared_error", optimizer = optimizer_sgd(lr = 0.005))
    #             serialize_model(model)
    #         },
    #         epochs = 1000, ## Number of maximum training epochs. The training is however stopped early if the loss on the validation set does not improve for 20 epochs. This early stopping is hardcoded in fitter_nn.
    #         validation_split = 0.2, ## Fraction of the training data used to monitor validation loss
    #         verbose = 0,
    #         batch_size = 128 ## Size of the minibatches for training.
    # ),

    # Passed to the third element, SVMs. See help(svm, "e1071") for possible arguments
    list(type = "nu-regression", cost = 8, nu=0.5, kernel="radial"),

    # Passed to the fourth element, fitter_glmnet. This should contain a mandatory argument `degree` which specifies the degree of the polynomial model (1 for linear, 2 for quadratic etc...). Here we use degree = 2 corresponding to our LASSO2 model Other arguments are passed to getS3method("cv.glmnet", "formula"),
    list(alpha = 1, nfolds=10, degree = 2),

    # Passed to the fourth element, fitter_linear. This only accepts a degree argument specifying the degree of the polynomial model. Here we use degree = 1 corresponding to a linear model.
    list(degree = 1)
)

## ---- pipeline execution, eval = TRUE-----------------------------------------
if(length(regression_functions) != length(extra_args_regression_params)){
    stop("Number of models and number of lists of hyperparameters mismatch")
}
imputed_data <- infinity_flow(
	regression_functions = regression_functions,
	extra_args_regression_params = extra_args_regression_params,
	path_to_fcs = path_to_fcs,
	path_to_output = path_to_output,
	path_to_intermediary_results = path_to_intermediary_results,
	backbone_selection_file = backbone_selection_file,
	annotation = targets,
	isotype = isotypes,
	input_events_downsampling = input_events_downsampling,
	prediction_events_downsampling = prediction_events_downsampling,
	verbose = TRUE,
	cores = cores
)

## ---- output------------------------------------------------------------------
   print(imputed_data$bgc[1:2, ])

## ---- cores-------------------------------------------------------------------
cores = 1L

## ---- eval = FALSE------------------------------------------------------------
#  optional_dependencies <- c("keras", "tensorflow")
#  unmet_dependencies <- setdiff(optional_dependencies, rownames(installed.packages()))
#  if(length(unmet_dependencies) > 0){
#       install.packages(unmet_dependencies)
#  }
#  for(pkg in optional_dependencies){
#      library(pkg, character.only = TRUE)
#  }
#  
#  invisible(eval(try(keras_model_sequential()))) ## avoids conflicts with flowCore...
#  
#  if(!is_keras_available()){
#       install_keras() ## Instal keras unsing the R interface - can take a while
#  }
#  
#  if (!requireNamespace("BiocManager", quietly = TRUE)){
#      install.packages("BiocManager")
#  }
#  BiocManager::install("infinityFlow")
#  
#  library(infinityFlow)
#  
#  data(steady_state_lung)
#  data(steady_state_lung_annotation)
#  data(steady_state_lung_backbone_specification)
#  
#  dir <- file.path(tempdir(), "infinity_flow_example")
#  input_dir <- file.path(dir, "fcs")
#  write.flowSet(steady_state_lung, outdir = input_dir)
#  
#  write.csv(steady_state_lung_backbone_specification, file = file.path(dir, "backbone_selection_file.csv"), row.names = FALSE)
#  
#  path_to_fcs <- file.path(dir, "fcs")
#  path_to_output <- file.path(dir, "output")
#  path_to_intermediary_results <- file.path(dir, "tmp")
#  backbone_selection_file <- file.path(dir, "backbone_selection_file.csv")
#  
#  targets <- steady_state_lung_annotation$Infinity_target
#  names(targets) <- rownames(steady_state_lung_annotation)
#  isotypes <- steady_state_lung_annotation$Infinity_isotype
#  names(isotypes) <- rownames(steady_state_lung_annotation)
#  
#  input_events_downsampling <- 1000
#  prediction_events_downsampling <- 500
#  
#  ## Passed to fitter_nn, e.g. neural networks through keras::fit. See https://keras.rstudio.com/articles/tutorial_basic_regression.html
#  regression_functions <- list(NN = fitter_nn)
#  
#  backbone_size <- table(read.csv(backbone_selection_file)[,"type"])["backbone"]
#  extra_args_regression_params <- list(
#          list(
#  		object = { ## Specifies the network's architecture, loss function and optimization method
#  		model = keras_model_sequential()
#  		model %>%
#  		layer_dense(units = backbone_size, activation = "relu", input_shape = backbone_size) %>%
#  		layer_dense(units = backbone_size, activation = "relu", input_shape = backbone_size) %>%
#  		layer_dense(units = 1, activation = "linear")
#  		model %>%
#  		compile(loss = "mean_squared_error", optimizer = optimizer_sgd(lr = 0.005))
#  		serialize_model(model)
#  		},
#  		epochs = 1000, ## Number of maximum training epochs. The training is however stopped early if the loss on the validation set does not improve for 20 epochs. This early stopping is hardcoded in fitter_nn.
#  		validation_split = 0.2, ## Fraction of the training data used to monitor validation loss
#  		verbose = 0,
#  		batch_size = 128 ## Size of the minibatches for training.
#  	)
#  )
#  
#  imputed_data <- infinity_flow(
#  	regression_functions = regression_functions,
#  	extra_args_regression_params = extra_args_regression_params,
#  	path_to_fcs = path_to_fcs,
#  	path_to_output = path_to_output,
#  	path_to_intermediary_results = path_to_intermediary_results,
#  	backbone_selection_file = backbone_selection_file,
#  	annotation = targets,
#  	isotype = isotypes,
#  	input_events_downsampling = input_events_downsampling,
#  	prediction_events_downsampling = prediction_events_downsampling,
#  	verbose = TRUE,
#  	cores = 1L
#  )

## -----------------------------------------------------------------------------
sessionInfo()

