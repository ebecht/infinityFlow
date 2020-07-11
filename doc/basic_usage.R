## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- installation, eval=FALSE------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE)){
#      install.packages("BiocManager")
#  }
#  BiocManager::install("infinityFlow")

## ---- load_package, eval = TRUE-----------------------------------------------
library(infinityFlow)
data(steady_state_lung)

## ---- load_example, eval=TRUE-------------------------------------------------
dir <- file.path(tempdir(), "infinity_flow_example")
print(dir)
input_dir <- file.path(dir, "fcs")
write.flowSet(steady_state_lung, outdir = input_dir) ## Omit this if you already have FCS files
list.files(input_dir)

## ---- load_annotation---------------------------------------------------------
data(steady_state_lung_annotation)
print(steady_state_lung_annotation)

## ---- eval = FALSE------------------------------------------------------------
#  backbone_specification <- select_backbone_and_exploratory_markers(list.files(input_dir, pattern = ".fcs", full.names = TRUE))

## ---- backbone specification input--------------------------------------------
data(steady_state_lung_backbone_specification)
print(head(steady_state_lung_backbone_specification))

## ---- backbone specification output-------------------------------------------
write.csv(steady_state_lung_backbone_specification, file = file.path(dir, "backbone_selection_file.csv"), row.names = FALSE)

## ---- inspect input directory-------------------------------------------------
list.files(dir)

## ---- input FCS files path argument-------------------------------------------
path_to_fcs <- file.path(dir, "fcs")
head(list.files(path_to_fcs, pattern = ".fcs"))

## ---- output path argument----------------------------------------------------
path_to_output <- file.path(dir, "output")

## ---- backbone selection file path argument-----------------------------------
list.files(dir)
backbone_selection_file <- file.path(dir, "backbone_selection_file.csv")
head(read.csv(backbone_selection_file))

## ---- targets and isotypes arguments------------------------------------------
targets <- steady_state_lung_annotation$Infinity_target
names(targets) <- rownames(steady_state_lung_annotation)
isotypes <- steady_state_lung_annotation$Infinity_isotype
names(isotypes) <- rownames(steady_state_lung_annotation)
head(targets)
head(isotypes)

## ---- input and output events downsampling argument---------------------------
input_events_downsampling <- 1000
prediction_events_downsampling <- 500
cores = 1L

## ---- input temporary directory path------------------------------------------
path_to_intermediary_results <- file.path(dir, "tmp")

## ---- pipeline execution, eval = TRUE-----------------------------------------
imputed_data <- infinity_flow(
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

## -----------------------------------------------------------------------------
head(list.files(path_to_fcs)) ## Input files
fcs_raw <- file.path(path_to_output, "FCS", "split")
head(list.files(fcs_raw)) ## Raw output FCS files
fcs_bgc <- file.path(path_to_output, "FCS_background_corrected", "split") ## Background-corrected output FCS files
head(list.files(fcs_bgc)) ## Background-corrected output FCS files

## -----------------------------------------------------------------------------
file.path(path_to_output, "umap_plot_annotated.pdf") ## Raw plot
file.path(path_to_output, "umap_plot_annotated_backgroundcorrected.pdf") ## Background-corrected plot

## ---- eval = FALSE, echo = FALSE----------------------------------------------
#  knitr::include_graphics(file.path(path_to_output, "umap_plot_annotated.pdf"))

## -----------------------------------------------------------------------------
sessionInfo()

## ---- debugging, eval = FALSE, echo = FALSE-----------------------------------
#  files = list.files(dir, recursive = TRUE)
#  sapply(
#  	files,
#  	function(x){
#  		file.copy(from = file.path(dir, x), to = "~/Desktop/test/", recursive = TRUE)
#  	}
#  )

