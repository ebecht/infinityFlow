#' Generates description for a flowFrame
#'
#' @param ff A flowframe
#' @importFrom Biobase pData
#' @return A flowframe with description consistent with pData(ff@parameters)
#' @noRd

generate_description<-function(ff){
    res <- ff
    res.desc <- pData(res@parameters)
    new.names <- sprintf("$P%sS",seq_len(nrow(res.desc)))
    new.values <- setNames(res.desc[,"desc"],new.names)
    new.values[is.na(new.values)] <- " "
    new.values <- as.list(new.values)
    res@description[names(res@description)%in%names(new.values)] <- NULL
    res@description <- c(res@description,new.values)

    new.names <- sprintf("$P%sR",seq_len(nrow(res.desc)))
    new.values <- setNames(ceiling(res.desc[,"maxRange"]-res.desc[,"minRange"]),new.names)
    new.values[is.na(new.values)] <- " "
    new.values <- as.list(new.values)
    res@description[names(res@description)%in%names(new.values)] <- NULL
    res@description <- c(res@description,new.values)

    new.names <- sprintf("$P%sE",seq_len(nrow(res.desc)))
    new.values <- setNames(rep("0,0",length(new.names)),new.names)
    res@description[names(res@description)%in%names(new.values)] <- NULL
    res@description <- c(res@description,new.values)

    new.names <- sprintf("$P%sB",seq_len(nrow(res.desc)))
    new.values <- setNames(rep("32",length(new.names)),new.names)
    res@description[names(res@description)%in%names(new.values)] <- NULL
    res@description <- c(res@description,new.values)

    new.names <- sprintf("$P%sG",seq_len(nrow(res.desc)))
    new.values <- setNames(rep("1",length(new.names)),new.names)
    res@description[names(res@description)%in%names(new.values)] <- NULL
    res@description <- c(res@description,new.values)
    res
}

#' Colors points of a biplot (2d-tSNE, 2d-PCA...) by the intensity of channels for each backbone and exploratory marker
#'
#' @param matrix A matrix
#' @param x_axis A column name of matrix used in biplots as the x axis
#' @param y_axis A column name of matrix used in biplots as the y axis
#' @param global_across_channels Boolean specificying whether the color key should be calibrated across all channels or individually per channel.
#' @param palette A vector of colors that'll be passed to colorRampPalette
#' @param resolution The resolution of the files that'll be imported in the pdf. Default is 72, increase for higher resolution. You may need to enlarge rasters accordingly.
#' @param data_transformation_reverse The colors will be linearly scaled across the range of the data. If the data was transformed, you may however want the labels which appear in your color-keys to reflect the raw intensities. In this case, this should be the inverse of the function you used to transform your data
#' @param file_name String of the location of the output file (should end in .pdf)
#' @param raster.width Width of each embedded raster. Defaults to 480
#' @param raster.height height of each embedded raster. Defaults to 480
#' @param ... passed to plot (suggested: pch=16, cex=0.5 or less)
#' @return NULL
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom raster tmpDir
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics plot.new
#' @importFrom graphics rect
#' @importFrom graphics text
#' @importFrom grid grid.raster
#' @importFrom grid grid.newpage
#' @importFrom grid grid.text
#' @importFrom grid unit
#' @importFrom grid gpar
#' @importFrom png readPNG
#' @importFrom utils tail
#' @importFrom parallel makeCluster clusterExport clusterEvalQ stopCluster
#' @importFrom pbapply pblapply
#' @note Since pdf files are vectorized, they can get really big if a lot of data point are plotted. This function thus used bitmap images that are stored in a temporary directory (tmpDir()) and then import them in a single pdf. If you're interested in using the bitmap images, you can fetch them in tmpDir()
#' @noRd

color_biplot_by_channels <- function(
                                     h5file,
                                     umap_group = "/umap/backbone/",
                                     backbone_data_group = "/input/expression_transformed/",
                                     predictions_group = "/predictions/raw/",
                                     annot,
                                     palette=c("blue","green","red"),
                                     resolution=72,
                                     transforms = transforms,
                                     chans = chans,
                                     yvar = yvar,
                                     ## data_transformation_reverse=transforms[[yvar]]$backward,
                                     file_name="biplot.pdf",
                                     raster.height=480,
                                     raster.width=480,
                                     chop_quantiles,
                                     cores,
                                     ... #pass to plot for e.g. tSNE biplot
                                     ){
    umap <- list()
    for(i in seq_len(nrow(annot))){
        umap[[i]] <- h5read(h5file, name = paste0(umap_group, i))
    }
    umap <- do.call(rbind, umap)
    colnames(umap) = paste0("UMAP", c(1,2))
    spl <- sample(seq_len(nrow(umap)))
    umap <- umap[spl, ]

    x <- umap[,1]
    y <- umap[,2]
    
    regular_channels <- h5readAttributes(h5file, paste0(predictions_group, "1"))$colnames
    regular_channels = c(chans, regular_channels)
    rasters_list <- list()

    color.scale <- unique(colorRampPalette(palette)(100))
    n <- length(color.scale)
    
    n_backbone <- length(chans)

    tmpdir <- eval(tmpDir())
    
    cl <- makeCluster(cores)
    env <- environment()
    clusterExport(
        cl,
        c("n_backbone", "chans", "h5file", "annot", "chop_quantiles", "color.scale", "tmpdir", "resolution", "raster.height", "raster.width", "x", "y", "regular_channels", "transforms", "spl")
      , envir = env
    )
    
    ##for(i in seq_len(n_backbone + nrow(annot))){
    rasters_list <- pblapply(
        seq_len(n_backbone + nrow(annot)),
        cl = cl,
        function(i){
            preds <- list()
            if(i < n_backbone + 1){
                col_index <- match(chans[i], h5readAttributes(h5file, "/input/expression/1")$colnames)
                for(j in seq_len(nrow(annot))){
                    w <- h5read(h5file, name = paste0("/sampling/predictions/", j)) == 1L
                    preds[[j]] <- h5read(h5file, name = paste0(backbone_data_group, j))[w, col_index]
                }
            } else {
                for(j in seq_len(nrow(annot))){
                    preds[[j]] <- h5read(h5file, name = paste0(predictions_group, j), index = list(NULL, i - n_backbone))
                }
            }
            preds <- do.call(c, preds)
            preds <- preds[spl]

            q <- quantile(preds,c(chop_quantiles,1-chop_quantiles))
            preds[preds<=q[1]] <- q[1]
            preds[preds>=q[2]] <- q[2]

            range <- range(preds)
            breaks <- unique(seq(range[1],range[2],length.out=n+1))

            if(length(unique(breaks))>1){
                points.colors <- as.character(cut(preds,breaks=breaks,labels=color.scale))
            } else {
                points.colors <- rep("lightsteelblue",length(preds))
            }
            mainplot <- paste(tmpdir,"/mainplot_",i,".png",sep="")
            png(mainplot,res=resolution,height=raster.height*resolution/72,width=raster.width*resolution/72)
            par("bty"="l")
            plot(
                x,
                y,
                col=points.colors,
                xlab="UMAP1",
                ylab="UMAP2",
                main=regular_channels[i],
                ...
            )
            dev.off()

            colorscale <- paste(tmpdir,"/colorscale_",i,".png",sep="")
            png(colorscale,res=resolution,height=raster.height/2*resolution/72,width=raster.width*resolution/72)
            plot.new()
            par("mar"=c(2,1,2,1))

            xlims <- par("usr")[c(1, 2)]
            x_coords <- seq(xlims[1],xlims[2],length.out=n+1)
            ylims <- par("usr")[c(3, 4)]

            if(i < n_backbone + 1){
                breaks_native <- transforms[[chans[i]]]$backward(breaks)
            } else {
                breaks_native <- transforms[[yvar]]$backward(breaks)
            }
            labels <- signif(breaks_native,2)
            labels <- labels[round(seq(1,length(labels),length.out=5))]
            labels.x_coords <- seq(x_coords[1],x_coords[length(x_coords)],length.out=5)

            rect(border=NA,ybottom=ylims[1],ytop=ylims[2],xleft=x_coords[-length(x_coords)],xright=x_coords[-1],col=color.scale)
            text(xpd=TRUE,y=ylims[1],pos=1,labels=labels,x=labels.x_coords)
            text(xpd=TRUE,y=ylims[2],pos=3,labels=paste(annot[i, "target"],"intensity"),x=mean(xlims))
            dev.off()

            ## rasters_list <- c(rasters_list, list(list(mainplot = mainplot, colorscale = colorscale)))
            return(list(mainplot = mainplot, colorscale = colorscale))
        }
    )
    
    stopCluster(cl)
    pdf(file_name)
    for(i in seq_along(rasters_list)){
        par("mar"=c(0,0,0,0))
        grid.newpage()
        
        grid.raster(readPNG(rasters_list[[i]]$mainplot,native=TRUE),y=0.6,height=0.8)
        grid.raster(readPNG(rasters_list[[i]]$colorscale,native=TRUE),y=0.1,height=0.2)
        
        grid.text(x=unit(1,"npc"),y=unit(1,"npc"),label=ifelse(i < n_backbone + 1, chans[i], annot[i - n_backbone, "target"]),just=c(1,1),gp=gpar(col="white",cex=0.1))
    }
    dev.off()
    
    invisible()
}

#' @title For each parameter in the FCS files, interactively prompts whether it is part of the Backbone, the Infinity (exploratory) markers or should be ignored.
#' @description This function will load the first of the input FCS files and extract the measured parameters as well as their labels. For each of these, it will ask the user whether it is part of the backbone measurements (which will be used as a predictor variable in regressions models), Infinity (exploratory) measurements (usually PE-conjugated or APC-conjugated, used as dependent/target variable in regressions) or discarded (e.g. for parameter such as Time, Sample IDs, Event number IDs, ...).
#' @param files character vector of paths to FCS files
#' @export
#' @importFrom flowCore read.FCS parameters
#' @importFrom Biobase pData
#' @importFrom methods is
#' @return A data.frame
#' @examples
#' data(steady_state_lung)
#' dir  <-  tempdir()
#' fcs_tmp  <-  file.path(dir, "tmp.fcs")
#' library(flowCore)
#' write.FCS(steady_state_lung[[1]], file  <-  fcs_tmp)
#' if(interactive()){
#'     select_backbone_and_exploratory_markers(fcs_tmp)
#' }

select_backbone_and_exploratory_markers <- function(files){
    if(!interactive()){
        stop("Interactive backbone selection requires an interactive R session")
    } else {
        requireNamespace("flowCore")
        representative_file <- read.FCS(files[1],truncate_max_range=FALSE,ignore.text.offset=TRUE)
        data_channels <- pData(parameters(representative_file)[,c("name","desc")])
        cat("For each data channel, enter either: backbone, exploratory or discard (can be abbreviated)\n")
        choices <- c("backbone","exploratory","discard")
        result <- cbind(data_channels,type=NA,stringsAsFactors=FALSE)
        for(i in seq_len(nrow(data_channels))){
            user_choice <- NA
            while(is.na(user_choice)){
                x <- data_channels[i,,drop=TRUE]
                user_choice <- readline(paste0(ifelse(!is.na(x["desc"]),x["desc"],x["name"])," (",x["name"],"):"))
                user_choice <- pmatch(user_choice,choices)
                if(is.na(user_choice))
                    cat("Incorrect selection, enter either: backbone, exploratory or discard\n")
            }
            result[i,"type"] <- choices[user_choice]
        }

        ## user_choice <- NA
        ## co_factor <- 1000
        ## result = cbind(result, cofactor = ifelse(!is.na(result$desc), co_factor, NA))
        ## while(is.na(user_choice)){
        ##     correct = readline(paste0("Does a co-factor of ", co_factor, " in the asinh transformation of the data suit you? (yes/no)"))
        ##     choices <- c("yes","no")
        ##     user_choice <- choices[pmatch(correct,choices)]
        ##     if(!is.na(user_choice)){
        ##         if(user_choice == "no"){
        ##             for(i in which(!is.na(data_channels$desc))){
        ##                 user_choice <- NA
        ##                 while(is.na(user_choice)){
        ##                     x <- data_channels[i,,drop=TRUE]
        ##                     user_choice <- suppressWarnings(as.numeric(readline(paste0("Enter cofactor value for ", ifelse(!is.na(x["desc"]),x["desc"],x["name"])," (",x["name"],"):"))))
        ##                     if(is.na(user_choice)){
        ##                         cat("Incorrect selection, enter a number\n")
        ##                     }
                            
        ##                 }
        ##                 result[i, "cofactor"] <- user_choice
        ##             }
        ##         }
        ##     }
        ## }
        
        user_choice <- NA
        while(is.na(user_choice)){
            print(result)
            correct <- readline("Is selection correct? (yes/no): ")
            choices <- c("yes","no")
            user_choice <- choices[pmatch(correct,choices)]
            if(is.na(user_choice)){
                print(result)
            } else if(user_choice=="no"){
                return(select_backbone_and_exploratory_markers(files=files))
            } else if(user_choice=="yes"){
                if(sum(result$type=="exploratory")!=1){
                    cat("At least one measurement must be exploratory\n")
                    return(select_backbone_and_exploratory_markers(files=files))
                }

                result$transformation = "asinh"
                result$cofactor = 1000
                w = grepl("^[FS]SC-?[AWH]?", result$name) | tolower(result$name) == "time"
                result[w, "transformation"] = "identity"
                result[w, "cofactor"] = 1

                cat("\n Auto-detecting channels for transformation. You can edit identity ('linear') / arcsinh manually\n")
                print(result)
                return(result)
            }
        }
    }
}

#' Split a matrix into a list of chunks. Faster than using split on a data.frame
#' @param mat A matrix
#' @param vector A vector of length nrow(mat) if byrow=TRUE, ncol(mat) if byrow=FALSE
#' @param byrow if TRUE split rows, if FALSE split columns
#' @noRd
split_matrix <- function (mat, vector, byrow = TRUE) 
{
    if (byrow & nrow(mat) != length(vector)) {
        stop("if byrow=TRUE, vector's length should have length nrow(mat)")
    }
    else if (!byrow & ncol(mat) != length(vector)) {
        !byrow & ncol(mat) != length(vector)
        stop("if byrow=FALSE, vector's length should have length ncol(mat)")
    }
    if (byrow) {
        levels  <-  split(seq_len(nrow(mat)), vector)
        res  <-  lapply(levels, function(x) {
            mat[x, , drop = FALSE]
        })
    }
    else {
        levels  <-  split(seq_len(ncol(mat)), vector)
        res  <-  lapply(levels, function(x) {
            mat[, x, drop = FALSE]
        })
    }
    res
}

#' Linear scale with chosen boundaries
#' @noRd
minmax_scale <- function(matrix,min=1,max=1000,na.rm=TRUE){
    apply(matrix,2,function(x){
        a <- (max-min)/(max(x,na.rm=na.rm)-min(x,na.rm=na.rm))
        b <- min-a*min(x, na.rm=na.rm)
        a*x+b
    })
}

#' Plot a 2D heatmap of frequencies using rectangular bins, mapping log10(counts per bin) to colors
#' @param x numeric vector for x-axis
#' @param y numeric vector for y-axis
#' @param breaks integer(1) number of horizontal and vertical breaks
#' @param na.rm Should we discard pairs of (x[i], y[i]) if either x[i] or y[i] is NA
#' @param palette Color palette encoding the 2D log10-density
#' @param add_white Should bins with no count be plotted in white (default : TRUE)
#' @param ... passed to image()
#' @noRd
freqplot <- function(
                     x,
                     y,
                     breaks=200,
                     na.rm=TRUE,
                     palette = rev(c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695")),
                     add_white=TRUE,
                     ...
                     ){

    w<-is.na(x)|is.na(y)|is.nan(x)|is.nan(y)|is.infinite(x)|is.infinite(y)
    if(any(w)){
        if(na.rm){
            x<-x[!w]
            y<-y[!w]
        } else {
            stop("NA values found and na.rm is FALSE")
        }
    }
    
    w.x<-length(unique(x))>1
    w.y<-length(unique(y))>1
    
    if(w.x){
        breaks.x<-seq(min(x),max(x),length.out=breaks)
        labels.x<-breaks.x[-length(breaks.x)]
        X<-cut(x,breaks=breaks.x,labels=labels.x,include.lowest=TRUE)
    } else {
        X<-x
    }
    if(w.y){
        breaks.y<-seq(min(y),max(y),length.out=breaks)
        labels.y<-breaks.y[-length(breaks.y)]
        Y<-cut(y,breaks=breaks.y,labels=labels.y,include.lowest=TRUE)
    } else {
        Y<-y
    }

    tab<-log10(1+table(X,Y))

    if(length(x)<1|length(y)<1){
        plot.new()
        return(tab)
    }
    
    if(w.x&w.y){
        if(add_white){
            null_color = "white"
        } else {
            null_color = NULL
        }
        image(tab,col=c(null_color,colorRampPalette(palette)(100)),x=breaks.x,y=breaks.y,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",...)
        ticks<-seq(0,1,by=0.25)
        if(par("xaxt") != "n"){
            axis(side=1,at=quantile(breaks.x,ticks),labels=signif(quantile(breaks.x,ticks),2),line=0.5)
        }
        if(par("yaxt") != "n"){
            axis(side=2,at=quantile(breaks.y,ticks),labels=signif(quantile(breaks.y,ticks),2),line=0.5)
        }
    } else {
        if(!w.x){
            X<-runif(length(x))
        } else {
            X<-x
        }
        if(!w.y){
            Y<-runif(length(y))
        } else {
            Y<-y
        }
        freqplot(X,Y,breaks=breaks,na.rm=na.rm,...)
    }
    tab
}

## Documentation of data objects

#' Subset of a massively parallel cytometry experiment of mouse lung single cells
#' @docType data
#' @format a flowSet containing 10 flowFrames (thus corresponding to 10 FCS files)
#' @source \url{https://flowrepository.org/id/FR-FCM-Z2LP}
#' @usage data(steady_state_lung)
"steady_state_lung"

#' Target and isotypes annotation for the data object infinityFlow::steady_state_lnug
#' @docType data
#' @format a data.frame specifying the Infinity antibody targets and isotypes for each flowFrame of the steady_state_lung flowSet
#' @usage data(steady_state_lung_annotation)
"steady_state_lung_annotation"

#' Backbone and Infinity antibodies specification for the data object infinityFlow::steady_state_lnug
#' @docType data
#' @format a data.frame specifying the Infinity antibody targets and isotypes for each flowFrame of the steady_state_lung flowSet
#' @usage data(steady_state_lung_backbone_specification) 
"steady_state_lung_backbone_specification"
