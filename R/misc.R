#' Generates description for a flowFrame
#'
#' @param ff A flowframe
#' @importFrom Biobase pData
#' @return A flowframe with description consistent with pData(ff@parameters)
#' @noRd

generate_description<-function(ff){
    res=ff
    res.desc=pData(res@parameters)
    new.names=sprintf("$P%sS",1:nrow(res.desc))
    new.values=setNames(res.desc[,"desc"],new.names)
    new.values[is.na(new.values)]=" "
    new.values=as.list(new.values)
    res@description[names(res@description)%in%names(new.values)]=NULL
    res@description=c(res@description,new.values)

    new.names=sprintf("$P%sR",1:nrow(res.desc))
    new.values=setNames(ceiling(res.desc[,"maxRange"]-res.desc[,"minRange"]),new.names)
    new.values[is.na(new.values)]=" "
    new.values=as.list(new.values)
    res@description[names(res@description)%in%names(new.values)]=NULL
    res@description=c(res@description,new.values)

    new.names=sprintf("$P%sE",1:nrow(res.desc))
    new.values=setNames(rep("0,0",length(new.names)),new.names)
    res@description[names(res@description)%in%names(new.values)]=NULL
    res@description=c(res@description,new.values)

    new.names=sprintf("$P%sB",1:nrow(res.desc))
    new.values=setNames(rep("32",length(new.names)),new.names)
    res@description[names(res@description)%in%names(new.values)]=NULL
    res@description=c(res@description,new.values)

    new.names=sprintf("$P%sG",1:nrow(res.desc))
    new.values=setNames(rep("1",length(new.names)),new.names)
    res@description[names(res@description)%in%names(new.values)]=NULL
    res@description=c(res@description,new.values)
    res
}

#' Colors points of a biplot (2d-tSNE, 2d-PCA...) by the intensity of channels for each flowframe in the flowset
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
#' @note Since pdf files are vectorized, they can get really big if a lot of data point are plotted. This function thus used bitmap images that are stored in a temporary directory (tmpDir()) and then import them in a single pdf. If you're interested in using the bitmap images, you can fetch them in tmpDir()
#' @noRd

color_biplot_by_channels <- function(
                                     matrix,
                                     x_axis,
                                     y_axis,
                                     global_across_channels=T,
                                     palette=c("blue","green","red"),
                                     resolution=72,
                                     data_transformation_reverse=identity,
                                     file_name="biplot.pdf",
                                     raster.height=480,
                                     raster.width=480,
                                     ... #pass to plot for e.g. tSNE biplot
                                     )
{
    regular_channels=setdiff(colnames(matrix),c(x_axis,y_axis))

    if(global_across_channels){
        data_range=range(matrix[,regular_channels],na.rm=T)
        data_range=matrix(rep(data_range,length(regular_channels)),ncol=length(regular_channels),nrow=2,byrow=F,dimnames=list(c("min","max"),regular_channels))
    } else {
        data_range=apply(matrix[,regular_channels],na.rm=T,2,range,na.rm=T)
        rownames(data_range)=c("min","max")
    }

    x = matrix[,x_axis]
    y = matrix[,y_axis]
    xp = matrix[,regular_channels,drop=FALSE]

    if(any(!(is.na(x)&is.na(y)))){
        rasters=sapply(regular_channels,function(pname,xp,data.range,x,y)
        {
            color.scale=unique(colorRampPalette(palette)(1000))
            n=length(color.scale)

            breaks=unique(seq(data.range["min",pname],data.range["max",pname],length.out=n+1))
            if(length(unique(breaks))>1){
                points.colors=as.character(cut(xp[,pname],breaks=breaks,labels=color.scale))
            } else {
                points.colors=rep("lightsteelblue",length(xp[,pname]))
            }
            mainplot=paste(tmpDir(),"/mainplot_",pname,".png",sep="")
            png(mainplot,res=resolution,height=raster.height*resolution/72,width=raster.width*resolution/72)
            par("bty"="l")
            plot(
                x,
                y,
                col=points.colors,
                xlab=x_axis,
                ylab=y_axis,
                main=pname,
                ...
            )
            dev.off()

            colorscale=paste(tmpDir(),"/colorscale_",pname,".png",sep="")
            png(colorscale,res=resolution,height=raster.height/2*resolution/72,width=raster.width*resolution/72)
            plot.new()
            par("mar"=c(2,1,2,1))


            xlims=par("usr")[1:2]
            x_coords=seq(xlims[1],xlims[2],length.out=n+1)
            ylims=par("usr")[3:4]

            labels=signif(data_transformation_reverse(breaks),2)
            labels=labels[round(seq(1,length(labels),length.out=5))]
            labels.x_coords=seq(x_coords[1],x_coords[length(x_coords)],length.out=5)

            rect(border=NA,ybottom=ylims[1],ytop=ylims[2],xleft=x_coords[-length(x_coords)],xright=x_coords[-1],col=color.scale)
            text(xpd=T,y=ylims[1],pos=1,labels=labels,x=labels.x_coords)
            text(xpd=T,y=ylims[2],pos=3,labels=paste(pname,"intensity"),x=mean(xlims))
            dev.off()

            return(list(main.file=mainplot,scale.file=colorscale))
        },simplify=F,xp=xp,data.range=data_range,x=x,y=y)

        file=file_name

        pdf(file)
        if(global_across_channels){
            plot.new()
            grid.raster(readPNG(rasters[[1]]$scale.file,native=T))
        }
        sapply(rasters,function(x){
            par("mar"=c(0,0,0,0))
            grid.newpage()
            label=sub(".png","",sub("mainplot_","",tail(strsplit(x$main.file,"/")[[1]],1),fixed=TRUE),fixed=TRUE)
            
            if(!global_across_channels){
                grid.raster(readPNG(x$main.file,native=T),y=0.6,height=0.8)
                grid.raster(readPNG(x$scale.file,native=T),y=0.1,height=0.2)
            }
            if(global_across_channels){
                grid.raster(readPNG(x$main.file,native=T))
            }
            grid.text(x=unit(1,"npc"),y=unit(1,"npc"),label=label,just=c(1,1),gp=gpar(col="white",cex=0.1))
            return(NULL)
        })
        dev.off()
    }
    NULL
}

#' @title For each parameter in the FCS files, interactively specify whether it is part of the backbone, the exploratory markers or should be ignored.
#' @description This This function will load the first of the input FCS files and extract the measured parameters as well as their labels. For each of these, it will ask the user whether it is part of the backbone measurements (and will be used as a predictor variable in regressions), exploratory measurements (used as dependent variable in regressions) or discarded (e.g. for Time, Sample IDs, Event number IDs, ...).
#' @param files character vector of paths to FCS files
#' @export
#' @importFrom flowCore read.FCS parameters
#' @importFrom Biobase pData
select_backbone_and_exploratory_markers=function(files){
    if(!interactive()){
        stop("Interactive backbone selection requires an interactive R session")
    } else {
        requireNamespace("flowCore")
        representative_file=read.FCS(files[1],truncate_max_range=FALSE,ignore.text.offset=TRUE)
        data_channels=pData(parameters(representative_file)[,c("name","desc")])
        cat("For each data channel, enter either: backbone, exploratory or discard (can be abbreviated)\n")
        choices=c("backbone","exploratory","discard")
        result=cbind(data_channels,type=NA,stringsAsFactors=FALSE)
        for(i in 1:nrow(data_channels)){
            user_choice=NA
            while(is.na(user_choice)){
                x=data_channels[i,,drop=TRUE]
                user_choice=readline(paste0(ifelse(!is.na(x["desc"]),x["desc"],x["name"])," (",x["name"],"):"))
                user_choice=pmatch(user_choice,choices)
                if(is.na(user_choice))
                    cat("Incorrect selection, enter either: backbone, exploratory or discard\n")
            }
            result[i,"type"]=choices[user_choice]
        }

        user_choice=NA
        print(result)
        
        while(is.na(user_choice)){
            correct=readline("Is selection correct? (yes/no): ")
            choices=c("yes","no")
            user_choice=choices[pmatch(correct,choices)]
            if(is.na(user_choice)){
                print(result)
            } else if(user_choice=="no"){
                return(select_backbone_and_exploratory_markers(files=files))
            } else if(user_choice=="yes"){
                if(!any(result$type=="exploratory")){
                    cat("At least one measurement must be exploratory\n")
                    return(select_backbone_and_exploratory_markers(files=files))
                }
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
split_matrix=function (mat, vector, byrow = TRUE) 
{
    if (byrow & nrow(mat) != length(vector)) {
        stop("if byrow=TRUE, vector's length should have length nrow(mat)")
    }
    else if (!byrow & ncol(mat) != length(vector)) {
        !byrow & ncol(mat) != length(vector)
        stop("if byrow=FALSE, vector's length should have length ncol(mat)")
    }
    if (byrow) {
        levels = split(1:nrow(mat), vector)
        res = lapply(levels, function(x) {
            mat[x, , drop = FALSE]
        })
    }
    else {
        levels = split(1:ncol(mat), vector)
        res = lapply(levels, function(x) {
            mat[, x, drop = FALSE]
        })
    }
    res
}

#' Linear scale with chosen boundaries
#' @noRd
minmax_scale=function(matrix,min=1,max=1000,na.rm=TRUE){
    apply(matrix,2,function(x){
        a=(max-min)/(max(x,na.rm=na.rm)-min(x,na.rm=na.rm))
        b=min-a*min(x, na.rm=na.rm)
        a*x+b
    })
}
