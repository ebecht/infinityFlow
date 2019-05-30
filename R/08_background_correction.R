require(flowCore)

#' @export
correct_background=function(
                            paths,
                            FCS_export,
                            CSV_export,
                            chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                            transforms_chan=readRDS(file.path(paths["rds"],"transforms_chan.Rds")),
                            transforms_pe=readRDS(file.path(paths["rds"],"transforms_pe.Rds")),
                            xp=readRDS(file.path(paths["rds"],"xp.Rds")),
                            xp_scaled=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds")),
                            umap=readRDS(file.path(paths["rds"],"umap.Rds")),
                            events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
                            sampling=readRDS(file.path(paths["rds"],"sampling_preds.Rds")),
                            preds=readRDS(file.path(paths["rds"],"predictions.Rds")),
                            prediction_colnames=readRDS(file.path(paths["rds"],"prediction_colnames.Rds")),
                            a=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE),
                            verbose=TRUE
                            ){

    if(verbose){
        message("Background correcting")
    }
    
    ## chans=readRDS(file.path(paths["rds"],"chans.Rds"));
    ## transforms_chan=readRDS(file.path(paths["rds"],"transforms_chan.Rds"));
    ## transforms_pe=readRDS(file.path(paths["rds"],"transforms_pe.Rds"));
    ## xp=readRDS(file.path(paths["rds"],"xp.Rds"));
    ## umap=readRDS(file.path(paths["rds"],"umap.Rds"));
    ## events.code=readRDS(file.path(paths["rds"],"pe.Rds"));
    ## sampling=readRDS(file.path(paths["rds"],"sampling_preds.Rds"));
    ## preds=readRDS(file.path(paths["rds"],"predictions.Rds"));
    ## xp_scaled=readRDS(file.path(paths["rds"],"xp_transformed_scaled.Rds"))
    ## CSV_export=TRUE;
    ## FCS_export=c("split","concatenated")

    a=read.csv(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE)
    
    rownames(a)=a[,"file"]
    a[,"target"]=gsub("/","-",a[,"target"])
    a[is.na(a$target),"target"]="Blank"
    if(any(a$target=="Blank")){
        a[a$target=="Blank","target"]=paste0("Blank",1:sum(a$target=="Blank"))
    }
    
    preds_raw=preds
    
    preds_rawbgc=lapply(
        preds_raw,
        function(x){
            x[,rownames(subset(a,isotype!="Auto"&!grepl("Isotype",target)))]
        }
    )

    for(i in seq_along(preds_rawbgc)){
        for(file in colnames(preds_rawbgc[[i]])){
            iso=a[a$target==paste0("Isotype_",a[file,"isotype"]),"file"]
            x=preds_raw[[i]][,iso]
            y=preds_raw[[i]][,file]
            lm=lm(y~x)$coefficients
            intercept=lm[1]
            slope=lm[2]
            orthogonal_residuals=(-slope*x+y-intercept)/sqrt(slope^2+1)
            preds_rawbgc[[i]][,file]=orthogonal_residuals
        }
    }

    saveRDS(preds_rawbgc,file=file.path(paths["rds"],"predictions_backgroundcorrected.Rds"))
    
    ## ##################
    ## Exporting Phenograph / BGC / UMAPs
    ## ##################

    if(verbose){
        message("\tTransforming background-corrected predictions. (Use logarithm to visualize)")
    }
    ## Preparing data for display using log transform
    preds_rawbgc_linear=lapply(
        preds_rawbgc,function(x){
            apply(x,2,function(x){
                10^(x-min(x))
            })
        }
    )
    
    ## Adding uncorrected data (isotypes and autofluorescence)
    for(i in seq_along(preds_rawbgc_linear)){
        preds_rawbgc_linear[[i]]=cbind(
            preds_rawbgc_linear[[i]],
            preds_raw[[i]][,setdiff(colnames(preds[[i]]),c(colnames(preds_rawbgc_linear[[i]]),chans)),drop=FALSE]
        )
        preds_rawbgc[[i]]=cbind(
            preds_rawbgc[[i]],
            preds_raw[[i]][,setdiff(colnames(preds[[i]]),c(colnames(preds_rawbgc[[i]]),chans)),drop=FALSE]
        )
    }
    
    preds_rawbgc_linear=lapply(
        preds_rawbgc_linear,
        function(x){
            x[,rownames(a)]
        }
    )
    preds_rawbgc=lapply(
        preds_rawbgc_linear,
        function(x){
            x[,rownames(a)]
        }
    )
    for(x in names(preds_rawbgc_linear)){
        colnames(preds_rawbgc_linear[[x]])=paste0(a[colnames(preds_rawbgc_linear[[x]]),"target"],".",x,"_bgc")
        colnames(preds_rawbgc[[x]])=paste0(a[colnames(preds_rawbgc[[x]]),"target"],".",x,"_bgc")
    }
    
    preds_rawbgc_linear=cbind(xp[sampling,],do.call(cbind,preds_rawbgc_linear),minmax_scale(umap))
    preds_rawbgc=cbind(xp_scaled[sampling,],do.call(cbind,preds_rawbgc),minmax_scale(umap))
    
    unique_pes=unique(events.code)
    PE_id=sapply(events.code[sampling],match,table=unique_pes)

    preds_rawbgc_linear=cbind(preds_rawbgc_linear,PE_id=PE_id)
    preds_rawbgc=cbind(preds_rawbgc,PE_id=PE_id)

    if(CSV_export){
        if(verbose){
            message("\t","Exporting as CSV")
        }
        dir.create(file.path(paths["output"],"FCS_background_corrected/"),showWarnings=FALSE,recursive=TRUE)
        require(data.table)
        fwrite(as.data.frame(preds_rawbgc_linear),row.names=FALSE,file=file.path(paths["output"],"predicted_data_background_corrected.csv"))
    }

    if(any(FCS_export=="split")){
        if(verbose){
            message("\t","Exporting FCS files (1 per well)")
        }
        dir.create(file.path(paths["output"],"FCS_background_corrected/","split"),recursive=TRUE,showWarnings=FALSE)
        FCS_list=split_matrix(preds_rawbgc_linear,events.code[sampling])
        FCS_list=lapply(FCS_list,function(x){
            x=flowFrame(x)
            x@parameters$desc=as.character(x@parameters$desc)
            x@parameters$name=as.character(x@parameters$name)
            generate_description(x)
        })
        lapply(
            names(FCS_list),
            function(x){
                invisible(write.FCS(FCS_list[[x]],file=file.path(paths["output"],"FCS_background_corrected/","split",paste0(sub(".fcs","",x),"_target_",a[x,"target"],".fcs"))))
            }
        )   
    }

    if(any(FCS_export=="concatenated")){
        if(verbose){
            message("\t","Exporting concatenated FCS file")
        }
        FCS=flowFrame(preds_rawbgc_linear)
        FCS@parameters$desc=as.character(FCS@parameters$desc)
        FCS@parameters$name=as.character(FCS@parameters$name)
        FCS=generate_description(FCS)
        dir.create(file.path(paths["output"],"FCS_background_corrected/","concatenated"),recursive=TRUE,showWarnings=FALSE)
        invisible(write.FCS(FCS,file=file.path(paths["output"],"FCS_background_corrected","concatenated","concatenated_results_background_corrected.fcs")))
    }

    saveRDS(preds_rawbgc,file.path(paths["rds"],"predictions_bgc_cbound.Rds"))
    preds_rawbgc
}

## ####################
## ## Three different UMAPs. Shows that we only lose resolution for one single cell subsets from the predictions.
## ## Also we kind of recovered one subset !
## ####################

## library(uwot)
## if(FALSE){
##     umap_bgc=umap(preds_bgc[,],n_neighbors=15L,min_dist=0.2,metric="euclidean",verbose=TRUE,n_epochs=1000L,n_threads=4L)
##     umap_preds=umap(apply(preds[,setdiff(colnames(preds),chans)],2,function(x){eb.autolgcl(x)(x)}),n_neighbors=15L,min_dist=0.2,metric="euclidean",verbose=TRUE,n_epochs=1000L,n_threads=4L)
##     colnames(umap_bgc)=paste0("UMAP_bgc",1:2)
##     colnames(umap_preds)=paste0("UMAP_preds",1:2)
##     save(umap_bgc,file="./rdata/umap_bgc.RData")
##     save(umap_preds,file="./rdata/umap_preds.RData")
## } else {
##     load("./rdata/umap_bgc.RData")
##     load("./rdata/umap_preds.RData")
## }

## png("./graphs/Comparison_3UMAPs.png",width=6000,height=2000,res=300)
## par(mfrow=c(1,3),bty="n")
## plot(umap,pch=16,cex=0.3,col=toColors_2d(umap,rank=FALSE))
## plot(umap_preds,pch=16,cex=0.3,col=toColors_2d(umap,rank=FALSE))
## plot(umap_bgc,pch=16,cex=0.3,col=toColors_2d(umap,rank=FALSE))
## dev.off()

## ####################
## ## Running Phenograph on backbone and BGC
## ####################

## library(Rphenograph)
## if(FALSE){
##     set.seed(123)
##     clust=membership(Rphenograph(preds[,chans])[[2]])
##     save(clust,file="./rdata/Phenograph.RData")

##     set.seed(123)
##     clust_bgc=membership(Rphenograph(preds_bgc[,setdiff(colnames(preds_bgc),chans)])[[2]]) ## This wasn't saved properly
##     save(clust_bgc,file="./rdata/Phenograph_bgc.RData")
## } else {
##     load("./rdata/Phenograph.RData")
##     load("./rdata/Phenograph_bgc.RData")
## }

## png("./graphs/Phenograph.png",width=6000,height=2000,res=300)
## par(bty="n")
## par(mfrow=c(1,3))
## plot(umap,pch=16,col=toColors_discrete(clust),cex=0.2,main="Backbone")
## coords=aggregate_df(umap,clust,margin=1,fun=median)
## text(labels=rownames(coords),x=coords[,1],y=coords[,2])

## plot(umap_preds,pch=16,col=toColors_discrete(clust),cex=0.2,main="Raw predictions")
## coords=aggregate_df(umap_preds,clust,margin=1,fun=median)
## text(labels=rownames(coords),x=coords[,1],y=coords[,2])

## plot(umap_bgc,pch=16,col=toColors_discrete(clust),cex=0.2,main="Background-corrected predictions")
## coords=aggregate_df(umap_bgc,clust,margin=1,fun=median)
## text(labels=rownames(coords),x=coords[,1],y=coords[,2])
## dev.off()

## ####################
## ## Summary figure for grant
## ####################

## png("./graphs/Figure_for_illustration.png",width=8000,height=3000,res=300)
## par(mar=c(4,4,2,2))
## layout(mat=matrix(nrow=1,data=c(1,2)),widths=c(3,5))
## plot(umap,pch=16,col=toColors_discrete(clust),cex=0.2)
## coords=aggregate_df(umap,clust,margin=1,fun=median)
## text(labels=rownames(coords),x=coords[,1],y=coords[,2])
## mfis=aggregate_df(preds_bgc,clust,median)
## colnames(mfis)=a[colnames(mfis),"target"]
## mfis=scale(mfis,center=TRUE,scale=FALSE)
## labCol=rep("",ncol(mfis))
## labCol[seq(1,ncol(mfis),by=20)]=paste0(" Expl. marker ",seq(1,ncol(mfis),by=20))
## labCol[seq(11,ncol(mfis),by=20)]="."
## labCol[seq(10,ncol(mfis),by=20)]="."
## labCol[seq(12,ncol(mfis),by=20)]="."
## hc=hclust(dist(t(mfis)))$order
## hcr=hclust(dist(mfis))$order
## mfis=mfis[hcr,hc]
## breaks=fastBreaks(mfis,101)
## library(RColorBrewer)
## ## mfis=t(mfis)
## ##labCol=substr(colnames(mfis),1,10)
## ##labCol=paste0(c("","           ","                      ","                                 ","                                            ","                                                       "),labCol)

## par(mar=c(4,4,8,2))
## imageWrapper(mfis,breaks=breaks,col=colorRampPalette(c("darkblue","blue","white","red","darkred"))(100),labCol=labCol,labRow=rownames(mfis))
## dev.off()

## ####################
## ## Dip tests
## ####################

## library(diptest)
## preds_split=split_matrix(preds_bgc,clust)
## dips=sapply(
##     preds_split,
##     function(x){
##         apply(
##             x,
##             2,
##             function(y){
##                 dip.test(y)$p.value
##             }
##         )
##     }
## )
## discoveries=which(dips<=1,arr.ind=TRUE) ## eheheh
## rownames(discoveries)=paste0("d",1:nrow(discoveries))

## pdf("~/test.pdf")
## apply(discoveries,1,function(x){
##     i=clust==x[2]
##     j=rownames(dips)[x[1]]
##     hist(preds_bgc[i,j],main=paste0(a[j,"target"]," in cluster ",x[2]),breaks=40,xlab="Predicted intensity",col="black")
##     p=dips[x[1],x[2]]
##     legend(x="topright",legend=paste0("p (Dip)=",signif(ifelse(p>0,p,.Machine$double.eps),2)),bty="n")
## })
## dev.off()

## ## Check which dip.test<0.05 also strongly correlate with measured data (make sure this is real signal)
## ## dips_cors=apply(discoveries,1,function(x){
## ##     i=clust==x[2]
## ##     j=rownames(dips)[x[1]]   
## ##     w=events.code[sampling]==j
## ##     ## hist(transforms_pe[[file]](xp[w,"FJComp-PE(yg)-A"]))
## ##     ## hist(transforms_pe[[file]](xp[,"FJComp-PE(yg)-A"][sampling[w]][clust[w]==x[2]]),breaks=20)
## ##     cor(
## ##         transforms_pe[[file]](xp[,"FJComp-PE(yg)-A"][sampling[w]][clust[w]==x[2]]),
## ##         preds_bgc[w,file][clust[w]==x[2]]
## ##     )
## ## })

## dips_data=apply(
##     discoveries,
##     1,
##     function(x){
##         i=clust==x[2]
##         j=rownames(dips)[x[1]]   
##         w=events.code[sampling]==j
##         ## hist(transforms_pe[[file]](xp[w,"FJComp-PE(yg)-A"]))
##         ## hist(transforms_pe[[file]](xp[,"FJComp-PE(yg)-A"][sampling[w]][clust[w]==x[2]]),breaks=20)
##         cbind(measured=transforms_pe[[file]](xp[,"FJComp-PE(yg)-A"][sampling[w]][clust[w]==x[2]]),predicted=preds_bgc[w,file][clust[w]==x[2]])
##     }
## )

## dips_cors=lapply(
##     dips_data,
##     cor
## )

## ## signif_dips=discoveries[order(abs(dips_cors)),]
## signif_dips=data.frame(discoveries,r=dips_cors,target=a[colnames(preds_bgc),"target"][signif_dips[,1]])
## signif_dips=signif_dips[order(signif_dips$r),]
## signif_dips=signif_dips[order(signif_dips[,"r"]),]

## png(width=2000,height=1000,res=300,"./graphs/Illustration_bimodality.png")
## par(mfrow=c(1,2))
## w=clust==21
## m1="CD144"
## m2="CD39"
## freqplot(breaks=50,preds_bgc[w,subset(a,target==m1)$file],preds_bgc[w,subset(a,target==m2)$file])
## title(xlab=m1,ylab=m2,main="Cluster 21 (CD4 T cells)")

## w=clust==17
## m1="CD199 (CCR9)"
## m2="CD138"
## freqplot(breaks=50,preds_bgc[w,subset(a,target==m1)$file],preds_bgc[w,subset(a,target==m2)$file])
## title(xlab=m1,ylab=m2,main="Cluster 17 (Neutrophils)")
## dev.off()

## png(width=4000,height=2000,res=600,"./graphs/Illustration_bimodality_v2.png")
## par(mfrow=c(1,3))
## w=clust%in%c(4,17)
## m1="CD193 (CCR3)"
## m2="Ly-6G"
## plot(preds_bgc[w,subset(a,target==m1)$file],preds_bgc[w,subset(a,target==m2)$file],pch=16,cex=0.1,col=ifelse(clust[w]==4,"red","blue"),xlab=m1,ylab=m2,main="Eosinophils and Neutrophils")
## legend(x="topright",legend=c("Eosinophils","Neutrophils"),col=c("red","blue"),pch=16,cex=0.5)

## w=clust==17
## m1="MD-1"
## m2="CD49a"
## freqplot(breaks=50,preds_bgc[w,subset(a,target==m1)$file],preds_bgc[w,subset(a,target==m2)$file])
## title(xlab=m1,ylab=m2,main="Neutrophils\n(Cluster 17)")

## w=clust==19
## m1="CD200R (OX2R)"
## m2="IL-21R"
## freqplot(breaks=50,preds_bgc[w,subset(a,target==m1)$file],preds_bgc[w,subset(a,target==m2)$file])
## title(xlab=m1,ylab=m2,main="CD11b+ cDCs\n(Cluster 19)")
## dev.off()

## ## Confirmations
## png("./graphs/Illustration_bimodality_v2_confirmatory_data.png",width=2000,height=1000,res=300)
## par(mfrow=c(1,2))
## plot(dips_data[["d88"]],xlab="Measured",ylab="Predicted",main="MD-1 on Neutrophils",pch=16)
## plot(dips_data[["d106"]],xlab="Measured",ylab="Predicted",main="CD49a on Neutrophils",pch=16)
## dev.off()

## ## Checking measured expression of CCR7 in endothelial cells

## ## Phenograph -> AUC

## ## preds_sorted=apply(preds,2,order,decreasing=FALSE)
## ## aucs_uni=sapply(sort(unique(clust)),function(x){
## ##     apply(preds_sorted,2,function(y){
## ##         aucFromRules(clust[y],x)
## ##     })
## ## })
## ## rownames(aucs_uni)[rownames(aucs_uni)%in%rownames(a)]=paste0(a[rownames(aucs_uni)[rownames(aucs_uni)%in%rownames(a)],"target"],".predicted")
## ## pairs=combn(sort(unique(clust)),2,simplify=FALSE)
## ## aucs_pairs=sapply(
## ##     pairs,
## ##     function(x){
## ##         print(x)
## ##         preds_sorted=apply(preds[clust%in%x,],2,order,decreasing=FALSE)
## ##         apply(
## ##             preds_sorted,2,function(y){
## ##                 aucFromRules(
## ##                     y,
## ##                     x[1]
## ##                 )
## ##             }
## ##         )
## ##     }
## ## )
## ## colnames(aucs_pairs)=paste0(pairs,collapse="_vs_")
## ## rownames(aucs_uni_bgc)[rownames(aucs_uni_bgc)%in%rownames(a)]=paste0(a[rownames(aucs_uni_bgc)[rownames(aucs_uni_bgc)%in%rownames(a)],"target"],".predicted")

## preds_bgc_sorted=apply(preds_bgc,2,order,decreasing=FALSE)
## aucs_uni_bgc=sapply(sort(unique(clust)),function(x){
##     apply(preds_bgc_sorted,2,function(y){
##         aucFromRules(clust[y],x)
##     })
## })
## rownames(aucs_uni_bgc)[rownames(aucs_uni_bgc)%in%rownames(a)]=paste0(a[rownames(aucs_uni_bgc)[rownames(aucs_uni_bgc)%in%rownames(a)],"target"],".predicted")
## pairs=combn(sort(unique(clust)),2,simplify=FALSE)
## aucs_pairs_bgc=sapply(
##     pairs,
##     function(x){
##         print(x)
##         preds_bgc_sorted=apply(preds_bgc[clust%in%x,],2,order,decreasing=FALSE)
##         apply(
##             preds_bgc_sorted,2,function(y){
##                 aucFromRules(
##                     y,
##                     x[1]
##                 )
##             }
##         )
##     }
## )
## colnames(aucs_pairs_bgc)=sapply(pairs,paste0,collapse="_vs_)"
## rownames(aucs_pairs_bgc)[rownames(aucs_pairs_bgc)%in%rownames(a)]=paste0(a[rownames(aucs_pairs_bgc)[rownames(aucs_pairs_bgc)%in%rownames(a)],"target"],".predicted")

## library(data.table)
## fwrite(aucs_uni_bgc,file=paste0("/fh/fast/headley_m/InfinityFlowShare/NovemberLegendScreens/Steady-State/exports_background_corrected/AUC_BGC_phenograph.csv"))
## fwrite(aucs_pairs_bgc,file=paste0("/fh/fast/headley_m/InfinityFlowShare/NovemberLegendScreens/Steady-State/exports_background_corrected/AUC_BGC_pairs_phenograph.csv"))
## ####################
## ## Meaningplot export
## ####################

## preds_plot=cbind(preds_raw[,chans],preds_bgc)
## scrbl=sample(1:nrow(preds_plot))
## umap_plot=umap[scrbl,]
## w=colnames(preds_plot)%in%rownames(a)
## colnames(preds_plot)[w]=paste0(a[colnames(preds_plot)[w],"target"],".pred_bgc")
## colnames(preds_plot)=gsub("/","_",colnames(preds_plot))

## chop_quantiles=0.005
## for(col in colnames(preds_plot)){
##     q=quantile(preds_plot[,col],c(chop_quantiles,1-chop_quantiles))
##     preds_plot[,col][preds_plot[,col]<=q[1]]=q[1]
##     preds_plot[,col][preds_plot[,col]>=q[2]]=q[2]
##     preds_plot[,col]=preds_plot[,col][scrbl]
## }

## channels.code=setNames(colnames(preds_plot),colnames(preds_plot))
## preds_plot=cbind(umap_plot,preds_plot[,!w],preds_plot[,sort(colnames(preds_plot)[w])])

## library(matlab)
## color_biplot_by_channels(
##     preds_plot,
##     x_axis="UMAP1",
##     y_axis="UMAP2",
##     global_across_channels=FALSE,
##     file_name=paste0("/fh/fast/headley_m/InfinityFlowShare/NovemberLegendScreens/Steady-State/exports_background_corrected/","predicted_data_background_corrected.pdf"),
##     palette=jet.colors(100),
##     pch=16,
##     cex=0.2,
##     res=72,
##     raster.height=360*4,
##     raster.width=360*4
## )

## ####################
## ## More graph, not particularly useful thus far
## ####################

## cordata=cor(preds_bgc)
## pdf("./graphs/Correlation_matrix_heatmap.pdf",height=42,width=42)
## library(RColorBrewer)
## library(pheatmap)
## library(viridis)
## ## ar=model.matrix(~isotype,data=a[s,])
## ## ar=ar[,-match("(Intercept)",colnames(ar))]
## ## colnames(ar)=sub("isotype","",colnames(ar))
## ## ar=cbind(ar,1-rowSums(ar))
## ## colnames(ar)[ncol(ar)]=setdiff(a$isotype,colnames(ar))
## ## mode(ar)="character"
## ## isotypes=colnames(ar)
## ## ar=data.frame(ar,q95=apply(preds[,rownames(a)],2,quantile,0.95),stringsAsFactors=FALSE)
## ## ar[,"q95"]=asinh(ar[,"q95"])
## ## ar[ar[,"q95"]<0,"q95"]=min(ar[,"q95"][ar[,"q95"]>0])
## ## ar=cbind(ar,r2preds=r2preds[rownames(ar)],r2iso=r2iso[rownames(ar)],r2auto=r2auto[rownames(ar)])
## ## rownames(ar)=a[s,"target"]
## hm=pheatmap(
##     cordata,
##     cluster_rows=hclust(as.dist(1-cordata)/2),
##     cluster_cols=hclust(as.dist(1-cordata)/2),
##     breaks=seq(-1,1,by=0.01),
##     col=colorRampPalette(c("blue","white","red"))(200)
##     ## annotation_row=ar,
##     ## annotation_colors=c(list(r2preds=viridis(50),q95=viridis(50),r2iso=rev(viridis(50)),r2auto=rev(viridis(50))),s
##     ## apply(isotypes,function(x){c("0"="white","1"="black")},simplify=FALSE))
## )
## dev.off()

## mfis=t(aggregate_df(preds_bgc,clust,fun=median))    

## pheatmap(
##     mfis,
##     cluster_row=hm$tree_row
## )

## pca=prcomp(preds_bgc)
## png("./graphs/PCAs on UMAP_bgc.png",height=2000*4,width=2000*4,res=300)
## par(mfrow=c(4,4))
## for(i in 1:16){
##     plot(
##         test2,
##         pch=16,
##         cex=0.2,
##         col=toColors_continuous(pca$x[,i])
##     )
## }
## dev.off()

## mfis=t(aggregate_df(pca$x[,1:50],clust,fun=median))
## pheatmap(
##     mfis,
##     cluster_row=FALSE,
##     breaks=fastBreaks(pca$x[,1:50],201,sym=T),
##     col=colorRampPalette(c("blue","white","red"))(200)
## )

## x=subset(a,grepl("CD3",target))[1,1]
## y=subset(a,grepl("CD4",target))[1,1]
## z=subset(a,grepl("CD8a",target))[1,1]

## par(mfrow=c(1,2))
## freqplot(preds_raw[,y],preds_raw[,z])
## title(xlab="CD4",ylab="CD8a",main="raw")
## freqplot(preds_bgc[,y],preds_bgc[,z])
## title(xlab="CD4",ylab="CD8a",main="bg corrected")

## ## n=4
## ## i=0
## ## j=1
## ## for(file in rownames(a)){
## ##     if(i%%(n^2)==0){
## ##         png(paste0("./graphs/correlations_to_isotype_",sprintf("%02d",j),".png"),height=n*1000,width=n*1000,res=300)
## ##         par(mfrow=c(n,n))
## ##     }
    
## ##     print(match(file,rownames(a)))
## ##     isotype=a[file,"isotype"]
## ##     if(isotype!="Auto"){
## ##         iso=a[a$target==paste0("Isotype_",isotype),"file"]
## ##     } else {
## ##         iso=a[a$target=="Blank1","file"]
## ##     }
    
## ##     x=preds[,iso]
## ##     y=preds[,file]
    
## ##     lgcl=eb.autolgcl(c(x,y))
## ##     x=lgcl(x)
## ##     y=lgcl(y)
    
## ##     freqplot(x,y)
## ##     title(xlab="isotype",ylab="predicted",main=paste(a[file,"target"],a[file,"isotype"],signif(ar[file,"r2preds"],3),signif(cor(x,y),2)))
## ##     abline(lm(y~x),col="black",lty=2)

## ##     abline(a=0,b=1)

## ##     ## Total least squares
## ##     v=prcomp(cbind(x,y))$rotation
## ##     beta=v[2,1]/v[1,1]
## ##     intercept=mean(y)-beta*mean(x)
## ##     abline(a=intercept,b=beta,col="black",lty=3)
    
## ##     i=i+1
## ##     if(i%%(n^2)==0|match(file,rownames(a))==nrow(a)){
## ##         dev.off()
## ##         j=j+1
## ##     }
## ## }
