tolower<-stringr::str_to_lower
#' function to quantile cut
#' @export
quantileCut<-function (x = NULL, lo = 0.025, hi = 0.975, maxIf0 = TRUE)
{
  q <- quantile(x, probs = c(lo, hi))
  if (q[2] == 0) {
    if (maxIf0) {
      q[2] <- max(x)
    }
  }
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}


#' function to plot Groups
#' @param seurat an seurat object.
#' @param assay which assay to use.
#' @param groupBy which column  in `meta.data` use as group.
#' @param name which feature in matrix
#' @param plotAS ridges or viloin plot.
#' @export
#'
MyplotGroups <- function(
  seurat = NULL,
  groupBy = "Sample",
  assay = "RNA",
  name = "Krt14",
  maxCells = 1000,
  quantCut = c(0.002, 0.998),
  pal = NULL,
  discreteSet = "stallion",
  ylim = NULL,
  size = 0.5,
  baseSize = 6,
  ratioYX = NULL,
  ridgeScale = 2,
  plotAs = "violin",
  ...
){
  require(Seurat)
  require(ggplot2)
  require(stringr)
  #Make Sure ColorBy is valid!
  if(length(assay) > 1){
    stop("assay must be of length 1!")
  }
  assayNames=names(seurat@assays)
  stopifnot(assay%in%assayNames)
  stopifnot(plotAs%in%c("ridges","violin"))

  metadata=seurat@meta.data
  groups <- metadata[,groupBy,drop=FALSE]
  groups[[groupBy]]=as.character(groups[[groupBy]])
  groupNames <- groups[,1]
  names(groupNames) <- rownames(groups)
  groupNames2 <- gtools::mixedsort(unique(groupNames))


  plotParams <- list(...)
  log2Norm <- TRUE

  Mat=GetAssayData(seurat,slot="data",assay=assay)
  colorMat=Mat[name,,drop=FALSE]

  colorList <- lapply(seq_len(nrow(colorMat)), function(x){
    colorParams <- list()
    colorParams$color <- colorMat[x, ]
    if(!is.null(discreteSet)){
      colorParams$pal <- suppressMessages(paletteDiscrete(values = groupNames2, set = discreteSet))
    }
    if(!is.null(pal)){
      colorParams$pal <- pal
    }
    colorParams
  })


  if(!is.null(maxCells)){
    splitGroup <- split(names(groupNames), groupNames)
    useCells <- lapply(splitGroup, function(x){
      if(length(x) > maxCells){
        sample(x, maxCells)
      }else{
        x
      }
    }) %>% unlist %>% as.vector
    idx <- match(useCells, names(groupNames))
  }else{
    idx <- seq_along(groupNames)
  }

  pl <- lapply(seq_along(colorList), function(x){

    message(paste0(x, " "), appendLF = FALSE)

    if(is.null(ylim)){
      ylim <- range(colorList[[x]]$color,na.rm=TRUE) %>% extendrange(f = 0.05)
    }

    plotParamsx <- plotParams
    plotParamsx$x <- groupNames[idx]
    if(!is.null(quantCut)){
      plotParamsx$y <- quantileCut(colorList[[x]]$color[idx], min(quantCut), max(quantCut))
    }else{
      plotParamsx$y <- colorList[[x]]$color[idx]
    }
    plotParamsx$xlabel <- groupBy
    plotParamsx$ylabel <- name[x]
    plotParamsx$baseSize <- baseSize
    plotParamsx$ridgeScale <- ridgeScale
    plotParamsx$ratioYX <- ratioYX
    plotParamsx$size <- size
    plotParamsx$plotAs <- plotAs
    plotParamsx$pal <- colorList[[x]]$pal

    p <- do.call(ggGroup, plotParamsx)

    p

  })

  names(pl) <- name
  message("")

  if(length(name)==1){
    pl[[1]]
  }else{
    pl
  }

}

#' function to plot feature Groups
#' @param seurat an seurat object.
#' @param assay which assay to use.
#' @param groupBy which column  in `meta.data` use as group.
#' @param features  features in matrix
#' @param combine whether combine plots
#' @export
#'
GgPlot<-function(seurat=NULL,
                 assay="RNA",
                 features=NULL,
                 outDir=NULL,
                 splitBy=NULL,
                 groupBy=NULL,
                 combine=FALSE,
                 width=12,
                 height=10){
  require(Seurat)
  require(ggplot2)
  if(is.null(outDir)){
    outDir <- "Plots"
  }
  if(!dir.exists(outDir)){
    dir.create(outDir,recursive = TRUE,showWarnings = TRUE)
  }

  if(!is.null(splitBy)){
    seurat_list=SplitObject(seurat,split.by=splitBy)
  }else{
    seurat_list=list("SeuratObject"=seurat)  # only one object list
  }

  for(i in seq_along(seurat_list)){
    sr=seurat_list[[i]]
    name=names(seurat_list)[i]
    cat(sprintf("INFO : [ %d of %d ] --- [ %s ]\n",i,length(seurat_list),name))
    plotList=MyplotGroups(seurat=sr,
                          assay=assay,
                          groupBy=groupBy,
                          name=features,
                          plotAs="violin")

    if(combine){
      MyplotPDF(plotList,name=name,outpath=outDir,width=width,height=height)
    }else{
      for(f in names(plotList)){
        cat(sprintf("INFO : Save --- [ %s ] \n",name))
        p=plotList[[f]] #+ ggpubr::stat_compare_means()
        #ggsave(file.path(outDir,paste0(name,"-",f,".pdf")),plot=p,width=12,height=10)
        MyplotPDF(p,name=paste0(name,"-",f),outpath=outDir,width=width,height=height)
      }
    }
  }
}
