####################################################################
# scRNA Trajectory Analysis Methods(Like ArchR)
####################################################################


.getQuantiles <- function (v = NULL, len = length(v))
{
  if (length(v) < len) {
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }
  else {
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if (length(v) < len) {
    p <- p[seq_along(v)]
  }
  return(p)
}


#' Add a Supervised Trajectory to an Seurat Object
#'
#' This function will fit a supervised trajectory in a lower dimensional space that
#' can then be used for downstream analyses.
#'
#' @param seurat An `Seurat` object.
#' @param name A string indicating the name of the fitted trajectory to be added in `meta.data`.
#' @param trajectory The order of cell groups to be used for constraining the initial supervised fitting procedure.
#' For example, to get a trajectory from Cluster1 to Cluster2 to Cluster3, input should be c("Cluster1", "Cluster2", "Cluster3").
#' Cells will then be used from these 3 groups to constrain an initial fit in the group order.
#' @param groupBy A string indicating the column name from `meta.data` that contains the cell group definitions used in
#' `trajectory` to constrain the initial supervised fitting procedure.
#' @param embedding A string indicating the name of the `embedding` object from the `seurat` that should be used for distance computation.
#' @param preFilterQuantile Prior to the initial supervised trajectory fitting, cells whose euclidean distance from the cell-grouping
#' center is above the provided quantile will be excluded.
#' @param postFilterQuantile After initial supervised trajectory fitting, cells whose euclidean distance from the cell-grouping center
#' is above the provided quantile will be excluded.
#' @param useAll A boolean describing whether to use cells outside of trajectory groups for post-fitting procedure.
#' @param dof The number of degrees of freedom to be used in the spline fit. See `stats::smooth.spline()` for more information.
#' @param spar The sparsity to be used in the spline fit. See `stats::smooth.spline()` for more information.
#' @param force A boolean value indicating whether to force the trajactory indicated by `name` to be overwritten if it already exists in the given `ArchRProject`.
#' @param seed A number to be used as the seed for random number generation for trajectory creation.
#' @export

addSeuratTrajectory <- function(
  object = NULL,
  name = "Trajectory",
  trajectory = NULL,
  groupBy = "Clusters",
  #reducedDims = "IterativeLSI",
  embedding = NULL,
  preFilterQuantile = 0.9,
  postFilterQuantile = 0.9,
  useAll = FALSE,
  dof = 250,
  spar = 1,
  force = FALSE,
  seed = 1
){
  #require(dplyr)
  #require(S4Vectors)

  stopifnot(class(object)=="Seurat")
  if(!is.null(seed)) set.seed(seed)
  meta.data <- object@meta.data
  stopifnot(groupBy%in%colnames(meta.data))
  groupDF <- meta.data[,groupBy,drop=FALSE]
  groupDF <- groupDF[groupDF[,1] %in% trajectory,,drop=FALSE]

  if(sum(unique(groupDF[,1]) %in% trajectory)==0){
    stop("trajectory does not span any groups in groupBy! Are you sure your input is correct?")
  }

  if(sum(unique(groupDF[,1]) %in% trajectory) < 3){
    if(!force){
      stop("trajectory must span at least 3 groups in groupBy!")
    }
  }

  if(is.null(embedding)){
    mat <- Embeddings(object,reduction ="pca")
  }else{
    mat <- Embeddings(object, reduction = embedding)
  }
  mat <- mat[rownames(groupDF),,drop = FALSE]

  ######################################################
  #Filter Outliers
  ######################################################
  filterObj <- lapply(seq_along(trajectory), function(x){

    #Subset
    groupsx <- rownames(groupDF)[groupDF[,1]==trajectory[x]]
    matx <- mat[groupsx,,drop = FALSE]

    #Filter Distance
    matMeanx <- colMeans(matx)
    diffx <- sqrt(colSums((t(matx) - matMeanx)^2))
    idxKeep <- which(diffx <= quantile(diffx, preFilterQuantile))

    #Filter
    list(mat = matx[idxKeep,,drop=FALSE], groups = groupsx[idxKeep])

  })

  matFilter <- lapply(seq_along(filterObj), function(x) filterObj[[x]]$mat) %>% Reduce("rbind", .)
  groupsFilter <- groupDF[lapply(seq_along(filterObj), function(x) filterObj[[x]]$groups) %>% Reduce("c", .),,drop=FALSE]

  ######################################################
  #Now Initial Alignment
  ######################################################
  message("INFO : Initial Alignment Before Spline Fit")
  initialTime <- lapply(seq_along(trajectory), function(x){

    groupsx <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x]]
    matx <- matFilter[groupsx,,drop = FALSE]

    #Get Differences
    if(x != length(trajectory)){
      groupsxp1 <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x + 1]]
      meanx <- colMeans(matFilter[groupsxp1,,drop = FALSE])
      diffx <- sqrt(colSums((t(matx) - meanx)^2))
      timex <- (1 - .getQuantiles(diffx)) + x
    }else{
      groupsxm1 <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x - 1]]
      meanx <- colMeans(matFilter[groupsxm1,,drop = FALSE])
      diffx <- sqrt(colSums((t(matx) - meanx)^2))
      timex <- .getQuantiles(diffx) + x
    }

    timex

  }) %>% unlist

  ######################################################
  #Fit Cubic Splines
  ######################################################
  message("INFO :Spline Fit")
  matSpline <- lapply(seq_len(ncol(matFilter)), function(x){
    tryCatch({
      stats::smooth.spline(
        x = initialTime,
        y = matFilter[names(initialTime), x],
        df = dof,
        spar = spar
      )[[2]]
    }, error = function(e){
      errorList <- list(
        it = x,
        x = initialTime,
        y = matFilter[names(initialTime), x],
        df = dof,
        spar = spar
      )
      #.logError(e, fn = "smooth.spline", info = "", errorList = errorList, logFile = logFile)
    })
  }) %>% Reduce("cbind",.) %>% data.frame()

  ######################################################
  # 1. KNN Fit vs Actual
  ######################################################
  message("INFO : KNN to Spline")
  knnObj <- nabor::knn(
    data =  matSpline,
    query = mat,
    k = 3
  )

  #Estimate place along trajectory
  knnIdx <- knnObj[[1]]
  knnDist <- knnObj[[2]]
  knnDiff <- ifelse(knnIdx[,2] > knnIdx[,3], 1, -1)
  knnDistQ <- .getQuantiles(knnDist[,1])

  #Filter Outlier Cells to Trajectory for High Resolution
  idxKeep <- which(knnDist[,1] <= quantile(knnDist[,1], postFilterQuantile))
  dfTrajectory <- DataFrame(
    row.names = rownames(mat),
    Distance = knnDist[, 1],
    DistanceIdx = knnIdx[, 1] + knnDiff * knnDistQ
  )[idxKeep, , drop = FALSE]

  ######################################################
  # 2. Fit cells not in trajectory clusters
  ######################################################
  if(useAll){
    message("INFO : Aligning cells not in trajectory")

    if(is.null(embedding)){
      mat2 <- Embeddings(object, reduction = "pca")
    }else{
      mat2 <- Embeddings(object, reduction = embedding)
    }
    meta.data <- object@meta.data
    groupDF <- meta.data[,groupBy,]
    groupDF <- groupDF[groupDF[,1] %ni% trajectory,,drop=FALSE]
    mat2 <- mat2[rownames(groupDF),,drop = FALSE]

    #Nearest Neighbors
    knnObj2 <- nabor::knn(
      data =  matSpline,
      query = mat2,
      k = 3
    )

    #Estimate place along trajectory
    knnIdx2 <- knnObj2[[1]]
    knnDist2 <- knnObj2[[2]]
    knnDiff2 <- ifelse(knnIdx2[,2] > knnIdx2[,3], 1, -1)
    knnDistQ2 <- .getQuantiles(knnDist2[,1])

    #Keep Cells that are within the maximum distance of a cluster
    idxKeep <- which(knnDist2[,1] < max(dfTrajectory[,1]))
    dfTrajectory2 <- DataFrame(
      row.names = rownames(mat2),
      Distance = knnDist2[, 1],
      DistanceIdx = knnIdx2[, 1] + knnDiff2 * knnDistQ2
    )[idxKeep, , drop = FALSE]

    #Final Output
    dfTrajectory3 <- rbind(dfTrajectory, dfTrajectory2)
  }else{
    dfTrajectory3 <- dfTrajectory
  }

  dfTrajectory3$Trajectory <- 100 * .getQuantiles(dfTrajectory3[,2])
  tr <- as.data.frame(dfTrajectory3)
  tr <- tr[,"Trajectory",drop=FALSE]
  object <- AddMetaData(object,tr,col.name=name)
  object
}


#' Plot Trajectory after function `getSeuratTrajectory`
#'
#' This function will plot Trajectory in seurat object
#'
#' @param seurat An `Seurat` object.
#' @param name A string indicating the name  in `meta.data`.
#' @param trajectory The order of cell groups to be used for constraining the initial supervised fitting procedure.
#' For example, to get a trajectory from Cluster1 to Cluster2 to Cluster3, input should be c("Cluster1", "Cluster2", "Cluster3").
#' Cells will then be used from these 3 groups to constrain an initial fit in the group order.
#' @param colorBy A string indicating the column name from `meta.data` or feature in slot data,eg(data,counts,scale.data)
#' @param embedding A string indicating the name of the `embedding` object from the `seurat` that should be used for embedding plot.
#' @param log2Norm A boolean value indicating whether a log2 transformation should be performed on the values from `colorBy`.
#' @param imputeWeights The weights to be used for imputing numerical values for each cell as a linear combination of other cells'
#' values. See `addImputationWeights()` and `getImutationWeights()` for more information.
#' @param pal The name of a custom palette from `ArchRPalettes` to use for coloring cells.
#' @param size A number indicating the size of the points to plot if `plotAs` is set to "points".
#' @param rastr A boolean value that indicates whether the plot should be rasterized. This does not rasterize lines and labels,
#' just the internal portions of the plot.
#' @param quantCut If this is not `NULL`, a quantile cut is performed to threshold the top and bottom of the distribution of numerical values.
#' This prevents skewed color scales caused by strong outliers. The format of this should be c(x,y) where x is the lower threshold and y is
#' the upper threshold. For example, quantileCut = c(0.025,0.975) will take the 2.5th percentile and 97.5 percentile of values and
#' set values below/above to the value of the 2.5th and 97.5th percentile values respectively.
#' @param quantHex The numeric xth quantile of all dots within each individual hexagon will determine the numerical value for
#' coloring to be displayed. This occurs when (i) `plotAs` is set to "hex" or (ii) `plotAs` is set to `NULL` and the values of `colorBy` are numeric.
#' @param discreteSet The name of a discrete palette from `ArchRPalettes` for visualizing `colorBy` in the embedding if a discrete color set is desired.
#' @param continuousSet The name of a continuous palette from `ArchRPalettes` for visualizing `colorBy` in the embedding if a continuous color set is desired.
#' @param randomize A boolean value that indicates whether to randomize points prior to plotting to prevent cells from one cluster
#' being present at the front of the plot.
#' @param keepAxis A boolean value that indicates whether the x and y axis ticks and labels should be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param addArrow A boolean value that indicates whether to add a smoothed arrow in the embedding based on the aligned trajectory.
#' @param plotAs A string that indicates whether points ("points") should be plotted or a hexplot ("hex") should be plotted. By default
#' if `colorBy` is numeric, then `plotAs` is set to "hex".
#' @param smoothWindow An integer value indicating the smoothing window for creating inferred Arrow overlay on to embedding.
#' @param ... Additional parameters to pass to `ggPoint()` or `ggHex()`.
#' @export

plotSeuratTrajectory <- function(
  object = NULL,
  embedding = "umap",
  trajectory = "Trajectory",
  name = "Trajectory",
  colorBy="cellcoldata",
  log2Norm = TRUE,
  pal = NULL,
  size = 0.2,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  quantHex = 0.5,
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  plotAs = NULL,
  smoothWindow = 5,
  addFit=NULL,
  addArrow=TRUE,
  ...
){

  #require("ggplot2")
  #require("ArchR")
  #require("Seurat")
  if(is.null(quantCut)){
    quantCut <- c(0, 1)
  }


  ##############################
  # Plot Helpers
  ##############################
  .summarizeHex <- function(x = NULL){
    quantile(x, quantHex, na.rm = TRUE)
  }

  ##############################
  # Get Trajectory
  ##############################
  meta.data <- object@meta.data
  dfT <- meta.data[,trajectory,drop=FALSE]
  idxRemove <- which(is.na(dfT[,1]))

  ##############################
  # Get Embedding
  ##############################
  df <- as.data.frame(Embeddings(object, reduction = embedding))
  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("x", "y", "PseudoTime")

  #Parameters
  plotParams <- list(...)
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of  Trajectory")
  plotParams$baseSize <- baseSize
  if(is.null(addFit)){
    plotParams$addFit="loess"
  }else{
    plotParams$addFit=addFit
  }

  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
    plotParams$color <- as.vector(meta.data[,name,drop=FALSE][rownames(df), 1])
    plotParams$discrete <- isDiscrete(plotParams$color)
    plotParams$continuousSet <- "horizonExtra"
    plotParams$discreteSet <- "stallion"
    plotParams$title <- paste(plotParams$title, " colored by\ncolData : ", name)
    if(is.null(plotAs)){
      plotAs <- "hexplot"
    }
  }else{
    plotParams$continuousSet <- "solarExtra"
    if(log2Norm){
      slot <- "data"
    }else{
      slot <- "counts"
    }
    plotParams$continuousSet <- "horizonExtra"
    MatSc <- GetAssayData(object,slot)
    color <- as.matrix(MatSc[name,,drop=FALSE])
    plotParams$color <- color[, rownames(df), drop = FALSE]
    plotParams$discrete <- FALSE
    plotParams$title <- sprintf("%s colored by\n%s : %s", plotParams$title, colorBy, name)
    if(is.null(plotAs)){
      plotAs <- "hexplot"
    }
    if(plotAs=="hexplot"){
      plotParams$fun <- .summarizeHex
    }

  }


  #Additional Params!
  plotParams$xlabel <- colnames(df)[1]
  plotParams$ylabel <- colnames(df)[2]

  if(!is.null(continuousSet)){
    plotParams$continuousSet <- continuousSet
  }
  if(!is.null(continuousSet)){
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize

  if(plotParams$discrete){
    plotParams$color <- paste0(plotParams$color)
  }

  if(!plotParams$discrete){

    plotParams$color <- as.vector(plotParams$color)

    if(name != trajectory){
      plotParams$color <- quantileCut(plotParams$color, min(quantCut), max(quantCut))
    }

    if(!is.null(log2Norm)){
      if(log2Norm){
        plotParams$color <- log2(plotParams$color + 1)
        plotParams$colorTitle <- paste0("Log2Norm")
      }else{
        plotParams$colorTitle <- "Counts"
      }
    }

    plotParams$color[idxRemove] <- NA
    plotParams$pal <- paletteContinuous(set = plotParams$continuousSet)

    if(tolower(plotAs) == "hex" | tolower(plotAs) == "hexplot"){
      plotParams$addPoints <- TRUE
      if(is.null(plotParams$bins)){
        plotParams$bins <- 150
      }
      message("Plotting")
      out <- do.call(ggHex, plotParams)
    }else{
      message("Plotting")
      out <- do.call(ggPoint, plotParams)
    }

  }else{
    message("Plotting")
    out <- do.call(ggPoint, plotParams)
  }

  if(!keepAxis){
    out <- out + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }
  if(!plotParams$discrete){
    message("Plotting Trajectory")
    #Prep Trajectory Vector
    dfT$value <- plotParams$color
    dfT <- dfT[order(dfT$PseudoTime), ]
    dfT <- dfT[!is.na(dfT$PseudoTime), ]

    #Plot Pseudo-Time
    out2 <- ggPoint(
      x = dfT$PseudoTime,
      y = dfT$value,
      color = dfT$PseudoTime,
      discrete = FALSE,
      xlabel = "PseudoTime",
      ylabel = name,
      pal = plotParams$pal,
      ratioYX = 0.5,
      rastr = TRUE
    ) + geom_smooth(color = "black")
    attr(out2, "ratioYX") <- 0.5
  }else{
    out2<-NULL
  }

  if(addArrow){

    dfArrow <-  split(dfT, floor(dfT$PseudoTime / 1.01)) %>%
      lapply(colMeans) %>% Reduce("rbind",.) %>% data.frame
    dfArrow$x <- centerRollMean(dfArrow$x, smoothWindow)
    dfArrow$y <- centerRollMean(dfArrow$y, smoothWindow)
    dfArrow <- rbind(dfArrow, dfT[nrow(dfT), ,drop = FALSE])

    out <- out + geom_path(
      data = dfArrow, aes(x, y, color=NULL), size= 1,
      arrow = arrow(type = "open", length = unit(0.1, "inches"))
    )
  }
  list(out, out2)

}


#seurat=readRDS("/Path/to/seurat")
#trajectory=c("Memory BC","Naive BC","Plasma BC")
#groupBy="label_fine"
#seurat=getSeuratTrajectory(seurat,groupBy=groupBy,trajectory=trajectory,embedding="umap")
#p1=plotSeuratTrajectory(seurat,embedding="umap",name="PAX6",colorBy="matrix")
#ggsave("plot2.pdf",plot=p1[[2]])
#ggsave("plot1.pdf",plot=p1[[1]])


getGroupMatrix<-function (seurat=NULL, features = NULL, groupList = NULL,
                          threads = 4,  verbose = TRUE,slot="data",splitBy=NULL,
                          asSparse = FALSE,useFM=FALSE)
{
  #require(ArchR)
  #require(stringr)

  Donors <- unlist(lapply(Cells(seurat),function(x)return(str_split(x,"-")[[1]][2])))
  seurat$Donors <- Donors
  meta <- seurat@meta.data
  Donors <- unique(Donors)
  #seqnames <- unique(featureDF$seqnames)
  if(!useFM){
    if(is.null(features) & !is.null(splitBy)){
      features <- rownames(seurat)
      Groups <- table(meta[[splitBy]])
      features <- lapply(names(Groups),function(i)return(sample(features,size=Groups[[i]],replace=TRUE)))
    }
  }else{
    if(is.null(splitBy)){stop("groupBy must not be NULL!")}
    require(future)
    message("Find Markers ...")
    plan("multiprocess", workers = threads)
    Idents(seurat)=meta[[splitBy]]
    markers <- FindAllMarkers(object=seurat,
                              logfc.threshold = 0.25,
                              test.use = "wilcox",
                              only.pos=FALSE,
                              slot=slot)
    Groups <- unique(markers$cluster)
    features <- lapply(Groups,function(i){
      f <- markers$gene[markers$cluster%in%i]
      return(f)})
  }
  #rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  allMat <- GetAssayData(seurat,slot=slot)
  cellNames <- unlist(groupList, use.names = FALSE)
  allCellsList <- lapply(seq_along(Donors), function(x) {
    allCells <- rownames(meta)[meta$Donors%in%x]
    allCells <- allCells[allCells %in% cellNames]
    if (length(allCells) != 0) {
      allCells
    }
    else {
      NULL
    }
  })
  mat <- .safelapply(seq_along(features), function(x) {
    featureDFx <- features[[x]]
    matChr <- matrix(0, nrow =length(featureDFx), ncol = length(groupList))
    colnames(matChr) <- names(groupList)
    rownames(matChr) <- featureDFx
    for (y in seq_along(Donors)) {
      allCells <- allCellsList[[y]]
      if (!is.null(allCells)) {
        maty <- allMat[featureDFx,allCells,drop=FALSE]
        for (z in seq_along(groupList)) {
          cellsGroupz <- groupList[[z]]
          idx <- BiocGenerics::which(colnames(maty) %in%
                                       cellsGroupz)
          if (length(idx) > 0) {
            matChr[, z] <- matChr[, z] + Matrix::rowSums(maty[,
                                                              idx, drop = FALSE])
          }
        }
        rm(maty)
      }
      if (y%%20 == 0 | y%%length(Donors) == 0) {
        gc()
      }
    }
    if (asSparse) {
      matChr <- as(matChr, "dgCMatrix")
    }
    cat(sprintf("Finished Group Matrix %s of %s\n",x, length(features)))
    matChr
  }, threads = threads) %>% Reduce("rbind", .)
  mat <- mat[unique(unlist(features)), , drop = FALSE]
  message("Successfully Created Group Matrix")
  gc()
  return(mat)
}

getGroupList <- function(
  seurat = NULL,
  name = "Trajectory",
  groupEvery = 1
){

  meta <-  seurat@meta.data
  trajectory <- meta[,name,drop=FALSE]
  trajectory <- trajectory[!is.na(trajectory[,1]),,drop=FALSE]
  breaks <- seq(0, 100, groupEvery)
  if(!all(is.numeric(trajectory[,1]))){
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if(!all(trajectory[,1] >= 0 & trajectory[,1] <= 100)){
    stop("Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?")
  }

  groupList <- lapply(seq_along(breaks), function(x){
    if(x == 1){
      NULL
    }else{
      rownames(trajectory)[which(trajectory[,1] > breaks[x - 1] & trajectory[,1] <= breaks[x])]
    }
  })[-1]
  names(groupList) <- paste0("T.", breaks[-length(breaks)], "_", breaks[-1])
  return(groupList)
}

#' function to get Trajectory Object.
#' @param seurat an `seurat` object.
#' @param splitBy column in `meta.data` use as ident for find marker or spliting for features.
#' @param useFM bool value,whether to use `FindMarkers`.
#' @param name the trajectory'name in `meta.data` after caculate specify trajectory.
#' @param slot slot name in seurat.
#' @param features list of features,each list include many genes.
#' @export
getSeuratTrajectory <- function(
  seurat = NULL,
  splitBy=NULL,
  useFM=FALSE, # use FindAllMarkers
  name = "Trajectory",
  groupEvery = 1,
  log2Norm = TRUE,
  slot="data",
  features=NULL,
  scaleTo = 10000,
  smoothWindow = 11,
  threads = 16
){
  message("Creating Trajectory Group List ..")
  groupList <- getGroupList(seurat,name=name,groupEvery=groupEvery)
  message("Creating Trajectory Group Matrix..")
  groupMat <- getGroupMatrix(seurat,
                             features=features,
                             useFM=useFM,
                             slot=slot,
                             groupList=groupList,
                             splitBy=splitBy,
                             threads=threads)
  if(!is.null(scaleTo)){
    if(any(groupMat < 0)){
      message("Some values are below 0, this could be a DeviationsMatrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
    }else{
      groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
    }
  }

  if(log2Norm){
    if(any(groupMat < 0)){
      message("Some values are below 0, this could be a DeviationsMatrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
    }else{
      groupMat <- log2(groupMat + 1)
    }
  }

  if(!is.null(smoothWindow)){
    message("Smoothing...")
    smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) .centerRollMean(x, k = smoothWindow))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <- SummarizedExperiment(
      assays = SimpleList(
        smoothMat = as.matrix(smoothGroupMat),
        mat = as.matrix(groupMat)))
  }else{
    colnames(groupMat) <- paste0(colnames(groupMat))
    #Create SE
    seTrajectory <- SummarizedExperiment(assays = SimpleList(mat = as.matrix(groupMat)))
  }

  metadata(seTrajectory)$Params <- list(matrixClass="matrix",
                                        scaleTo = scaleTo,
                                        log2Norm = log2Norm,
                                        smoothWindow = smoothWindow,
                                        date = Sys.Date())
  seTrajectory

}


.centerRollMean <- function (v = NULL, k = NULL)
{
  o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
  if (k%%2 == 0) {
    o2 <- c(rep(o1[k], floor(k/2) - 1), o1[-seq_len(k - 1)],
            rep(o1[length(o1)], floor(k/2)))
  }
  else if (k%%2 == 1) {
    o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k - 1)],
            rep(o1[length(o1)], floor(k/2)))
  }
  else {
    stop("Error!")
  }
  o2
}



