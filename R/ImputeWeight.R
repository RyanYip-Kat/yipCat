######################
### caculate impute weight in scRNA
######################


#' function to caculate impute weights
#' @param object an seurat object.
#' @param reducedDims Embedding name in seurat object.
#' @param dimsToUse,integer number of dims to use.
#' @param sampleCells number of sample cells used
#' @export
ImputeWeights <- function(
  object = NULL,
  reducedDims = "IterativeLSI",
  dimsToUse = NULL,
  corCutOff = 0.75,
  td = 3,
  ka = 4,
  sampleCells = 5000,
  nRep = 2,
  k = 15,
  epsilon = 1,
  randomSuffix = FALSE,
  threads = 8,
  seed = 1,
  verbose = TRUE,
  ...
){

  require(Seurat)
  require(SummarizedExperiment)
  require(dplyr)

  message("INFO : addImputeWeights Input-Parameters")

  #Adapted From
  #https://github.com/dpeerlab/magic/blob/master/R/R/run_magic.R

  set.seed(seed)

  tstart <- Sys.time()

  #Get Reduced Dims
  stopifnot(reducedDims%in%Reductions(object))
  matDR <- Embeddings(object,reduction=reducedDims)
  if(!is.null(dimsToUse)){
    matDR <- matDR[,1:dimsToUse]
  }
  N <- nrow(matDR)
  rn <- rownames(matDR)

  if(!is.null(sampleCells)){
    if(sampleCells > nrow(matDR)){
      sampleCells <- NULL
    }
  }
  if(is.null(sampleCells)){
    binSize <- N
    nRep <- 1
  }else{
    cutoffs <- lapply(seq_len(1000), function(x){
      N / x
    }) %>% unlist
    binSize <- min(cutoffs[order(abs(cutoffs - sampleCells))[1]] + 1, N)
  }
  if(!dir.exists("ImputeWeights")){
    dir.create("ImputeWeights", showWarnings = FALSE)
  }
  if(randomSuffix){
    weightFiles <- .tempfile("Impute-Weights", tmpdir ="ImputeWeights")
    weightFiles <- paste0(weightFiles, "-Rep-", seq_len(nRep))
  }else{
    weightFiles <- file.path("ImputeWeights", paste0("Impute-Weights-Rep-", seq_len(nRep)))
  }

  o <- suppressWarnings(file.remove(weightFiles))

  weightList <- .safelapply(seq_len(nRep), function(y){

    cat(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s)\n", y, nRep))

    if(!is.null(sampleCells)){
      idx <- sample(seq_len(nrow(matDR)), nrow(matDR))
      blocks <- split(rownames(matDR)[idx], ceiling(seq_along(idx)/binSize))
    }else{
      blocks <- list(rownames(matDR))
    }

    weightFile <- weightFiles[y]


    blockList <- lapply(seq_along(blocks), function(x){

      if(x %% 10 == 0){
        cat(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s, Iteration %s of %s)", y, nRep, x, length(blocks)))
      }

      ix <- blocks[[x]]
      Nx <- length(ix)

      #Compute KNN
      knnObj <- nabor::knn(data = matDR[ix,], query = matDR[ix, ], k = k)
      knnIdx <- knnObj$nn.idx
      knnDist <- knnObj$nn.dists
      rm(knnObj)

      if(ka > 0){
        knnDist <- knnDist / knnDist[,ka]
      }

      if(epsilon > 0){
        W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x=c(knnDist), dims = c(Nx, Nx))
      } else {
        W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x=1, dims = c(Nx, Nx)) # unweighted kNN graph
      }
      W <- W + Matrix::t(W)

      #Compute Kernel
      if(epsilon > 0){
        W@x <- exp(-(W@x / epsilon^2))
      }

      #Markov normalization
      W <- W / Matrix::rowSums(W)

      #Initialize Matrix
      Wt <- W

      #Computing Diffusion Matrix
      for(i in seq_len(td)){
        Wt <- Wt %*% W
      }
      #rownames(Wt) <- rownames(matDR)[ix]
      #colnames(Wt) <- rownames(matDR)[ix]
      rownames(Wt) <- ix
      colnames(Wt) <- ix
      #rm(knnIdx)
      #rm(knnDist)
      #rm(W)
      #gc()
      return(Wt)
    }) %>% SimpleList

    names(blockList) <- paste0("b",seq_along(blockList))
    blockList
  }, threads = threads) %>% SimpleList
  names(weightList) <- paste0("w",seq_along(weightList))

  cat(sprintf("Completed Getting Magic Weights!", round(object.size(weightList) / 10^9, 3)))

  imputeWeights <- SimpleList(
    Weights = weightList,
    Params =
      list(
        reducedDims = reducedDims,
        td = td,
        k = k,
        ka = ka,
        epsilon = epsilon
      )
  )
  object@misc[["imputeWeights"]]<-imputeWeights  # store in misc
  return(object)
  #imputeWeights

}



#' impute matrix use impute weights
#' @param mat an matrix or sparse matrix.
#' @param imputeWeights impute weight object from function `ImputeWeights`.
#' @export
MyimputeMatrix <- function(
  mat = NULL,
  imputeWeights = NULL,
  threads = 16,
  verbose = FALSE
){


  if(!inherits(imputeWeights$Weights, "SimpleList") & !inherits(imputeWeights$Weights, "list")){
    message("Weights are not a list, Please re-run ImputeWeights (update)!")
    stop("Weights are not a list, Please re-run addImputeWeights (update)!")
  }

  message("INFO : imputeMatrix Input-Parameters")

  weightList <- imputeWeights$Weights

  tstart <- Sys.time()

  imputeMat <- lapply(seq_along(weightList), function(x){

    cat(sprintf("Imputing Matrix (%s of %s)\n", x, length(weightList)))

    message("Using weights in memory")

    matx <- .safelapply(seq_along(weightList[[x]]), function(y){

      if(verbose) message(y, " ", appendLF = FALSE)
      message(paste0(y, " of ", length(weightList[[x]])))
      Matrix::t(as.matrix(weightList[[x]][[y]]) %*% Matrix::t(mat[, paste0(colnames(weightList[[x]][[y]])), drop = FALSE]))

    }, threads = threads) %>% Reduce("cbind", .)

    if(verbose) message("")

    matx[, colnames(mat)] #Return Ordered

  }) %>% Reduce("+", .)

  #Compute Average
  imputeMat <- imputeMat / length(weightList)

  message("Finished Imputing Matrix")

  imputeMat

}

