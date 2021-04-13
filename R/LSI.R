#' compute LSI(pass)
#' @export
.computeLSI <- function(
  mat = NULL,
  LSIMethod = 1,
  scaleTo = 10^4,
  nDimensions = 50,
  binarize = TRUE,
  outlierQuantiles = c(0.02, 0.98),
  seed = 1,
  verbose = FALSE
){
  #require(S4Vectors)
  #require(Matrix)
  #require(dplyr)
  out2 <- tryCatch({

    set.seed(seed)

    #TF IDF LSI adapted from flyATAC
    if(binarize){
      mat@x[mat@x > 0] <- 1
    }

    #Compute Col Sums
    message("Computr col Sums")
    colSm <- Matrix::colSums(mat)
    if(any(colSm == 0)){
      exclude <- which(colSm==0)
      mat <- mat[,-exclude, drop = FALSE]
      colSm <- colSm[-exclude]
    }else{
      exclude <- c()
    }

    cn <- colnames(mat)
    filterOutliers <- 0
    if(!is.null(outlierQuantiles)){
      qCS <- quantile(colSm, probs = c(min(outlierQuantiles), max(outlierQuantiles)))
      idxOutlier <- which(colSm <= qCS[1] | colSm >= qCS[2])
      if(length(idxOutlier) > 0){
        #.safeSaveRDS(mat, "temp.rds", compress = FALSE)
        matO <- mat[, idxOutlier, drop = FALSE]
        mat <- mat[, -idxOutlier, drop = FALSE]
        mat2 <- mat[, head(seq_len(ncol(mat)), 10), drop = FALSE] # A 2nd Matrix to Check Projection is Working
        colSm <- colSm[-idxOutlier]
        filterOutliers <- 1
      }
    }

    #Clean up zero rows
    message("Clean up zero rows")
    rowSm <- Matrix::rowSums(mat)
    idx <- which(rowSm > 0)
    mat <- mat[idx, ]
    rowSm <- rowSm[idx]

    #TF - Normalize
    message("TF - Normalize")
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))

    if(LSIMethod == 1 | tolower(LSIMethod) == "tf-logidf"){

      #Adapted from Casanovich et al.

      #LogIDF
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")

      #TF-LogIDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

    }else if(LSIMethod == 2 | tolower(LSIMethod) == "log(tf-idf)"){

      #Adapted from Stuart et al.

      #IDF
      idf   <- as(ncol(mat) / rowSm, "sparseVector")

      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

      #Log transform TF-IDF
      mat@x <- log(mat@x * scaleTo + 1)

    }else if(LSIMethod == 3 | tolower(LSIMethod) == "logtf-logidf"){

      #LogTF
      mat@x <- log(mat@x + 1)

      #LogIDF
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")

      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    }else{

      stop("LSIMethod unrecognized please select valid method!")

    }

    gc()

    #Calc SVD then LSI
    message("Calc SVD then LSI")
    svd <- irlba::irlba(mat, nDimensions, nDimensions)
    svdDiag <- matrix(0, nrow=nDimensions, ncol=nDimensions)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))

    #Return Object
    out <- SimpleList(
      matSVD = matSVD,
      rowSm = rowSm,
      nCol = length(colSm),
      exclude = exclude,
      idx = idx,
      svd = svd,
      binarize = binarize,
      scaleTo = scaleTo,
      nDimensions = nDimensions,
      LSIMethod = LSIMethod,
      outliers = NA,
      date = Sys.Date(),
      seed = seed
    )

    if(filterOutliers == 1){
      #Quick Check LSI-Projection Works
      pCheck <- .projectLSI(mat = mat2, LSI = out, verbose = verbose)
      pCheck2 <- out[[1]][rownames(pCheck), ]
      pCheck3 <- lapply(seq_len(ncol(pCheck)), function(x){
        cor(pCheck[,x], pCheck2[,x])
      }) %>% unlist
      if(min(pCheck3) < 0.95){
        stop("Error with LSI-projection! Cor less than 0.95 of re-projection. Please report bug to github!")
      }
      #Project LSI Outliers
      out$outliers <- colnames(matO)
      outlierLSI <- .projectLSI(mat = matO, LSI = out, verbose = verbose)
      allLSI <- rbind(out[[1]], outlierLSI)
      allLSI <- allLSI[cn, , drop = FALSE] #Re-Order Correctly to original
      out[[1]] <- allLSI
    }
    rm(mat)
    gc()

    out

  }, error = function(e){
    message("ProjectLSI Error!")
  })

  out2

}

#' project LSI(pass)
#' @export
.projectLSI <- function(
  mat = NULL,
  LSI = NULL,
  returnModel = FALSE,
  verbose = FALSE
){

  #require(S4Vectors)
  out2 <- tryCatch({

    require(Matrix)
    set.seed(LSI$seed)

    #Get Same Features
    mat <- mat[LSI$idx,]

    #Binarize Matrix
    if(LSI$binarize){
      mat@x[mat@x > 0] <- 1
    }

    #TF
    colSm <- Matrix::colSums(mat)
    if(any(colSm == 0)){
      exclude <- which(colSm==0)
      mat <- mat[,-exclude]
      colSm <- colSm[-exclude]
    }
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))

    if(LSI$LSIMethod == 1 | tolower(LSI$LSIMethod) == "tf-logidf"){

      #Adapted from Casanovich et al.

      #LogIDF
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")

      #TF-LogIDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

    }else if(LSI$LSIMethod == 2 | tolower(LSI$LSIMethod) == "log(tf-idf)"){

      #Adapted from Stuart et al.

      #IDF
      idf   <- as(LSI$nCol / LSI$rowSm, "sparseVector")

      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

      #Log transform TF-IDF
      mat@x <- log(mat@x * LSI$scaleTo + 1)

    }else if(LSI$LSIMethod == 3 | tolower(LSI$LSIMethod) == "logtf-logidf"){

      #LogTF
      mat@x <- log(mat@x + 1)

      #LogIDF
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")

      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    }else{

      stop("LSIMethod unrecognized please select valid method!")

    }

    gc()

    #Clean Up Matrix
    idxNA <- Matrix::which(is.na(mat),arr.ind=TRUE)
    if(length(idxNA) > 0){
      mat[idxNA] <- 0
    }

    #Calc V
    V <- Matrix::t(mat) %*% LSI$svd$u %*% Matrix::diag(1/LSI$svd$d)

    #LSI Diagonal
    svdDiag <- matrix(0, nrow=LSI$nDimensions, ncol=LSI$nDimensions)
    diag(svdDiag) <- LSI$svd$d
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
    matSVD <- as.matrix(matSVD)
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))

    if(returnModel){
      X <- LSI$svd$u %*% diag(LSI$svd$d) %*% t(V)
      out <- list(matSVD = matSVD, V = V, X = X)
    }else{
      out <- matSVD
    }

    out


  }, error = function(e){
    message("ProjectLSI Error!")

  })

  out2

}
