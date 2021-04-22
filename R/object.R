#' create seurat object from count matrix
#' @param file 10X count or aggr filtered_feature_bc_matrix or filtered_feature_bc_matrix.h5
#' @param species human or others,if human,genes names upper,otherwise,title
#' @param addPercent whether add percent metrics
#' @param addID whther add ID name in barcode
#' @param preclust whether pre load analysis cluster
#' @export
create10xObj<-function(file=NULL,
                       species="human",
                       addPercent=TRUE,
                       addID=NULL,
                       preclust=TRUE){
  count<-tryCatch({
    message("loading matrix..")
    v<-str_split(basename(file),"\\.")[[1]]
    if(v[length(v)]=="h5"){
      count<-Read10X_h5(file)
    }else{
      count<-Read10X(file)
    }
    if(!is.null(addID)){
      cells<-colnames(count)
      cells<-unlist(lapply(cells,function(cell){
        v<-str_split(cell,"-")[[1]]
        paste0(v[1],"-",addID)
      }))
      colnames(count)<-cells
    }
    genes<-rownames(count)
    if(str_to_lower(species)=="human"){
      rownames(count)<-str_to_upper(genes)
    }else{
      rownames(count)<-str_to_title(genes)
    }
    count
  },error=function(e){
    message("loading count matrix error,please check file!")
    NULL
  })

  obj<-NULL
  if(!is.null(count)){
    message("create object...")
    obj <- CreateSeuratObject(counts =count,
                              min.cells = 0,
                              min.features =0)
  }

  if(preclust){
    metaList<-list()
    if(!"analysis"%in%list.files(dirname(file))){
      stop("Please check analysis path whether in cellranger aggr or count output if you want to preload!!!")
    }
    clusterDir<-file.path(dirname(file),"analysis/clustering")
    clusterNames<-list.files(clusterDir)
    if(!is.null(clusterNames)){
      for(Name in clusterNames){
        f<-file.path(clusterDir,Name,"clusters.csv")
        m<-tryCatch({
          d<-read.csv(f,stringsAsFactors=FALSE)
          rownames(d)<-as.character(d$Barcode)
          d<-subset(d,select=Cluster)
          colnames(d)<-Name
          d
        },error=function(e){
          cat(sprintf("Check file -- %s whether is invalid!!!\n",f))
          NULL
        })
        if(!is.null(m)){
          metaList[[Name]]<-m
        }
      }
    }

    if(length(metaList)>0){
      if(!is.null(obj)){
        message("add metadata ..")
        for(i in 1:length(metaList)){
          meta<-metaList[[i]]
          Name<-names(metaList)[i]
          cat(sprintf("add %s metadata\n",Name))
          obj<-AddMetaData(obj,meta,col.name =Name )
        }
      }
    }
  }
  if(addPercent){
    if(!is.null(obj)){
      message("add metrics..")
      if(str_to_lower(species)=="human"){
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
        obj[["percent.rpl"]] <- PercentageFeatureSet(obj, pattern = "^RPL")
        obj[["percent.rps"]] <- PercentageFeatureSet(obj, pattern = "^RPS")

        genes<-rownames(obj)
        keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]
        keep_genes<-keep_genes[!str_detect(keep_genes,"\\.")]
        obj@misc[["pureGene"]]<-keep_genes

      }else{
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^Mt-")
        obj[["percent.rpl"]] <- PercentageFeatureSet(obj, pattern = "^Rpl")
        obj[["percent.rps"]] <- PercentageFeatureSet(obj, pattern = "^Rps")

        genes<-rownames(obj)
        keep_genes<-genes[!str_detect(genes,"^Mt-|^Rpl|^Rps")]
        keep_genes<-keep_genes[!str_detect(keep_genes,"\\.")]
        obj@misc[["pureGene"]]<-keep_genes
      }
    }
  }
  obj
}

#' merge single seurat objects into multiple
#' @param objs seurat objects or seurat rds files
#' @param samples sample names relative objs
#' @param status statu names relative obj
#' @param format seurat or file
#' @export
makeMultObject<-function(objs=NULL,
                         samples=NULL,
                         status=NULL,
                         format="seurat"){
  stopifnot(length(objs)==length(samples))
  stopifnot(length(objs)==length(status))
  stopifnot(tolower(format)%in%c("seurat","file"))
  countList<-list()
  metaList<-list()
  ##############
  for(i in seq_along(objs)){
    ob<-objs[[i]]
    sam<-samples[i]
    st<-status[i]
    cat(sprintf("loading object -- [ sample : %s ] -- [ statu : %s ] \n",sam,st))
    if(tolower(format)=="seurat"){
      obj<-ob
    }else{
      obj<-readRDS(ob)
    }
    obj@meta.data$Sample<-sam
    obj@meta.data$Status<-st
    obj@meta.data$Barcode<-Cells(obj)
    countList[[sam]]<-GetAssayData(obj,"counts")
    metaList[[sam]]<-obj@meta.data[,c("Barcode","Sample","Status")]
  }

  message("merge object lists..")
  count<-do.call(cbind,countList)
  metaData<-do.call(rbind,metaList)
  rownames(metaData)<-metaData$Barcode
  obj<-CreateSeuratObject(counts = count,min.cells=1,min.features=1)
  obj<-AddMetaData(obj,metaData)

}

#' main pipeline for seurat object
#' you need to install SeuratWrappers,hamony first
#' @param obj seurat object
#' @param nHVG deafult 2000,number of HVG
#' @param batch_correct whether batch correct
#' @param batch_key if batch correct,which key as batch
#' @param batch_method if batch correct,which method use
#' @param CellCycle whether do CellCycle to remove CellCycle effect
#' @import SeuratWrappers
#' @export
objectPipeline<-function(obj=NULL,
                         nHVG=2000,
                         batch_correct=FALSE,
                         batch_key=NULL,
                         batch_method="harmony",
                         CellCycle=TRUE){
  #CellCycleScoring
  message("NormalizeData..")
  obj<-NormalizeData(obj)
  message("FindVariableFeatures..")
  obj<- FindVariableFeatures(obj, selection.method = "vst",
                                 nfeatures =nHVG,verbose = FALSE)

  if(CellCycle){
    message("CellCycle Scoring..")
    s.genes <- cc.genes[["s.genes"]]
    g2m.genes <- cc.genes[["g2m.genes"]]
    obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    var.to.regress<-c("S.Score", "G2M.Score")
  }else{
    var.to.regress<-c("nFeature_RNA")
  }
  message("Run ScaleData")
  obj<-ScaleData(obj,features=VariableFeatures(obj),model.use = "linear",
                    vars.to.regress = var.to.regress,verbose =FALSE)

  message("Run PCA")
  obj<-RunPCA(obj,npcs = 50,verbose = FALSE)
  use_rep<-"pca"

  if(batch_correct){
    require("SeuratWrappers")
    require("harmony")
    vars=ifelse(!is.null(batch_key),batch_key,"Sample")
    stopifnot(vars%in%colnames(obj@meta.data))
    if(batch_method=="harmony"){
      use_rep<-"harmony"
      cat(sprintf("Run Harmony with %s\n",vars))
      obj<-RunHarmony(obj, group.by.vars =vars,
                        reduction = "pca",dims.use=1:30)
    }else if(batch_method=="fastmnn"){
      cat(sprintf("Run fastmnn with %s\n",vars))
      obj.list <-SplitObject(obj, split.by = vars)
      obj<-RunFastMNN(obj.list,features=3000)
      use_rep<-"mnn"
    }else{
      stop("Invalid batch correct method !")
    }
  }


  message("Run UMAP..")
  obj <- RunUMAP(obj, reduction =use_rep, dims = 1:30)
  message("Run TSNE..")
  obj<-RunTSNE(obj,reduction=use_rep)

  message("Find Clusters")
  obj <- FindNeighbors(obj, reduction =use_rep, dims = 1:30)
  obj <- FindClusters(obj,resolution=0.8)
  obj
}

