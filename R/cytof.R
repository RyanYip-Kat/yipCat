##############################
###### function to transform cytof data into Seurat object
##############################

#' function to convert cytof data into Seurat
#'
#' @param sample_csv at least three columns(Path,Sample,Condition)
#' @param config_csv at least two columns(channel(Cd111Di,Nd143Di...),marker(CD8A,CD11C..))
#' @param N sample size,default sample 20000 each fcs file
#' @param useFeatures feature names provide to run scaleData and runpca.
#' @param batch_correct whether run batch correct to remove batch effect
#' @param batch_key when remove batch correct,which key use  to run.
#' @param path2barcode10X 10X barcode file,can be found in cellranger path
#' @param transformation whether transfrom data
#' @export
#'
#'

Cytof2Seurat<-function(sample_csv=NULL,
                       config_csv=NULL,
                       N=20000,
                       cofactor=5,
                       useFeatures=NULL,
                       batch_correct=FALSE,
                       batch_key="sample",
                       transformation=FALSE,
                       path2barcode10X="3M-february-2018.txt",
                       outdir="cytof"
){

  #require(stringr)
  #require(Seurat)
  #require(flowCore)
  #library(DropletUtils)

  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = TRUE,showWarnings = FALSE)
  }
  message("read 10x barcode..")
  barcode_df=data.table::fread(path2barcode10X,stringsAsFactors=FALSE,header=FALSE)
  barcodes=as.character(barcode_df$V1)

  message("read sample file..")
  samples_DF=read.csv(sample_csv,stringsAsFactors=FALSE,sep=",")
  stopifnot("Condition"%in%colnames(samples_DF))
  stopifnot("Path"%in%colnames(samples_DF))  # fcs file
  stopifnot("Sample"%in%colnames(samples_DF)) # fcs sample name
  cat(sprintf("There are : %d sample fcs file\n",nrow(samples_DF)))

  message("loading config file..")
  config=read.csv(config_csv,stringsAsFactors=FALSE,sep=",")
  stopifnot("marker"%in%colnames(config))
  stopifnot("channel"%in%colnames(config))

  channel=config$channel
  markers=config$marker
  pattern=stringr::str_extract(channel,"\\d+")  # extract channel number(unique)

  markers_df=data.frame("markers"=markers,"pattern"=pattern)
  rownames(markers_df)=markers_df$pattern

  pattern=paste(pattern,collapse="|")
  names(markers)=channel

  message("read flowSet")
  fset=list()
  for(i in 1:length(samples_DF$Path)){
    cat(sprintf("INFO : read [ %d ] fcs file\n",i))
    fr=read.FCS(samples_DF$Path[i],which.lines=N,column.pattern=pattern,transformation =transformation)
    column=flowCore::colnames(fr)

    beat=stringr::str_extract(column,"\\d+")
    marker=markers_df[beat,]$markers
    flowCore::colnames(fr)=marker

    fset[[samples_DF$Sample[i]]]=fr
  }

  message("rename fset exprs")
  exprs_list=list()
  for(i in 1:nrow(samples_DF)){
    cat(sprintf("INFO : rename [ %d ] fset object\n",i))
    sample=samples_DF$Sample[i]
    condition=samples_DF$Condition[i]
    rowName=paste0(condition,"#",paste(sample,1:nrow(exprs(fset[[sample]])),sep="#"))
    barcode=sample(barcodes,size=nrow(exprs(fset[[sample]])),replace=FALSE)
    cells=paste(barcode,1,sep="-")
    rowName=paste0(rowName,"#",cells)

    rownames(exprs(fset[[sample]]))=rowName
    expr=exprs(fset[[sample]])
    if(!transformation){
      expr <- asinh(expr / cofactor) # transform data
    }
    exprs(fset[[sample]])=expr
    exprs_list[[i]]=expr
    barcodes=setdiff(barcodes,barcode)
  }
  message("Save flowSet")
  #saveRDS(fset,file.path(outdir,"flowSet.rds"))

  message("Concat dataset..")
  exprs_matrix=do.call(rbind,exprs_list)
  cat(sprintf("Theraa are : %s cells\n",nrow(exprs_matrix)))

  Names=rownames(exprs_matrix)
  status=unlist(lapply(Names,function(x){return(str_split(x,"#")[[1]][1])}))
  samples=unlist(lapply(Names,function(name){return(str_split(name,"#")[[1]][2])}))
  cells=unlist(lapply(Names,function(name){return(str_split(name,"#")[[1]][4])}))
  rownames(exprs_matrix)=cells
  metadata=data.frame("cell_id"=cells,"sample"=samples,"status"=status,"condition"=status)
  Names=rownames(exprs_matrix)


  message("Create Seurat object..")
  seurat<-CreateSeuratObject(counts=t(exprs_matrix),
                             project ="cytof",
                             min.cells=0,
                             min.features=1)
  rownames(metadata)=Names
  seurat<-AddMetaData(seurat,metadata)

  if(is.null(useFeatures)){
    useFeatures<-rownames(seurat)
  }
  if(transformation){
    seurat <- NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)
  }

  use_rep <-"pca"
  ndims <- 20

  message("Run PCA")
  seurat<-ScaleData(seurat,features=useFeatures)
  seurat<-RunPCA(seurat,npcs=30,features=useFeatures,verbose=FALSE)

  if(batch_correct){
     require(harmony)
     message("Run Harmony")
     if(is.null(batch_key)){
       batch_key<-"sample"
     }
     harmony_embeddings=HarmonyMatrix(exprs_matrix, metadata, batch_key) #v3
     rownames(harmony_embeddings)=Names
     colnames(harmony_embeddings)=paste("Harmony",1:ncol(harmony_embeddings),sep="_")

     use_rep<-"harmony"
     ndims<-20

     message("add harmony")
     mat=as.matrix(harmony_embeddings)
     colnames(mat)<-paste("Harmony_",1:ncol(mat),sep = "")
     seurat[["harmony"]]<-CreateDimReducObject(embeddings =mat,
                                              key = "Harmony_",
                                              assay = DefaultAssay(seurat))
  }


  message("Run UMAP")
  seurat <- RunUMAP(object = seurat, reduction =use_rep, dims = 1:ndims)
  message("Run TSNE")
  seurat <- RunTSNE(object = seurat, reduction = use_rep, dims = 1:ndims)
  message("Run FindNeighbors")
  seurat <- FindNeighbors(object = seurat,reduction =use_rep, dims = 1:ndims)

  message("Run Find Clusters")
  seurat <- FindClusters(object = seurat,resolution=0.8, verbose =TRUE)

  message("write counts matrix into 10x")
  counts<-GetAssayData(seurat,"counts")
  rownames(counts)=str_to_upper(rownames(counts))

  filename=file.path(outdir,"filtered_feature_cytof_matrix.h5")
  cat(sprintf("Write counts matrix into : %s\n",filename))
  write10xCounts(x =counts, path=filename,type="HDF5")

  filename=file.path(outdir,"filtered_feature_cytof_matrix")
  cat(sprintf("Write counts matrix into : %s\n",filename))
  write10xCounts(x =counts, path=filename)

  td <- system.file()
  #cmd<-sprintf("cd %s  && awk '{print $1"\t"$2"\tGene Expression"}' genes.tsv > features.tsv && rm genes.tsv && sed -i 's/\./\-/' barcodes.tsv && gzip *\"",filename)

  return(seurat)
}

