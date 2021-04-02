#' a function to caculate aging score
#' @param seurat a seurat object.
#' @param geneSet a vector of gene names
#' @param identColumn which column in metadata as idents
#' @param plot whether show aging score plot
#' @export
ageScore<-function(seurat=NULL,
                   geneSet=NULL,
                   identColumn=NULL,
                   plot=FALSE){
   require(Seurat)
  .gsScore<-function(seurat=NULL,
                    geneSet=NULL,
                    identColumn="celltype",
                    CellID="pDC"){
    stopifnot(class(seurat)=="Seurat")
    metadata <- seurat@meta.data
    Idents(seurat)<-metadata[[identColumn]]
    cells <- WhichCells(seurat,idents=CellID)
    gs_seurat <- subset(seurat,cells = cells,features = geneSet)
    gsij_umi <- gs_seurat@meta.data$nCount_RNA

    cj_seurat <- subset(seurat,idents = CellID)
    cij_umi <- cj_seurat@meta.data$nCount_RNA

    n_cells <- length(cells)
    Cj <- gsij_umi/cij_umi*100 #score of Cj cell for GSx gene set
    D <- data.frame("GS"=Cj,"Cluster"=CellID)
    return(D)
  }
  stopifnot(class(seurat)=="Seurat")
  metadata <- seurat@meta.data
  DATA <- lapply(unique(metadata[[identColumn]]),.gsScore,seurat=seurat,identColumn=identColumn,geneSet=geneSet)
  score <- do.call(rbind,DATA)
  if(plot){
    require(ggplot2)
    p <- ggplot(score,aes(x=Cluster,y=GS,fill=Cluster))+
      #geom_violin()+geom_jitter()+theme_bw()+
      geom_violin()+scale_y_sqrt()+
      theme_bw()+
      theme(legend.position = "none",
            axis.text.x = element_text(size = 15,angle=90))
    print(p)

  }
  return(score)
}
