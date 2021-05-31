#' using ggstatsplot to visual gene expression or continous name in metadata
#' @param obj  Seurat Object
#' @param colorBy matrix or metadata
#' @param groupBy which group in metadata as groupBy factor
#' @param name which variable will be plotting
#' @param plot.type which type to plot ,boxviolin,box or violin
#' @param type A character specifying the type of statistical approach. Four possible options:parametric,nonparametric,robust,bayes
#' @param slot which slot in Seurat Object to use,data ot count
#' @param p.adjust.method fdr(default),Adjustment method for p-values for multiple comparisons. Possible methods are: holm,hochberg, hommel, bonferroni, BH, BY,fdr, none.
#' @param ...  detail see ggstatsplot docs
#' @return a ggplot object
#' @export
statsBoxViolin<-function(obj,
                         colorBy="matrix",
                         groupBy="seurat_clusters",
                         name="CD8A",
                         plot.type="boxviolin",
                         type = "parametric",
                         slot="data",
                         p.adjust.method = "fdr",
                         ...){
  stopifnot(plot.type%in%c("box","violin","boxviolin"))
  stopifnot(colorBy%in%c("matrix","metadata"))

  metaData<-obj@meta.data
  if(colorBy=="matrix"){
    M<-GetAssayData(obj,slot=slot)
    DATA<-as.data.frame(t(as.matrix(M[name,,drop=FALSE])))
    metaDF<-metaData[,groupBy,drop=FALSE]
    DF<-cbind(metaDF,DATA)
    colnames(DF)<-c("Group","Expr")
  }else{
    if(isDiscrete(metaData[[name]])){
      stop(sprintf("%s can not be Discrete!!!",name))
    }
    M<-as.data.frame(M[,name,drop=FALSE])
    metaDF<-metaData[,groupBy,drop=FALSE]
    DF<-cbind(metaDF,DATA)
    colnames(DF)<-c("Group","Expr")
  }

  set.seed(123)
  p<-ggbetweenstats(data=DF,
                 x=Group,
                 y=Expr,
                 plot.type = plot.type,
                 type=type,
                 xlab = groupBy,
                 ylab = name,
                 p.adjust.method = p.adjust.method,
                 ...)
  return(p)
}
