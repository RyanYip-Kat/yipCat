#' Volcano Plot of Difference Analysis Markers
#' @param deDF markers DataFrame
#' @param origin the markers from which platform,only support seurat and scanpy
#' @param showGenes the genew want to show in Volcano
#' @param cutoff the value to cutoff log2fc
#' @param x_limit the limits value for x_axis
#' @param y_limit the limits value for y_axis
#' @param down which value as down and show in left
#' @param up which value as up and show in right
#' @param colorManual color platte for scale_color_manual
#' @param outdir path to save plot
#' @export
DAPlot<-function(deDF=NULL,
                 origin="seurat",
                 showGenes=NULL,
                 cutoff=0.5,
                 x_limit=5,
                 y_limit=1e-500,
                 down="HC",
                 up="VKH",
                 colorManual=c('skyblue', 'gray', 'red'),
                 outdir="Plots"){
  colNames<-c("pvals","log2fc","pvals_adj","cluster","gene")
  if(tolower(origin)=="seurat"){
    selectCol<-c("p_val","avg_log2FC","p_val_adj","cluster","gene")
    markers<-deDF[,selectCol]
    colnames(markers)<-colNames
  }else if(tolower(origin)=="scanpy"){
    selectCol<-c("pvals","logfoldchanges","pvals_adj","cluster","names")
    markers<-deDF[,selectCol]
    colnames(markers)<-colNames
  }else{
    stop("Please check the markers's origin!!!")
  }

  if(!is.null(x_limit)){
    markers$log2fc[markers$log2fc > x_limit] = x_limit
    markers$log2fc[markers$log2fc < -x_limit] = -x_limit
    markers$pvals[markers$pvals < y_limit]= y_limit
  }
  if(!is.null(y_limit)){
    markers$pvals[markers$pvals < y_limit]= y_limit
  }

  uniqCluster<-unique(markers$cluster)
  stopifnot(length(uniqCluster)==2)# length must be 2
  if(is.null(down)|is.null(up)){
    down<-uniqCluster[1]
    up<-uniqCluster[2]
  }
  message("select down an up markers dataframe..")
  down_markers<-subset(markers,cluster == down & log2fc > 0)
  down_markers$log2fc <- -down_markers$log2fc
  up_markers<-subset(markers,cluster == up & log2fc > 0)
  markers<-rbind(up_markers,down_markers)

  message("Plot..")
  if(!is.null(showGenes)){
    message("only show text specify markers..")
    #direction<-with(markers,ifelse(gene%in%showGenes,"Down","Up"))
    direction<-with(markers,ifelse(log2fc > cutoff,"Up",
                                   ifelse(log2fc < -cutoff,"Down","NS")))
    markers$direction<-direction
    p<-ggplot(data = markers,
           aes(x = log2fc,
               y = -log10(pvals))) +
      geom_point(size = 3,aes(color = direction),show.legend = T) +
      scale_color_manual(values = colorManual) +
      labs(x = 'Log2(fold change)', y = '-log10(p-value)')+
      ggstatsplot::theme_ggstatsplot() +
      cowplot::theme_half_open() +
      labs(x = 'Log2(fold change)',
           y = '-log10(p-value)')
    if(!is.null(x_limit)){
      p<-p+xlim(-x_limit,x_limit)
    }
    p<-p+ggrepel::geom_text_repel(data = markers%>%
                        filter(gene%in%showGenes),
                      aes(x = log2fc,y = -log10(pvals),
                          label = gene),
                      size = 5,box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.5, "lines"),
                      segment.color = "grey50",
                      show.legend = FALSE,
                      colour = "black")
  }else{
    direction<-with(markers,ifelse(log2fc > cutoff,"Up",
                                    ifelse(log2fc < -cutoff,"Down","NS")))
    markers$direction<-direction
    label<-with(markers,ifelse(abs(log2fc) > cutoff,gene,""))
    markers$label<-label

    p <- ggplot(data = markers,
                 aes(x = log2fc,
                     y = -log10(pvals))) +
      geom_point(size =3.5,aes(color = direction),show.legend = T) +
      scale_color_manual(values=colorManual)+
      labs(x = 'Log2(fold change)', y = '-log10(p-value)') +
      ggstatsplot::theme_ggstatsplot() +
      cowplot::theme_half_open() +
      labs(x = 'Log2(fold change)',
           y = '-log10(p-value)')+
      geom_hline(yintercept = -log10(0.05),
                   linetype = 'dotdash',
                   color = 'grey30') +
      geom_vline(xintercept = c(-cutoff,cutoff),
                 linetype = 'dotdash',
                 color = 'grey30')
    if(!is.null(x_limit)){
      p<-p+xlim(-x_limit,x_limit)
    }

    p<-p+ggrepel::geom_text_repel(data = markers, aes(x = log2fc,
                                               y = -log10(pvals),
                                               label = label),
                           size =2.5,box.padding = unit(0.35, "lines"),
                           point.padding = unit(0.5, "lines"),
                           segment.color = "grey50",
                           show.legend = TRUE,
                           colour = "black")+theme_bw()
  }
  makedir(outdir)
  ggsave(file.path(outdir,"DAplot.pdf"),plot=p,width = 16,height = 12)
  return(p)
}
