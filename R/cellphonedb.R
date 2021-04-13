.suppressAll<-function (expr = NULL)
{
  suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
}


.checkPath<-function (u = NULL, path = NULL, throwError = TRUE)
{
  require(dplyr)
  if (is.null(u)) {
    out <- TRUE
  }
  out <- lapply(u, function(x, error = TRUE) {
    if (Sys.which(x) == "") {
      if (!is.null(path) && file.exists(file.path(path,
                                                  x))) {
        o <- TRUE
      }
      else {
        if (throwError) {
          stop(x, " not found in path, please add ",
               x, " to path!")
        }
        else {
          o <- FALSE
        }
      }
    }
    else {
      o <- TRUE
    }
    return(o)
  }) %>% unlist %>% all
  return(out)
}

#' function to find cellphonedb
#'
#' @export
findCellPhoneDB<-function ()
{
  message("Searching For cellphonedb..")
  if (.suppressAll(.checkPath("cellphonedb", throwError = FALSE))) {
    message("Found with $path!")
    return("cellphonedb")
  }
  message("Not Found in $PATH")
  search2 <- suppressWarnings(tryCatch({
    system2("pip", "show cellphonedb", stdout = TRUE, stderr = NULL)
  }, error = function(x) {
    "ERROR"
  }))
  search3 <- suppressWarnings(tryCatch({
    system2("pip3", "show cellphonedb", stdout = TRUE, stderr = NULL)
  }, error = function(x) {
    "ERROR"
  }))
  if (length(search2) > 0) {
    if (search2[1] != "ERROR") {
      path2Install <- gsub("Location: ", "", search2[grep("Location",
                                                          search2, ignore.case = TRUE)])
      path2Bin <- gsub("lib/python/site-packages", "bin/cellphonedb",
                       path2Install)
      if (.suppressAll(.checkPath(path2Bin, throwError = error))) {
        message("Found with pip!")
        return(path2Bin)
      }
    }
  }
  message("Not Found with pip")
  if (length(search3) > 0) {
    if (search3[1] != "ERROR") {
      path2Install <- gsub("Location: ", "", search3[grep("Location",
                                                          search3, ignore.case = TRUE)])
      path2Bin <- gsub("lib/python/site-packages", "bin/cellphonedb",
                       path2Install)
      if (.suppressAll(.checkPath(path2Bin, throwError = error))) {
        message("Found with pip3!")
        return(path2Bin)
      }
    }
  }
  message("Not Found with pip3")
  stop("Could Not Find Macs2! Please install w/ pip, add to your $PATH, or just supply the macs2 path directly and avoid this function!")
}

#' function to install cellphonedb
#'
#' @export
installCellphonedb<-function(){
  message("Install cellphonedb")
  search2 <- suppressWarnings(tryCatch({
    system2("pip", "install cellphonedb", stdout = TRUE, stderr = NULL)
  }, error = function(x) {
    "Install ERROR"
  }))
}


#' a function to export count and metadata for cellphonedb
#'
#' @param seurat a seurat object
#' @param cells cell barcode vector.
#' @param cells feature names vector.
#' @param slot which slot to use.
#' @param selectCol which column to export as cell_type.
#' @param outdir path to save tables
#' @param cellphonedbPath  cellphonedb path
#' @param runCPDB whether to run cellphonedb
#' @export
exportCellPhoneDB<-function(seurat=NULL,
                            cells=NULL,
                            features=NULL,
                            slot="data",
                            selectCol="seurat_clusters",
                            species="human",
                            outdir="cellphoneDB",
                            path2CPDB=findCellPhoneDB(),
                            runCPDB=FALSE){
  require(Seurat)
  stopifnot(class(seurat)=="Seurat")
  stopifnot(slot%in%c("counts","data","scale.data"))
  if(species=="human"){
    require(org.Hs.eg.db)
    orgDB <- org.Hs.eg.db
  }else{
    require(org.Mm.eg.db)
    orgDB <- org.Mm.eg.db
  }
  if(!is.null(cells)){
    seurat <- subset(seurat,cells=cells)
  }
  if(!is.null(features)){
    seurat <- subset(seurat,features=features)
  }
  message("get counts..")
  counts <- GetAssayData(seurat,slot = slot)
  counts <- as.data.frame(as.matrix(counts))

  symbols <- mapIds(x=orgDB,keys=rownames(counts),keytype="SYMBOL",column="ENSEMBL")
  Genes <- data.frame(Gene=as.character(symbols))

  counts <- cbind(Genes,counts)
  counts <- na.omit(counts)

  message("get metadata..")
  metaData <- seurat@meta.data
  stopifnot(selectCol%in%colnames(metaData))
  celltypes <- metaData[[selectCol]]
  cells <- colnames(seurat)
  meta.data<-data.frame(Cell=cells,cell_type=celltypes,stringsAsFactors=FALSE)

  message("export..")

  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = TRUE,showWarnings = TRUE)
  }
  write.table(counts,file.path(outdir,"cell_counts.txt"),sep="\t",row.names = FALSE,quote=F)
  write.table(meta.data,file.path(outdir,"cell_meta.txt"),sep="\t",row.names = FALSE,quote=F)

  if(runCPDB){
    if(is.null(path2CPDB)){
      installCellphonedb()
    }
    message("Run cellphonedb..")
    countFile <- file.path(outdir,"cell_counts.txt")
    metaFile <- file.path(outdir,"cell_meta.txt")
    cmd <- sprintf("method statistical_analysis %s %s  --project-name out --output-path %s --threads 8 --pvalue 0.05",
                   metaFile,countFile,outdir)
    system2(path2CPDB,cmd, stdout = TRUE, stderr = NULL)
  }
}



##########################
#' function to prepare data from cellphone output for dot_plot
#'
#' @param selected_rows  vector of row names in pvalue.txt' interacting_pair
#' @param selected_columns vector of col names in pvalue.txt
#' @param means_path the cellphone method analysis out's means.txt file
#' @param pvalues_path the cellphone method analysis out's pvalues.txt file
#' @export
CPDBDotData<-function(selected_rows = NULL,
           selected_columns = NULL,
           means_path = './means.txt',
           pvalues_path = './pvalues.txt',
           means_separator = '\t',
           pvalues_separator = '\t'
  ){

    all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
    all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)

    intr_pairs = all_pval$interacting_pair
    all_pval = all_pval[,-c(1:11)]
    all_means = all_means[,-c(1:11)]

    if(is.null(selected_rows)){
      selected_rows = intr_pairs
    }

    if(is.null(selected_columns)){
      selected_columns = colnames(all_pval)
    }

    sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
    sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

    df_names = expand.grid(selected_rows, selected_columns)
    pval = unlist(sel_pval)
    pval[pval==0] = 0.0009
    plot.data = cbind(df_names,pval)
    pr = unlist(as.data.frame(sel_means))
    pr[pr==0] = 1
    plot.data = cbind(plot.data,log2(pr))
    colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
    return(plot.data)
  }


#' function to prepare data from cellphone output for heatmap_plot
#'
#' @param meta_file meta data file,cell_meta.txt
#' @param pvalues_file the cellphone method analysis out's pvalues.txt file
#' @param pvaue pvalue Threshold.
#' @export
CPDBHeatmapData<-function(meta_file=NULL,
                          pvalues_file=NULL,
                          meta_sep='\t',
                          pvalues_sep='\t',
                          pvalue=0.05){
  #######   Network
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)

  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]


  split_sep = '\\|'
  join_sep = '|'

  pairs1_all = unique(meta[,2])

  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))

  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]

    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))

    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }

  all_count = all_count[-1,]
  #######   count interactions

  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]

    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))

  }

  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])
    return(count_matrix)
  }
  else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}


#' dot plot for cellphone output
#'
#' @param selected_rows  vector of row names in pvalue.txt' interacting_pair
#' @param selected_columns vector of col names in pvalue.txt
#' @param means_path the cellphone method analysis out's means.txt file
#' @param pvalues_path the cellphone method analysis out's pvalues.txt file
#' @param filename   plot filename
#' @export
#'
CPDBDotplot <- function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.pdf',
                    width = 8,
                    height = 10,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){
  require(ggplot2)
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)

  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }

  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }

  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

  p<-ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

  if (output_extension == '.pdf') {
    ggsave(filename, plot=p,width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
    ggsave(filename, plot=p,width = width, height = height, limitsize=F)
  }
}




#' heatmap plot for cellphone output
#'
#' @param meta_file meta data file,cell_meta.txt
#' @param pvalues_file the cellphone method analysis out's pvalues.txt file
#' @param pvaue pvalue Threshold.
#' @export
#'
CPDBHeatmaps <- function(meta_file=NULL,
                         pvalues_file=NULL,
                         count_filename="count_heatmap.pdf",
                         log_filename="logcount_heatmap.pdf",
                         count_network_filename="count_network.txt",
                         interaction_count_filename="interaction_count.txt",
                         count_network_separator="\t",
                         interaction_count_separator="\t",
                         show_rownames = T,
                         show_colnames = T,
                         scale="none",
                         cluster_cols = T,
                         border_color='white',
                         cluster_rows = T,
                         fontsize_row=11,
                         fontsize_col = 11,
                         main = '',
                         treeheight_row=0,
                         family='Arial',
                         treeheight_col = 0,
                         col1 = "dodgerblue4",
                         col2 = 'peachpuff',
                         col3 = 'deeppink4',
                         meta_sep='\t',
                         pvalues_sep='\t',
                         pvalue=0.05){
  #######   Network
  require(heatmap)
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)

  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]


  split_sep = '\\|'
  join_sep = '|'

  pairs1_all = unique(meta[,2])

  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))

  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]

    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))

    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }

  all_count = all_count[-1,]
  write.table(all_count, count_network_filename, sep=count_network_separator, quote=F, row.names = F)

  #######   count interactions

  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]

    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))

  }

  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])

    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    write.table(all_sum, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=F)

    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )

    pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)

    pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}
