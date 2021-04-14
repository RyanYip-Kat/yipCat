#' combine 10X contig files
#' @param metaDF meta csv file,must include Path,Sample,Batch columns.Path is the path of filtered_contig_annotations.csv.
#' @param pattern contig is PCR or TCR
#' @examples  data("pCRexample")
#' Combine10XContig(pCRexample[["meta"]])
#' str(pCRexample[["combined"]])
#' @export
Combine10XContig<-function(metaDF=NULL,
                        pattern="TCR"){
  metaDF$ID<-1:nrow(metaDF)
  stopifnot("Path"%in%colnames(metaDF))
  stopifnot("Sample"%in%colnames(metaDF))
  stopifnot("Batch"%in%colnames(metaDF))
  stopifnot(pattern%in%c("TCR","BCR"))

  message("reformat contig file...")
  contig_list<-lapply(1:length(metaDF$Path),function(i){
    cat(sprintf("reformat NO. [ %d/%d ]\n",i,nrow(metaDF)))

    file<-metaDF$Path[i]
    contig<-read.csv(file,stringsAsFactors=FALSE)
    barcode<-unlist(lapply(contig$barcode,function(b)return(str_split(b,"-")[[1]][1])))
    barcode<-paste(barcode,i,sep="-")
    contig$barcode<-barcode

    contig_id<-unlist(lapply(contig$contig_id,function(b)return(str_split(b,"_")[[1]][2])))
    contig_id<-paste(barcode,contig_id,sep="_")
    contig$contig_id<-contig_id
    return(contig)
  })
  message("combined contig list...")
  if(pattern=="TCR"){
    combined <- combineTCR(contig_list, samples = metaDF$Sample, ID =metaDF$Batch , cells ="T-AB")
  }else{
    combined<-combineBCR(contig_list,samples = metaDF$Sample, ID =metaDF$Batch)
  }
  message("add Variable..")
  combined<- addVariable(combined, name = "batch", variables = metaDF$Batch)
  return(combined)
}

#' Adding clonotype information to a seurat or SCE object
#' This function adds the immune receptor information to the seurat or SCE object to the meta data. By defualt this function also
#' calculates the frequencies of the clonotypes by sequencing run
#' (groupBy = "none").  To change how the frequencies are calculated,
#' select a column header for the groupBy variable. Importantly,
#' before using combineExpression() ensure the barcodes of the seurat
#' or SCE object match the barcodes in the output of the
#' combinedContig() call. Check changeNames() to change the prefix of
#' the seurat object. If the dominant clonotypes have a greater
#' frequency than 500, adjust the cloneTypes variable.
#' @param seurat seurat object
#' @param combined object come from Combine10XContig
#' @param groupby the column label in the combined contig object in which clonotype frequency will be calculated.
#' @param cloneCall  How to call the clonotype - CDR3 gene (gene), CDR3 nucleotide (nt) CDR3 amino acid (aa), or CDR3 gene+nucleotide (gene+nt).
#' @export
contigWrapperSeurat<-function(seurat=NULL,
                              combined=NULL,
                              groupby="sample",
                              cloneCall="gene+nt"){

  stopifnot(cloneCall%in%c("gene+nt","aa","gene","nt"))
  message("combineExpression with seurat...")
  seurat <- combineExpression(combined, seurat, cloneCall=cloneCall, groupBy = groupby)

  slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType,
                                                levels = c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)",
                                                           "Medium (5 < X <= 20)", "Small (1 < X <= 5)",
                                                           "Single (0 < X <= 1)", NA))
  return(seurat)
}
