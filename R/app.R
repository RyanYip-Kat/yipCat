#' Rename one ident to a new name
#' @param obj seurat object
#' @param oldname old ident name you want to rename
#' @param newname new name you want to use
#' @examples
#' x<-IdentRename(obj=pbmc_small,oldname = "0",newname = "C0")
#' table(Idents(x))
#' @export
IdentRename<-function(obj=NULL,oldname=NULL,newname=NULL){
  cluster.ids <-levels(obj)
  cluster.ids[which(cluster.ids  == oldname)] <- newname
  names(cluster.ids) <- levels(obj)
  obj <- RenameIdents(obj, cluster.ids)
  obj
}

#' re use column as Idents
#' @param obj seurat object
#' @param useCol which column as Idents
#' @export
ReIdent<-function(obj=NULL,useCol="seurat_clusters"){
  metaData<-obj@meta.data
  if(!useCol%in%colnames(metaData)){
    stop("Invalid column use!!!")
  }
  col<-metaData[[useCol]]
  Idents(obj)<-col
  obj
}


isColor <- function(x = NULL) {
  unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)),
                                        error = function(e) FALSE)))
}


isDiscrete<-function (x = NULL)
{
  is.factor(x) || is.character(x) || is.logical(x)
}

#' save image
#' @param plotList plot object list
#' @param filename filename to save the plot list
#' @param width image width
#' @param height image height
#' @export
saveImage<-function(plotList,filename,width=12,height=16){
  o <- tryCatch({
    height=12
    width=16
    pdf(filename, width = width, height = height, useDingbats = FALSE)
    for (i in seq_along(plotList)) {
      if (inherits(plotList[[i]], "gg")) {
        message("Plotting Ggplot!")
        if (!is.null(attr(plotList[[i]], "ratioYX"))) {
          fixPlotSize(plotList[[i]], plotWidth = width,
                      plotHeight = height, height = attr(plotList[[i]],
                                                         "ratioYX"), newPage = FALSE)
        }
        else {
          fixPlotSize(plotList[[i]], plotWidth = width,
                      plotHeight = height, newPage = FALSE)
        }
        if (i != length(plotList)) {
          grid::grid.newpage()
        }
      }
      else if (inherits(plotList[[i]], "gtable")) {
        message("Plotting Gtable!")
        print(grid::grid.draw(plotList[[i]]))
        if (i != length(plotList)) {
          grid::grid.newpage()
        }
      }
      else if (inherits(plotList[[i]], "HeatmapList") |
               inherits(plotList[[i]], "Heatmap")) {
        message("Plotting ComplexHeatmap!")
        padding <- 15
        draw(plotList[[i]], padding = unit(c(padding,
                                             padding, padding, padding), "mm"), heatmap_legend_side = "bot",
             annotation_legend_side = "bot")
      }
      else {
        message("Plotting Other")
        print(plotList[[i]])
      }
    }
    dev.off()
  }, error = function(x) {
    message(x)
  })
}

#' a simple Shiny App for Visualization
#' @param obj seurat object
#' @param browserTheme shiny browser theme
#' @param host host ip address
#' @param port ip port
#' @param launch.browser whether launch browser
#' @example
#' setAppServer(pbmc_small,host = NULL,launch.browser = TRUE)
#' @export
setAppServer <- function(
  obj = NULL,
  browserTheme = "cosmo",
  host="10.100.44.33",
  port=5013,
  launch.browser=FALSE
){


  require("rhandsontable")
  require("DT")
  #Determine Grouping Methods
  ccd <- obj@meta.data
  discreteCols <- lapply(seq_len(ncol(ccd)), function(x){
    isDiscrete(ccd[, x])
  }) %>% unlist %>% {colnames(ccd)[.]}
  if("Sample" %in% discreteCols){
    selectCols <- "Sample"
  }else{
    selectCols <- "seurat_clusters"
  }


  #####################
  #Shiny App UI
  theme <- shinythemes::shinytheme(browserTheme)
  ui <- fluidPage(
    theme = theme,
    titlePanel(
      h1(div(HTML("<b>Shiny Seurat Report Browser </b>")), align = "left")
    ),
    sidebarLayout(
      sidebarPanel(
        column(8,h2(div(HTML(paste0("<b>nCell : ",ncol(obj),"</b>"))),align="left")),
        column(8,h2(div(HTML(paste0("<b>nGene : ",nrow(obj),"</b>"))),align="left")),
      ),
      mainPanel(
        tabsetPanel(
          imageOutput("p1", width = "100%", height = "400px", click = NULL,
                      dblclick = NULL, hover = NULL,
                      brush = NULL,
                      inline = FALSE),
          tabPanel("FeaturePlot",
                   sidebarLayout(sidebarPanel(
                     column(8,selectInput("nameF",
                                          label = "Gene Symbol",
                                          choices = rownames(obj),
                                          multiple = FALSE,
                                          selected = rownames(obj)[1]
                     )),
                     column(8,selectInput("embF",
                                          label = "Embedding name",
                                          choices = Reductions(obj),
                                          multiple = FALSE,
                                          selected = Reductions(obj)[1])),
                     column(8,selectInput("continuousSetF",
                                          label="continuous color set",
                                          choices=names(ArchRPalettes),
                                          multiple=FALSE,
                                          selected="solarExtra")),
                     column(8,sliderInput("sizeF", "Plot point size:", min = 0.1, max = 2.0, value = 0.2,step = 0.1)),
                     column(8,downloadButton("downF", "Download",
                                             style="color: #fff; background-color: green; border-color: Black;"))),
                     mainPanel(plotOutput("featurePlot"),width = "1200px",height = "1200px"))),
          tabPanel("ReductionPlot",
                   sidebarLayout(sidebarPanel(
                     column(8,selectInput("groupR",
                                          label = "groupBy",
                                          choices = discreteCols,
                                          multiple = FALSE,
                                          selected = selectCols)),
                     column(8,selectInput("embR",
                                          label = "Embedding name",
                                          choices = Reductions(obj),
                                          multiple = FALSE,
                                          selected = Reductions(obj)[1])),
                     column(8,selectInput("discreteSetR",
                                          label="discrete  color set",
                                          choices = names(ArchRPalettes),
                                          multiple = FALSE,
                                          selected="stallion")),
                     column(8,sliderInput("sizeR", "Plot point size:", min = 0.1, max = 2.0, value = 0.2,step = 0.1)),
                     column(8,downloadButton("downR", "Download",
                                             style="color: #fff; background-color: green; border-color: Black;"))),
                   mainPanel(plotOutput("embPlot"),width = "1200px",height = "1200px"))),
          tabPanel("GroupPlot",
                   sidebarLayout(sidebarPanel(
                     column(8,selectInput("groupG",
                                          label = "groupBy",
                                          choices = discreteCols,
                                          multiple = FALSE,
                                          selected = selectCols)),
                     column(8,selectInput("nameG",
                                          label = "Gene Symbol",
                                          choices = rownames(obj),
                                          multiple = FALSE,
                                          selected = rownames(obj)[1]
                     )),
                     column(8,selectInput("discreteSetG",
                                          label="discrete  color set",
                                          choices = names(ArchRPalettes),
                                          multiple = FALSE,
                                          selected="stallion")),
                     column(8,selectInput("plotAS",
                                          label = "plotAS",
                                          choices = c("ridges","violin"),
                                          multiple = FALSE,
                                          selected = "violin")),
                     column(8,downloadButton("downG", "Download",
                                             style="color: #fff; background-color: green; border-color: Black;"))
                   ),
                   mainPanel(plotOutput("groupPlot"),width = "1200px",height = "1200px"))),
          tabPanel("Heatmap",
                   sidebarLayout(sidebarPanel(
                     #h4("Please Select Gene List CSV file"),
                     #fileInput("file1","Choose Gene List File",
                     #           accept = c("text/csv/tsv"))),
                     column(8,selectInput("groupH",
                                          label = "groupBy",
                                          choices = discreteCols,
                                          multiple = FALSE,
                                          selected = selectCols)),
                     column(8,selectInput("slot",label="slot",choices = c("data","scale.data"),selected="scale.data")),
                     column(8,checkboxInput("scale","scale",value = FALSE)),
                     column(8,numericInput("N","nTop Genes",value = 30,min=5,max=100,step = 1)),
                     column(8,h4("Limits value")),
                     column(8,numericInput("Up","Limit Up",value = 2.5,min = 1.0,max = 10,step = 0.5)),
                     column(8,numericInput("Down","Limit Down",value = -2.5,min = -10,max = -1,step = 0.5)),
                     column(8,h4("Whether Do Heatmap")),
                     column(8,checkboxInput("Doheatmap","Do Heatmap",value = FALSE)),
                     column(8,h4("Donwload Heatmap Plot")),
                     column(8,downloadButton("downH", "Download",
                                    style="color: #fff; background-color: green; border-color: Black;"))),

                     mainPanel(plotOutput("heatmap")))),

          tabPanel('table',
                   h3("Object MetaData"),
                   dataTableOutput('metadata'),
                   h3("Crossdata"),
                   sidebarLayout(sidebarPanel(selectInput("Column1",
                                             label="label1",
                                             choices = discreteCols,
                                             multiple = FALSE,
                                             selected = discreteCols[1]),

                                 selectInput("Column2",
                                             label="label2",
                                             choices = discreteCols,
                                             multiple = FALSE,
                                             selected = discreteCols[2]),
                                 downloadButton("downT", "Download")),
                                 #column(2, style = "margin-top: 200px;",
                                 #       downloadButton("down3", "Download",
                                 #                       style="color: #fff; background-color: red; border-color: Black;")),
                                 #downloadButton("down3", "Download"),
                                 mainPanel(
                                   fluidRow(dataTableOutput('crosstable')))))
        )
      )
    )
  )

  plotList<-list()
  #Shiny app server
  server<-function(input,output,session){
    output$embPlot <- renderPlot({

      p<-MyplotEmbedding(obj,
                      assay=DefaultAssay(obj),
                      embedding = input$embR,
                      colorBy = "metadata",
                      name=input$groupR,
                      size = input$sizeR,
                      discreteSet=input$discreteSetR,
                      continuousSet = NULL)
      print(p)
      plotList<-list(p)
      output$downR <- downloadHandler(

        filename <- function(){
          doc <- gsub(":", "-", stringr::str_split(Sys.time(),
                                                   pattern = " ", simplify = TRUE)[1, 2])
          paste0(input$embR,"-",input$groupR,"-EmbeddingPlot", "_Date-",Sys.Date(), "_Time-", doc, ".pdf")
        },
        content = function(file) {
          saveImage(plotList,filename=file,width = 12,height = 16)
        })
    },height = 800, width =800)


    output$groupPlot<-renderPlot({
      p<-MyplotGroups(seurat=obj,
                   assay=DefaultAssay(obj),
                   groupBy = input$groupG,
                   name=input$nameG,
                   discreteSet=input$discreteSetG,
                   plotAs=input$plotAS)
      print(p)
      plotList<-list(p)
      output$downG <- downloadHandler(

        filename <- function(){
          doc <- gsub(":", "-", stringr::str_split(Sys.time(),
                                                   pattern = " ", simplify = TRUE)[1, 2])
          paste0(input$groupG,"-",input$nameG,"-GroupPlot", "_Date-",Sys.Date(), "_Time-", doc, ".pdf")
        },
        content = function(file) {
          saveImage(plotList,filename=file,width = 12,height = 10)
        })
    })

    output$featurePlot<-renderPlot({

      p<-MyplotEmbedding(obj,
                      assay=DefaultAssay(obj),
                      embedding = input$embF,
                      colorBy = "matrix",
                      name=input$nameF,
                      size = input$sizeF,
                      discreteSet=NULL,
                      continuousSet = input$continuousSetF)
      print(p)
      plotList<-list(p)
      output$downF <- downloadHandler(

        filename <- function(){
          doc <- gsub(":", "-", stringr::str_split(Sys.time(),
                                                   pattern = " ", simplify = TRUE)[1, 2])
          paste0(input$embF,"-",input$nameF,"-featurePlot", "_Date-",Sys.Date(), "_Time-", doc, ".pdf")
        },
        content = function(file) {
          saveImage(plotList,filename=file,width = 12,height = 16)
        })

    },height = 800, width =800)

    #output$contents <-  DT::renderDataTable({
    #  inFile <- input$file1
    #  if (is.null(inFile)){
    #    return(NULL)
    #  }else{
    #    pb<- read.csv(inFile$datapath)
    #    pb
    #    req(pb)

    #  }
    #})
    output$heatmap<-renderPlot({

      #inFile <- input$file1
      #if(is.null(inFile)){
      #  geneSet<-NULL
      #}else{
      #  geneSet<-as.character(read.csv(inFile$datapath,sep=",",stringsAsFactors = FALSE)$V1)
      #}
      if(input$Doheatmap){

        limits<-c(input$Down,input$Up)
        scale <- input$scale
        slot <- input$slot
        p<-DoLikeArchRHeatmap(seurat=obj,groupBy = input$groupH,slot=slot,scale=scale,limits = limits,N = input$N)
        print(p)

        plotList<-list(p)
        output$downH <- downloadHandler(

          filename <- function(){
            doc <- gsub(":", "-", stringr::str_split(Sys.time(),
                                                     pattern = " ", simplify = TRUE)[1, 2])
            paste0(input$groupH,"-",input$N,"-Heatmap", "_Date-",Sys.Date(), "_Time-", doc, ".pdf")
          },
          content = function(file) {
            saveImage(plotList,filename=file,width = 12,height = 16)
          })
      }
    })

    output$metadata<-renderDataTable({
      obj@meta.data
      #as.data.frame(as.matrix(table(obj@meta.data[[input$Column1]],obj@meta.data[[input$Column2]])))
    },options = list(pageLength = 10))

    output$crosstable<-renderDataTable({
      m<-as.matrix(table(obj@meta.data[[input$Column1]],obj@meta.data[[input$Column2]]))
      mm<-matrix(array(m),nrow = nrow(m),ncol=ncol(m))
      rownames(mm)<-rownames(m)
      colnames(mm)<-colnames(m)
      as.data.frame(mm)
    })
    output$downT <- downloadHandler(

      filename <- function(){
        doc <- gsub(":", "-", stringr::str_split(Sys.time(),
                                                 pattern = " ", simplify = TRUE)[1, 2])
        paste0(input$Column1,"-",input$Column2,"-CrossTable", "_Date-",Sys.Date(), "_Time-", doc, ".csv")
      },
      content = function(file) {
        m<-as.matrix(table(obj@meta.data[[input$Column1]],obj@meta.data[[input$Column2]]))
        mm<-matrix(array(m),nrow = nrow(m),ncol=ncol(m))
        rownames(mm)<-rownames(m)
        colnames(mm)<-colnames(m)
        d<-as.data.frame(mm)
        write.table(d,file=file,sep=",",quote=FALSE)
      })
  }

  shiny::runApp(list(ui = ui, server = server),host=host,port=port, launch.browser = launch.browser)
}

