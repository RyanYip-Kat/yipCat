#' function to save plotlist
#' @param plotlist ggplot object list or complex heatmap object list
#' @param name filename to save plot
#' @param width image width
#' @param height image height
#' @param outpath out dir to save image
#' @param addDOC whether to add date message
#' @export
MyplotPDF<-function (..., name = "Plot", width = 6, height = 6, outpath = NULL,
                     addDOC = TRUE, useDingbats = FALSE, plotList = NULL)
{
  #require("ggplot2")
  #require("gridExtra")
  #require("grid")
  #require("ComplexHeatmap")

  if (is.null(plotList)) {
    plotList <- list(...)
    plotList2 <- list()
    for (i in seq_along(plotList)) {
      if (inherits(plotList[[i]], "list")) {
        for (j in seq_along(plotList[[i]])) {
          plotList2[[length(plotList2) + 1]] <- plotList[[i]][[j]]
        }
      }
      else {
        plotList2[[length(plotList2) + 1]] <- plotList[[i]]
      }
    }
    plotList <- plotList2
    rm(plotList2)
    gc()
  }
  else {
    plotList2 <- list()
    for (i in seq_along(plotList)) {
      if (inherits(plotList[[i]], "list")) {
        for (j in seq_along(plotList[[i]])) {
          plotList2[[length(plotList2) + 1]] <- plotList[[i]][[j]]
        }
      }
      else {
        plotList2[[length(plotList2) + 1]] <- plotList[[i]]
      }
    }
    plotList <- plotList2
    rm(plotList2)
    gc()
  }
  name <- gsub("\\.pdf", "", name)
  if (is.null(outpath)) {
    outDir <- "Plots"
  }
  else {
    outDir <- file.path(outpath, "Plots")
  }
  if(!dir.exists(outDir)){
    dir.create(outDir,recursive=TRUE)
  }

  if (addDOC) {
    doc <- gsub(":", "-", stringr::str_split(Sys.time(),
                                             pattern = " ", simplify = TRUE)[1, 2])
    filename <- file.path(outDir, paste0(name, "_Date-",
                                         Sys.Date(), "_Time-", doc, ".pdf"))
  }
  else {
    filename <- file.path(outDir, paste0(name, ".pdf"))
  }
  o <- tryCatch({
    pdf(filename, width = width, height = height, useDingbats = useDingbats)
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
  return(invisible(0))
}

#' set plot size
#' @export
fixPlotSize <- function(
  p = NULL,
  plotWidth = unit(6, "in"),
  plotHeight = unit(6, "in"),
  margin = 0.25,
  height = 1,
  it = 0.05,
  newPage = FALSE
){

  require(grid)
  require(gridExtra)
  if(!inherits(plotWidth, "unit")){
    plotWidth <- unit(plotWidth, "in")
  }

  if(!inherits(plotHeight, "unit")){
    plotHeight <- unit(plotHeight, "in")
  }

  #adapted from https://github.com/jwdink/egg/blob/master/R/set_panel_size.r
  g <- ggplotGrob(p)

  legend <- grep("guide-box", g$layout$name)
  if(length(legend)!=0){
    gl <- g$grobs[[legend]]
    g <- ggplotGrob(p + theme(legend.position = "none"))
  }else{
    gl <- NULL
    g <- ggplotGrob(p)
  }

  panels <- grep("panel", g$layout$name)
  panel_index_w <- unique(g$layout$l[panels])
  panel_index_h <- unique(g$layout$t[panels])

  nw <- length(panel_index_w)
  nh <- length(panel_index_h)

  pw <- grid::convertWidth(plotWidth, unitTo = "in", valueOnly = TRUE)
  ph <- grid::convertWidth(plotHeight, unitTo = "in", valueOnly = TRUE)

  pw <- pw * 0.95
  ph <- ph * 0.95

  x <- 0
  width <- 1
  sm <- FALSE

  while(!sm){

    x <- x + it

    w <- unit(x * width, "in")
    h <- unit(x * height / width, "in")
    m <- unit(x * margin / width, "in")

    g$widths[panel_index_w] <-  rep(w, nw)
    g$heights[panel_index_h] <- rep(h, nh)

    sw <- grid::convertWidth(
      x = sum(g$widths)  + m,
      unitTo = "in",
      valueOnly = TRUE
    )

    sh <- grid::convertHeight(
      x = sum(g$heights) + m,
      unitTo = "in",
      valueOnly = TRUE
    )

    sm <- sw > pw | sh > ph

  }

  if(length(legend)!=0){

    sgh <- grid::convertHeight(
      x = sum(g$heights),
      unitTo = "in",
      valueOnly = TRUE
    )

    sgw <- grid::convertWidth(
      x = sum(g$widths),
      unitTo = "in",
      valueOnly = TRUE
    )

    slh <- grid::convertHeight(
      x = sum(gl$heights),
      unitTo = "in",
      valueOnly = TRUE
    )

    slw <- grid::convertWidth(
      x = sum(gl$widths),
      unitTo = "in",
      valueOnly = TRUE
    )

    size <- 15
    wh <- 0.1
    it <- 0

    while(slh > 0.2 * ph | slw > pw){

      it <- it + 1

      if(it > 3){
        break
      }

      size <- size * 0.8
      wh <- wh * 0.8

      gl <- ggplotGrob(
        p + theme(
          legend.key.width = unit(wh, "cm"),
          legend.key.height = unit(wh, "cm"),
          legend.spacing.x = unit(0, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.text = element_text(size = max(size, 2))
        ) + guides(fill = guide_legend(ncol = 4), color =  guide_legend(ncol = 4))
      )$grobs[[legend]]

      slh <- convertHeight(
        x = sum(gl$heights),
        unitTo = "in",
        valueOnly = TRUE
      )

      slw <- grid::convertWidth(
        x = sum(gl$widths),
        unitTo = "in",
        valueOnly = TRUE
      )

    }

    p <- grid.arrange(g, gl, ncol=1, nrow=2,
                      heights = unit.c(unit(sgh,"in"), unit(min(slh, 0.2 * pw), "in")),
                      newpage = newPage
    )

  }else{

    p <- grid.arrange(g, newpage = newPage)

  }


  invisible(p)

}


#########################
#' archr theme
#' a function
#' @export
theme_ArchR<-function (color = "black", textFamily = "sans", baseSize = 10,
                       baseLineSize = 0.5, baseRectSize = 0.5, plotMarginCm = 1,
                       legendPosition = "bottom", legendTextSize = 5, axisTickCm = 0.1,
                       xText90 = FALSE, yText90 = FALSE)
{
  theme <- theme_bw() + theme(text = element_text(family = textFamily),
                              axis.text = element_text(color = color, size = baseSize),
                              axis.title = element_text(color = color, size = baseSize),
                              title = element_text(color = color, size = baseSize),
                              plot.margin = unit(c(plotMarginCm, plotMarginCm, plotMarginCm,
                                                   plotMarginCm), "cm"), panel.background = element_rect(fill = "transparent",
                                                                                                         colour = NA), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA,
                                                                                              color = color, size = (4/3) * baseRectSize * as.numeric(grid::convertX(grid::unit(1,
                                                                                                                                                                                "points"), "mm"))), axis.ticks.length = unit(axisTickCm,
                                                                                                                                                                                                                             "cm"), axis.ticks = element_line(color = color, size = baseLineSize *
                                                                                                                                                                                                                                                                (4/3) * as.numeric(grid::convertX(grid::unit(1, "points"),
                                                                                                                                                                                                                                                                                                  "mm"))), legend.key = element_rect(fill = "transparent",
                                                                                                                                                                                                                                                                                                                                     colour = NA), legend.text = element_text(color = color,
                                                                                                                                                                                                                                                                                                                                                                              size = legendTextSize), legend.box.background = element_rect(color = NA),
                              legend.position = legendPosition, strip.text = element_text(size = baseSize,
                                                                                          color = "black"))
  if (xText90) {
    theme <- theme %+replace% theme(axis.text.x = element_text(angle = 90,
                                                               hjust = 1))
  }
  if (yText90) {
    theme <- theme %+replace% theme(axis.text.y = element_text(angle = 90,
                                                               vjust = 1))
  }
  return(theme)
}


###########################
#' List of color palettes that can be used in plots
#'
#' A collection of some original and some borrowed color palettes to provide appealing color aesthetics for plots in ArchR
#'
#' @export
ArchRPalettes <- list(

  #DISCLOSURE: This is a collection of palettes that includes some original palettes and some palettes originally
  #implemented by others in other packages.
  #They are included here for convenience because they help improve plot aesthetics.

  #PS: all palettes included in the "Primarily Continuous Palettes" section should also work for discrete usage but not vice versa.
  #Each continuous palette has been ordered by color to generate a visually appealing discrete palette.

  #---------------------------------------------------------------
  # Primarily Discrete Palettes
  #---------------------------------------------------------------

  #20-colors
  stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D"),

  stallion2 = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767"),

  calm = c("1"="#7DD06F", "2"="#844081", "3"="#688EC1", "4"="#C17E73", "5"="#484125", "6"="#6CD3A7", "7"="#597873","8"="#7B6FD0", "9"="#CF4A31", "10"="#D0CD47",
           "11"="#722A2D", "12"="#CBC594", "13"="#D19EC4", "14"="#5A7E36", "15"="#D4477D", "16"="#403552", "17"="#76D73C", "18"="#96CED5", "19"="#CE54D1", "20"="#C48736"),

  kelly = c("1"="#FFB300", "2"="#803E75", "3"="#FF6800", "4"="#A6BDD7", "5"="#C10020", "6"="#CEA262", "7"="#817066", "8"="#007D34", "9"="#F6768E", "10"="#00538A",
            "11"="#FF7A5C", "12"="#53377A", "13"="#FF8E00", "14"="#B32851", "15"="#F4C800", "16"="#7F180D", "17"="#93AA00", "18"="#593315", "19"="#F13A13", "20"="#232C16"),

  #16-colors
  bear = c("1"="#faa818", "2"="#41a30d","3"="#fbdf72", "4"="#367d7d",  "5"="#d33502", "6"="#6ebcbc", "7"="#37526d",
           "8"="#916848", "9"="#f5b390", "10"="#342739", "11"="#bed678","12"="#a6d9ee", "13"="#0d74b6",
           "14"="#60824f","15"="#725ca5", "16"="#e0598b"),

  #15-colors
  ironMan = c("9"='#371377',"3"='#7700FF',"2"='#9E0142',"10"='#FF0080', "14"='#DC494C',"12"="#F88D51","1"="#FAD510","8"="#FFFF5F","4"='#88CFA4',
              "13"='#238B45',"5"="#02401B", "7"="#0AD7D3","11"="#046C9A", "6"="#A2A475", "15"='grey35'),

  circus = c("1"="#D52126", "2"="#88CCEE", "3"="#FEE52C", "4"="#117733", "5"="#CC61B0", "6"="#99C945", "7"="#2F8AC4", "8"="#332288",
             "9"="#E68316", "10"="#661101", "11"="#F97B72", "12"="#DDCC77", "13"="#11A579", "14"="#89288F", "15"="#E73F74"),

  #12-colors
  paired = c("9"="#A6CDE2","1"="#1E78B4","3"="#74C476","12"="#34A047","11"="#F59899","2"="#E11E26",
             "10"="#FCBF6E","4"="#F47E1F","5"="#CAB2D6","8"="#6A3E98","6"="#FAF39B","7"="#B15928"),

  #11-colors
  grove = c("11"="#1a1334","9"="#01545a","1"="#017351","6"="#03c383","8"="#aad962","2"="#fbbf45","10"="#ef6a32","3"="#ed0345","7"="#a12a5e","5"="#710162","4"="#3B9AB2"),

  #7-colors
  summerNight = c("1"="#2a7185", "2"="#a64027", "3"="#fbdf72","4"="#60824f","5"="#9cdff0","6"="#022336","7"="#725ca5"),

  #5-colors
  zissou = c("1"="#3B9AB2", "4"="#78B7C5", "3"="#EBCC2A", "5"="#E1AF00", "2"="#F21A00"), #wesanderson
  darjeeling = c("1"="#FF0000", "2"="#00A08A", "3"="#F2AD00", "4"="#F98400", "5"="#5BBCD6"), #wesanderson
  rushmore = c("1"="#E1BD6D", "5"="#EABE94", "2"="#0B775E", "4"="#35274A" , "3"="#F2300F"), #wesanderson
  captain = c("1"="grey","2"="#A1CDE1","3"="#12477C","4"="#EC9274","5"="#67001E"),

  #---------------------------------------------------------------
  # Primarily Continuous Palettes
  #---------------------------------------------------------------

  #10-colors
  horizon = c("1"='#000075',"4"='#2E00FF', "6"='#9408F7', "10"='#C729D6', "8"='#FA4AB5', "3"='#FF6A95', "7"='#FF8B74', "5"='#FFAC53', "9"='#FFCD32', "2"='#FFFF60'),

  #9-colors
  horizonExtra =c("1"="#000436","4"="#021EA9","6"="#1632FB","8"="#6E34FC","3"="#C732D5","9"="#FD619D","7"="#FF9965","5"="#FFD32B","2"="#FFFC5A"),
  blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D"),
  sambaNight = c("6"='#1873CC',"2"='#1798E5',"8"='#00BFFF',"5"='#4AC596',"1"='#00CC00',"4"='#A2E700',"9"='#FFFF00',"7"='#FFD200',"3"='#FFA500'), #buencolors
  solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D'),  #buencolors
  whitePurple = c("9"='#f7fcfd',"6"='#e0ecf4',"8"='#bfd3e6',"5"='#9ebcda',"2"='#8c96c6',"4"='#8c6bb1',"7"='#88419d',"3"='#810f7c',"1"='#4d004b'),
  whiteBlue = c("9"='#fff7fb',"6"='#ece7f2',"8"='#d0d1e6',"5"='#a6bddb',"2"='#74a9cf',"4"='#3690c0',"7"='#0570b0',"3"='#045a8d',"1"='#023858'),
  whiteRed = c("1"="white", "2"="red"),
  comet = c("1"="#E6E7E8","2"="#3A97FF","3"="#8816A7","4"="black"),

  #7-colors
  greenBlue = c("4"='#e0f3db',"7"='#ccebc5',"2"='#a8ddb5',"5"='#4eb3d3',"3"='#2b8cbe',"6"='#0868ac',"1"='#084081'),

  #6-colors
  beach = c("4"="#87D2DB","1"="#5BB1CB","6"="#4F66AF","3"="#F15F30","5"="#F7962E","2"="#FCEE2B"),

  #5-colors
  coolwarm = c("1"="#4858A7", "4"="#788FC8", "5"="#D6DAE1", "3"="#F49B7C", "2"="#B51F29"),
  fireworks = c("5"="white","2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  greyMagma = c("2"="grey", "4"="#FB8861FF", "5"="#B63679FF", "3"="#51127CFF", "1"="#000004FF"),
  fireworks2 = c("5"="black", "2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  purpleOrange = c("5"="#581845", "2"="#900C3F", "4"="#C70039", "3"="#FF5744", "1"="#FFC30F")
)

paletteDiscrete <- function(
  values = NULL,
  set = "stallion",
  reverse = FALSE
){


  values <- gtools::mixedsort(values)
  n <- length(unique(values))
  pal <- ArchRPalettes[[set]]
  palOrdered <- pal[gtools::mixedsort(names(pal))] #mixed sort gets 1,2,3,4..10,11,12

  if(n > length(palOrdered)){
    message("Length of unique values greater than palette, interpolating..")
    palOut <- colorRampPalette(pal)(n)
  }else{
    palOut <- palOrdered[seq_len(n)]
  }

  if(reverse){
    palOut <- rev(palOut)
  }

  names(palOut) <- unique(values)

  return(palOut)

}

#'continuous paltette color
#' @param set The name of a color palette provided in the `ArchRPalettes` list object.
#' @param n The number of unique colors to generate as part of this continuous color palette.
#' @param reverse A boolean variable that indicates whether to return the palette colors in reverse order.
#' @export
paletteContinuous <- function(
  set = "solarExtra",
  n = 256,
  reverse = FALSE
){

  pal <- ArchRPalettes[[set]]
  palOut <- colorRampPalette(pal)(n)

  if(reverse){
    palOut <- rev(palOut)
  }

  return(palOut)

}

