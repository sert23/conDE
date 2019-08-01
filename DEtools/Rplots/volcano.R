library(rjson)
library(manhattanly)
library(magrittr)


#rewrite volcanoly function
 if(TRUE){
volcanoly2 <- function(x,
                      col = c("#252525"),
                      point_size = 5,
                      effect_size_line = c(-1,1),
                      effect_size_line_color = "grey",
                      effect_size_line_width = 0.5,
                      effect_size_line_type = 2,
                      genomewideline = -log10(1e-5),
                      genomewideline_color = "grey",
                      genomewideline_width = 0.5,
                      genomewideline_type = 2,
                      highlight = NULL,
                      highlight_color = "red",
                      xlab = NULL,
                      ylab = "-log10(p)",
                      title = "Volcano Plot", ...) {
  
  UseMethod("volcanoly")
  
}

#' @export
volcanoly.default <- function(x,
                              col = c("#252525"),
                              point_size = 5,
                              effect_size_line = c(-1,1),
                              effect_size_line_color = "grey",
                              effect_size_line_width = 0.5,
                              effect_size_line_type = 2,
                              genomewideline = -log10(1e-5),
                              genomewideline_color = "grey",
                              genomewideline_width = 0.5,
                              genomewideline_type = 2,
                              highlight = NULL,
                              highlight_color = "red",
                              xlab = NULL,
                              ylab = "-log10(p)",
                              title = "Volcano Plot", ...) {
  
  mh <- volcanor(x, ...)
  volcanoly.volcanor(mh,
                     col = col,
                     point_size = point_size,
                     effect_size_line = effect_size_line,
                     effect_size_line_color = effect_size_line_color,
                     effect_size_line_width = effect_size_line_width,
                     effect_size_line_type = effect_size_line_type,
                     genomewideline = genomewideline,
                     genomewideline_color = genomewideline_color,
                     genomewideline_width = genomewideline_width,
                     genomewideline_type = genomewideline_type,
                     highlight = highlight,
                     highlight_color = highlight_color,
                     xlab = xlab,
                     ylab = ylab,
                     title = title)
}


#' @export
volcanoly.volcanor <- function(x,
                               col = c("#252525"),
                               point_size = 5,
                               effect_size_line = c(-1,1),
                               effect_size_line_color = "grey",
                               effect_size_line_width = 0.5,
                               effect_size_line_type = 2,
                               genomewideline = -log10(1e-5),
                               genomewideline_color = "grey",
                               genomewideline_width = 0.5,
                               genomewideline_type = 2,
                               highlight = NULL,
                               highlight_color = "red",
                               xlab = NULL,
                               ylab = "-log10(p)",
                               title = "Volcano Plot",
                               ...) {
  
  
  d <- x$data
  head(d)
  pName <- x$pName
  log10pName <- "LOG10P"
  effectName <- x$effectName
  snpName <- x$snpName
  geneName <- x$geneName
  annotation1Name <- x$annotation1Name
  annotation2Name <- x$annotation2Name
  labs <- x$labs
  xlabel <- x$xlabel
  # print(effect_size_line)
  
  if (!is.null(highlight) & is.na(snpName)) stop("You're trying to highlight snps, but havent provided a snp column")
  if (!is.logical(effect_size_line)) {
    if (length(effect_size_line) < 2) stop("'effect_size_line' must be a numeric vector of length 2")
    if (length(effect_size_line) > 2) message("More than two values provided to 'effect_size_line'. Only the first two elements will be used")
    if (effect_size_line[1] > effect_size_line[2]) stop("First element of 'effect_size_line' must be smaller than second element")
  }
  if (is.logical(effect_size_line)) {
    if (effect_size_line) stop("If effect_size_line is a logical, it must be set to FALSE")
  }
  if (is.logical(genomewideline)) {
    if (genomewideline) stop("If genomewideline is a logical, it must be set to FALSE")
  }
  if (is.null(highlight) & is.logical(effect_size_line) & is.logical(genomewideline)) 
    message("Since both effect_size_line and genomewideline are set to FALSE, no points will be highlighted")
  
  if (!is.null(highlight) && is.logical(highlight) && highlight) stop("'highlight' argument must be set to either NULL, FALSE, or a character vector of SNPs to highlight")
  # Initialize plot
  # xmax = ceiling(max(d$EFFECTSIZE) * 1.03)
  # xmin = floor(max(d$EFFECTSIZE) * -0.03)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Initalize plotly
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # rm(p)
  
  TEXT <- paste(paste0("test1"," ","meme"))
  
  p <- ggplot2::ggplot(d, ggplot2::aes_string(x = effectName, y = log10pName)) + 
    ggplot2::geom_point() + 
    ggplot2::theme_classic() + 
    ggplot2::labs(x = if(!is.null(xlab)) xlab else xlabel,
                  y = ylab,
                  title = title)
  
  if (!is.logical(effect_size_line) & !is.logical(genomewideline)) {
    
    eline1 <- ggplot2::geom_vline(xintercept = effect_size_line[1], 
                                  linetype = effect_size_line_type, 
                                  size = effect_size_line_width, 
                                  color = effect_size_line_color)
    eline2 <- ggplot2::geom_vline(xintercept = effect_size_line[2],
                                  linetype = effect_size_line_type,
                                  size = effect_size_line_width,
                                  color = effect_size_line_color)
    pline <- ggplot2::geom_hline(yintercept = genomewideline[1],
                                 linetype = genomewideline_type,
                                 size = genomewideline_width,
                                 color = genomewideline_color)
    
    p <- p + eline1 + eline2 + pline
  }
  
  if (is.logical(effect_size_line) & !is.logical(genomewideline)) {
    
    pline <- ggplot2::geom_hline(yintercept = genomewideline[1],
                                 linetype = genomewideline_type,
                                 size = genomewideline_width,
                                 color = genomewideline_color)
    
    p <- p + pline
  }
  
  
  if (!is.logical(effect_size_line) & is.logical(genomewideline)) {
    
    eline1 <- ggplot2::geom_vline(xintercept = effect_size_line[1], 
                                  linetype = effect_size_line_type, 
                                  size = effect_size_line_width, 
                                  color = effect_size_line_color)
    eline2 <- ggplot2::geom_vline(xintercept = effect_size_line[2],
                                  linetype = effect_size_line_type,
                                  size = effect_size_line_width,
                                  color = effect_size_line_color)
    
    p <- p + eline1 + eline2
  }
  

  
  p <- plotly::ggplotly(p)
  # p %<>%  
  #   plotly::ggplotly() %<>%  
  #   plotly::add_markers(
  #     marker = list(
  #       color = col,
  #       size = point_size,
  #       text = TEXT)) 
  
  
  if (is.null(highlight)) {
    if (!is.na(snpName)) {
      
      # Highlight snps automatically to be those greater than genomewideline and effect_size_line
      if ((is.null(highlight) & !is.logical(effect_size_line)) | (is.null(highlight) & !is.logical(genomewideline))) {
        
        # if both lines are provided
        if (!is.logical(effect_size_line) & !is.logical(genomewideline)) {
          
          highlight_index <- c(which((d$EFFECTSIZE < effect_size_line[1]) & (d$LOG10P > genomewideline)), 
                               which((d$EFFECTSIZE > effect_size_line[2]) & (d$LOG10P > genomewideline)))
          
        } else if (!is.logical(effect_size_line) & is.logical(genomewideline)) {
          
          # if only effect_size_line is provided
          highlight_index <- c(which(d$EFFECTSIZE < effect_size_line[1]), 
                               which(d$EFFECTSIZE > effect_size_line[2]))
          
        } else if (is.logical(effect_size_line) & !is.logical(genomewideline)) {
          
          # if only genomewideline is provided
          highlight_index <- which(d$LOG10P > genomewideline)
        }  
        
        if (length(highlight_index)==0) message("No points are beyond the effect_size_line or genomewideline, therefore no points will be highlighted")
        if (length(highlight_index)>0) {
          
          d.highlight <- d[highlight_index, ] 
          
          TEXT <- paste0("log2(FoldChange):",d$EFFECTSIZE,"<br>", 
                         "-log10(pval):",d$LOG10P, "<br>","name:", d$GENE , sep = " ")
          
          p %<>% plotly::add_trace(x = d$EFFECTSIZE, y = d$LOG10P,
                                   type = "scatter",
                                   mode = "markers",
                                   text = TEXT,
                                   marker = list(color = "black",
                                                 size = point_size),
                                   name = "")
          TEXT <- paste0("log2(FoldChange):",d.highlight$EFFECTSIZE,"<br>", 
                         "-log10(pval):",d.highlight$LOG10P, "<br>","name:", d.highlight$GENE,  sep = " ")
          
          #test1
          p %<>% plotly::add_trace(x = d.highlight$EFFECTSIZE, y = d.highlight$LOG10P,
                                   type = "scatter",
                                   mode = "markers",
                                   text = TEXT,
                                   marker = list(color = highlight_color,
                                                 size = point_size),
                                   name = "")
        }
        
      }
    }
  }  
  
  p
}
 }


args <- commandArgs(TRUE)
# json_data <- fromJSON(file="/Users/ernesto/PycharmProjects/conDE/upload/AA99/plot_config.json")
json_data <- fromJSON(file=args[1])
matfile <- read.delim(json_data[["input_matrix"]], header=TRUE, row.names=1)   # Input the input delimited text file containing the count matrix
groups <- unlist(strsplit( json_data[["matrixDesc"]], ","))  # Sample description
sampletypevalues <- rev(unique(groups))  # Getting the group levels
outdir <- json_data[["folder"]]
FC <- as.numeric(json_data[["FC"]])
descending <- json_data[["descending"]]
basename <- ""
top_n <- as.numeric(json_data[["top_n"]])  # Cutoff for plotting top genes
volcano_title <- json_data[["title"]]
pval <- as.numeric(json_data[["pval"]])

if (is.na(FC)){
  FC<-1
}

if (is.na(pval)){
  pval<-0.001
}

FC_vector <- c(-abs(FC),abs(FC))
log10pval <- -log10(pval)

colnames(matfile) = gsub("log2FoldChange","EFFECTSIZE",  colnames(matfile))
matfile$name<- rownames(matfile)
matfile$gene<- rownames(matfile)
matfile$GENE<- rownames(matfile)
volcano_obj<-volcanor(matfile,p="pvalue", snp= "name", gene="GENE")
print("FC_vector")
print(FC_vector)
p<-volcanoly2(volcano_obj, title=volcano_title, gene="GENE", 
             effect_size_line=FC_vector, genomewideline=log10pval, xlab = "log2(Fold Change)")

htmlwidgets::saveWidget(p, paste(outdir,"volcano.html",sep="/"))
