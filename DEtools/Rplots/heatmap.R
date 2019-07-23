library("heatmaply")
library("RColorBrewer")
library("edgeR")

### Notes: Script arguments: 1) <input matrix>.mat 2) <cs list of the groups> 3) path/to/project_folder 4) baseName 5) percentage/top_n_genes 6) title of heatmap within quotation marks
### Example Command: Rscript visualisation.R matfile.mat cell,cell,cell,exosomes,exosomes,exosomes /projects/project_name mature_sense_minExpr5_RCadj 0.1 "Top Genes"
args <- commandArgs(TRUE)
json_data <- fromJSON(file="/Users/ernesto/PycharmProjects/conDE/upload/AA88/plot_config.json")
# json_data <- fromJSON(file=args[1])
matfile <- read.delim(json_data[["input_matrix"]], header=TRUE, row.names=1)   # Input the input delimited text file containing the count matrix
groups <- unlist(strsplit( json_data[["matrixDesc"]], ","))  # Sample description
sampletypevalues <- rev(unique(groups))  # Getting the group levels
outdir <- json_data[["folder"]]
FC <- as.numeric(json_data[["FC"]])
descending <- json_data[["descending"]]
basename <- ""
top_n <- as.numeric(json_data[["top_n"]])  # Cutoff for plotting top genes
heatmap_title <- json_data[["title"]]
setwd(outdir)

sort_genes_av <- function(m, desc=TRUE){
  m[c("FoldChange","log2FoldChange","pvalue","padj")]<-NULL
  if (desc){
    return(m[ order(-rowMeans(m)), ])
  }else{
    return(m[ order(rowMeans(m)), ])
  }
}

sort_genes_cv <- function(m, desc=TRUE){
  m[c("FoldChange","log2FoldChange","pvalue","padj")]<-NULL
  if (desc){
    m = m[order(-apply(m, 1, function(x) sd(x)/mean(x))),]
  }else{
    m = m[order(-apply(m, 1, function(x) sd(x)/mean(x))),]
  }
}

if(json_data[["folder"]] == "Average"){
  m<-sort_genes_av(matfile, descending)
}else if(json_data[["folder"]] == "CV"){
  m<-sort_genes_cv(matfile, descending)
}else{
  m<-matfile
  m[c("FoldChange","log2FoldChange","pvalue","padj")]<-NULL
}

data <- DGEList(counts=m, group=factor(groups))  # Summarise the input data



# Getting log2 of cpm data, 2 is being added to raw data
log2counts <- cpm(data$counts, prior.count=2, log=TRUE)
if(!is.null(top_n)){
  log2counts<-head(log2counts,top_n)
}
top_genes <- log2counts

# Colouring the different conditions
if (length(sampletypevalues) == 2) {
  col_condition <- c("#7FC97F", "#BEAED4")[data$samples$group]
} else {
  col_condition <- c(brewer.pal(n = length(sampletypevalues), name = "Accent"))[data$samples$group]
}



# Creating the heatmaps. If samples are more than 100, then sample labels are being omitted
if (length(groups) <= 100) {
	# Heatmap of top selected normalised genes
	heatmaply(top_genes, file = paste(outdir,"heatmap.html",sep ="/"),
          limits = NULL, colors = brewer.pal(11,"Spectral"), scale = "row", main = heatmap_title, 
          key.title=NULL, col_side_colors = data.frame(groups), hide_colorbar = FALSE, 
          column_text_angle=60, fontsize_col = 9, fontsize_row = 8,  showticklabels=c(TRUE,TRUE))
} else {

# Heatmap of top selected genes, No sample-names 
heatmaply(top_genes, file = paste(outdir,"heatmap.html",sep ="/"),
          limits = NULL, colors = brewer.pal(11, "Spectral"), scale = "row", main = heatmap_title, 
          key.title=NULL, col_side_colors = data.frame(groups), hide_colorbar = FALSE, 
          showticklabels = c(FALSE, TRUE), fontsize_row = 8)
}

