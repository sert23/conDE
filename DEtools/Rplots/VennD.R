library(venn)
library(rjson)
library(Cairo)

args <- commandArgs(TRUE)
x11 = function (...) grDevices::x11(...,type='cairo')


# json_data <- fromJSON(file="/Users/ernesto/PycharmProjects/conDE/upload/AA88/config.json")
json_data <- fromJSON(file=args[1])
methods <- json_data[["methods"]]
folder <- json_data[["folder"]]
FC <- as.numeric(json_data[["FC"]])
pval <- as.numeric(json_data[["pval"]])
oset <- json_data[["set"]]


select_genes <- function(ifile, FC, pval, outset ){
  
  df <- read.delim(ifile, stringsAsFactors=FALSE)
  if (FC<1){
    FC<- 1/FC
  }
  sbs <- subset(df, (df$pvalue < pval ) & ((df$FoldChange>FC)|(df$FoldChange<1/FC)))
  if (outset=="Under"){
    sbs<-subset(sbs, (df$FoldChange < 1 ) )
  }
  if (outset=="Over"){
    sbs<-subset(sbs, (df$FoldChange > 1 ) )
  }
  return(sbs$name)
  # return(genes)
}

l = list()
v = c()
i <- 1
for (method in methods) {
  current_string<-method
  input_file <- paste(folder,"de",method, "allGenes.csv", sep="/") 
  current_vector <-select_genes(input_file,FC,pval,oset)
  # print(length(current_vector))
  print(current_vector)
  fake_list<- list()
  fake_list[[method]]<- current_vector
  l<- append(l,fake_list)
  # l[[ val]] <- current_vector
  v <- c(v,toString(method))
  # print(l[[val]])
  i <- i + 1
}

out_path <- paste(folder,"plots","Venn.jpg", sep="/")
CairoJPEG(out_path, width = 8, height = 5, units = 'in', res = 300)
# jpeg(out_path, width = 8, height = 5, units = 'in', res = 300)
venn(l,ilab=TRUE, zcolor = "style")
zeroset <- matrix(1000*c(0,1,1,0,0,0,0,1,1,0), ncol = 2)
lines(zeroset, col='white', lwd=5)
dev.off()