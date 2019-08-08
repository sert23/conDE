### Stavros Giannoukakos ###

### Notes: Script arguments: 1) <input matrix>.mat 2) <cs list of the groups> 3) path/to/project_folder 4) baseName 
### Example Command: Rscript de_deseq.R matfile.mat cell,cell,cell,exosomes,exosomes,exosomes /projects/project_name mature_sense_minExpr5_RCadj
#/opt/local/R-3.5.3/bin/Rscript /opt/sRNAtoolboxDB/r/de_deseq.R mature_sense_minExpr1_RCadj.mat exosome,exosome,exosome,exosome,exosome,exosome,cell,cell,cell,cell,cell,cell /opt/sRNAtoolbox_prod/sRNAtoolboxweb/upload/32E6QU5QY1P9I2V/de/de_deseq2 mature_sense_minExpr1_RCadj


# Input arguments and error control
args <- commandArgs(TRUE)
if (length(args) == 4) {
  if (!file.exists(args[1])) {
    cat("ERROR - The input matrix does NOT exist...\nEXITING!\n")
    quit()
  }
  # Input the input delimited text file containing the count matrix
  matfile <- read.delim(args[1], header=TRUE, row.names=1)   
  groups <- unlist(strsplit(args[2], ","))  # Sample description
  sampletypevalues <- unique(groups)  # Getting the group levels
  if (!dir.exists(args[3])) {
    cat("ERROR - The output directory does NOT exist...\nEXITING!\n")
    quit()
  }
  outdir <- args[3]  # Output directory
  basename <- args[4]  # Base name
  pvalue <- 1  # Default differential expression cutoff
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# matfile <- read.delim("/Users/stavris/R/projects/example_project/example_matrix.mat", header=TRUE, row.names=1)
# groups <- unlist(strsplit("cell,cell,cell,exosomes,exosomes,exosomes", ","))
# sampletypevalues <- unique(groups)
# outdir <- "/Users/stavris/R/projects/example_project"
# basename <- "mature_sense_minExpr5_RCadj"
# pvalue <- 0.05


# Initiating packrat environment and
# loading all necessary libraries
library("DESeq")
library("vsn")
library("ggplot2")
setwd(outdir)

# Redirecting all output to a log file
sink(file=paste(outdir, "/deseq_", basename, ".log", sep=""), append = TRUE)

print("DESeq Differential Expression analysis of RNA-seq expression profiles is now running...")

# Collecting the necessary information for the DESeq analysis
cds <- newCountDataSet(matfile, groups)
# Estimate the effective library size. The following function
# estimates the size factors from the count data
cds <- estimateSizeFactors(cds)
# Obtaining the 'normalised' matrix (counts/geometric mean)
normalised_counts <- counts(cds, normalized=TRUE)

# Exporting the 'normalised' table
# print(paste("Exporting DESeq normalised table containing all genes to: ", outdir,"/",basename,"_deseq_normTable.csv", sep=""))
# write.table(data.frame("name"=rownames(normalised_counts), normalised_counts), file=paste(outdir,"/",basename,"_deseq_RLEnormTable.csv", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Check whether there are replicates for the two sample types and estimating the dispersions
# This function estimates the typical relationship between the data’s variance and their mean (data’s dispersion and their mean)
if(sum(duplicated(groups))<2) {
  cds = estimateDispersions(cds, method='blind', sharingMode="fit-only", fitType="local")
} else {
  cds = estimateDispersions(cds, fitType="local")
}


for(i in 1:(length(sampletypevalues)-1)) {
    for(j in (i+1):length(sampletypevalues)){
      
        # Differential gene expression analysis based on the negative binomial distribution
        print(paste("Performing a test for differential expression between:", sampletypevalues[i], "and", sampletypevalues[j]))
        curresult <- nbinomTest(cds, sampletypevalues[i], sampletypevalues[j])
        
        if(sum(is.na(curresult))>0) {
          print("WARNING: missing values (NA) appear in DESeq result, probably due to 0 counts in both samples")
          print("Rows with NA values will be omitted...")
        }
        
        # Selecting the samples
        selected_samples <- (which(groups==sampletypevalues[i] | groups==sampletypevalues[j]))
        # Obtaining the list of genes with absolute log2FoldChange greater than 1 and p adjusted value lower than the user-input value
        selected <- (curresult$id[(abs(curresult$log2FoldChange)>=1) & (curresult$padj<=pvalue)])
        # Removing rows with missing values on columns
        selected <- na.omit(selected)
        # Combing the normalised data along with statistical analysis results ("foldChange", "log2FoldChange", "pval", "padj")
        ncounts_selected <- as.data.frame(cbind(normalised_counts[ ,selected_samples], curresult$foldChange, curresult$log2FoldChange, curresult$pval, curresult$padj))
        # Naming the new columns
        colnames(ncounts_selected) <- c(head(colnames(ncounts_selected), n=-4),"FoldChange","log2FoldChange","pvalue","padj")
        # Obtaining the final matrix of selected genes
        result <- ncounts_selected[selected, ]

        # Exporting the normalised results table containing ALL Genes
        print(paste("Exporting DESeq normalised table containing all genes to: ", outdir,"/allGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(ncounts_selected), ncounts_selected), file=paste(outdir,"/allGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        # Generating a matrix containing the mean normalised results per group per (ALL) gene
        mean_ncounts_selected <- as.data.frame(cbind(mean_Agroup=rowMeans(ncounts_selected[ ,which(groups==sampletypevalues[i])]), 
                                                     mean_Bgroup=rowMeans(ncounts_selected[ ,which(groups==sampletypevalues[j])]),
                                                     ncounts_selected[,c("FoldChange","log2FoldChange","pvalue","padj")]))
        colnames(mean_ncounts_selected)[1:2] <- c(paste("mean_",sampletypevalues[i],sep=""), paste("mean_",sampletypevalues[j],sep=""))
        # Order dataframe based on the padj
        mean_ncounts_selected <- mean_ncounts_selected[order(mean_ncounts_selected$padj),]
        # Exporting the generated matrix containing the mean normalised results per group per (ALL) gene
        print(paste("Exporting DESeq mean normalised results per group per (all) gene to: ", outdir,"/",basename,"_",sampletypevalues[i],"VS",sampletypevalues[j],"_deseq_meanGroupsAllGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(mean_ncounts_selected), mean_ncounts_selected), file=paste(outdir,"/",basename,"_",sampletypevalues[i],"VS",sampletypevalues[j],"_deseq_meanGroupsAllGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
}

## DESeq DE analysis finished ##
print("DESeq Differential Expression analysis has finished successfully!")

### DESeq analysis is now being terminated...
print(" ---------------------------------------------- DESeq ANALYSIS FINISHED ---------------------------------------------- ")

# Deactivating sink(). Revert output back to the console
sink()
