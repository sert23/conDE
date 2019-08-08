### Stavros Giannoukakos ###

### Notes: Script arguments: 1) <input matrix>.mat 2) <cs list of the groups> 3) path/to/project_folder 4) baseName
### Example Command: Rscript edgeR.R matfile.mat cell,cell,cell,exosomes,exosomes,exosomes /projects/project_name mature_sense_minExpr5_RCadj


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
# names(matfile) <- gsub(x = names(matfile), pattern = "\\.", replacement = " ")
# groups <- unlist(strsplit("cell,cell,cell,exosomes,exosomes,exosomes", ","))
# sampletypevalues <- unique(groups)
# outdir <- "/Users/stavris/R/projects/example_project"
# basename <- "mature_sense_minExpr5_RCadj"
# pvalue <- 0.05

# Initiating packrat environment and 
# loading all necessary libraries
library("edgeR")
setwd(outdir)

# Redirecting all output to a log file
sink(paste(outdir, "/edger_", basename, ".log", sep=""), append = TRUE)

print("EdgeR Differential Expression analysis of RNA-Seq expression profiles is now running...")
# Storing the raw read counts table in a simple list-based data object called a DGEList.
edgeR_table <- DGEList(counts=matfile, group=factor(groups))  # From now on, all the necessary info will be stored in this variable
# Normalisation for RNA composition by finding a set of scaling factors for the library sizes 
# that minimize the log-fold changes between the samples for most genes. The default method
# for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
edgeR_table <- calcNormFactors(edgeR_table)
# The first major step in the analysis of DGE  data using the NB model is to estimate the dispersion
# parameter for each tag, a measure of the degree of inter-library variation for that tag. 
# Estimating the common dispersion gives an idea of overall variability across the genome for this dataset.
edgeR_table <- estimateCommonDisp(edgeR_table)
# For routine differential expression analysis, we use empirical Bayes tagwise dispersions. 
# Note that common dispersion needs to be estimated before estimating tagwise dispersions.
edgeR_table <- estimateTagwiseDisp(edgeR_table)

for(i in 1:(length(sampletypevalues)-1)) {
    for(j in (i+1):length(sampletypevalues)){
        
        print(paste("Performing an exact test between:", sampletypevalues[i], "and", sampletypevalues[j]))
        # Performing an exact test for the difference in expression between each pair of conditions. 
        # More explicitly, computing gene-wise exact tests for differences in the means between two 
        # groups of negative-binomially distributed counts.
        dgeTest <- exactTest(edgeR_table, pair=c(sampletypevalues[i], sampletypevalues[j]))
        # The test results for the most significant tags are conveniently in the following function.
        topTags <- topTags(dgeTest, n=Inf)
        # Selecting the samples
        selected_samples <- (which(groups==sampletypevalues[i] | groups==sampletypevalues[j]))
        # Obtaining the pseudo-counts for each sample
        z <- as.data.frame(edgeR_table$pseudo.counts)[ ,selected_samples]
        # Obtaining the statistical results from the edge test
        t <- as.data.frame(topTags)
        # Merging the above info into one table
        data <- as.data.frame(merge(z, t, by="row.names"))
        row.names(data) <- data$Row.names  # row names manipulation
        data$Row.names <- NULL  # row names manipulation
        # Combing the normalised data along with statistical analysis results ("log2FoldChange", "pval", "padj")
        data <- as.data.frame(cbind(data[ , selected_samples], data$logFC, data$PValue, data$FDR))
        # Naming the new columns
        colnames(data) <- c(head(colnames(data), n=-3), "log2FoldChange", "pvalue", "padj")
        # Selecting only genes with padj lower than the input p-value
        selected <- which(data$padj<=pvalue) 
        
        # Generating a matrix containing the mean normalised results per group per (ALL) gene
        mean_ncounts_selected <- as.data.frame(cbind(mean_Agroup=rowMeans(data[ ,which(groups==sampletypevalues[i])]), 
                                                     mean_Bgroup=rowMeans(data[ ,which(groups==sampletypevalues[j])]),
                                                     data[,c("log2FoldChange", "pvalue", "padj")]))
        colnames(mean_ncounts_selected)[1:2] <- c(paste("mean_",sampletypevalues[i],sep=""), paste("mean_",sampletypevalues[j],sep=""))
        # Inserting FoldChange calculations
        mean_ncounts_selected <- as.data.frame(cbind(mean_ncounts_selected[,1:2],
                                                     transform(mean_ncounts_selected[0], FoldChange = (mean_ncounts_selected[ ,2]+1) / (mean_ncounts_selected[ ,1]+1)),
                                                     mean_ncounts_selected[,3:5]))
        # Inserting FoldChange calculations
        data <- as.data.frame(merge(data[ ,selected_samples], mean_ncounts_selected[ ,c("FoldChange", "log2FoldChange", "pvalue", "padj")],  by="row.names"))
        row.names(data) <- data$Row.names  # row names manipulation
        data$Row.names <- NULL  # row names manipulation
        
        # Order dataframe based on the padj
        mean_ncounts_selected <- mean_ncounts_selected[order(mean_ncounts_selected$padj), ]
        
        # Obtaining the final matrix of selected genes
        result <- data[selected, ] 
        
        # Exporting the normalised results table containing ALL genes
        print(paste("Exporting edgeR normalised table containing all genes to: ", outdir,"/allGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(data), data), file=paste(outdir,"/allGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        # Exporting the generated matrix containing the mean normalised results per group per (ALL) gene
        print(paste("Exporting edgeR mean normalised results per group per (all) gene to: ", outdir,"/",basename,"_",sampletypevalues[i],"VS",sampletypevalues[j],"_edger_meanGroupsAllGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(mean_ncounts_selected), mean_ncounts_selected), file=paste(outdir,"/",basename,"_",sampletypevalues[i],"VS",sampletypevalues[j],"_edger_meanGroupsAllGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
}

## edgeR DE analysis finished ##
print("edgeR Differential Expression analysis has finished successfully!")

### edgeR analysis is now being terminated...
print(" ---------------------------------------------- edgeR ANALYSIS FINISHED ---------------------------------------------- ")

# Deactivating sink(). Revert output back to the console
sink()