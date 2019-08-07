### Stavros Giannoukakos ###
## Limma is a package for the analysis of gene expression data arising from microarray or RNA-Seq technologies.
## In the limma approach to RNA-seq, read counts are converted to log2-counts-per-million (logCPM) and the mean-variance
## relationship is modelled either with precision weights or with an empirical Bayes prior trend. The
## precision weights approach is called "voom" and the prior trend approach is called "limma-trend"

### Notes: Script arguments: 1) <input matrix>.mat 2) <cs list of the groups> 3) path/to/project_folder 4) baseName 
### Example Command: Rscript de_limma-voom.R matfile.mat cell,cell,cell,exosomes,exosomes,exosomes /projects/project_name mature_sense_minExpr5_RCadj

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
  sampletypevalues <- rev(unique(groups))  # Getting the group levels
  if (!dir.exists(args[3])) {
    cat("ERROR - The output directory does NOT exist...\nEXITING!\n")
    quit()
  }
  outdir <- args[3]  # Output directory
  basename <- args[4]  # Base name
  pvalue <- 1  # Default differential expression cutoff
}  else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# # 3 groups 
# matfile <- read.delim("/Users/stavris/R/projects/example_project/example_matrix_3.mat", header=TRUE, row.names=1)
# groups <- unlist(strsplit("cell,cell,cell,exosomes,exosomes,exosomes,platelets,platelets,platelets", ","))
# # 2 groups
# matfile <- read.delim("/Users/stavris/R/projects/example_project/example_matrix.mat", header=TRUE, row.names=1)
# groups <- unlist(strsplit("cell,cell,cell,exosomes,exosomes,exosomes", ","))
# # Rest of the parameters
# sampletypevalues <- unique(groups)
# outdir <- "/Users/stavris/R/projects/example_project"
# basename <- "mature_sense_minExpr5_RCadj"
# pvalue <- 1


# Initiating packrat environment and
# loading all necessary libraries
library("limma")
library("edgeR")
setwd(outdir)

# Redirecting all output to a log file
sink(paste(outdir, "/limma-trend_", basename, ".log", sep=""), append = TRUE)

print("Limma-voom Differential Expression analysis of RNA-Seq expression profiles is now running...")

# Storing the raw read counts table in a simple list-based data object called a DGEList.
limma_table <- DGEList(counts=matfile, group=factor(groups))  # From now on, all the necessary info will be stored in this variable

# Define the design matrix based on the experimental design
design <- model.matrix(~ 0 + factor(groups))
colnames(design) <- sampletypevalues

# The filterByExpr function keeps rows that have worthwhile counts in a minumum number of samples
keep <- filterByExpr(limma_table, design)  # removing rows that consistently have zero or very low counts
limma_table <- limma_table[keep, ,keep.lib.sizes=F] # Discarting lowly expressed genes

# Normalisation for RNA composition by finding a set of scaling factors for the library sizes 
# that minimize the log-fold changes between the samples for most genes. The default method
# for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
limma_table <- calcNormFactors(limma_table)



### DIFFERENCIAL EXPRESSION: limma-voom
# When the library sizes are quite variable between samples,
# then the voom approach is theoretically more powerful 
voom_transf <- voom(limma_table, design, plot=F) # the voom transformation is applied to thenormalized and filtered DGEList

z <- voom_transf$E # Obtaining the voom normalised log2 counts

# lmFit fits a linear model using weighted least squares for each gene
voom_fit <- lmFit(voom_transf, design)


# Empirical Bayes smoothing of standard errors 
# (shrinks standard errors that are much larger or smaller
# than those from other genes towards the average standard error)
voom_fit <- eBayes(voom_fit)

# Extracting the differentially expressed genes
t <- topTable(voom_fit, coef=ncol(design), p.value=pvalue, lfc=0, sort.by = "P", n = Inf) # Obtaining the statistical results

# Merging the above info into one table
data <- as.data.frame(merge(z, t, by="row.names"))
row.names(data) <- data$Row.names  # row names manipulation
data$Row.names <- NULL  # row names manipulation


# Combing the normalised data along with statistical analysis results ("log2FoldChange", "pval", "padj")
selected_samples <- (which(groups==sampletypevalues[1] | groups==sampletypevalues[2]))
data <- as.data.frame(cbind(data[ , selected_samples], data$logFC, data$P.Value, data$adj.P.Val))
# Naming the new columns
colnames(data) <- c(head(colnames(data), n=-3), "log2FoldChange", "pvalue", "padj")

# Generating a matrix containing the mean normalised results per group per (ALL) gene
mean_ncounts_selected <- as.data.frame(cbind(mean_Agroup=rowMeans(data[ ,which(groups==sampletypevalues[1])]), 
                                             mean_Bgroup=rowMeans(data[ ,which(groups==sampletypevalues[2])]),
                                             data[,c("log2FoldChange", "pvalue", "padj")]))
colnames(mean_ncounts_selected)[1:2] <- c(paste("mean_",sampletypevalues[1], sep=""), paste("mean_",sampletypevalues[2], sep=""))
# Inserting FoldChange calculations
mean_ncounts_selected <- as.data.frame(cbind(mean_ncounts_selected[,1:2],
                                             transform(mean_ncounts_selected[0], FoldChange = (mean_ncounts_selected[ ,2]+1) / (mean_ncounts_selected[ ,1]+1)),
                                             mean_ncounts_selected[,3:5]))
# Inserting FoldChange calculations
data <- as.data.frame(merge(data[ ,selected_samples], mean_ncounts_selected[ ,c("FoldChange", "log2FoldChange", "pvalue", "padj")],  by="row.names"))
row.names(data) <- data$Row.names  # row names manipulation
data$Row.names <- NULL  # row names manipulation

# Exporting the normalised results table containing ALL Genesde
print(paste("Exporting Limma-Voom normalised table containing all genes to: ", outdir,"/",basename,"_",sampletypevalues[1],"VS",sampletypevalues[2],"_limmaVoom_allGenes.csv", sep=""))
write.table(data.frame("name"=rownames(data), data), file=paste(outdir,"/",basename,"_",sampletypevalues[1],"VS",sampletypevalues[2],"_limmaVoom_allGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Order dataframe based on the padj
mean_ncounts_selected <- mean_ncounts_selected[order(mean_ncounts_selected$padj), ]

# Exporting the generated matrix containing the mean normalised results per group per (ALL) gene
print(paste("Exporting Limma-Voom mean normalised results per group per (all) gene to: ", outdir,"/",basename,"_",sampletypevalues[1],"VS",sampletypevalues[2],"_limmaVoom_meanGroupsAllGenes.csv", sep=""))
write.table(data.frame("name"=rownames(mean_ncounts_selected), mean_ncounts_selected), file=paste(outdir,"/",basename,"_",sampletypevalues[1],"VS",sampletypevalues[2],"_limmaVoom_meanGroupsAllGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

##  DE analysis finished ##
print("Limma-voom Differential Expression analysis has finished successfully!")

### Limma-voom analysis is now being terminated...
print(" ---------------------------------------------- Limma-voom ANALYSIS FINISHED ---------------------------------------------- ")

# Deactivating sink(). Revert output back to the console
sink()



