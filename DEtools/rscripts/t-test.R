### Stavros Giannoukakos ###

### Notes: Script arguments: 1) <input matrix>.mat 2) <cs list of the groups> 3) path/to/project_folder 4) baseName 
### Example Command: Rscript t-test.R matfile.mat cell,cell,cell,exosomes,exosomes,exosomes /projects/project_name mature_sense_minExpr5_RCadj

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
}  else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}


# matfile <- read.delim("/Users/stavris/R/projects/example_project/example_matrix.mat", header=TRUE, row.names=1)
# groups <- unlist(strsplit("cell,cell,cell,exosomes,exosomes,exosomes", ","))
# sampletypevalues <- unique(groups)
# outdir <- "/Users/stavris/R/projects/example_project"
# basename <- "mature_sense_minExpr5_RCadj"


# Initiating packrat environment and
# loading all necessary libraries
set.seed(12345)
library("edgeR")
setwd(outdir)

# Redirecting all output to a log file
sink(paste(outdir, "/ttest_", basename, ".log", sep=""), append = TRUE)

print("Student t-test is now running...")

# Performing CPM transformation on the raw data
transf_table <- cpm(matfile, log=F, prior.count=0)

# Applying the paired t-test 
stat_value <- as.data.frame(apply(transf_table, 1, function(x)
                            t.test(x[which(groups==sampletypevalues[1])], x[which(groups==sampletypevalues[2])], paired=T)$p.value))

stat_value <- cbind(stat_value, p.adjust(stat_value[ ,1], method="BH"))

colnames(stat_value) <- c("pvalue", "padj")  # Naming the new columns

# Generating a matrix containing the mean normalised results per group per (ALL) gene
mean_ncounts_selected <- as.data.frame(cbind(mean_Agroup=rowMeans(transf_table[ ,which(groups==sampletypevalues[1])]), 
                                             mean_Bgroup=rowMeans(transf_table[ ,which(groups==sampletypevalues[2])])))
                                             
colnames(mean_ncounts_selected)[1:2] <- c(paste("mean_",sampletypevalues[1],sep=""), paste("mean_",sampletypevalues[2],sep=""))

# Inserting pvalue and padjusted
mean_ncounts_selected <- as.data.frame(merge(mean_ncounts_selected, stat_value,  by="row.names"))
row.names(mean_ncounts_selected) <- mean_ncounts_selected$Row.names  # row names manipulation
mean_ncounts_selected$Row.names <- NULL  # row names manipulation

# Insert FoldChange and log2FoldChange  calculations
mean_ncounts_selected <- as.data.frame(cbind(mean_ncounts_selected[ ,1:2], 
                                             FoldChange = ((mean_ncounts_selected[ ,2]+1) / (mean_ncounts_selected[ ,1]+1)),
                                             log2FoldChange = log2(((mean_ncounts_selected[ ,2]+1) / (mean_ncounts_selected[ ,1]+1))),
                                             mean_ncounts_selected[ ,3:4]))

# Generating the normalised data along with statistical analysis results ("FoldChange", "log2FoldChange", "pval", "padj")
data <- as.data.frame(merge(transf_table, mean_ncounts_selected[ ,3:6],  by="row.names"))
row.names(data) <- data$Row.names  # row names manipulation
data$Row.names <- NULL  # row names manipulation


# Exporting the normalised results table containing ALL Genesde
print(paste("Exporting t-test normalised table containing all genes to: ", outdir,"/",basename,"_",sampletypevalues[1],"VS",sampletypevalues[2],"_ttest_allGenes.csv", sep=""))
write.table(data.frame("name"=rownames(data), data), file=paste(outdir,"/",basename,"_",sampletypevalues[1],"VS",sampletypevalues[2],"_ttest_allGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Order dataframe based on the padj
mean_ncounts_selected <- mean_ncounts_selected[order(mean_ncounts_selected$padj), ]

# Exporting the generated matrix containing the mean normalised results per group per (ALL) gene
print(paste("Exporting t-test mean normalised results per group per (all) gene to: ", outdir,"/",basename,"_",sampletypevalues[1],"VS",sampletypevalues[2],"_ttest_meanGroupsAllGenes.csv", sep=""))
write.table(data.frame("name"=rownames(mean_ncounts_selected), mean_ncounts_selected), file=paste(outdir,"/",basename,"_",sampletypevalues[1],"VS",sampletypevalues[2],"_ttest_meanGroupsAllGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


## Student t-testing analysis finished
print("Student t-testing analysis has finished successfully!")

### Student t-testing analysis is now being terminated...
print(" ---------------------------------------------- Student t-testing ANALYSIS FINISHED ---------------------------------------------- ")

# Deactivating sink(). Revert output back to the console
sink()