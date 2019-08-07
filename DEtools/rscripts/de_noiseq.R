### Stavros Giannoukakos ###

### Notes: Script arguments: 1) <input matrix>.mat 2) <cs list of the groups> 3) path/to/project_folder 4) baseName
### Example Command: Rscript NOISeq.R matfile.mat cell,cell,cell,exosomes,exosomes,exosomes /projects/project_name


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
  noiseq_prob <- 0  # Default differential expression cutoff
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# matfile <- read.delim("/Users/stavris/R/projects/example_project/example_matrix.mat", header=TRUE, row.names=1)
# groups <- unlist(strsplit("cell,cell,cell,exosomes,exosomes,exosomes", ","))
# sampletypevalues <- rev(unique(groups))
# outdir <- "/Users/stavris/R/projects/example_project"
# basename <- "mature_sense_minExpr5_RCadj"
# noiseq_prob <- 0

# Initiating packrat environment and
# loading all necessary libraries
library("NOISeq")
setwd(outdir)

# Redirecting all output to a log file
sink(paste(outdir, "/noiseq_", basename, ".log", sep=""), append = TRUE)

print("NOISeq Differential Expression analysis of RNA-Seq expression profiles is now running...")
  
# Designing the data's factors which indicate the experimental group for each sample
samplefactors <- data.frame(row.names=colnames(matfile), condition = factor(groups, levels=sampletypevalues))
# Importing all necessary information into a NOISeq object
countdata <- readData(data = matfile, factors = samplefactors)
# Trimmed Mean of M values (TMM) normalisation to correct the sequencing depth bias
# No length is taken into account
TMMvalues = tmm(assayData(countdata)$exprs, long = 1000, lc = 0)

for(i in 1:(length(sampletypevalues)-1)) {
    for(j in (i+1):length(sampletypevalues)){ 
      
        # Computing the differential expression between experimental conditions from the filtered read count data
        NOISeq <- noiseq(countdata, k=0.5, lc=0, norm="tmm", factor="condition", conditions=c(sampletypevalues[i], sampletypevalues[j]))

        # Extract the table containing the test results for each gene of the original count table
        curresult <- NOISeq@results[[1]]
        
        if(sum(is.na(curresult))>0){
          print("WARNING: NA values appear in NOISeq result, probably due to 0 counts in both samples")
          print("Rows with NA values will be omitted...")}
        
        # Selecting the samples
        selected_samples <- (which(groups==sampletypevalues[i] | groups==sampletypevalues[j]))
        # Obtaining the list of genes with probability higher than the user-input value (default 0.8)
        selected <- rownames(degenes(NOISeq, q = noiseq_prob, M = NULL))
          
        selected <- na.omit(selected)
        # Combing the normalised data along with statistical analysis results ("M", "prob", "padj")
        TMMvalues_selected <- as.data.frame(cbind(TMMvalues[, selected_samples], curresult$M, curresult$prob, (1-curresult$prob)))
        # Naming the new columns
        colnames(TMMvalues_selected) <- c(head(colnames(TMMvalues_selected),n=-3), "log2FoldChange", "prob", "padj")
        
        # Generating a matrix containing the mean normalised results per group per (ALL) gene
        mean_ncounts_selected <- as.data.frame(cbind(mean_Agroup=rowMeans(TMMvalues_selected[ ,which(groups==sampletypevalues[j])]), 
                                                     mean_Bgroup=rowMeans(TMMvalues_selected[ ,which(groups==sampletypevalues[i])]),
                                                     TMMvalues_selected[,c("log2FoldChange", "prob", "padj")]))
        colnames(mean_ncounts_selected)[1:2] <- c(paste("mean_",sampletypevalues[j],sep=""), paste("mean_",sampletypevalues[i],sep=""))
        # Insert FoldChange calculations
        mean_ncounts_selected <- as.data.frame(cbind(mean_ncounts_selected[,1:2],
                                                     transform(mean_ncounts_selected[0], FoldChange = (mean_ncounts_selected[ ,2]+1) / (mean_ncounts_selected[ ,1]+1)),
                                                     mean_ncounts_selected[,3:5]))
        
        
        
        # Inserting FoldChange calculations
        TMMvalues_selected <- as.data.frame(merge(TMMvalues_selected[ ,selected_samples], mean_ncounts_selected[ ,c("FoldChange", "log2FoldChange", "prob", "padj")],  by="row.names"))
        row.names(TMMvalues_selected) <- TMMvalues_selected$Row.names  # row names manipulation
        TMMvalues_selected$Row.names <- NULL  # row names manipulation
        
        # Order dataframe based on the FDR
        mean_ncounts_selected <- mean_ncounts_selected[order(mean_ncounts_selected$padj), ]
        
        # Obtaining the final matrix of selected genes
        result <- TMMvalues_selected[selected, ]
        
        # Exporting the normalised results table containing ALL genes
        print(paste("Exporting NOISeq normalised table containing all genes to: ", outdir,"/allGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(TMMvalues_selected), TMMvalues_selected), file=paste(outdir,"/allGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
        # Exporting the generated matrix containing the mean normalised results per group per (ALL) gene
        print(paste("Exporting NOISeq mean normalised results per group per (all) gene to: ", outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_meanGroupsAllGenes.csv", sep=""))
        write.table(data.frame("name"=rownames(mean_ncounts_selected), mean_ncounts_selected), file=paste(outdir,"/",basename,"_",sampletypevalues[j],"VS",sampletypevalues[i],"_noiseq_meanGroupsAllGenes.csv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
} 

## NOISeq DE analysis finished ##
print("NOISeq Differential Expression analysis has finished successfully!")

### NOISeq analysis is now being terminated...
print(" ---------------------------------------------- NOISeq ANALYSIS FINISHED ---------------------------------------------- ")

# Deactivating sink(). Revert output back to the console
sink()
