#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GO.db")
#BiocManager::install("DESeq2")
#BiocManager::install("goseq")
#BiocManager::install("org.Hs.eg.db")
#install.packages("dplyr")

library(dplyr)
library(ggplot2)
library(BiocManager)
library(GO.db)
library(DESeq2)
library(goseq)
library(org.Hs.eg.db)


#read in an example RNAseq file 
example_childhood_AML <- read.table("TARGET-NCI-AML-RNAseq-data/TARGET-20-PABGKN-09A-01R.gene.quantification.txt",
                                    sep="\t",
                                    header=T)

#import the list of RNAseq files
RNAseq_files <- list.files(path="TARGET-NCI-AML-RNAseq-data",
                           pattern="*.txt", 
                           full.names=TRUE, 
                           recursive=FALSE) 

#Note: switched to manually pre-processing data files, but run this if not pre-processed
#Only obtain files that are RNAseq gene quant.
#for (f in files){                          
#  if (grepl("gene", f, fixed=TRUE) == F){
#    file.remove(f)
#  }
#}

#initialize a master dataframe of all RNAseq raw counts by gene
RNAseq_df <- data.frame(example_childhood_AML$gene)
RNAseq_df <- RNAseq_df %>% rename(example_childhood_AML.gene = "gene")

#Obtain all RNAseq raw counts from each patient sample and add it to the master dataframe
for (file in RNAseq_files){
  file_df = read.table(file, sep = "\t",header=T)
  
  #Process the file names by removing the extraneous elements that do not match the names in the metadata
  #Remove the file path and .gene.txt
  file_name <- unlist(strsplit(file, split=".", fixed=TRUE))[1]
  file_name <- unlist(strsplit(file_name, split="/", fixed=TRUE))[2]
  #Remove -0#R from the string
  file_name <- unlist (strsplit(file_name, split="-", fixed=TRUE))
  file_name <- head(file_name, -1)
  file_name <- paste(file_name, collapse='-')
  
  #Add the raw counts column from each patient matching to the patient/sample name, only if it is not duplicated
  if (!(file_name %in% colnames(RNAseq_df))){
    RNAseq_df$raw_counts <- file_df$raw_counts
    RNAseq_df <- RNAseq_df %>% rename(raw_counts = file_name)
  }
}

#import the metadata file
AML_metadata <- read.table("TARGET_AML_mRNA-seq_metadata.txt",sep="\t",header=T)

#drop the columns we do not require and remove duplicate entries
AML_metadata <- dplyr::select(AML_metadata, 
                              -Source.Name, 
                              -Provider,
                              -Material.Type, 
                              -Characteristics.Organism., 
                              -Performer, 
                              -Protocol.REF)
AML_metadata <- unique(AML_metadata) 

#Keep only the metadata that we have RNAseq data for
AML_metadata <- AML_metadata[(AML_metadata$Sample.Name %in% colnames(RNAseq_df)), ]


#Creates new RNAseq dataframe with genes as rownames, and metadata dataframe with sample names as rownames
RNAseq_df_RN <- data.frame(RNAseq_df, row.names=1)
AML_metadata <- data.frame(AML_metadata, row.names=4)

#For some reason, when converting gene names in RNAseq_RN to row names, the "-" in sample names is replaced with "."
#This creates errors when making the deseq2 object. So the following code replaces the "." with "-" in the column names
new_names = list()
for (name in colnames(RNAseq_df_RN)){
  name <- gsub("\\.","-", name)
  new_names <- append(new_names, name)
}
colnames(RNAseq_df_RN) <- c(new_names)

#Subset the metadata to sample tissue: Bone Marrow vs. Peripheral Blood
AML_metadata_BM <- AML_metadata[AML_metadata$Characteristics.OrganismPart. == 'Bone Marrow',]
AML_metadata_PB <- AML_metadata[AML_metadata$Characteristics.OrganismPart. == 'Peripheral Blood',]

#Subset the metadata to disease state: Primary AML vs. Recurrent AML
AML_metadata_Primary <- AML_metadata[AML_metadata$Characteristics.DiseaseState. == 'Childhood Acute Myeloid Leukemia',]
AML_metadata_Recurrent <- AML_metadata[AML_metadata$Characteristics.DiseaseState. == 'Recurrent Childhood Acute Myeloid Leukemia',]

#Subset the RNAseq dataset according to the metadata subsets
RNAseq_df_BM <- RNAseq_df_RN[, rownames(AML_metadata_BM)]
RNAseq_df_PB <- RNAseq_df_RN[, rownames(AML_metadata_PB)]
RNAseq_df_Primary <- RNAseq_df_RN[, rownames(AML_metadata_Primary)]
RNAseq_df_Recurrent <- RNAseq_df_RN[, rownames(AML_metadata_Recurrent)]

#Create the deseq2 objects for each group: 
#using only BM, compare primary vs recurrent childhood AML 
ddseq_obj_BM <- DESeqDataSetFromMatrix(countData=RNAseq_df_BM, 
                                    colData=AML_metadata_BM, 
                                    design=~Characteristics.DiseaseState.)
#using only PB, compare primary vs recurrent childhood AML 
ddseq_obj_PB <- DESeqDataSetFromMatrix(countData=RNAseq_df_PB, 
                                    colData=AML_metadata_PB, 
                                    design=~Characteristics.DiseaseState.)
#using only Primary, compare BM vs PB
ddseq_obj_Primary <- DESeqDataSetFromMatrix(countData=RNAseq_df_Primary, 
                                    colData=AML_metadata_Primary, 
                                    design=~Characteristics.OrganismPart.)
#using only Recurrent, compare BM vs PB
ddseq_obj_Recurrent <- DESeqDataSetFromMatrix(countData=RNAseq_df_Recurrent, 
                                    colData=AML_metadata_Recurrent, 
                                    design=~Characteristics.OrganismPart.)


#The user can set their own rowsum and fdr threshold.
rowsum.threshold <- 1
fdr.threshold <- 0.1

#filters out genes that are not expressed at all across samples
ddseq_go_BM <- ddseq_obj_BM[rowSums(counts(ddseq_obj_BM)) > rowsum.threshold,]
ddseq_go_PB <- ddseq_obj_PB[rowSums(counts(ddseq_obj_PB)) > rowsum.threshold,]
ddseq_go_Primary <- ddseq_obj_Primary[rowSums(counts(ddseq_obj_Primary))> rowsum.threshold,]
ddseq_go_Recurrent <- ddseq_obj_Recurrent[rowSums(counts(ddseq_obj_Recurrent)) > rowsum.threshold,]

#DESeq fills in the DEseq object with information including 
#normalization, dispersion estimates, differential expression 
#results, etc. This code can take a little while to run
ddseq_go_BM <- DESeq(ddseq_go_BM)
ddseq_go_PB <- DESeq(ddseq_go_PB)
ddseq_go_Primary <- DESeq(ddseq_go_Primary)
ddseq_go_Recurrent <- DESeq(ddseq_go_Recurrent)

# Helper function! perform GO enrichment analysis, generate plots
plot_GO_enrichment <- function(ddseq_obj, contrast_values, fdr_threshold, title) {
  # Obtain DESeq results and filter by count threshold
  res <- results(ddseq_obj,
                 contrast = contrast_values,
                 independentFiltering = FALSE)
  
  # Extract differentially expressed genes based on FDR threshold
  assayed_genes <- rownames(res)
  de_genes <- rownames(res)[which(res$padj < fdr_threshold)]
  
  # Create a binary vector indicating differential expression
  gene_vector <- as.integer(assayed_genes %in% de_genes)
  names(gene_vector) <- assayed_genes
  
  # Create the null distribution to account for transcript length bias
  pwf <- nullp(gene_vector, "hg19", "ensGene")
  
  # Perform Gene Ontology enrichment analysis
  GO_BP <- goseq(pwf, 'hg19', 'ensGene', test.cats = c("GO:BP"))
  
  # Plot the top over-represented GO groups
  GO_BP %>%
    top_n(20, wt = -over_represented_pvalue) %>%
    mutate(hitsPerc = numDEInCat * 100 / numInCat) %>%
    ggplot(aes(x = hitsPerc,
               y = term,
               colour = over_represented_pvalue,
               size = numDEInCat)) +
    geom_point() +
    expand_limits(x = 0) +
    labs(x = "Hits (%)", y = "GO term", colour = "p value", size = "Count") +
    ggtitle(paste("GO terms:", title))
  
  return(GO_BP)
}

contrast_by_DiseaseState = c("Characteristics.DiseaseState.",
                               "Recurrent Childhood Acute Myeloid Leukemia", 
                               "Childhood Acute Myeloid Leukemia")
contrast_by_OrganismPart = c("Characteristics.OrganismPart.",
                              "Bone Marrow", 
                              "Peripheral Blood")

GO_BP_BM <- plot_GO_enrichment(ddseq_go_BM, contrast_by_DiseaseState, fdr.threshold, 
                               "BM - Primary vs Recurrent AML")
GO_BP_PB <- plot_GO_enrichment(ddseq_go_PB, contrast_by_DiseaseState, fdr.threshold, 
                               "PB - Primary vs Recurrent AML")
GO_BP_Primary <- plot_GO_enrichment(ddseq_go_Primary, contrast_by_OrganismPart, fdr.threshold, 
                                    "Primary Childhood AML - BM vs PB")
GO_BP_Recurrent <- plot_GO_enrichment(ddseq_go_Recurrent, contrast_by_OrganismPart, fdr.threshold,
                                      "Recurrent Childhood AML - BM vs PB")



#Helper function! RNAseq analysis:
perform_rnaseq_analysis <- function(ddseq_obj, contrast_values, title) {
  # Visualize the ddseq object in a dataframe
  res <- as_tibble(results(ddseq_obj,
                           contrast = contrast_values,
                           tidy = TRUE))
  
  # Add a column, sig, that is TRUE if padj is <0.05, false otherwise
  res <- res %>% mutate(sig = padj < 0.05)
  
  # Make a volcano plot to show differential gene expression for contrast
  volcano_plot <- res %>% ggplot(aes(log2FoldChange, -1 * log10(pvalue), col = sig)) + 
    geom_point() + 
    ggtitle("Volcano plot of gene expression for ", title) + 
    xlab("Log2 Fold Change") + 
    ylab("-1 * Log10(p-value)")
  
  # Create a PCA, coloring by disease state
  rld <- vst(ddseq_obj)
  pca_plot <- plotPCA(rld, intgroup = contrast_values[1]) +
                        labs(title = paste("PCA plot -", title))
  
  # Return the volcano plot and PCA plot
  return(list(volcano_plot = volcano_plot, pca_plot = pca_plot))
}

perform_rnaseq_analysis(ddseq_go_BM, contrast_by_DiseaseState, "Bone Marrow - Primary vs Recurrent Childhood AML")
perform_rnaseq_analysis(ddseq_go_PB, contrast_by_DiseaseState, "Peripheral Blood - Primary vs Recurrent Childhood AML")
perform_rnaseq_analysis(ddseq_go_Primary, contrast_by_OrganismPart, "Primary Childhood AML - BM vs PB")
perform_rnaseq_analysis(ddseq_go_Recurrent, contrast_by_OrganismPart, "Recurrent Childhood AML - BM vs PB")
