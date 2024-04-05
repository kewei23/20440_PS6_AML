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
example_childhood_AML <- read.table("TARGET-NCI-AML-RNAseq-data/TARGET-20-PABGKN-09A-01R.gene.quantification.txt",sep="\t",header=T)

#import the list of RNAseq files
RNAseq_files <- list.files(path="TARGET-NCI-AML-RNAseq-data", pattern="*.txt", full.names=TRUE, recursive=FALSE) 

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
AML_metadata <- dplyr::select(AML_metadata, -Source.Name, -Provider,-Material.Type, -Characteristics.Organism., -Performer, -Protocol.REF)
AML_metadata <- unique(AML_metadata) 

#Keep only the metadata that we have RNAseq data for
AML_metadata <- AML_metadata[(AML_metadata$Sample.Name %in% colnames(RNAseq_df)), ]

#We are only interested in bone marrow samples, so subset the metadata to only include BM samples.
AML_metadata <- AML_metadata[AML_metadata$Characteristics.OrganismPart. == 'Bone Marrow',]


#Creates new RNAseq dataframe with genes as rownames, and metadata dataframe with sample names as rownames
RNAseq_df_RN <- data.frame(RNAseq_df, row.names=1)
AML_metadata_RN <- data.frame(AML_metadata, row.names=4)

#For some reason, when converting gene names in RNAseq_RN to row names, the "-" in sample names is replaced with "."
#This creates errors when making the deseq2 object. So the following code replaces the "." with "-" in the column names
new_names = list()
for (name in colnames(RNAseq_df_RN)){
  name <- gsub("\\.","-", name)
  new_names <- append(new_names, name)
}
colnames(RNAseq_df_RN) <- c(new_names)

#We need to only keep the RNAseq data for BM samples.
RNAseq_df_RN<- RNAseq_df_RN[,(c(BM_metadata$Sample.Name))]

#Create the deseq2 object, design separating by childhood AML vs recurrent AML ()
ddseq_obj <- DESeqDataSetFromMatrix(countData=RNAseq_df_RN, colData=AML_metadata_RN, design=~Characteristics.DiseaseState.)


#The user can set their own rowsum and fdr threshold.
rowsum.threshold <- 1
fdr.threshold <- 0.1

rs <- rowSums(counts(ddseq_obj))
#filters out genes that are not expressed at all across samples
ddseq_go <- ddseq_obj[rs > rowsum.threshold,]
#calculate a new deseq_object
ddseq_go <- DESeq(ddseq_go)
#obtain deseq results and filter by count threshold
#The states in contrast can be flipped.
res_go <- results(ddseq_go, contrast=c("Characteristics.DiseaseState.", "Recurrent Childhood Acute Myeloid Leukemia","Childhood Acute Myeloid Leukemia"), independentFiltering=FALSE)
assayed.genes <- rownames(res_go)
de.genes <- rownames(res_go)[which(res_go$padj < fdr.threshold)]
#obtain a vector of differentially-expressed genes, with each element being
#either 1 or 0, with 1 as differentially-expressed, and 0 as not.
gene.vector=as.integer(assayed.genes %in%de.genes)
names(gene.vector)=assayed.genes

#Create the null distribution to account for transcript length bias affecting De
#supportedOrganisms()[supportedOrganisms()$Genome=="hg19",]
pwf = nullp(gene.vector, "hg19", "ensGene")

#Obtain the list of bio-processes that are over/under represented in childhood
#versus recurrent AML.
GO.BP=goseq(pwf,'hg19','ensGene',test.cats=c("GO:BP"))

#Plot the top over-represented GO groups.
GO.BP %>% 
  top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

#GO.BP %>% 
#  top_n(20, wt=-under_represented_pvalue) %>% 
#  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
#  ggplot(aes(x=hitsPerc, 
#             y=term, 
#             colour=over_represented_pvalue, 
#             size=numDEInCat)) +
#  geom_point() +
#  expand_limits(x=0) +
#  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

#Can also obtain all GO groups.
#GO.wall=goseq(pwf,"hg19","ensGene")
#head(GO.wall)

#DESeq fills in the DEseq object with information including 
#normalization, dispersion estimates, differential expression 
#results, etc. This code can take a little while to run
ddseq_obj <- DESeq(ddseq_obj)

#visualize the ddseq object in a dataframe
res <- as_tibble(results(ddseq_obj,contrast=c("Characteristics.DiseaseState.", "Childhood Acute Myeloid Leukemia","Recurrent Childhood Acute Myeloid Leukemia"),tidy=TRUE))
#Add a column, sig, that is TRUE if padj is <0.05, false otherwise.
res <- res %>% mutate(sig=padj<0.05)

#res %>% ggplot(aes(baseMean, log2FoldChange, col=sig)) + geom_point() + scale_x_log10() + ggtitle("MA plot")

#Make a volcano plot to show differential gene expression for childhood vs recurrent AML.
#Future edits to this figure (not yet implemented) will examine peripheral vs bone marrow samples and label genes in microenvironment sets.
res %>% ggplot(aes(log2FoldChange, -1*log10(pvalue), col=sig)) + geom_point() + ggtitle("Volcano plot of gene expression for childhood vs recurrent AML")

#Create a PCA, coloring by disease state.
rld <- vst(ddseq_obj)
plotPCA(rld, intgroup="Characteristics.DiseaseState.")



