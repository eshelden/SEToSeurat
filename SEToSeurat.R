library(SummarizedExperiment)
library(tidyr)
library(biomaRt)
library(plyr)
library(Seurat)
#########################################################################################
#this function takes a StandardizedExperiment object and returns a SeuratObject.
#If running on a windows computer, it will also delete intermediate files
#if running on another system, change the lines containing 'del' to a system appropriate
#command
#########################################################################################
SEToSeurat<-function(se,sName, verbose = FALSE){
  sGSE<-as.data.frame(assay(se))
  mGeneIDs<-rownames(sGSE) #ensembl formatted IDs
  if (!grepl("ENSG",mGeneIDs[1]) && verbose){
    print("Warning: rownames of SummarizedExperiment do not have ensemble gene ID format.")
  }
  
  sGSE<-cbind(mGeneIDs,sGSE)
  colnames(sGSE)[1]<-"GeneID"
  longT<-gather(sGSE, key = "Sample", value = "mVal", -GeneID) #gather needs tidyr
  names<-as.factor(longT[,1])
  cells<-as.factor(longT[,2])
  sampleLevels<-levels(cells)
  nameLevels<-levels(names)
  
  if (verbose) print("Obtaining external gene names from biomaRt.")
  #we want the "real" names of the genes
  mart <- useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
  allIDs<-getBM(attributes = c("external_gene_name","entrezgene", "ensembl_gene_id"), values = mGeneIDs, mart = mart, uniqueRows = TRUE)
  my_genes_ann = allIDs[match(mGeneIDs, allIDs$ensembl_gene_id),]
  nameLevels<-cbind(nameLevels,as.character(my_genes_ann$external_gene_name))
  newData<-cbind(names,cells,longT[,3])
  colnames(newData)<-c("IDs","samples","values")
  
  if (verbose) print("Writing intermediate files")

  write.table(nameLevels, file = "genes.tsv", append = FALSE, sep = "\t", quote = FALSE, dec = ".",row.names = FALSE, col.names = FALSE)
  write.table(sampleLevels, file = "barcodes.tsv", append = FALSE, sep = "\t", quote = FALSE, dec = ".",row.names = FALSE, col.names = FALSE)
  write.table("%%MatrixMarket matrix coordinate integer general", file = "matrix.mtx", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table("%", file = "matrix.mtx", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  oString<-paste(length(nameLevels)/2,length(sampleLevels),dim(longT)[1])
  write.table(oString, file = "matrix.mtx", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(newData, file = "matrix.mtx", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  if (verbose) print("Reading 10X data")
  mySeuratData <- Read10X(data.dir = getwd())
  if (verbose) print("Creating Seurat object and cleaning up.")
  mSeurat<-CreateSeuratObject(raw.data = mSeuratData, min.cells = 3, min.genes = 200, project = sName)
  system("del genes.tsv")
  system("del barcodes.tsv")
  system("del matrix.mtx")
  return(mSeurat)
}