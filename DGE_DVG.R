library(ggplot2, warn.conflicts = FALSE) #for plotting
library(dplyr, warn.conflicts = FALSE) #for manipulating data
library(stringr, warn.conflicts = FALSE) #for manipulating character data
library(openxlsx) #to send tables to excell
library(Seurat) #for seurat analysis
library(Matrix) #for creating matrix
library(SeuratObject) #for seurat object
library(matrixStats)
library(tidyverse)
library(tibble)
library(DESeq2)


#load into seurat object
expression_matrix<- ReadMtx(
  mtx = "/gpfs/fs2/scratch/sspandau/Yan_lab/Sample_2dpi/filtered_feature_bc_matrix/matrix1.mtx.gz",features = "/gpfs/fs2/scratch/sspandau/Yan_lab/Sample_2dpi/filtered_feature_bc_matrix/features1.tsv.gz",
  cells = "/gpfs/fs2/scratch/sspandau/Yan_lab/Sample_2dpi/filtered_feature_bc_matrix/barcodes1.tsv.gz", )

s_2dpi<- CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 200, max.features = 5000) 
s_2dpi[["percent.mt"]] <- PercentageFeatureSet(s_2dpi, pattern = "^MT-")
s_2dpi <- subset(s_2dpi, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#subset to only infected
Q1<-subset(x = s_2dpi, subset = MT020881.1 > 0) #creating object with at least one covid count
infectioncounts<- Q1@assays$RNA@counts
infectioncounts<-as.data.frame(infectioncounts)

#load in dvg status data frame
dvg_status<-read.csv("dvg_status_2dpi.csv")
colnames(dvg_status)

#align gene matrix barcodes with dvg barcodes
infectioncounts<-as.data.frame(t(infectioncounts))
infectioncounts<-tibble::rownames_to_column(infectioncounts, "row_names")
infectioncounts<-merge(infectioncounts, dvg_status, by = "row_names")
rownames(infectioncounts)<-infectioncounts$row_names

dvg_status<-merge(dvg_status, infectioncounts, by = "row_names")
dvg_status<-dvg_status[,c(1,2)]
colnames(dvg_status)<-c("row_names", "dvg_status") #fix colnames after adjusting number of cells to matrix cells

rownames(dvg_status)<-dvg_status$row_names #makes cell barcodes row names

infectioncounts<-infectioncounts[,-c(1, 29575)] #remove row_names from gene matrix
colnames(infectioncounts)
# reinput infectioncounts into seurat
infectioncounts<-t(infectioncounts)
s_2dpi_infected<- CreateSeuratObject(counts = infectioncounts)
#add meta data
s_2dpi_infected<-AddMetaData(s_2dpi_infected, dvg_status$dvg_status, col.name = 'dvg_status')
#set meta data to ident
s_2dpi_infected<-SetIdent(s_2dpi_infected, value = s_2dpi_infected@meta.data$dvg_status)

#standardization and normalization
s_2dpi_infected<-SCTransform(s_2dpi_infected, vars.to.regress = "MT020881.1", )

#Find marker genes
#DeSeq2
l_2dpi_deseq2<-FindMarkers(s_2dpi_infected, ident.1 = "Y", ident.2 = "N", test.use = "DeSeq2")
write.csv(l_2dpi_deseq2, "2dpi_DGE_list_DeSeq2.csv")
l_2dpi_deseq2<-as.data.frame(l_2dpi_deseq2)
l_2dpi_deseq2<-tibble::rownames_to_column(l_2dpi_deseq2, "row_names")

#MAST
l_2dpi_mast<-FindMarkers(s_2dpi_infected, ident.1 = "Y", ident.2 = "N", test.use = "MAST")
write.csv(l_2dpi_mast, "2dpi_DGE_list_MAST.csv")
l_2dpi_mast<-as.data.frame(l_2dpi_mast)
l_2dpi_mast<-tibble::rownames_to_column(l_2dpi_mast, "row_names")

# average expression for the two idents
s_2dpi_infected<-NormalizeData(object = s_2dpi_infected)
avg_E <-AverageExpression(s_2dpi_infected, features = l_2dpi_mast$row_names)
avg_E<-as.data.frame(avg_E)
avg_E<-tibble::rownames_to_column(avg_E, "row_names")
l_2dpi_mast_avgE<-merge(l_2dpi_mast, avg_E, by = "row_names")
write.csv(l_2dpi_mast_avgE, '2dpi_DGE_list_avgE_MAST.csv')