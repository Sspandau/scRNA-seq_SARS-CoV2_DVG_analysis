library(ggplot2, warn.conflicts = FALSE) #for plotting
library(dplyr, warn.conflicts = FALSE) #for manipulating data
library(stringr, warn.conflicts = FALSE) #for manipulating character data
library(Seurat) #for seurat analysis
library(Matrix) #for creating matrix
library(SeuratObject) #for seurat object
library(matrixStats)
library(tidyverse)
library(sctransform)
library(uwot) ## for umap

#reading in matrix from cellranger if one lane
expression_matrix<- ReadMtx(
  mtx = "/gpfs/fs2/scratch/ysun81_lab/Simone_Terry/PHLE_sample/PHLE_cellranger/outs/filtered_feature_bc_matrix/matrix.mtx.gz",features = "/gpfs/fs2/scratch/ysun81_lab/Simone_Terry/PHLE_sample/PHLE_cellranger/outs/filtered_feature_bc_matrix/features.tsv.gz",
  cells = "/gpfs/fs2/scratch/ysun81_lab/Simone_Terry/PHLE_sample/PHLE_cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", )
Expression_table <- as.data.frame(as.matrix(expression_matrix))

#if multiple lanes
expression_matrix_l1<- ReadMtx(
  mtx = "/gpfs/fs2/scratch/ysun81_lab/Simone_Terry/Sample_231_i48/231_l1/filtered_feature_bc_matrix/matrix.mtx.gz",features = "/gpfs/fs2/scratch/ysun81_lab/Simone_Terry/Sample_231_i48/231_l1/filtered_feature_bc_matrix/features.tsv.gz",
  cells = "/gpfs/fs2/scratch/ysun81_lab/Simone_Terry/Sample_231_i48/231_l1/filtered_feature_bc_matrix/barcodes.tsv.gz", )
expression_matrix_l2<- ReadMtx(
  mtx = "/gpfs/fs2/scratch/ysun81_lab/Simone_Terry/Sample_231_i48/231_l2/filtered_feature_bc_matrix/matrix.mtx.gz",features = "/gpfs/fs2/scratch/ysun81_lab/Simone_Terry/Sample_231_i48/231_l2/filtered_feature_bc_matrix/features.tsv.gz",
  cells = "/gpfs/fs2/scratch/ysun81_lab/Simone_Terry/Sample_231_i48/231_l2/filtered_feature_bc_matrix/barcodes.tsv.gz", )
#convert expression matrix to data frame
Expression_table_l1 <- as.data.frame(as.matrix(expression_matrix_l1))
Expression_table_l2 <- as.data.frame(as.matrix(expression_matrix_l2))
Expression_table <- cbind(Expression_table_l1, Expression_table_l2)
Expression_table<-aggregate(. ~ barcode, data=Expression_table, FUN=sum) #fixes the dual lane issue
rm("expression_matrix_l1")
rm("expression_matrix_l2")
rm("Expression_table_l1")
rm("Expression_table_l2")

#adding DVG UMI
#adding dvg matrix from R filter after scRNAseq ViReMa
DVG_UMI<- read.csv("/gpfs/fs2/scratch/ysun81_lab/Simone_Terry/Sample_231_i48/231_dvgmatrix.csv")
DVG_Umi<-as.data.frame(t(DVG_UMI))
DVG_Umi$barcode<-NULL
colnames(DVG_Umi)=colnames(Expression_table)
Expression_DVG<-rbind(Expression_table, DVG_Umi)
rm("DVG_UMI")
rm("DVG_Umi")
rm("Expression_table")
Expression_DVG<-as.matrix(Expression_DVG)
library(data.table)
#naming DVG row as DVG, the row number below
#may be different depending on the number of features
#can check with this code : rownames(Expression_DVG)
rownames(Expression_DVG)[60667]<-"DVG"
rownames(Expression_DVG) 

#Cell Type
Expression_DVG<-as.data.frame(Expression_DVG)
# to check if data frame is in right format 
head(rownames(Expression_DVG))
head(colnames(Expression_DVG))

#load into seurat
seurat_Object<-CreateSeuratObject(counts = Expression_DVG, min.cells = 3, min.features = 200)

#filtering mitochondrial genes
seurat_Object[["percent.mt"]] <- PercentageFeatureSet(seurat_Object, pattern = "^MT-")
seurat_Object <- subset(seurat_Object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalize
#variable features
#scale data
seurat_Object <- SCTransform(seurat_Object, vars.to.regress = c("MT020881.1", "DVG"), verbose = FALSE)

#PCA
seurat_Object <- RunPCA(seurat_Object, features = VariableFeatures(object = seurat_Object))

#Cluster
seurat_Object <- FindNeighbors(seurat_Object, dims = 1:10)
seurat_Object<- FindClusters(seurat_Object, resolution = 0.2) #resolution varies by number of cells in sample

#gene markers for clusters
all_markers <-FindAllMarkers(seurat_Object, pval.type = "all", direction = "all")
top10_markers <- as.data.frame(all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
top10_markers

#renaming clusters based on cell types
new.cluster.ids.0.2res <- c("Basal", "Secretory", "Ciliated")
names(new.cluster.ids.0.2res) <- levels(seurat_Object)
seurat_Object<-RenameIdents(seurat_Object, new.cluster.ids.0.2res)

#adding cell type to data frame
Celltype<-seurat_Object@assays$RNA@counts
Celltype<-as.data.frame(Celltype)
cell<-data.frame(seurat_Object@active.ident)
cell<-t(cell)
colnames(cell)=colnames(Celltype)
Celltype<-rbind(Celltype, cell)

# first graphs
library(ggplot2)
library(scales)
library(ggnewscale)

#TSNE
seurat_Object<-RunTSNE(seurat_Object)
DimPlot(seurat_Object, reduction = "tsne") + theme(legend.title = element_text(size = 3), 
               legend.text = element_text(size = 3))
FeaturePlot(seurat_Object, reduction = "tsne", features = "MT020881.1")
FeaturePlot(seurat_Object, features = "DVG", reduction = "tsne")

Covid_boxplota<- ggplot(Celltype, aes(x = seurat_Object.active.ident, y = as.numeric(MT020881.1)))+ geom_jitter()

Covid_boxplotb<-boxplot(MT020881.1~seurat_Object.active.ident, data = Celltype_2dpi)

DVG_dotplot<-ggplot(Celltype, aes(x = seurat_Object.active.ident, y = as.numeric(DVG)))+ geom_jitter()

#Umap
seurat_Object_2dpi<-RunUMAP(seurat_Object, umap.method = uwot)
DimPlot(seurat_Object, reduction = "umap") + theme(legend.title = element_text(size = 3), 
               legend.text = element_text(size = 3))
Covid_Umap<-FeaturePlot(seurat_Object, reduction = "umap", feature = "MT020881.1")

DVG_Umap<-FeaturePlot(seurat_Object, reduction = "umap", feature = "DVG")