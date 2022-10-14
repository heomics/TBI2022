library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

#######################################standard data processing and integration#############################################
#load data
Pe1=Read10X(data.dir = "sample1_pe/filtered_feature_bc_matrix")
Co1=Read10X(data.dir = "sample1_co/filtered_feature_bc_matrix")

Pe2=Read10X(data.dir = "sample2_pe/filtered_feature_bc_matrix")
Co2=Read10X(data.dir = "sample2_co/filtered_feature_bc_matrix")

Pe3=Read10X(data.dir = "sample3_pe/filtered_feature_bc_matrix")
Co3=Read10X(data.dir = "sample3_co/filtered_feature_bc_matrix")

#make names specific
colnames(Pe1)=paste0("Pe1_", colnames(Pe1))
colnames(Pe2)=paste0("Pe2_", colnames(Pe2))
colnames(Pe3)=paste0("Pe3_", colnames(Pe3))
colnames(Co1)=paste0("Co1_", colnames(Co1))
colnames(Co2)=paste0("Co2_", colnames(Co2))
colnames(Co3)=paste0("Co3_", colnames(Co3))

#creat seurat object 
Pe1.seurat <- CreateSeuratObject(counts = Pe1, min.cells = 3, min.features = 200)
Pe2.seurat <- CreateSeuratObject(counts = Pe2, min.cells = 3, min.features = 200)
Pe3.seurat <- CreateSeuratObject(counts = Pe3, min.cells = 3, min.features = 200)
Co1.seurat <- CreateSeuratObject(counts = Co1, min.cells = 3, min.features = 200)
Co2.seurat <- CreateSeuratObject(counts = Co2, min.cells = 3, min.features = 200)
Co3.seurat <- CreateSeuratObject(counts = Co3, min.cells = 3, min.features = 200)

dat.seurat <- merge(x = Pe1.seurat, y = c(Pe2.seurat,  Pe3.seurat, Co1.seurat, Co2.seurat, Co3.seurat) )
table(sub("_.*", "", colnames(dat.seurat)))
##  Co1  Co2  Co3  Pe1  Pe2  Pe3 
## 6991 2211  643  893 1070  995

dat.seurat@meta.data[, "sample"]=paste0("P",sub("Pe|Co", "", dat.seurat@meta.data$orig.ident))

dat.seurat[["percent.mt"]] <- PercentageFeatureSet(dat.seurat, pattern = "^MT-")
dat.seurat <- subset(dat.seurat, subset =  percent.mt < 20)
dim(dat.seurat)
## [1] 21142 12438

table(sub("_.*", "", colnames(dat.seurat)))
##  Co1  Co2  Co3  Pe1  Pe2  Pe3 
## 6778 2190  637  836 1049  948

dat.seurat <- NormalizeData(dat.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
dat.seurat <- FindVariableFeatures(dat.seurat, selection.method = "vst", nfeatures = 2000)

#scale data
all.genes <- rownames(dat.seurat)
dat.seurat <- ScaleData(dat.seurat, features = all.genes)

#run PCA
dat.seurat <- RunPCA(dat.seurat, features = VariableFeatures(object = dat.seurat))

#determine dimension
dat.seurat <- JackStraw(dat.seurat, num.replicate = 100)
ElbowPlot(dat.seurat, ndims = 50)

#Clustering
dat.seurat <- FindNeighbors(dat.seurat, reduction = "pca", dims = 1:30)
dat.seurat <- FindClusters(dat.seurat, resolution = 0.2)

#UMAP after find clusters
dat.seurat <- RunUMAP(dat.seurat, reduction = "pca", dims = 1:30)

p1= DimPlot(dat.seurat, reduction = "umap", group.by = "orig.ident" )
p2 <- DimPlot(dat.seurat, reduction = "umap", label = TRUE)
plot_grid(p1,p2)

#integration 6samples by Seurat rPCA
# set up future for parallelization
library(future)
library(future.apply)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 20000 * 1024^2)

dat.seurat.list <- SplitObject(dat.seurat, split.by = "sample")
dat.seurat.list <- future_lapply(X = dat.seurat.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = dat.seurat.list)
dat.seurat.list <- future_lapply(X = dat.seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = dat.seurat.list, reduction = "rpca", 
                                  dims = 1:30)

dat.seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

dat.seurat.integrated <- ScaleData(dat.seurat.integrated, verbose = FALSE)
#run PCA
dat.seurat.integrated <- RunPCA(dat.seurat.integrated, verbose = FALSE)
# UMAP and Clustering
dat.seurat.integrated <- RunUMAP(dat.seurat.integrated, dims = 1:30)
dat.seurat.integrated <- FindNeighbors(dat.seurat.integrated, reduction = "pca", dims = 1:30)
dat.seurat.integrated <- FindClusters(dat.seurat.integrated, resolution = 0.2)

p1 = DimPlot(dat.seurat.integrated, reduction = "umap", group.by = "orig.ident" )
p2 <- DimPlot(dat.seurat.integrated, reduction = "umap", label = TRUE)
plot_grid(p1,p2)


##############################################Doublet filtering#################################################
library(DoubletFinder)

#work on individual dataset
dat.seurat.list <- SplitObject(dat.seurat.integrated, split.by = "orig.ident")

###
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 20000 * 1024^2)

dat.seurat.list <- future_lapply(X = dat.seurat.list, FUN = function(x) {
  DefaultAssay(x)="RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x)
  x <- RunPCA(x)
  x <- RunUMAP(x, dims = 1:10)
  
  ## pK Identification (no ground-truth) 
  sweep.res.list <- paramSweep_v3(x, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_value <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  
  ## Homotypic Doublet Proportion Estimate 
  annotations <- x@meta.data$integrated_snn_res.0.2
  homotypic.prop <- modelHomotypic(annotations)  
  
  ## Assuming 0.8%/2 doublet formation rate
  DoubletRate<-ncol(x)*4*1e-6
  
  nExp_poi <- round(DoubletRate*length(colnames(x)))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder 
  pN_value <- 0.25
  pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
  x <- doubletFinder_v3(x, PCs = 1:10, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
})

table(dat.seurat.list[[1]]@meta.data[,10])
## Doublet Singlet 
##       3     833
table(dat.seurat.list[[2]]@meta.data[,10])
## Doublet Singlet 
##     184    6594
table(dat.seurat.list[[3]]@meta.data[,10])
## Doublet Singlet 
##       4    1045
table(dat.seurat.list[[4]]@meta.data[,10])
## Doublet Singlet 
##      19    2171
table(dat.seurat.list[[5]]@meta.data[,10])
## Doublet Singlet 
##       4     944
table(dat.seurat.list[[6]]@meta.data[,10])
## Doublet Singlet 
##       2     635

dat.seurat.doublefinder.result=lapply(dat.seurat.list, function(x){
  data.frame(cell.id=rownames(x@meta.data), doublet.finder=x@meta.data[,10])
})

dat.seurat.doublefinder.result=do.call(rbind, dat.seurat.doublefinder.result)
dim(dat.seurat.doublefinder.result)
## [1] 12438     2

rownames(dat.seurat.doublefinder.result)=dat.seurat.doublefinder.result$cell.id

cell.id.single=with(dat.seurat.doublefinder.result, as.vector(cell.id)[doublet.finder=="Singlet"])

dat.seurat.integrated.f=subset(dat.seurat.integrated, cells = cell.id.single)
dim(dat.seurat.integrated.f)
## [1]  2000 12222

table(Idents(dat.seurat.integrated.f))
## 0    1    2    3    4    5    6    7 
## 3069 2783 1694 1240 1135 1034  787  480 



#sessionInfo
#R version 3.6.1 (2019-07-05)
#other attached packages: Seurat_3.1.1; DoubletFinder_2.0.3; ggplot2_3.3.3; dplyr_0.8.3; cowplot_1.0.0   









