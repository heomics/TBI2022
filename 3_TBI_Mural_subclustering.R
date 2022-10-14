library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)


#############################################Mural cell subcluter###############################################
Mural.seurat <- subset(dat.seurat.integrated.f, idents = 7 )
dim(Mural.seurat)
## [1] 21142   480

#Integrated cells by Seurat SCT transform
Mural.seurat.list <- SplitObject(Mural.seurat, split.by = "sample")
Mural.seurat.list <- future_lapply(X = Mural.seurat.list[1:2], FUN = function(x) {
  x <- SCTransform(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = Mural.seurat.list, nfeatures = 3000)
Mural.seurat.list <- PrepSCTIntegration(object.list = Mural.seurat.list, anchor.features = features,verbose=FALSE)

k.filter <- min(200, min(sapply(Mural.seurat.list, ncol)))

anchors <- FindIntegrationAnchors(object.list = Mural.seurat.list, normalization.method = "SCT", 
                                  anchor.features = features,verbose=FALSE, k.filter = k.filter)


Mural.seurat.integrated<- IntegrateData(anchorset = anchors, normalization.method = "SCT",verbose=T)

Mural.seurat.integrated <- RunPCA(Mural.seurat.integrated,verbose = FALSE)
Mural.seurat.integrated <- RunUMAP(Mural.seurat.integrated, reduction = "pca", dims = 1:30)
Mural.seurat.integrated <- FindNeighbors(Mural.seurat.integrated, reduction = "pca", dims = 1:30)
#used 1
Mural.seurat.integrated<- FindClusters(Mural.seurat.integrated,resolution = 1)

p3 = DimPlot(Mural.seurat.integrated, reduction = "umap", group.by = "orig.ident" )
p4 <- DimPlot(Mural.seurat.integrated, reduction = "umap", label = TRUE)
plot_grid(p3, p4)

DefaultAssay(Mural.seurat.integrated)<-"RNA"
VlnPlot(Mural.seurat.integrated, features = c("nCount_RNA","nFeature_RNA"),  pt.size = 0.01, ncol = 2)








