library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

##############################################Fibroblast subcluter#############################################
DefaultAssay(dat.seurat.integrated.f)<-"RNA"

Idents(dat.seurat.integrated.f)<-"seurat_clusters"
FB.seurat <- subset(dat.seurat.integrated.f, idents = 6 )
dim(FB.seurat)
## [1] 21142   787

FB.seurat <- NormalizeData(FB.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
FB.seurat <- FindVariableFeatures(FB.seurat, selection.method = "vst", nfeatures = 2000)

#scale data
FB.seurat <- ScaleData(FB.seurat)

##run PCA
FB.seurat <- RunPCA(FB.seurat, features = VariableFeatures(object = FB.seurat))

##determine dimension
FB.seurat <- JackStraw(FB.seurat, num.replicate = 100)
ElbowPlot(dat.seurat, ndims = 50)


FB.seurat <- RunUMAP(FB.seurat, dims = 1:30)
FB.seurat <- FindNeighbors(FB.seurat, dims = 1:30)
FB.seurat <- FindClusters(FB.seurat, resolution = 0.2)

DimPlot(FB.seurat, reduction = "umap", label = TRUE)

table(Idents(FB.seurat))
##   0   1   2 
## 401 206 180

VlnPlot(FB.seurat,features = c("nCount_RNA","nFeature_RNA"))
#Find fibroblast subcluster markers
FB.seurat.markers <- FindAllMarkers(FB.seurat, only.pos = TRUE)
FB.seurat.markers.f<-subset(FB.seurat.markers,FB.seurat.markers$avg_logFC>0.58&FB.seurat.markers$p_val_adj<0.05)

cluster.k.order<-0:2
##
g.select=lapply(cluster.k.order, function(i){
  t.k=FB.seurat.markers[FB.seurat.markers$cluster==i,]
  #t.k.genes=with(t.k, gene[avg_logFC >1 & p_val_adj <0.05 ])
  t.k.genes=with(t.k, gene[p_val_adj <0.05 ])
  t.k.genes=head(t.k.genes, 20)
})
g.select=unlist(g.select)

library(ggplot2)

DoHeatmap(FB.seurat, slot="data", features = g.select, label = F ) + scale_fill_gradientn(colors = c("blue", "white", "red"))









