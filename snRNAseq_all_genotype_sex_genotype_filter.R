library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)
library(cowplot)
library(harmony)

PS19_1<- Read10X("/data/mouse-AD-project/PS19_1/outs/filtered_feature_bc_matrix")
PS19_2<- Read10X("/data/mouse-AD-project/PS19_2/outs/filtered_feature_bc_matrix")
PS19_3<- Read10X("/data/mouse-AD-project/PS19_3/outs/filtered_feature_bc_matrix")

WT_1<- Read10X("/data/mouse-AD-project/WT_1/outs/filtered_feature_bc_matrix")
WT_2<- Read10X("/data/mouse-AD-project/WT_2/outs/filtered_feature_bc_matrix")
WT_3<- Read10X("/data/mouse-AD-project/WT_3/outs/filtered_feature_bc_matrix")

CL_1<- Read10X("/data/mouse-AD-project/CL_1/outs/filtered_feature_bc_matrix")
CL_2<- Read10X("/data/mouse-AD-project/CL_2/outs/filtered_feature_bc_matrix")
CL_3<- Read10X("/data/mouse-AD-project/CL_3/outs/filtered_feature_bc_matrix")
CL_4<- Read10X("/data/mouse-AD-project/CL_4/outs/filtered_feature_bc_matrix")

CL_PS19_1<- Read10X("/data/mouse-AD-project/CL_PS19_1/outs/filtered_feature_bc_matrix")
CL_PS19_2<- Read10X("/data/mouse-AD-project/CL_PS19_2/outs/filtered_feature_bc_matrix")
CL_PS19_3<- Read10X("/data/mouse-AD-project/CL_PS19_3/outs/filtered_feature_bc_matrix")
CL_PS19_4<- Read10X("/data/mouse-AD-project/CL_PS19_4/outs/filtered_feature_bc_matrix")

PS19_1<- CreateSeuratObject(counts = PS19_1, project = "PS19_1_batch1_male", min.cells = 3, min.features = 200)
PS19_2<- CreateSeuratObject(counts = PS19_2, project = "PS19_2_batch2_male", min.cells = 3, min.features = 200)
PS19_3<- CreateSeuratObject(counts = PS19_3, project = "PS19_3_batch3_female", min.cells = 3, min.features = 200)

WT_1<- CreateSeuratObject(counts = WT_1, project = "WT_1_batch1_male", min.cells = 3, min.features = 200)
WT_2<- CreateSeuratObject(counts = WT_2, project = "WT_2_batch2_male", min.cells = 3, min.features = 200)
WT_3<- CreateSeuratObject(counts = WT_3, project = "WT_3_batch3_female", min.cells = 3, min.features = 200)

CL_1<- CreateSeuratObject(counts = CL_1, project = "CL_1_batch1_male", min.cells = 3, min.features = 200)
CL_2<- CreateSeuratObject(counts = CL_2, project = "CL_2_batch1_male", min.cells = 3, min.features = 200)
CL_3<- CreateSeuratObject(counts = CL_3, project = "CL_3_batch2_male", min.cells = 3, min.features = 200)
CL_4<- CreateSeuratObject(counts = CL_4, project = "CL_4_batch3_female", min.cells = 3, min.features = 200)

CL_PS19_1<- CreateSeuratObject(counts = CL_PS19_1, project = "CL_PS19_1_batch1_male", min.cells = 3, min.features = 200)
CL_PS19_2<- CreateSeuratObject(counts = CL_PS19_2, project = "CL_PS19_2_batch1_male", min.cells = 3, min.features = 200)
CL_PS19_3<- CreateSeuratObject(counts = CL_PS19_3, project = "CL_PS19_3_batch2_male", min.cells = 3, min.features = 200)
CL_PS19_4<- CreateSeuratObject(counts = CL_PS19_4, project = "CL_PS19_4_batch3_female", min.cells = 3, min.features = 200)


PS19_1[["percent.mt"]] <- PercentageFeatureSet(PS19_1, pattern = "^mt-")
PS19_2[["percent.mt"]] <- PercentageFeatureSet(PS19_2, pattern = "^mt-")
PS19_3[["percent.mt"]] <- PercentageFeatureSet(PS19_3, pattern = "^mt-")

WT_1[["percent.mt"]] <- PercentageFeatureSet(WT_1, pattern = "^mt-")
WT_2[["percent.mt"]] <- PercentageFeatureSet(WT_2, pattern = "^mt-")
WT_3[["percent.mt"]] <- PercentageFeatureSet(WT_3, pattern = "^mt-")

CL_1[["percent.mt"]] <- PercentageFeatureSet(CL_1, pattern = "^mt-")
CL_2[["percent.mt"]] <- PercentageFeatureSet(CL_2, pattern = "^mt-")
CL_3[["percent.mt"]] <- PercentageFeatureSet(CL_3, pattern = "^mt-")
CL_4[["percent.mt"]] <- PercentageFeatureSet(CL_4, pattern = "^mt-")

CL_PS19_1[["percent.mt"]] <- PercentageFeatureSet(CL_PS19_1, pattern = "^mt-")
CL_PS19_2[["percent.mt"]] <- PercentageFeatureSet(CL_PS19_2, pattern = "^mt-")
CL_PS19_3[["percent.mt"]] <- PercentageFeatureSet(CL_PS19_3, pattern = "^mt-")
CL_PS19_4[["percent.mt"]] <- PercentageFeatureSet(CL_PS19_4, pattern = "^mt-")


WT_1 <- subset(WT_1, subset = nFeature_RNA > 200 & percent.mt < 5)
WT_2 <- subset(WT_2, subset = nFeature_RNA > 200 & percent.mt < 5)
WT_3 <- subset(WT_3, subset = nFeature_RNA > 200 & percent.mt < 5)

PS19_1 <- subset(PS19_1, subset = nFeature_RNA > 200 & percent.mt < 5)
PS19_2 <- subset(PS19_2, subset = nFeature_RNA > 200 & percent.mt < 5)
PS19_3 <- subset(PS19_3, subset = nFeature_RNA > 200 & percent.mt < 5)

CL_1 <- subset(CL_1, subset = nFeature_RNA > 200 & percent.mt < 5)
CL_2 <- subset(CL_2, subset = nFeature_RNA > 200 & percent.mt < 5)
CL_3 <- subset(CL_3, subset = nFeature_RNA > 200 & percent.mt < 5)
CL_4 <- subset(CL_4, subset = nFeature_RNA > 200 & percent.mt < 5)

CL_PS19_1 <- subset(CL_PS19_1, subset = nFeature_RNA > 200 & percent.mt < 5)
CL_PS19_2 <- subset(CL_PS19_2, subset = nFeature_RNA > 200 & percent.mt < 5)
CL_PS19_3 <- subset(CL_PS19_3, subset = nFeature_RNA > 200 & percent.mt < 5)
CL_PS19_4 <- subset(CL_PS19_4, subset = nFeature_RNA > 200 & percent.mt < 5)

WT_1<- NormalizeData(WT_1)
WT_1<- FindVariableFeatures(WT_1, selection.method = "vst", nfeatures = 2000)
WT_1<- ScaleData(WT_1)
WT_1<- RunPCA(WT_1)
WT_1<- RunUMAP(WT_1, dims = 1:30)
nExp <- round(ncol(WT_1) * 0.05)
WT_1 <- doubletFinder_v3(WT_1, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
WT_1 <- subset(WT_1, cells=rownames(WT_1@meta.data)[which(WT_1@meta.data$DF.classification == "Singlet")])


WT_2<- NormalizeData(WT_2)
WT_2<- FindVariableFeatures(WT_2, selection.method = "vst", nfeatures = 2000)
WT_2<- ScaleData(WT_2)
WT_2<- RunPCA(WT_2)
WT_2<- RunUMAP(WT_2, dims = 1:30)
nExp <- round(ncol(WT_2) * 0.05)
WT_2 <- doubletFinder_v3(WT_2, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
WT_2 <- subset(WT_2, cells=rownames(WT_2@meta.data)[which(WT_2@meta.data$DF.classification == "Singlet")])


WT_3<- NormalizeData(WT_3)
WT_3<- FindVariableFeatures(WT_3, selection.method = "vst", nfeatures = 2000)
WT_3<- ScaleData(WT_3)
WT_3<- RunPCA(WT_3)
WT_3<- RunUMAP(WT_3, dims = 1:30)
nExp <- round(ncol(WT_3) * 0.05)
WT_3 <- doubletFinder_v3(WT_3, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
WT_3 <- subset(WT_3, cells=rownames(WT_3@meta.data)[which(WT_3@meta.data$DF.classification == "Singlet")])


PS19_1<- NormalizeData(PS19_1)
PS19_1<- FindVariableFeatures(PS19_1, selection.method = "vst", nfeatures = 2000)
PS19_1<- ScaleData(PS19_1)
PS19_1<- RunPCA(PS19_1)
PS19_1<- RunUMAP(PS19_1, dims = 1:30)
nExp <- round(ncol(PS19_1) * 0.05)
PS19_1 <- doubletFinder_v3(PS19_1, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
PS19_1 <- subset(PS19_1, cells=rownames(PS19_1@meta.data)[which(PS19_1@meta.data$DF.classification == "Singlet")])

PS19_2<- NormalizeData(PS19_2)
PS19_2<- FindVariableFeatures(PS19_2, selection.method = "vst", nfeatures = 2000)
PS19_2<- ScaleData(PS19_2)
PS19_2<- RunPCA(PS19_2)
PS19_2<- RunUMAP(PS19_2, dims = 1:30)
nExp <- round(ncol(PS19_2) * 0.05)
PS19_2 <- doubletFinder_v3(PS19_2, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
PS19_2 <- subset(PS19_2, cells=rownames(PS19_2@meta.data)[which(PS19_2@meta.data$DF.classification == "Singlet")])

PS19_3<- NormalizeData(PS19_3)
PS19_3<- FindVariableFeatures(PS19_3, selection.method = "vst", nfeatures = 2000)
PS19_3<- ScaleData(PS19_3)
PS19_3<- RunPCA(PS19_3)
PS19_3<- RunUMAP(PS19_3, dims = 1:30)
nExp <- round(ncol(PS19_3) * 0.05)
PS19_3<- doubletFinder_v3(PS19_3, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
PS19_3 <- subset(PS19_3, cells=rownames(PS19_3@meta.data)[which(PS19_3@meta.data$DF.classification == "Singlet")])


CL_1<- NormalizeData(CL_1)
CL_1<- FindVariableFeatures(CL_1, selection.method = "vst", nfeatures = 2000)
CL_1<- ScaleData(CL_1)
CL_1<- RunPCA(CL_1)
CL_1<- RunUMAP(CL_1, dims = 1:30)
nExp <- round(ncol(CL_1) * 0.05)
CL_1 <- doubletFinder_v3(CL_1, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
CL_1 <- subset(CL_1, cells=rownames(CL_1@meta.data)[which(CL_1@meta.data$DF.classification == "Singlet")])

CL_2<- NormalizeData(CL_2)
CL_2<- FindVariableFeatures(CL_2, selection.method = "vst", nfeatures = 2000)
CL_2<- ScaleData(CL_2)
CL_2<- RunPCA(CL_2)
CL_2<- RunUMAP(CL_2, dims = 1:30)
nExp <- round(ncol(CL_2) * 0.05)
CL_2 <- doubletFinder_v3(CL_2, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
CL_2 <- subset(CL_2, cells=rownames(CL_2@meta.data)[which(CL_2@meta.data$DF.classification == "Singlet")])

CL_3<- NormalizeData(CL_3)
CL_3<- FindVariableFeatures(CL_3, selection.method = "vst", nfeatures = 2000)
CL_3<- ScaleData(CL_3)
CL_3<- RunPCA(CL_3)
CL_3<- RunUMAP(CL_3, dims = 1:30)
nExp <- round(ncol(CL_3) * 0.05)
CL_3 <- doubletFinder_v3(CL_3, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
CL_3 <- subset(CL_3, cells=rownames(CL_3@meta.data)[which(CL_3@meta.data$DF.classification == "Singlet")])

CL_4<- NormalizeData(CL_4)
CL_4<- FindVariableFeatures(CL_4, selection.method = "vst", nfeatures = 2000)
CL_4<- ScaleData(CL_4)
CL_4<- RunPCA(CL_4)
CL_4<- RunUMAP(CL_4, dims = 1:30)
nExp <- round(ncol(CL_4) * 0.05)
CL_4 <- doubletFinder_v3(CL_4, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
CL_4 <- subset(CL_4, cells=rownames(CL_4@meta.data)[which(CL_4@meta.data$DF.classification == "Singlet")])


CL_PS19_1<- NormalizeData(CL_PS19_1)
CL_PS19_1<- FindVariableFeatures(CL_PS19_1, selection.method = "vst", nfeatures = 2000)
CL_PS19_1<- ScaleData(CL_PS19_1)
CL_PS19_1<- RunPCA(CL_PS19_1)
CL_PS19_1<- RunUMAP(CL_PS19_1, dims = 1:30)
nExp <- round(ncol(CL_PS19_1) * 0.05)
CL_PS19_1 <- doubletFinder_v3(CL_PS19_1, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
CL_PS19_1 <- subset(CL_PS19_1, cells=rownames(CL_PS19_1@meta.data)[which(CL_PS19_1@meta.data$DF.classification == "Singlet")])

CL_PS19_2<- NormalizeData(CL_PS19_2)
CL_PS19_2<- FindVariableFeatures(CL_PS19_2, selection.method = "vst", nfeatures = 2000)
CL_PS19_2<- ScaleData(CL_PS19_2)
CL_PS19_2<- RunPCA(CL_PS19_2)
CL_PS19_2<- RunUMAP(CL_PS19_2, dims = 1:30)
nExp <- round(ncol(CL_PS19_2) * 0.05)
CL_PS19_2 <- doubletFinder_v3(CL_PS19_2, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
CL_PS19_2 <- subset(CL_PS19_2, cells=rownames(CL_PS19_2@meta.data)[which(CL_PS19_2@meta.data$DF.classification == "Singlet")])

CL_PS19_3<- NormalizeData(CL_PS19_3)
CL_PS19_3<- FindVariableFeatures(CL_PS19_3, selection.method = "vst", nfeatures = 2000)
CL_PS19_3<- ScaleData(CL_PS19_3)
CL_PS19_3<- RunPCA(CL_PS19_3)
CL_PS19_3<- RunUMAP(CL_PS19_3, dims = 1:30)
nExp <- round(ncol(CL_PS19_3) * 0.05)
CL_PS19_3<- doubletFinder_v3(CL_PS19_3, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
CL_PS19_3 <- subset(CL_PS19_3, cells=rownames(CL_PS19_3@meta.data)[which(CL_PS19_3@meta.data$DF.classification == "Singlet")])

CL_PS19_4<- NormalizeData(CL_PS19_4)
CL_PS19_4<- FindVariableFeatures(CL_PS19_4, selection.method = "vst", nfeatures = 2000)
CL_PS19_4<- ScaleData(CL_PS19_4)
CL_PS19_4<- RunPCA(CL_PS19_4)
CL_PS19_4<- RunUMAP(CL_PS19_4, dims = 1:30)
nExp <- round(ncol(CL_PS19_4) * 0.05)
CL_PS19_4<- doubletFinder_v3(CL_PS19_4, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
CL_PS19_4 <- subset(CL_PS19_4, cells=rownames(CL_PS19_4@meta.data)[which(CL_PS19_4@meta.data$DF.classification == "Singlet")])


merged.alldata<- merge(CL_1, y = c(CL_2, CL_3, CL_4, CL_PS19_1, CL_PS19_2,CL_PS19_3, CL_PS19_4, PS19_1, PS19_2, PS19_3, WT_1, WT_2,WT_3), add.cell.ids = c("CL_1_batch1_male","CL_2_batch1_male","CL_3_batch2_male","CL_4_batch3_female","CL_PS19_1_batch1_male","CL_PS19_2_batch1_male","CL_PS19_3_batch2_male","CL_PS19_4_batch3_female","PS19_1_batch1_male","PS19_2_batch2_male","PS19_3_batch3_female","WT_1_batch1_male","WT_2_batch2_male","WT_3_batch3_female"), project = "all.genotype")


merged.alldata@meta.data$genotype <- c(rep("CL", ncol(CL_1)), rep("CL", ncol(CL_2)), rep("CL", ncol(CL_3)),rep("CL", ncol(CL_4)), rep("CL_PS19", ncol(CL_PS19_1)), rep("CL_PS19", ncol(CL_PS19_2)), rep("CL_PS19", ncol(CL_PS19_3)), rep("CL_PS19", ncol(CL_PS19_4)), rep("PS19", ncol(PS19_1)), rep("PS19", ncol(PS19_2)), rep("PS19", ncol(PS19_3)),rep("WT", ncol(WT_1)), rep("WT", ncol(WT_2)), rep("WT", ncol(WT_3)))


merged.alldata@meta.data$batch <- c(rep("batch1", ncol(CL_1)), rep("batch1", ncol(CL_2)), rep("batch2", ncol(CL_3)),rep("batch3", ncol(CL_4)), rep("batch1", ncol(CL_PS19_1)), rep("batch1", ncol(CL_PS19_2)), rep("batch2", ncol(CL_PS19_3)), rep("batch3", ncol(CL_PS19_4)), rep("batch1", ncol(PS19_1)), rep("batch2", ncol(PS19_2)), rep("batch3", ncol(PS19_3)),rep("batch1", ncol(WT_1)), rep("batch2", ncol(WT_2)), rep("batch3", ncol(WT_3)))


merged.alldata@meta.data$sex<- c(rep("male", ncol(CL_1)), rep("male", ncol(CL_2)), rep("male", ncol(CL_3)),rep("female", ncol(CL_4)), rep("male", ncol(CL_PS19_1)), rep("male", ncol(CL_PS19_2)), rep("male", ncol(CL_PS19_3)), rep("female", ncol(CL_PS19_4)), rep("male", ncol(PS19_1)), rep("male", ncol(PS19_2)), rep("female", ncol(PS19_3)),rep("male", ncol(WT_1)), rep("male", ncol(WT_2)), rep("female", ncol(WT_3)))


merged.alldata<-NormalizeData(merged.alldata,normalization.method="LogNormalize",scale.factor=10000)
merged.alldata<-FindVariableFeatures(merged.alldata,selection.method="vst",nfeatures=2000)
merged.alldata<-ScaleData(merged.alldata,verbose=FALSE)
merged.alldata<-RunPCA(merged.alldata,npcs=50,verbose=FALSE)
merged.alldata<-RunUMAP(merged.alldata, reduction="pca",dims=1:30)

merged.alldata.harmony<-RunHarmony(merged.alldata,group.by.vars=c("batch","sex"), plot_convergence = TRUE)


test <- merged.alldata.harmony %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
    FindClusters(resolution = 0.6) %>% 
identity()

test$genotype <- factor(x = test$genotype, levels = c("WT", "PS19","CL", "CL_PS19"))

p<-DimPlot(test, reduction = "umap", label=TRUE,repel = TRUE,pt.size=0.1,raster=TRUE)
ggsave(p, file="all_genotype_umap_cluster.pdf", scale=1,useDingbats=FALSE,height=8,width=12)

p1<-DimPlot(test, reduction = "umap", repel = TRUE,pt.size=0.1,raster=TRUE,group.by="batch")
p2<-DimPlot(test, reduction = "umap", repel = TRUE,pt.size=0.1,raster=TRUE,group.by="sex")
p3<-DimPlot(test, reduction = "umap", repel = TRUE,pt.size=0.1,raster=TRUE,group.by="genotype")
p4<-p1+p2+p3
ggsave(p4, file="all_genotype_umap_cluster_batch_sex_genotype.pdf", scale=1,useDingbats=FALSE,height=8,width=20)

p5<-FeaturePlot(test, features = c("Trem2","Hexb","Il3ra","Cx3cr1","Csf1r","C1qa","P2ry12"), max.cutoff = 3,cols = c("grey", "red"))
ggsave(p5, file="all_genotype_umap_microglial.pdf", scale=1,useDingbats=FALSE,height=20,width=20)

p6<-FeaturePlot(test, features = c("Grin1","Syt1","Rbfox3","Snap25","Grin2a","Slc17a7","Satb2"), max.cutoff = 3,cols = c("grey", "red"))
ggsave(p6, file="all_genotype_umap_Excitatory_neurons.pdf", scale=1,useDingbats=FALSE,height=20,width=20)

p7<-FeaturePlot(test, features = c("Gad1","Gad2","Sst","Npy","Pde10a","Tac1","Penk"), max.cutoff = 3,cols = c("grey", "red"))
ggsave(p7, file="all_genotype_umap_Interneurons.pdf", scale=1,useDingbats=FALSE,height=20,width=20)

p8<-FeaturePlot(test, features = c("Plp1","Mbp","Cldn11","Mog"), max.cutoff = 3,cols = c("grey", "red"))
ggsave(p8, file="all_genotype_umap_Oligodendrocytes.pdf", scale=1,useDingbats=FALSE,height=20,width=20)

p9<-FeaturePlot(test, features = c("Slc1a2","Gja1","Aqp4"), max.cutoff = 3,cols = c("grey", "red"))
ggsave(p9, file="all_genotype_umap_Astrocytes.pdf", scale=1,useDingbats=FALSE,height=20,width=20)

p10<-FeaturePlot(test, features = c("Pdgfra","Vcan","Cspg4","Olig1"), max.cutoff = 3,cols = c("grey", "red"))
ggsave(p10, file="all_genotype_umap_OPCs.pdf", scale=1,useDingbats=FALSE,height=20,width=20)

p11<-FeaturePlot(test, features = c("Vtn","Flt1","Cldn5"), max.cutoff = 3,cols = c("grey", "red"))
ggsave(p11, file="all_genotype_umap_vasulaor.pdf", scale=1,useDingbats=FALSE,height=20,width=20)

p12<-FeaturePlot(test, features = c("Prox1"), max.cutoff = 3,cols = c("grey", "red"))
ggsave(p12, file="all_genotype_umap_Granule.pdf", scale=1,useDingbats=FALSE,height=20,width=20)

saveRDS(test,"all_genotype_merged.alldata1.doublet.harmony.filter.data.rds")


