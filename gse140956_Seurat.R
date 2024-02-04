library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)
library(clustree)
library(DEsingle)
setwd("~/gse140956")

a1 <- Read10X(data.dir = "~/gse140956/control")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "control.rds")

a1 <- Read10X(data.dir = "~/gse140956/empd")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "empd.rds")

control<-readRDS(file="control.rds")
empd<-readRDS(file="empd.rds")

control<-RenameCells(control,add.cell.id="control",for.merge=T)
control@meta.data$tech<-"control"
control@meta.data$celltype<-"control"

empd<-RenameCells(empd,add.cell.id="empd",for.merge=T)
empd@meta.data$tech<-"empd"
empd@meta.data$celltype<-"empd"

con_empd<-merge(control,empd)

saveRDS(con_empd, file="con_empd.rds")
hms<-con_empd

hms<-readRDS(file="con_empd.rds")
#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
#reference.list <- pancreas.list[c("Blood_P57","Tumor_P57","Blood_P58","Tumor_P58", 
#"Blood_P60","Tumor_P60", "Blood_P61", "Tumor_P61","Juxta_P60","Juxta_P61")]
reference.list <- pancreas.list[c("empd","control" )]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")
saveRDS(pancreas.integrated, file = "P57586061_after_integrated.rds")

hms_individual_integrated<-readRDS(file="hms_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1

#pbmc <- JackStraw(hms_individual_integrated, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
#ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test_1.2.rds")

hms_cluster<-readRDS(file="hms_cluster_test_1.2.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.05, 
                                logfc.threshold = 0.05
)
write.table(scRNA.markers,file="cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="top20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top20_marker_genes_1.2.csv",top20_table,row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

hms_cluster<-readRDS(file="hms_cluster_test_1.2.rds")
DimPlot(hms_cluster, reduction = "umap")

new.cluster.ids <- c( "Epithelial", "Epithelial","T", "Paget", "Epithelial","T","Progenitor","B", "Dendritic_cell","Epithelial","B",
                    "T", "T", "Epithelial","B", "MAST","Paget","B", "Melanocyte","Macrophage","Epithelial",
                    "Epithelial", "Epithelial","Melanocyte","Progenitor", "Epithelial","Fibroblast","Macrophage","Progenitor","B","Natural_killer") 
names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")

paget<-subset(hms_cluster_id, idents=c('Paget'))
DimPlot(paget, reduction = "umap")
saveRDS(paget, file="paget.rds")

#input each cluster
Paget<-readRDS("paget.rds")

#deg in Resident_Memory_CD8_T
paget<-readRDS("paget.rds")
a<-paget@meta.data
write.table(a,"a")
class(paget)
paget.sec<-as.SingleCellExperiment(paget)
group<-factor(c(rep(1,2023),rep(2,73)))

rds<-readRDS('paget.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"paget")
write.table(results.classified,"results.paget")

res <- read.csv("deg_paget", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-200,200),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

a<-read.table("biocarta_pathway", head=TRUE, sep="\t")
#fill=padj fill颜色填充，使用fdr
p <- ggplot(data=a,aes(x=Term,y=Count,fill=FDR))
#coord_flip()颠倒坐标轴
p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(color='gray'),
                 axis.text.y=element_text(color="black",face='bold',size=10),axis.text.x = element_text(color="black",face='bold',size=20))
#ylim(0,30) 更改横坐标的范围这里坐标轴颠倒了，虽然看起来是x轴，但其实是y轴
p3 <- p2 + ylim(0,55) + scale_fill_gradient(low="red",high="blue") 
p4 <- p3 + scale_x_discrete(limits=rev(a[,1])) +labs(x="",y="",face="bold")
p4

a<-read.table("kegg_pathway", head=TRUE, sep="\t")
#fill=padj fill颜色填充，使用fdr
p <- ggplot(data=a,aes(x=Term,y=Count,fill=FDR))
#coord_flip()颠倒坐标轴
p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(color='gray'),
                 axis.text.y=element_text(color="black",face='bold',size=10),axis.text.x = element_text(color="black",face='bold',size=20))
#ylim(0,30) 更改横坐标的范围这里坐标轴颠倒了，虽然看起来是x轴，但其实是y轴
p3 <- p2 + ylim(0,250) + scale_fill_gradient(low="red",high="blue") 
p4 <- p3 + scale_x_discrete(limits=rev(a[,1])) +labs(x="",y="",face="bold")
p4
