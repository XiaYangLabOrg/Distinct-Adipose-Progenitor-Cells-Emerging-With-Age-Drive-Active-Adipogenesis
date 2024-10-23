library(Seurat) 
library(Matrix)
library(dplyr)
library(metap)
library(stringr)
library(ggplot2)

directories = list.dirs(path = "./Aging", full.names = TRUE, recursive = TRUE)
directories =directories[-1]

#iterate through the samples
for(sample in directories){
  #sample=directories[2]
  print(sample)
  splitSample = unlist(strsplit(sample,"/"))
  sampleName = splitSample[3]
  mid=unlist(strsplit(sampleName,"-"))
  tissue="gWAT"
  geneticbackground=mid[1]
  treatment=mid[2]
  group=str_sub(sampleName,1,-2)
  dropEST.data <- Read10X(data.dir = sample) 
  dropEST.data@Dimnames[[2]] = paste0(sampleName,"_",dropEST.data@Dimnames[[2]])
  dropEST.seurat <- CreateSeuratObject(counts = dropEST.data, project = sampleName)
  print(sampleName)
  rm(dropEST.data)
  dropEST.seurat@meta.data$data.geneticbackground = geneticbackground
  dropEST.seurat@meta.data$data.sampleName = sampleName
  dropEST.seurat@meta.data$data.tissue = tissue
  dropEST.seurat@meta.data$data.treatment = treatment
  dropEST.seurat@meta.data$data.group = group
  if(!exists("dropEST.combined")){
    dropEST.combined <- dropEST.seurat
    firstSampleName = sampleName
    firstSample = TRUE
  } else{
    if(firstSample==TRUE){
      dropEST.combined <- merge(x = dropEST.combined, y = dropEST.seurat, project = "Aging")
      firstSample = FALSE
    } else{
      dropEST.combined <- merge(x = dropEST.combined, y = dropEST.seurat, project = "Aging")
    }
  }
  rm(dropEST.seurat)
}

saveRDS(dropEST.combined, file = "./gWATraw.rds")


#####QC
mito.features <- grep(pattern = "^mt-", x = rownames(x = dropEST.combined), value = TRUE) 
percent.mito <- Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")[mito.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")) #mito.features/全体

# get ribosome %
ribo.features <- grep(pattern = "^Rps", x = rownames(x = dropEST.combined), value = TRUE)
ribo.features <- c(ribo.features, grep(pattern = "^Rpl", x = rownames(x = dropEST.combined), value = TRUE))
percent.ribo <- Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")[ribo.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts"))

# get predicted genes % 
pred.features <- grep(pattern = "^Gm1", x = rownames(x = dropEST.combined), value = TRUE)
pred.features <- c(pred.features,grep(pattern = "^Gm2", x = rownames(x = dropEST.combined), value = TRUE)) #c添加内容
pred.features <- c(pred.features,grep(pattern = "^Gm3", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm4", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm5", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm6", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm7", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm8", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm9", x = rownames(x = dropEST.combined), value = TRUE))
percent.pred <- Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")[pred.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts"))

# Add % mito, ribo and pred to the meta data
dropEST.combined[["percent.mito"]] <- percent.mito
dropEST.combined[["percent.ribo"]] <- percent.ribo
dropEST.combined[["percent.pred"]] <- percent.pred

table(dropEST.combined @meta.data$ data.sampleName)
ifelse(!dir.exists(file.path("./Plots/")), dir.create(file.path("./Plots/")), FALSE)
ifelse(!dir.exists(file.path("./Plots/","./QCPlots")), dir.create(file.path("./Plots/","./QCPlots")), FALSE)
ifelse(!dir.exists(file.path("./Plots/QCPlots","./ViolinPlots")), dir.create(file.path("./Plots/QCPlots","./ViolinPlots")), FALSE)
ifelse(!dir.exists(file.path("./Plots/QCPlots/ViolinPlots","./PreFilter")), dir.create(file.path("./Plots/QCPlots/ViolinPlots","./PreFilter")), FALSE)

source("./HelperFunctions.R")
cellQualityPlot(seuratObject=dropEST.combined,fileName="./Plots/QCPlots/ViolinPlots/PreFilter/SamplesQuality.pdf",H=9,W=40,
                featuresPlot=c("nFeature_RNA","nCount_RNA","percent.mito","percent.ribo","percent.pred"),identPlot = "orig.ident",
                pointSize=0.5)
#Generate QCPlots
dropEST.combined = SetIdent(dropEST.combined, value = "data.sampleName")
pdf(file="./Plots/QCPlots/ViolinPlots/PreFilter/iWAT_PreFilter.pdf",height=9,width=40)
VlnPlot(object = dropEST.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "percent.pred"),pt.size=0)
dev.off()


setThresholds(MinCells=3,MinGenes=200,MinUMIs=700,MinPercentMT=-Inf,MaxGenes=6000,
              MaxPercentMT=0.25,MaxUMIs=22000,MaxPercentRibo=1)

dropEST.combined.filtered = subset(x = dropEST.combined, subset = nFeature_RNA > mingenes & nFeature_RNA < maxgenes & percent.mito < maxPercentMT 
                                   & percent.ribo < maxPercentRibo & nCount_RNA < maxUMIs)

table(dropEST.combined.filtered@meta.data$ data.sampleName)


ifelse(!dir.exists(file.path("./Plots/QCPlots/ViolinPlots","./PostFilter")), dir.create(file.path("./Plots/QCPlots/ViolinPlots","./PostFilter")), FALSE)

#Generate QCPlots
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "data.sampleName")
pdf(file="./Plots/QCPlots/ViolinPlots/PostFilter/iWAT_PostFilterchoice.pdf",height=9,width=40)
VlnPlot(object = dropEST.combined.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "percent.pred"),pt.size = 0)
dev.off()


dropEST.combined.filtered <- NormalizeData(object = dropEST.combined.filtered, normalization.method = "LogNormalize",scale.factor = 10000)

dropEST.combined.filtered <- FindVariableFeatures(object = dropEST.combined.filtered)
length(dropEST.combined.filtered@assays$RNA@var.features) 
dropEST.combined.filtered <- ScaleData(dropEST.combined.filtered, verbose = FALSE)

#Perform PCA on the scaled data (uses the highly var genes)
dropEST.combined.filtered <- RunPCA(object = dropEST.combined.filtered, verbose = T, npcs = 75, ndims.print = 1:5, nfeatures.print = 10)

#Plot PCA Plots
ifelse(!dir.exists(file.path("./Plots/QCPlots","./PCA")), dir.create(file.path("./Plots/QCPlots","./PCA")), FALSE)
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "condition")
pdf(file="./Plots/QCPlots/PCA/PCA_HighlyVarGenes.pdf")
DimPlot(object = dropEST.combined.filtered, dims = c(1,2), reduction = "pca")
dev.off()

#Plot Heatmaps of PCs
ifelse(!dir.exists(file.path("./Plots/QCPlots","./PCHeatmap")), dir.create(file.path("./Plots/QCPlots","./PCHeatmap75")), FALSE)
pdf(file="./Plots/QCPlots/PCHeatmap75/PCA_Heatmaps.pdf",height=20,width=20)
DimHeatmap(object = dropEST.combined.filtered, dims = 1:75, balanced = TRUE, cells = 100, reduction = "pca")
dev.off()


#Jackstraw permutation to determine the number of "Significant" PCs to use for tSNE 2D projection
dropEST.combined.filtered <- JackStraw(object = dropEST.combined.filtered, reduction = "pca", num.replicate = 50, 
                                       verbose = TRUE, dims = 50)
dropEST.combined.filtered = ScoreJackStraw(object = dropEST.combined.filtered, dims = 1:50, reduction = "pca")

#Visualize the Jackstraw permutations
ifelse(!dir.exists(file.path("./Plots/QCPlots","./JackStraw")), dir.create(file.path("./Plots/QCPlots","./JackStraw")), FALSE)

pdf(file="./Plots/QCPlots/JackStraw/ASC_JackStraw.pdf",height=45,width=10) #50 Sig PCs
JackStrawPlot(object = dropEST.combined.filtered, dims = 1:50)
dev.off()

dropEST.combined.filtered <- FindNeighbors(object = dropEST.combined.filtered, reduction = "pca", dims = 1:50, k.param = 25)
dropEST.combined.filtered <- FindClusters(object = dropEST.combined.filtered, resolution = c(0.1,0.2,0.3,0.5,0.7,0.9,1.0,1.2,1.4),verbose = T, reduction = "pca")

dropEST.combined.filtered <- RunUMAP(object = dropEST.combined.filtered, reduction = "pca", dims = 1:50)
dropEST.combined.filtered <- RunTSNE(object = dropEST.combined.filtered, reduction = "pca", dims = 1:50,check_duplicates = FALSE)
resolutions <- c("0.9","0.5","0.7","1.2","1.4")
resolutions <- c("0.1","0.2","0.3")
for(reso in resolutions){
  pdf(file=paste0("/Volumes/data/VlnPlot/umap、",reso,".pdf"))
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = paste0("RNA_snn_res.",reso))
  print(DimPlot(dropEST.combined.filtered, label = T, reduction = "umap"))
  dev.off()
}

for(reso in resolutions){
  pdf(file=paste0("./Plots/resoTSNE/TSNEres",reso,".pdf"))
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = paste0("RNA_snn_res.",reso))
  print(DimPlot(dropEST.combined.filtered, label = T, reduction = "tsne"))
  dev.off()
}

#####DEG
dropEST.combined.filtered=readRDS(file = "APC.rds")
cellTypeDEGsStandard = list()
table(dropEST.combined.filtered @meta.data$ cell.type.new)

dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "cell.type.new")


for(cellType in levels(dropEST.combined.filtered@active.ident)){
  print(cellType)
  cellTypeSubset = subset(dropEST.combined.filtered,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = "data.condition")
  StandardDEGs = FindMarkers(cellTypeSubset, ident.1 = "Age", ident.2 = "Young")
  cellTypeDEGsStandard[[cellType]] = StandardDEGs
}  

saveRDS(cellTypeDEGsStandard, file = "./APCDEGs.rds")

for(cellType in names(cellTypeDEGsStandard)){
  print(cellType)
  cluster_sig = cellTypeDEGsStandard[[cellType]]
#Also add other thresholding
  cluster_sig <- cluster_sig[which(cluster_sig$ p_val_adj < 0.05),]
  write.table(cluster_sig,file=paste0("./DEG/",cellType,".txt"),sep='\t') 
}

########Pathway cite Jessica utils here
all_DEG=combined
dir.create(file.path("/Volumes/data/",CT))

colnames(all_DEG)[colnames(all_DEG) == "\"GENE\""] <- 'GENE'
colnames(all_DEG)[colnames(all_DEG) == "\"p_val\""] <- 'p_val'
colnames(all_DEG)[colnames(all_DEG) == "\"avg_log2FC\""] <- 'avg_logFC'
colnames(all_DEG)[colnames(all_DEG) == "\"p_val_adj\""] <- 'p_val_adj'


types=unique(all_DEG$CellType)
types
head(all_DEG)


for(celltype_cluster in types){
  print(celltype_cluster)
  
  DEG_df=all_DEG[all_DEG$CellType == celltype_cluster,]
  rownames(DEG_df) <-DEG_df$GENE
  DEG_df$`Cell_type` = DEG_df$CellType
  DEG_df=addInfoAndTrim(DEG_df)
  identifier=celltype_cluster
  add_direction<-TRUE
  total=pathway_enrichment(deg_list=DEG_df,FDR_threshold=NULL, pval_threshold=0.01, celltype_cluster,identifier, pVal = TRUE)
  pathways=total
  pathways = pathways[pathways$nOverlap>3,]
  Direction=c()
  meanFC=c()
  if(add_direction){
    new_overlap = c()
    for(g in pathways$Overlap){
      
      genes = unlist(strsplit(g, split = ","))
      
      UP_genes = c("UP:")
      
      upcount=0
      DOWN_genes = c("DOWN:")
      downcount=0
      sumFC=0
      nFC=0
      
      for(gene in genes){
        
        FC = DEG_df$avg_logFC[DEG_df$HUMAN==gene]
        FC = FC[!is.na(FC)]
        if( length(FC)>0){
          if(length(FC)>1){
            FC=FC[1]
            #dup_genes = append(dup_genes, gene)
            #cat("dup\n")
          }
          
          #if( length(FC)>0){}
          if(FC>0){
            upcount=upcount+1
            UP_genes = append(UP_genes, paste0(gene,","))
            sumFC=sumFC+FC
            nFC=nFC+1
          }
          else{
            downcount=downcount+1
            DOWN_genes = append(DOWN_genes, paste0(gene,","))
            sumFC=sumFC+FC
            nFC=nFC+1
          }
        }
      }
      UP_genes = append(UP_genes, DOWN_genes)
      meanFC=c(meanFC,sumFC/nFC)
      detailed_overlap = concatenate(UP_genes, mysep = "")
      new_overlap = append(new_overlap, detailed_overlap)
      if(upcount>downcount){
        Direction=c(Direction,"UP")
      }
      else{
        Direction=c(Direction,"Down")
      }
    }
    pathways$Overlap = new_overlap
    pathways$Direction = Direction
    pathways$meanFC = meanFC
  }
  
  
  write.csv(pathways,file=paste0("/Volumes/data/",CT,"/",celltype_cluster,"path.csv")) 
}


directories = list.files(path = paste0("/Volumes/data/",CT,"/"), full.names = TRUE, recursive = TRUE)
directories = directories[grep("csv",directories)]
selected=data.frame()
for(sample in directories){
  print(sample)
  pathwaylist=(read.csv(sample)) 
  selected = rbind(selected,pathwaylist)
  rm(pathwaylist)
}
head(selected)
selected$CellType=CT
selected=selected[selected$FDR<0.05,]
write.csv(selected,file=paste0("/Volumes/data/",CT,"path.csv")) 

library(dplyr)
list_of_pathway_databases = list.files(path = "/Volumes/data//Resources/", pattern = "*.txt", full.names = TRUE)
pathB = read.table("/Volumes/data/Resources/Biocarta.txt", header=T, sep='\t', check.names=F, quote=NULL)
pathH= read.table("/Volumes/data/Resources/Hallmark.txt", header=T, sep='\t', check.names=F, quote=NULL)
pathK = read.table("/Volumes/data/Resources/Kegg.txt", header=T, sep='\t', check.names=F, quote=NULL)
pathR= read.table("/Volumes/data/Resources/Reactome.txt", header=T, sep='\t', check.names=F, quote=NULL)
path=rbind(pathB,pathH,pathK,pathR)

df <- path %>%
  group_by(module) %>%
  mutate(gene_count = n())

df_summary <- df %>%
  group_by(module) %>%
  summarize(gene_count = first(gene_count))


CB <- selected %>%
  left_join(df_summary, by = c("Pathway" = "module")) %>%
  rename(pathway_size = gene_count)



library(Seurat) 
library(Matrix)
library(dplyr)
library(metap)
library(stringr)
library(ggplot2)
setwd("/Volumes/data/SVF/")


#######Cell Type dim plot
p=DimPlot(dropEST.combined.filtered,label = T,pt.size = 0.5, reduction = "umap", cols = ASCcolors)+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("APC_celltypes_withname",".jpeg"), plot=p, device="jpeg",
       path="/", height=4,width=5.5, units="in", dpi=300)


ASCcolors=c("#FFE4B5","#90EE90","#FFA500","#9ACD32","#008B00")
p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.5, reduction = "umap", split.by = "data.condition", cols = ASCcolors,repel = T)+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("APC_celltypes_groupsplit_noname",".jpeg"), plot=p, device="jpeg",
       path="/", height=4,width=8.5, units="in", dpi=500)

dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "data.condition")
Age=subset(dropEST.combined.filtered,idents="Age")
Age = SetIdent(Age, value = "cell.type.new")
my_levels <- c("ASC","IAP","CP-1","CP-2","CP-A")
Idents(Age) <- factor(Idents(Age), levels= my_levels)
p=DimPlot(Age,label = F,pt.size = 0.5, reduction = "tsne", cols = ASCcolors)+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text = element_blank(),axis.ticks = element_blank())
ggsave(filename=paste0("AgeAPC_celltypes",".jpeg"), plot=p, device="jpeg",
       path="/", height=4,width=5.5, units="in", dpi=300)

Young=subset(dropEST.combined.filtered,idents="Young")
Young = SetIdent(Young, value = "cell.type.new")
my_levels <- c("ASC","IAP","CP-1","CP-2","CP-A")
Idents(Young) <- factor(Idents(Young), levels= my_levels)
p=DimPlot(Young,label = F,pt.size = 0.5, reduction = "tsne", cols = ASCcolors)+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text = element_blank(),axis.ticks = element_blank())
ggsave(filename=paste0("YoungAPC_celltypes",".jpeg"), plot=p, device="jpeg",
       path="/", height=4,width=5.5, units="in", dpi=300)



dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.condition")
p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.5, reduction = "tsne", cols = c("#4472C4","#ED7D31"))+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("APC_AgeYoubgTsne",".jpeg"), plot=p, device="jpeg",
       path="/", height=4,width=5.7, units="in", dpi=300)

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.sampleName")
p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.5, reduction = "tsne", cols = c("red","pink","orange","blue","cyan","lightblue"))+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("APC_sampleTsne",".jpeg"), plot=p, device="jpeg",
       path="/", height=4,width=5.7, units="in", dpi=300)

dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "cell.type.new")
my_levels <- c("CP-A","CP-2","CP-1","IAP","ASC")
Idents(dropEST.combined.filtered) <- factor(Idents(dropEST.combined.filtered), levels= my_levels)
p=DotPlot(dropEST.combined.filtered,features = c("Dpp4","Pi16","Cd55","Hspa1b","Scg3","Sned1","C7","Igf1","Apoe","Cilp","Mgp","Mfap4","Thbs1","Ccl11","Spry1"),cols = c("lightgrey","red"))+ RotatedAxis()+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("APC_markerdotplot",".jpeg"), plot=p, device="jpeg",
       path="/", height=4,width=9.5, units="in", dpi=300)


for(i in 1:length(cell_type_genes)){
  p = VlnPlot(dropEST.combined.filtered,features = cell_type_genes[i],pt.size=0,#cols = ASCcolors
  )+
    theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
          axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
  ggsave(filename=paste0("nonstem_VLNPlots_",cell_type_genes[i],".jpeg"), plot=p, device="jpeg",
         path="/", height=5,width=8, units="in", dpi=300)
}


p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.5, reduction = "tsne", split.by = "data.sampleName",repel = T)+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("all_samplesplit_noname",".jpeg"), plot=p, device="jpeg",
       path="/", height=4,width=20, units="in", dpi=500)


for(i in 1:length(cell_type_genes)){
  #jpeg(file=paste0("/Users/gaoyanli/Desktop/font/Wider2_Featureplot_",cell_type_genes,".jpeg"),width = 500, height = 400, res = 150)
  p=FeaturePlot(object =dropEST.combined.filtered, features = cell_type_genes[i], reduction = "tsne", min.cutoff = 0,max.cutoff = 1,pt.size = 0.001) +
    theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
          axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
  ggsave(filename=paste0("human_TSNEPlots_",cell_type_genes[i],".jpeg"), plot=p, device="jpeg",
         path="/", height=4, width=5, units="in", dpi=300)
}

p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.05, reduction = "tsne")+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("human_celltypes_number",".jpeg"), plot=p, device="jpeg",
       path="/", height=4,width=5, units="in", dpi=300)

for(i in 1:length(cell_type_genes)){
    p = VlnPlot(dropEST.combined.filtered,features = cell_type_genes[i],pt.size=0,#cols = ASCcolors
  )+
    theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
          axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
  ggsave(filename=paste0("human_VLNPlots_",cell_type_genes[i],".jpeg"), plot=p, device="jpeg",
         path="/", height=4,width=5, units="in", dpi=300)
}




##########trajectory
library(Seurat) 
library(SingleCellExperiment)
library(slingshot)

dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "data.condition")
dropEST.combined.filtered=subset(dropEST.combined.filtered,idents = c("Young"))
dropEST.combined.filtered <- NormalizeData(object = dropEST.combined.filtered, normalization.method = "LogNormalize",scale.factor = 10000)
dropEST.combined.filtered <- FindVariableFeatures(object = dropEST.combined.filtered)
length(dropEST.combined.filtered@assays$RNA@var.features) 
dropEST.combined.filtered <- RunPCA(object = dropEST.combined.filtered, verbose = T, npcs = 75, ndims.print = 1:5, nfeatures.print = 10)
dropEST.combined.filtered <- FindNeighbors(object = dropEST.combined.filtered, reduction = "pca", dims = 1:50, k.param = 25)

dropEST.combined.filtered <- FindClusters(object = dropEST.combined.filtered, resolution = c(0.7),verbose = T, reduction = "pca")
table(dropEST.combined.filtered @meta.data$ cell.type.new)
dropEST.combined.filtered <- RunTSNE(object = dropEST.combined.filtered, reduction = "pca", dims = 1:50,check_duplicates = FALSE)

dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "cell.type.new")
my_levels <- c("ASC","IAP","CP-1","CP-2","CP-A")
Idents(dropEST.combined.filtered) <- factor(Idents(dropEST.combined.filtered), levels= my_levels)

sce <- as.SingleCellExperiment(dropEST.combined.filtered)
table(sce$ident)
sce <- slingshot(sce, clusterLabels = dropEST.combined.filtered$cell.type.new, reducedDim = "PCA",
                 start.clus="ASC", stretch = 0)
saveRDS(sce, file = "./Plots/CellTypesUMAP/newallslingshot.rds")


library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

str(reducedDims(sce))

plot(reducedDims(sce)$PCA, col = ASCcolors, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

ASCcolors=c("#FFE4B5","#90EE90","#FFA500","#9ACD32","#008B00")

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(unique(sce$ident), brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(sce$seurat_clusters, hue_pal())

library(tradeSeq)

sim=readRDS(file = "./Plots/CellTypesUMAP/newallslingshot.rds")


sim <- fitGAM(sim)

ATres <- associationTest(sim)
topgenes <- c("Dpp4","Pi16","Cd55","Hspa1b","Scg3","Sned1","Thbs1","Ccl11","Spry1","C7","Igf1","Apoe","Cilp","Mgp","Mfap4")
topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sim$slingPseudotime_3, na.last = NA)
heatdata <- assays(sim)$counts[topgenes, pst.ord]
heatclus <- 5
ColSideColors = brewer.pal(9,"Set1")
heatclus <- sim$GMM[pst.ord]
pheatmap(log1p(heatdata))
heatmap(heatdata, Colv = NA,ColSideColors = brewer.pal(9,"Set1"))

library(pheatmap)
pdf("./Plots/CellTypesUMAP/newheatmaptrypseudo3.pdf",width = 5, height = 6)
pheatmap(log1p(heatdata))
dev.off()

cell_type_genes222 <- c("Dpp4","Pi16","Cd55","Hspa1b","Scg3","Sned1")

VlnPlot(object = dropEST.combined.filtered, features.plot = features.plot, x.lab.rot = TRUE)
RidgePlot(object = dropEST.combined.filtered, features.plot = features.plot, nCol = 2)

pdf(file="./Plots/conditiontraject/GeneralConditionTypesTSNE2.pdf",height = 4, width = 5.5)
DimPlot(dropEST.combined.filtered,label = T,pt.size = 0.5, reduction = "tsne")
dev.off()


ASCcolors=c("#FFE4B5","#90EE90","#FFA500","#9ACD32","#008B00")

pdf(file="./Plots/conditiontraject/newGeneralCellTypetsnesplit.pdf",height = 4, width = 9)
DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.5, reduction = "tsne", split.by = "data.condition",
        repel = T)
dev.off()

salls <- slingshot(Embeddings(dropEST.combined.filtered, "tsne"), 
                   clusterLabels = dropEST.combined.filtered$cell.type.condition, 
                   start.clus = "Young-ASC",stretch = 0)

pdf("./Plots/conditiontraject/tsneline.pdf")
plot(reducedDim(salls),col = c("pink"), pch = 16, cex = 0)
lines(salls, lwd = 2, type = 'lineages', col = 'black')
#lines(salls, lwd = 2, col = 'black')
dev.off()


##########compare with established cell types
setwd("/Volumes/data/newSVF/")
genes <- read.table("/Volumes/data/newSVF/markerstocompare.txt", header=T, sep='\t', check.names=F, quote=NULL,stringsAsFactors=FALSE)


dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$FAP1[1:50]), ctrl = 50, name = 'FAP1markers')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$FAP2[1:50]), ctrl = 50, name = 'FAP2markers')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$FAP3[1:50]), ctrl = 50, name = 'FAP3markers')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$FAP4[1:50]), ctrl = 50, name = 'FAP4markers')

dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$ASC1[1:50]), ctrl = 50, name = 'ASC1markers')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$ASC2[1:50]), ctrl = 50, name = 'ASC2markers')

dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$P1[1:50]), ctrl = 50, name = 'P1markers')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$P2[1:50]), ctrl = 50, name = 'P2markers')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$P3[1:50]), ctrl = 50, name = 'P3markers')

dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$APC[1:50]), ctrl = 50, name = 'APCmarkers')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$CP[1:50]), ctrl = 50, name = 'CPmarkers')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$FIP[1:50]), ctrl = 50, name = 'FIPmarkers')

dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$FIP[1:50]), ctrl = 50, name = 'FIPmarkers')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(genes$FIP[1:50]), ctrl = 50, name = 'FIPmarkers')

dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(cell_type_genes1), ctrl = 50, name = 'ARC')

names(x = dropEST.combined.filtered[[]])

cell_type_genes=names(x = dropEST.combined.filtered[[]])[XX:XX]
for(i in 1:length(cell_type_genes)){
  p=FeaturePlot(object = dropEST.combined.filtered, features = cell_type_genes[i],reduction = "tsne",
  min.cutoff=0, max.cutoff=1,pt.size = 0.01,cols =c("lightgrey", "red"))+
    theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
          axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
  ggsave(filename=paste0("eWAcompare",cell_type_genes[i],".jpeg"), plot=p, device="jpeg",
         path="/", height=4,width=5, units="in", dpi=300)
}

dropEST.combined.filtered[['cellmarkermodule']]<- CreateAssayObject(data = t(x = FetchData(object = dropEST.combined.filtered))) 
my_levels <- c("ASC","IAP","CP-1","CP-2","CP-A")
Idents(dropEST.combined.filtered) <- factor(Idents(dropEST.combined.filtered), levels= my_levels)
dropEST.combined.filtered[['FAPmodule']]<- CreateAssayObject(data = t(x = FetchData(object = dropEST.combined.filtered, 
                                                                                    vars = c('FAP1markers1','FAP2markers1','FAP3markers1','FAP4markers1'))))

p=DoHeatmap(object = dropEST.combined.filtered, features = c('FAP1markers1','FAP2markers1','FAP3markers1','FAP4markers1'), assay = 'FAPmodule', slot = 'data')+ scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-10,10),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("Heatmapcompare_Savari",".jpeg"), plot=p, device="jpeg",
       path="/", height=4.5,width=7, units="in", dpi=300)


dropEST.combined.filtered[['Burlmodule']]<- CreateAssayObject(data = t(x = FetchData(object = dropEST.combined.filtered, 
                                                                                     vars = c('ASC1markers1','ASC2markers1'))))
p=DoHeatmap(object = dropEST.combined.filtered, features = c('ASC1markers1','ASC2markers1'), assay = 'Burlmodule', slot = 'data')+ scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-10,10),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("Heatmapcompare_Burl",".jpeg"), plot=p, device="jpeg",
       path="/", height=4.5,width=7, units="in", dpi=300)

dropEST.combined.filtered[['Helpermodule']]<- CreateAssayObject(data = t(x = FetchData(object = dropEST.combined.filtered, 
                                                                                       vars = c('APCmarkers1','CPmarkers1','FIPmarkers1'))))
p=DoHeatmap(object = dropEST.combined.filtered, features = c('APCmarkers1','CPmarkers1','FIPmarkers1'), assay = 'Helpermodule', slot = 'data')+ scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-10,10),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("Heatmapcompare_Helper",".jpeg"), plot=p, device="jpeg",
       path="/", height=4.5,width=7, units="in", dpi=300)

dropEST.combined.filtered[['Schwaliemodule']]<- CreateAssayObject(data = t(x = FetchData(object = dropEST.combined.filtered, 
                                                                                         vars = c('P1markers1','P2markers1','P3markers1'))))
p=DoHeatmap(object = dropEST.combined.filtered, features = c('P1markers1','P2markers1','P3markers1'), assay = 'Schwaliemodule', slot = 'data')+ scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-10,10),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("eatmapcompare_Schwalie",".jpeg"), plot=p, device="jpeg",
       path="/", height=4.5,width=7, units="in", dpi=300)


dropEST.combined.filtered <- readRDS("intermediate_data/SeuratObj/mASPC_SeuratObj.rds")
Emont_mASPC_markers <- read_excel("resource/NewMarkersToCompare.xlsx", 
                                  sheet = "Emont_Nat2022_m_ASPC")
unique(Emont_mASPC_markers$Cell_type)
filtered_df <- Emont_mASPC_markers %>% 
  filter(p_val_adj < 0.05)
result_df <- filtered_df %>% 
  group_by(Cell_type) %>% 
  arrange(desc(avg_log2FC)) %>% 
  ungroup()
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$Sig_markers[result_df$Cell_type == "mASPC1"][1:20]), ctrl = 50, name = 'mASPC1')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$Sig_markers[result_df$Cell_type == "mASPC2"][1:20]), ctrl = 50, name = 'mASPC2')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$Sig_markers[result_df$Cell_type == "mASPC3"][1:20]), ctrl = 50, name = 'mASPC3')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$Sig_markers[result_df$Cell_type == "mASPC4"][1:20]), ctrl = 50, name = 'mASPC4')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$Sig_markers[result_df$Cell_type == "mASPC5"][1:20]), ctrl = 50, name = 'mASPC5')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$Sig_markers[result_df$Cell_type == "mASPC6"][1:20]), ctrl = 50, name = 'mASPC6')
dropEST.combined.filtered$mASPC1 <- dropEST.combined.filtered$mASPC11
dropEST.combined.filtered$mASPC2 <- dropEST.combined.filtered$mASPC21
dropEST.combined.filtered$mASPC3 <- dropEST.combined.filtered$mASPC31
dropEST.combined.filtered$mASPC4 <- dropEST.combined.filtered$mASPC41
dropEST.combined.filtered$mASPC5 <- dropEST.combined.filtered$mASPC51
dropEST.combined.filtered$mASPC6 <- dropEST.combined.filtered$mASPC61
dropEST.combined.filtered$mASPC11 <- NULL
dropEST.combined.filtered$mASPC21 <- NULL
dropEST.combined.filtered$mASPC31 <- NULL
dropEST.combined.filtered$mASPC41 <- NULL
dropEST.combined.filtered$mASPC51 <- NULL
dropEST.combined.filtered$mASPC61 <- NULL
dropEST.combined.filtered[['Emont.mModule']]<- CreateAssayObject(data = t(x = FetchData(object = dropEST.combined.filtered, 
                                                                             vars = c('mASPC1','mASPC2','mASPC3','mASPC4','mASPC5','mASPC6'))))
p = DoHeatmap(object = dropEST.combined.filtered, 
              features = c('mASPC1','mASPC2','mASPC3','mASPC4','mASPC5','mASPC6'), 
              assay = 'Emont.mModule',
              slot = 'data'
) + 
  scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +
  theme(text=element_text(size=22, family="Arial"),
        plot.title=element_text(size=24, face="italic"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-10,10),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18)) 
ggsave(filename = "Emont_mouse_top20markers_50ctrPerMarker.jpeg", plot = p, device = "jpeg", path = "output_plot/NewMarkerComparison/Emont_mASPC/top20markers_50ctrPerMarker/", height=4.5,width=7, units="in", dpi=300)
DefaultAssay(dropEST.combined.filtered) <- 'RNA'
names(x = dropEST.combined.filtered[[]])
for (subtype in c('mASPC1','mASPC2','mASPC3','mASPC4','mASPC5','mASPC6')){
  p <- FeaturePlot(dropEST.combined.filtered, 
                   features = subtype, 
                   reduction = "tsne",
                   min.cutoff = 0,
                   max.cutoff = 1,
                   pt.size = 0.01,
                   cols =c("lightgrey", "red"))  +
    scale_colour_gradient(limits = c(0.00,1.00), breaks = c(0.0, 0.25, 0.50, 0.75, 1.00), low = "lightgrey", high = "red") +
    theme(text=element_text(size=22,family="Arial"), 
          plot.title=element_text(size=24, face = "plain"),
          axis.text.x=element_text(size=18),
          axis.text.y=element_text(size=18)) + 
    labs(title = paste0(subtype," markers"))
  ggsave(filename=paste0("Emont_",subtype,"_min0max1_featureplot.jpeg"), 
         plot=p, 
         device="jpeg",
         path="output_plot/NewMarkerComparison/Emont_mASPC/top20markers_50ctrPerMarker/", height=4,width=5, units="in", dpi=300)
}

dropEST.combined.filtered <- readRDS("intermediate_data/SeuratObj/hASPC_SeuratObj_Gaoyan.rds")
Emont_hASPC_markers <- read_excel("resource/Emont_Nat2022_HumanSVFMarkers.xlsx", 
                                  sheet = "ASPC markers")
filtered_df <- Emont_hASPC_markers %>% 
  filter(p_val_adj < 0.05)
result_df <- filtered_df %>% 
  group_by(cluster) %>% 
  arrange(desc(avg_log2FC)) %>% 
  ungroup()
DefaultAssay(dropEST.combined.filtered) <- "RNA"
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$cluster == "hASPC1"][1:20]), ctrl = 50, name = 'hASPC1')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$cluster == "hASPC2"][1:20]), ctrl = 50, name = 'hASPC2')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$cluster == "hASPC3"][1:20]), ctrl = 50, name = 'hASPC3')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$cluster == "hASPC4"][1:20]), ctrl = 50, name = 'hASPC4')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$cluster == "hASPC5"][1:20]), ctrl = 50, name = 'hASPC5')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$cluster == "hASPC6"][1:20]), ctrl = 50, name = 'hASPC6')
dropEST.combined.filtered$hASPC1 <- dropEST.combined.filtered$hASPC11
dropEST.combined.filtered$hASPC2 <- dropEST.combined.filtered$hASPC21
dropEST.combined.filtered$hASPC3 <- dropEST.combined.filtered$hASPC31
dropEST.combined.filtered$hASPC4 <- dropEST.combined.filtered$hASPC41
dropEST.combined.filtered$hASPC5 <- dropEST.combined.filtered$hASPC51
dropEST.combined.filtered$hASPC6 <- dropEST.combined.filtered$hASPC61
dropEST.combined.filtered$hASPC11 <- NULL
dropEST.combined.filtered$hASPC21 <- NULL
dropEST.combined.filtered$hASPC31 <- NULL
dropEST.combined.filtered$hASPC41 <- NULL
dropEST.combined.filtered$hASPC51 <- NULL
dropEST.combined.filtered$hASPC61 <- NULL
dropEST.combined.filtered[['Emont.hModule']]<- CreateAssayObject(data = t(x = FetchData(object = dropEST.combined.filtered, 
                                                                              vars = c('hASPC1','hASPC2','hASPC3','hASPC4','hASPC5','hASPC6'))))
p = DoHeatmap(object = dropEST.combined.filtered, 
              features = c('hASPC1','hASPC2','hASPC3','hASPC4','hASPC5','hASPC6'), 
              assay = 'Emont.hModule',
              slot = 'data'
) + 
  scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +
  theme(text=element_text(size=22, family="Arial"),
        plot.title=element_text(size=24, face="italic"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-10,10),
        axis.text.y=element_text(size=18),
        axis.text.x= element_blank()
  ) 
ggsave(filename = "Emont_human_top20markers_50ctrPerMarker.jpeg", plot = p, device = "jpeg", path = "output_plot/NewMarkerComparison/Emont_hASPC/top20markers_50ctrPerMarker/", height=4.5,width=7, units="in", dpi=300)

DefaultAssay(dropEST.combined.filtered) <- 'RNA' 
for (subtype in c('hASPC1','hASPC2','hASPC3','hASPC4','hASPC5','hASPC6')){
  p <- FeaturePlot(dropEST.combined.filtered, 
                   features = subtype, 
                   reduction = "tsne", 
                   cols = c("grey","red"), 
                   slot = "data",
                   min.cutoff = 0,
                   max.cutoff = 1,
                   pt.size = 0.01) +
    scale_colour_gradient(limits = c(0.00,1.00), breaks = c(0.0, 0.25, 0.50, 0.75, 1.00), low = "lightgrey", high = "red") +
    theme(text = element_text(size=22,family="Arial"),
          plot.title=element_text(size=24,face = "plain"),
          axis.text.x=element_text(size=18),
          axis.text.y=element_text(size=18)) + 
    labs(title = paste0(subtype," markers"))
  ggsave(filename=paste0("Emont_",subtype,"_min0max1_featureplot.jpeg"), 
         plot=p, 
         device="jpeg",
         path="output_plot/NewMarkerComparison/Emont_hASPC/top20markers_50ctrPerMarker/", height=4,width=5, units="in", dpi=300)
}


Massier_ofC_markers <- read_excel("resource/NewMarkersToCompare.xlsx", 
                                  sheet = "Massier_NatCom2023_h_omFAP")
Massier_ofC_markers$Cell_type <- paste0("ofC",Massier_ofC_markers$Cell_type)
filtered_df <- Massier_ofC_markers %>% 
  filter(p_val_adj < 0.05)
result_df <- filtered_df %>% 
  group_by(Cell_type) %>% 
  arrange(desc(avg_log2FC)) %>% 
  ungroup()
DefaultAssay(dropEST.combined.filtered) <- "RNA"
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC0"][1:20]), ctrl = 50, name = 'ofC0')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC1"][1:20]), ctrl = 50, name = 'ofC1')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC2"][1:20]), ctrl = 50, name = 'ofC2')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC3"][1:20]), ctrl = 50, name = 'ofC3')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC4"][1:20]), ctrl = 50, name = 'ofC4')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC5"][1:20]), ctrl = 50, name = 'ofC5')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC6"][1:20]), ctrl = 50, name = 'ofC6')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC7"][1:20]), ctrl = 50, name = 'ofC7')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC8"][1:20]), ctrl = 50, name = 'ofC8')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC9"][1:20]), ctrl = 50, name = 'ofC9')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC10"][1:20]), ctrl = 50, name = 'ofC10')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC11"][1:20]), ctrl = 50, name = 'ofC11')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC12"][1:20]), ctrl = 50, name = 'ofC12')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC13"][1:20]), ctrl = 50, name = 'ofC13')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$Cell_type == "ofC14"][1:20]), ctrl = 50, name = 'ofC14')
dropEST.combined.filtered$ofC0 <- dropEST.combined.filtered$ofC01
dropEST.combined.filtered$ofC1 <- dropEST.combined.filtered$ofC11
dropEST.combined.filtered$ofC2 <- dropEST.combined.filtered$ofC21
dropEST.combined.filtered$ofC3 <- dropEST.combined.filtered$ofC31
dropEST.combined.filtered$ofC4 <- dropEST.combined.filtered$ofC41
dropEST.combined.filtered$ofC5 <- dropEST.combined.filtered$ofC51
dropEST.combined.filtered$ofC6 <- dropEST.combined.filtered$ofC61
dropEST.combined.filtered$ofC7 <- dropEST.combined.filtered$ofC71
dropEST.combined.filtered$ofC8<- dropEST.combined.filtered$ofC81
dropEST.combined.filtered$ofC9 <- dropEST.combined.filtered$ofC91
dropEST.combined.filtered$ofC01 <- NULL
dropEST.combined.filtered$ofC11 <- NULL
dropEST.combined.filtered$ofC21 <- NULL
dropEST.combined.filtered$ofC31 <- NULL
dropEST.combined.filtered$ofC41 <- NULL
dropEST.combined.filtered$ofC51 <- NULL
dropEST.combined.filtered$ofC61 <- NULL
dropEST.combined.filtered$ofC71 <- NULL
dropEST.combined.filtered$ofC81 <- NULL
dropEST.combined.filtered$ofC91 <- NULL
dropEST.combined.filtered$ofC10 <- dropEST.combined.filtered$ofC101
dropEST.combined.filtered$ofC11 <- dropEST.combined.filtered$ofC111
dropEST.combined.filtered$ofC12 <- dropEST.combined.filtered$ofC121
dropEST.combined.filtered$ofC13 <- dropEST.combined.filtered$ofC131
dropEST.combined.filtered$ofC14 <- dropEST.combined.filtered$ofC141
dropEST.combined.filtered$ofC101 <- NULL
dropEST.combined.filtered$ofC111 <- NULL
dropEST.combined.filtered$ofC121 <- NULL
dropEST.combined.filtered$ofC131 <- NULL
dropEST.combined.filtered$ofC141 <- NULL
dropEST.combined.filtered[['Massier.hModule']]<- CreateAssayObject(data = t(x = FetchData(object = dropEST.combined.filtered, 
                                                                                vars = c('ofC0','ofC1','ofC2','ofC3','ofC4','ofC5','ofC6','ofC7','ofC8','ofC9','ofC10','ofC11','ofC12','ofC13','ofC14'))))
p = DoHeatmap(object = dropEST.combined.filtered, 
              features = c('ofC0','ofC1','ofC2','ofC3','ofC4','ofC5','ofC6','ofC7','ofC8','ofC9','ofC10','ofC11','ofC12','ofC13','ofC14'), 
              assay = 'Massier.hModule',
              slot = 'data'
) + 
  scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +
  theme(text=element_text(size=22, family="Arial"),
        plot.title=element_text(size=24, face="italic"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-10,10),
        axis.text.x= element_blank(),
        axis.text.y=element_text(size=18)
  ) 
ggsave(filename = "Massier_human_top20markers_50ctrPerMarker.jpeg", plot = p, device = "jpeg", path = "output_plot/NewMarkerComparison/Massier_hASPC/top20markers_50ctrPerMarker/", height=4.5,width=7, units="in", dpi=300)

DefaultAssay(dropEST.combined.filtered) <- 'RNA' 
for (subtype in c('ofC0','ofC1','ofC2','ofC3','ofC4','ofC5','ofC6','ofC7','ofC8','ofC9','ofC10','ofC11','ofC12','ofC13','ofC14')){
  p <- FeaturePlot(dropEST.combined.filtered, 
                   features = subtype, 
                   reduction = "tsne", 
                   cols = c("grey","red"), 
                   slot = "data",
                   min.cutoff = 0,
                   max.cutoff = 1,
                   pt.size = 0.01) +
    theme(text = element_text(size=22,family="Arial"),
          plot.title=element_text(size=24,face = "plain"),
          axis.text.x=element_text(size=18),
          axis.text.y=element_text(size=18)) + 
    labs(title = paste0(subtype," markers"))
  
  png(paste0("output_plot/NewMarkerComparison/Massier_hASPC/top20markers_50ctrPerMarker/Massier_",subtype,"_min0max1_featureplot.png"),width = 1700, height = 1500, res = 300, type = "cairo")
  plot(p)
  dev.off()
}

eLife_Holman_SeuratObj <- readRDS("./resource/SeuratObj_forRefMarkers/Holman2024_SingleCell.rds")
View(eLife_Holman_SeuratObj)
unique(eLife_Holman_SeuratObj$Named_clusters)
eLife_Holman_ASPC_SeuratObj <- subset(eLife_Holman_SeuratObj, idents = c("Icam1+ Preadipocytes",
                                                                         "Dpp4+ Fibroblasts",
                                                                         "Cd142+ Fibroblasts",
                                                                         "Spp1+ Fibroblasts"))
eLife_Holman_ASPC_markers <- FindAllMarkers(eLife_Holman_ASPC_SeuratObj,
                                            only.pos = TRUE,
                                            assay = "RNA")
eLife_Holman_ASPC_markers.sig <- eLife_Holman_ASPC_markers[eLife_Holman_ASPC_markers$p_val_adj < 0.05,]
View(eLife_Holman_ASPC_markers.sig)
saveRDS(eLife_Holman_ASPC_markers.sig,"resource/Homan_eLife2024_ASPC_markers.rds")
dropEST.combined.filtered <- readRDS("intermediate_data/SeuratObj/mASPC_SeuratObj.rds")
eLife_Holman_ASPC_markers.sig <- readRDS("resource/Homan_eLife2024_ASPC_markers.rds")
unique(eLife_Holman_ASPC_markers.sig$cluster)
filtered_df <- eLife_Holman_ASPC_markers.sig %>% 
  filter(p_val_adj < 0.05)
result_df <- filtered_df %>% 
  group_by(cluster) %>% 
  arrange(desc(avg_log2FC)) %>% 
  ungroup()
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$cluster == "Dpp4+ Fibroblasts"][1:20]), ctrl = 50, name = 'Dpp4+ Fibroblasts')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$cluster == "Icam1+ Preadipocytes"][1:20]), ctrl = 50, name = 'Icam1+ Preadipocytes')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$cluster == "Cd142+ Fibroblasts"][1:20]), ctrl = 50, name = 'Cd142+ Fibroblasts')
dropEST.combined.filtered <- AddModuleScore(object = dropEST.combined.filtered, features = list(result_df$gene[result_df$cluster == "Spp1+ Fibroblasts"][1:20]), ctrl = 50, name = 'Spp1+ Fibroblasts')
dropEST.combined.filtered$`Dpp4+ Fibroblasts` <- dropEST.combined.filtered$`Dpp4+ Fibroblasts1`
dropEST.combined.filtered$`Icam1+ Preadipocytes` <- dropEST.combined.filtered$`Icam1+ Preadipocytes1`
dropEST.combined.filtered$`Cd142+ Fibroblasts` <- dropEST.combined.filtered$`Cd142+ Fibroblasts1`
dropEST.combined.filtered$`Spp1+ Fibroblasts` <- dropEST.combined.filtered$`Spp1+ Fibroblasts1`
dropEST.combined.filtered$`Dpp4+ Fibroblasts1` <- NULL
dropEST.combined.filtered$`Icam1+ Preadipocytes1` <- NULL
dropEST.combined.filtered$`Cd142+ Fibroblasts1` <- NULL
dropEST.combined.filtered$`Spp1+ Fibroblasts1` <- NULL
cell_type_genes=names(x = dropEST.combined.filtered[[]])[38:49]
for(i in 1:length(cell_type_genes)){
  p=FeaturePlot(object = dropEST.combined.filtered, 
                features = cell_type_genes[i],
                reduction = "tsne",
                min.cutoff = 0,
                max.cutoff = 1,
                pt.size = 0.01,
                cols =c("lightgrey", "red"))+
    theme(text=element_text(size=22,family="Arial"), 
          plot.title=element_text(size=24,face="italic"),
          axis.text.x=element_text(size=18),
          axis.text.y=element_text(size=18))
  ggsave(filename=paste0("eWAcompare",cell_type_genes[i],".jpeg"), plot=p, device="jpeg",
         path="/Users/gaoyanli/Desktop/Agingdraft/LIFRtry/compare", height=4,width=5, units="in", dpi=300)
}

dropEST.combined.filtered[['Holman.mModule']]<- CreateAssayObject(data = t(x = FetchData(object = dropEST.combined.filtered, 
                                                                              vars = c('Dpp4+ Fibroblasts',
                                                                                       'Icam1+ Preadipocytes',
                                                                                       'Cd142+ Fibroblasts',
                                                                                       'Spp1+ Fibroblasts'))))
p = DoHeatmap(object = dropEST.combined.filtered, 
              features = c('Dpp4+ Fibroblasts',
                           'Icam1+ Preadipocytes',
                           'Cd142+ Fibroblasts',
                           'Spp1+ Fibroblasts'), 
              assay = 'Holman.mModule',
              slot = 'data'
) + 
  scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 0, guide = "colourbar", aesthetics = "fill") +
  theme(text=element_text(size=22, family="Arial"),
        plot.title=element_text(size=24, face="italic"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-10,10),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=18)) 
ggsave(filename = "Holman_mouse_top20markers_50ctrPerMarker.jpeg", plot = p, device = "jpeg", path = "output_plot/NewMarkerComparison/Holman_mASPC/top20markers_50ctrPerMarker/", height=4.5,width=7, units="in", dpi=300)

DefaultAssay(dropEST.combined.filtered) <- 'RNA'
names(x = dropEST.combined.filtered[[]])
for (subtype in c('Dpp4+ Fibroblasts',
                  'Icam1+ Preadipocytes',
                  'Cd142+ Fibroblasts',
                  'Spp1+ Fibroblasts')){
  p <- FeaturePlot(dropEST.combined.filtered, 
                   features = subtype, 
                   reduction = "tsne",
                   min.cutoff = 0,
                   max.cutoff = 1,
                   pt.size = 0.01,
                   cols =c("lightgrey", "red"))  +
    scale_colour_gradient(limits = c(0.00,1.00), breaks = c(0.0, 0.25, 0.50, 0.75, 1.00), low = "lightgrey", high = "red") +
    theme(text=element_text(size=22,family="Arial"), 
          plot.title=element_text(size=24, face = "plain"),
          axis.text.x=element_text(size=18),
          axis.text.y=element_text(size=18)) + 
    labs(title = paste0(subtype," markers"))
  ggsave(filename=paste0("Holman_",subtype,"_min0max1_featureplot.jpeg"), 
         plot=p, 
         device="jpeg",
         path="output_plot/NewMarkerComparison/Holman_mASPC/top20markers_50ctrPerMarker/", height=4,width=5, units="in", dpi=300)
}


### adipocyte snRNAseq analysis
library(fgsea)
library(ggplot2)
library(patchwork)
library(gggsea)
m1_counts <- Read10X_h5("./raw_data/adipocyte snseq_mTmG/atac_53160_gex_53158/outs/filtered_feature_bc_matrix.h5")
m1_fragpath <- "./raw_data/adipocyte snseq_mTmG/atac_53160_gex_53158/outs/atac_fragments.tsv.gz"
m2_counts <- Read10X_h5("./raw_data/adipocyte snseq_mTmG/atac_53161_gex_53159/outs/filtered_feature_bc_matrix.h5")
m2_fragpath <- "./raw_data/adipocyte snseq_mTmG/atac_53161_gex_53159/outs/atac_fragments.tsv.gz"
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "Ensembl"
genome(annotation) <- "mm10"
SeuratObj_m1 <- CreateSeuratObject(
  counts = m1_counts$`Gene Expression`,
  assay = "RNA"
)
SeuratObj_m2 <- CreateSeuratObject(
  counts = m2_counts$`Gene Expression`,
  assay = "RNA"
)
m1_fragment <- CreateFragmentObject(path = m1_fragpath)
m2_fragment <- CreateFragmentObject(path = m2_fragpath)
m1_atac_counts = m1_counts$Peaks
m2_atac_counts = m2_counts$Peaks
SeuratObj_m1[["ATAC"]] <- CreateChromatinAssay(
  counts = m1_atac_counts,
  sep = c(":", "-"),
  fragments = m1_fragpath,
  annotation = annotation
)
SeuratObj_m2[["ATAC"]] <- CreateChromatinAssay(
  counts = m2_atac_counts,
  sep = c(":", "-"),
  fragments = m2_fragpath,
  annotation = annotation
)
DefaultAssay(SeuratObj_m1) <- "ATAC"
SeuratObj_m1 <- NucleosomeSignal(SeuratObj_m1)
SeuratObj_m1 <- TSSEnrichment(SeuratObj_m1)
DefaultAssay(SeuratObj_m2) <- "ATAC"
SeuratObj_m2 <- NucleosomeSignal(SeuratObj_m2)
SeuratObj_m2 <- TSSEnrichment(SeuratObj_m2)
SeuratObj_m1$orig.ident <- "m1" 
SeuratObj_m2$orig.ident <- "m2"
SeuratObj_merged <- merge(
  x = SeuratObj_m1,
  y = SeuratObj_m2,
  add.cell.ids = c("m1","m2")
)
SeuratObj_merged[["percent.mt"]] <- PercentageFeatureSet(SeuratObj_merged, pattern = "^mt-", assay = "RNA")
SeuratObj_filtered = subset(x = SeuratObj_merged,
                            subset = nFeature_RNA > 200 &
                              nFeature_RNA <  7500 & 
                              nCount_RNA >  500 &
                              nCount_RNA <   30000 & 
                              percent.mt < 0.15 & 
                              TSS.enrichment > 1 & 
                              nucleosome_signal < 2 )
DefaultAssay(SeuratObj_filtered) <- "RNA"
SeuratObj_filtered <- SCTransform(SeuratObj_filtered, 
                                  vars.to.regress = "percent.mt",
                                  variable.features.n = nrow(SeuratObj_filtered@assays$RNA),
                                  verbose = FALSE, vst.flavor = "v2")
SeuratObj_filtered <- RunPCA(SeuratObj_filtered, verbose = FALSE)
SeuratObj_filtered <- RunUMAP(SeuratObj_filtered, dims = 1:30, verbose = FALSE)
DefaultAssay(SeuratObj_filtered) <- "ATAC"
peaks <- CallPeaks(SeuratObj_filtered,
                   macs2.path = "./envs/PeakCalling_analysis/bin/macs2")
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
macs2_counts <- FeatureMatrix(
  fragments = Fragments(SeuratObj_filtered),
  features = peaks,
  cells = colnames(SeuratObj_filtered)
)
SeuratObj_filtered[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = list(m1_fragment,m2_fragment),
  annotation = annotation
)
DefaultAssay(SeuratObj_filtered) <- "peaks"
SeuratObj_filtered <- FindTopFeatures(SeuratObj_filtered, min.cutoff = NULL) 
SeuratObj_filtered <- RunTFIDF(SeuratObj_filtered)
SeuratObj_filtered <- RunSVD(SeuratObj_filtered)
DepthCor(SeuratObj_filtered)

DefaultAssay(SeuratObj_filtered) <- "SCT"
SeuratObj_filtered <- FindMultiModalNeighbors(
  object = SeuratObj_filtered,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
SeuratObj_filtered <- FindClusters(SeuratObj_filtered, graph.name = "wsnn", algorithm = 3)
SeuratObj_filtered@active.assay <- "RNA"
SeuratObj_filtered <- NormalizeData(SeuratObj_filtered)
SeuratObj_filtered <- FindVariableFeatures(SeuratObj_filtered, selection.method = "vst", nfeatures = nrow(SeuratObj_filtered))
SeuratObj_filtered <- ScaleData(SeuratObj_filtered, features = rownames(SeuratObj_filtered))
GFPpos <- WhichCells(SeuratObj_filtered, expression = GFP > 0)
SeuratObj_filtered$GFP_exp <- ifelse(colnames(SeuratObj_filtered) %in% GFPpos, "GFP.pos", "GFP.neg")
new.cluster.ids <- c("0" = "Adipocyte",
                     "1" = "Adipocyte",
                     "2" = "Adipocyte",
                     "3" = "Adipocyte",
                     "4" = "Adipocyte",
                     "5" = "Adipocyte",
                     "6" = "Adipocyte",
                     "7" = "Mac/Mono",
                     "8" = "Male_epi",
                     "9" = "Adipocyte",
                     "10" = "ASPC",
                     "11" = "Adipocyte",
                     "12" = "Adipocyte",
                     "13" = "Mac/Mono",
                     "14" = "ASPC",
                     "15" = "Endo",
                     "16" = "NK/T",
                     "17" = "Endo",
                     "18" = "Adipocyte",
                     "19" = "Meso",
                     "20" = "Pericyte",
                     "21" = "LEC")
SeuratObj_filtered <- RenameIdents(SeuratObj_filtered, new.cluster.ids)
SeuratObj_filtered$Cell_type <- SeuratObj_filtered@active.ident
SeuratObj_adipocyte <- subset(SeuratObj_filtered, idents = "Adipocyte")
Idents(SeuratObj_adipocyte) <- "GFP_exp"
DefaultAssay(SeuratObj_adipocyte) <- "RNA"
unique(Idents(SeuratObj_adipocyte))
DEGs <- FindMarkers(SeuratObj_adipocyte,
                    ident.1 = "GFP.neg",
                    ident.2 = "GFP.pos",
                    logfc.threshold = 0,
                    test.use = "wilcox",
                    min.pct = 0.1)
DEGs$gene <- rownames(DEGs)
DEGs$comparison <- "GFPneg.vs.GFPpos"
DEGs_ranked <- DEGs %>% arrange(desc(avg_log2FC)) 
DEGs_ranked_vector <- setNames(DEGs_ranked$avg_log2FC, DEGs_ranked$gene)
Saul_NatCom2022_SenMayoPanelGenes_mouse <- read_excel("resource/Saul_NatCom2022_SenMayoPanelGenes.xlsx", 
                                                      sheet = "mouse")
SenMayo_m <- list(Saul_NatCom2022_SenMayoPanelGenes_mouse$`Gene(murine)`)
names(SenMayo_m) <- "Senescence Markers From SenMayo (Saul et al. 2022)"
fgsea_results <- fgsea(pathways = SenMayo_m, stats = DEGs_ranked_vector, minSize = 15, maxSize = 500)
df <- gseaCurve(DEGs_ranked_vector, SenMayo_m, fgsea_results)
pdf("output_plot/Senencence_markers/SenMayo_mouse_GSEA_GFPneg.v.GFPpos_Adipocyte_ColAdjusted.pdf", width = 7, height = 5)
ggplot2::ggplot() + 
  geom_gsea(df,titlelength = 100) +
  theme_gsea(12) 
dev.off()
