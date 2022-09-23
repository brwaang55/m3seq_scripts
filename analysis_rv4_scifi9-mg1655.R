rm(list=ls())
set.seed(4913, kind = "L'Ecuyer-CMRG")
library(FNN)
library(data.table)
library(readr) 
library(zoo)
library(stringr)
library(RColorBrewer)
library(Seurat)
library(sctransform)
library(ggplot2)
library(gplots)
library(ggpmisc)
library(MASS)
library(fitdistrplus)
library(psych)
library(ggrepel)
figs.out = '~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/scifi9_figs/'
dir.create(figs.out)
#Colors
Breaks=seq(0,60,1)
#dir.create('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/slides/Figures/2020-07-20/')
Colors=rev(brewer.pal(11,"Spectral"))
colors=colorRampPalette(Colors)(120)
#Exponential + Stationary
#data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/petri_seq_with_anti.csv')
data = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/expression_filtered_scifi9_MG1655_25.csv')
cols = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/gene_index_scifi9_MG1655_25.csv')
annotation = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/cell_index_scifi9_MG1655_25.csv')
#colnames(data) = cols[,1]
#colnames(data) = cols[,1]
annotation$treatment = gsub('Lib1_HB_Tn5_MG1655_','',annotation$treatment)
data['Label'] = annotation$treatment
data['Identity'] = annotation$identity
#data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/Data/MG1655_genes_no_ribo_stationary.csv')
"Label" %in% names(data)
#data = data[,-1]
#### data_past ###
data_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/expression_filtered_scifi10_MG1655_25.csv')
cols_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/gene_index_scifi10_MG1655_25.csv')
annotation_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/cell_index_scifi10_MG1655_25.csv')
annotation_past$treatment = gsub('Lib2_HB_Tn5_MG1655_','',annotation_past$treatment)
data_past['Label'] = annotation_past$treatment
data_past['Identity'] = annotation_past$identity
#annotation_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/cell_index_scifi7_len40.csv')
#annotation_past$treatment = str_replace(annotation_past$treatment,'2HR_CEFOZOLIN','2hr_Cipro')
#annotation_past$treatment = str_replace(annotation_past$treatment,'2HR_CIPRO','2hr_Cef')
#annotation_past$treatment = str_replace(annotation_past$treatment,'4HR_CEFOZOLIN','4hr_Cef')
#colnames(data) = cols[,1]
annotation_past$treatment = gsub('KLENOW_Lib1_Tn5_','Past_',annotation_past$treatment)
data_past['Label'] = annotation_past$treatment
data_past['Identity'] = annotation_past$identity
#data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/Data/subtilis_genes_no_ribo_stationary.csv')
cols_past = names(data_past)
# Make minimum 5 observations

### Make this only  E. coli ####
#ec_names = names(data)[str_detect(names(data),'EC_')]
lambda_names = names(data)[str_detect(names(data),'Lambda_')]
nin_names = names(data2)[str_detect(names(data),'Nin')]
ec_names = names(data)[str_detect(names(data),'EC_')]
#data = data[ec_names]
#data_past = data_past[ec_names]
data['Label'] = annotation$treatment
data['Identity'] = annotation$identity
data_past['Label'] = annotation_past$treatment
data_past['Identity'] = annotation_past$identity
#keep_rows  = which(data_past$Identity == 'coli')
#keep_rows  = which(data$Identity == 'coli' & data$Label !="CRISPRI" & data$Label != "M9"& data$Label != "30MIN_FIX")
#data_past = data_past[keep_rows,]
data_total = rbind(data,data_past)
data_past = 9
#data_total = rbind(data)
#data_total = rbind(data_past)
data2 = data_total
#data2 = data
#keep_rows  = which(data2$Identity == 'MG1655')
#data2 = data[keep_rows,]
#keep_rows  = which(data$Identity == 'MG1655'  & data$Label != "PENTA_EXP"& data$Label != "30MIN_FIX")
#unique(data_past$Label)
#data_total = rbind(data2,data_past)
#data_total = rbind(data2)
#data_total = rbind(data2,data_past)
#write.csv(data_total,file='~/Documents/Data/subtilis_scifi7_8.csv')
#write.csv(data_total,file='~/Documents/Data/subtilis_coli7_8.csv')
unique(data_total$Label)
### Make this only E. MG1655 ####
data = data_total
#data = data[,-(dim(data)[2])]
cols = names(data)
# Make minimum 5 observations
k = colSums(data!=0)
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]

unique(data2$Label)
keep_rows  = which(data$Label !="Overnight")#| data$Label == "Past_STATIONARY" |data$Label == "Past_EXP")
data = data[keep_rows,]
keep_rows  = which(data$Label !="Lambda90")#| data$Label == "Past_STATIONARY" |data$Label == "Past_EXP")
data = data[keep_rows,]
keep_rows  = which(data$Label !="Lambda30")#| data$Label == "Past_STATIONARY" |data$Label == "Past_EXP")
data = data[keep_rows,]
unique(data$Label)
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','EC-ECU-04345'))])
#counts = data.frame(log(m, base = 10),data$Label)
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
ggplot(data = counts, aes(x = `count`, fill = `annotation`,color = `annotation`)) + 
  geom_histogram(alpha = 0.3 ,position = 'identity') + 
  #  geom_point(aes(color = factor(Label) ))+
  #  labs(x = "Total mRNA", y = "rRNA", color = "Label") +
  #  ggtitle(sprintf("Size Scaling of all cells")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 


# Filter out the operons with less than 10 obs total

#sc_precursor = data[,-which(names(data) %in% c('Label','Identity',lambda_names))]
sc_precursor = data[,-which(names(data) %in% c('Label','Identity'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
rm(sc_precursor)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc, normalization.method = 'LogNormalize',scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
#sct = SCTransform(sc)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sct), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sct)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)
rm(mat)
rm(sc_precursor)
rm(sc)
rm(data)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
#FindMarkers(sct, ident.1 = 'Exponential')

print(sct[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sct, dims = 1:2, reduction = "pca",nfeatures = 20, balanced = TRUE)
#ggsave(paste(figs.out,'dim_loadings_MG1655','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/pca_loadings.png')
#DimHeatmap(sct, dims = 1:3, cells = 500, balanced = TRUE)
ElbowPlot(sct,ndims = 15)

names(data)[str_detect(names(data),'stf')]

DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_MG1655_scifi7','.png'))
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/PCA_Relative_Count.png')

sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_MG1655_scifi_clusters','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scBactoSeq/Figures/2020-05-12/UMAP_Relative_Count.png')
#DimPlot(sct, label = TRUE, reduction = 'tsne')

#ggsave(paste(figs.out,'umap_clustered_MG1655_scifi3','.png'))
#DimPlot(sct, label = TRUE ,  reduction = 'pca') + ggtitle('Clusters (PCA)- Relative Counts') + theme_bw() +
#  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
g = Idents(sct)
#g[which(g == 2)]
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
#cluster1.markers <- FindMarkers(sct, ident.1 = 1, ident.2=0,min.pct = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2,  min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3,  min.pct = 0.05,logfc.threshold = 0.1)
#cluster3.markers <- FindMarkers(sct, ident.1 = 3, ident.2 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, ident.2 = 0, min.pct = 0.1,logfc.threshold = 0.1)
#cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, ident.2 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6,logfc.threshold = 0.1)
#cluster6.markers <- FindMarkers(sct, ident.1 = 6,ident.2=4, min.pct = 0.1)
#cluster7.markers <- FindMarkers(sct, ident.1 = 7, ident.2 = 1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7,ident.2 = 0,logfc.threshold = 0.1)
cluster8.markers <- FindMarkers(sct, ident.1 = 8, ident.2 = 0,min.pct = 0.1,logfc.threshold = 0.1)
cluster9.markers <- FindMarkers(sct, ident.1 = 9, ident.2 = 0,min.pct = 0.1,logfc.threshold = 0.1)
cluster10.markers <- FindMarkers(sct, ident.1 = 10,  ident.2 = 0,min.pct = 0.1, logfc.threshold = 0.1)
cluster11.markers <- FindMarkers(sct, ident.1 = 11,  ident.2 = 0, min.pct = 0.1, logfc.threshold = 0.1)
cluster14.markers <- FindMarkers(sct, ident.1 = 14,  min.pct = 0.1, logfc.threshold = 0.1)
cluster15.markers <- FindMarkers(sct, ident.1 = 15,  min.pct = 0.1, logfc.threshold = 0.1)
cluster16.markers <- FindMarkers(sct, ident.1 = 16,  min.pct = 0.1, logfc.threshold = 0.1)
cluster13.markers <- FindMarkers(sct, ident.1 = 13,  min.pct = 0.1, logfc.threshold = 0.1)
cluster12.markers <- FindMarkers(sct, ident.1 = 12,  ident.2 = 0,min.pct = 0.1, logfc.threshold = 0.1)
cluster17.markers <- FindMarkers(sct, ident.1 = 17,  min.pct = 0.1, logfc.threshold = 0.1)
cluster18.markers <- FindMarkers(sct, ident.1 = 18,  min.pct = 0.1, logfc.threshold = 0.1)
#pinR, tfaR, ydfK, rhsB/rhsA
cluster_markers = data.frame(FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1))
cluster_markers$cluster = 0
cluster_markers$gene = rownames(cluster_markers)
idents = sct@active.ident
for (ident in unique(idents)){
  print(ident)
  cluster = data.frame(FindMarkers(sct, ident.1 = ident, min.pct = 0.1,logfc.threshold = 0.1))
  cluster$cluster = ident
  cluster$gene = rownames(cluster)
  cluster_markers = rbind(cluster_markers, cluster)
  
  
}
write.csv(x = cluster_markers, file = paste(figs.out, 'static_cluster_markers','.csv', sep = ''))
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
ggsave(paste(figs.out,'umap_condition_MG1655_scifi9_static','.png'))
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

ggsave(paste(figs.out,'umap_clustered_MG1655_scifi9_static','.png'))
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_MG1655_scifi3','.png'))

#### Go Term ####

library(topGO)
ontology = 'BP'
# install.packages("org.Mm.eg.db") 
# Take from the markers
#cluster2.marker_subset = cluster2.markers[which(cluster2.markers$avg_log2FC > 0),]
#rownames(cluster2.marker_subset) = gsub(pattern = 'EC-',replacement = '',x=rownames(cluster2.marker_subset))
use_markers = cluster2.markers
expressed.genes <- rownames(use_markers[intersect(which(use_markers$avg_log2FC > 0), 
                                                       which(use_markers$p_val< 0.05)),])
down.genes <- rownames(use_markers[intersect(which(use_markers$avg_log2FC < 0), 
                                                  which(use_markers$p_val< 0.05)),])
expressed.genes = gsub(pattern = 'EC-',replacement = '',x = expressed.genes)
down.genes = gsub(pattern = 'EC-',replacement = '',x = down.genes)
all.genes = gsub(pattern = 'EC-',replacement = '', x = rownames(sct@assays$RNA))

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
geneList2 <- ifelse(all.genes %in% down.genes, 1, 0)
names(geneList) <- all.genes
names(geneList2) <- all.genes
GOdata <- new("topGOdata",
              ontology = ontology, # use biological process ontology
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.EcK12.eg.db", ID = "symbol")
#results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher", cutOff = 0.05)
goEnrichment = GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
#goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment$Fisher <- as.numeric(goEnrichment$Fisher)
goEnrichment <- goEnrichment[goEnrichment$Fisher<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","Fisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))


GOdata <- new("topGOdata",
              ontology = ontology, # use biological process ontology
              allGenes = geneList2,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.EcK12.eg.db", ID = "symbol")
#results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
resultFisher2 <- runTest(GOdata, algorithm = "elim", statistic = "fisher",cutOff = 0.05)
goEnrichment2 = GenTable(GOdata, Fisher = resultFisher2, topNodes = 20, numChar = 60)
#goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment2$Fisher <- as.numeric(goEnrichment2$Fisher)
goEnrichment2 <- goEnrichment2[goEnrichment2$Fisher<0.05,]
goEnrichment2 <- goEnrichment2[,c("GO.ID","Term","Fisher")]
goEnrichment2$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment2$Term)
goEnrichment2$Term <- gsub("\\.\\.\\.$", "", goEnrichment2$Term)
goEnrichment2$Term <- paste(goEnrichment2$GO.ID, goEnrichment2$Term, sep=", ")
goEnrichment2$Term <- factor(goEnrichment2$Term, levels=rev(goEnrichment2$Term))


goEnrichment$Fish = -log10(goEnrichment$Fisher)
goEnrichment2$Fish = log10(goEnrichment2$Fisher)
goEnrichment = rbind(goEnrichment, goEnrichment2)

goEnrichment$Term = gsub(pattern = 'GO:[0-9]+, ', replacement = '',x = goEnrichment$Term)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
library(stringr)
goEnrichment$Term = firstup(goEnrichment$Term)
goEnrichment = goEnrichment[order(goEnrichment$Fish),]
#goEnrichment$Term <- factor(goEnrichment$Term, levels=(goEnrichment$Term))
#goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
#goEnrichment$Fish <- factor(goEnrichment$Fish, levels=unique(sort(goEnrichment$Fish)))
ggplot(goEnrichment, aes(x=reorder(Term,Fish), y=Fish, fill = `Fish`)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p value)") +
  ggtitle("") +  
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Fisher)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    ##    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=8, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, face="bold", hjust=1.10, color = 'black'),
    axis.text.y=element_text(angle=0, size=10, face="bold", vjust=0.5, color = 'black'),
    axis.title=element_text(size=0.01, face="bold"),
    ##    legend.key=element_blank(),     #removes the border
    ##    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    #    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

#### UMAP Variability ###

#### Variablee Genes


umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out$label = gsub('Post_Exponential','Batch2_OD_0.6',umap_out$label)
#umap_out$label = gsub('Post_Stationary','Batch2_OD_2.8',umap_out$label)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)
umap_out$label = gsub('Post',' ',umap_out$label)
umap_out$label = gsub('3hr Cipro Rep1','90 min Cipro',umap_out$label)
umap_out$label = gsub('3hr Cef Rep1','90 min Cef',umap_out$label)
umap_out$label = gsub('90 min','T90',umap_out$label)
umap_out$label = gsub('8hr','T360',umap_out$label)

umap_out$label = gsub('Gent','Gentamycin',umap_out$label)
umap_out$label = gsub('Cipro','Ciprofloxacin',umap_out$label)
umap_out$label = gsub('Cef','Cefazolin',umap_out$label)
umap_out$label = gsub('Nal','Nalidixic acid',umap_out$label)
#umap_out$label = gsub('Cycloserineserine','Cycloserine',umap_out$label)
umap_out$label = gsub('Cyclo','Cycloserine',umap_out$label)
#umap_out$label = gsub('Erythromycinmycin','Erythromycin',umap_out$label)
umap_out$label = gsub('Erythro','Erythromycin',umap_out$label)
#umap_out$label = gsub('Tetracyclineracycline','Tetracycline',umap_out$label)
umap_out$label = gsub('Tet','Tetracycline',umap_out$label)
umap_out$label = gsub('Chlor','Chloramphenicol',umap_out$label)
unique(umap_out$label)
#umap_out$label = gsub('8hr','T360',umap_out$label)
#colors = data.frame(read.csv(paste(figs.out,'color_map.csv'), sep = ','))[,-1]
#umap_out = umap_out[umap_out$label=='STATIONARY',]
colors
#umap_out = merge(umap_out, colors,by.x = 'label',by.y = "umap_out.label")
library(ggsci)

library(scales)
colors = hue_pal()(length(unique(umap_out$ident)))
umap_out$ident <- factor(umap_out$ident, levels=as.character(sort(as.numeric(unique(umap_out$ident)))), labels=as.character(sort(as.numeric(unique(umap_out$ident)))))

aggregate(umi_counts~ident,umap_out,mean)
sizes = aggregate(umi_counts~ident,umap_out,length)
total = sum(aggregate(umi_counts~ident,umap_out,length)['umi_counts'])
sizes['umi_counts']/total
umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.5,alpha = 1.5) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("E. coli")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
  scale_color_manual(values =unique(umap_out$colour))
#ggsave(paste(figs.out, 'umap_MG1655_clean.png'))
#
umap_out = umap_out[umap_out$label != 'Overnight',]
#umap_out = umap_out[umap_out$label != 'Lambda90',]
#umap_out = umap_out[umap_out$label != 'Lambda30',]


umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.2,alpha = 0.8, stroke = 0.15) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=6)))+
  #  scale_color_simpsons()s
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,12)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2.png'))
#scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,10,5)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette1.png'))
ggsave(paste(figs.out, 'umap_MG1655_clean_palette_40pcs_all_antibiotics.png'))


levs = as.character(0:40)
umap_out$ident = factor(umap_out$ident, levels = levs)
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.15,alpha = 0.8, stroke = 0.2) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_clusters.png'))
#scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,10,5)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette1.png'))
#ggsave(paste(figs.out, 'umap_MG1655_clusters_40pcs.png'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.05,alpha = 1.5) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
 # ggtitle(sprintf("E. coli")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))#+
  #  scale_color_simpsons()s
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,5)])
ggsave(paste(figs.out, 'umap_MG1655_tet_chlor_clusters.png'))
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2_30_pcs.png'))



#### Gad on the whole thing 

#gad_vec = c('EC_gadA','EC_gadB','EC_gadC')
gad_vec = c('EC_gadA','EC_gadB')
gad_df = data.frame(gad_vec, rep(1,length(gad_vec)))
names(gad_df) = c('gene_id','gad_stage')
gad_df$gene_id = gsub('_','-',gad_df$gene_id)
gad_score = matrix(data = 0,ncol=1,nrow = dim(sct@assays$RNA)[1])
gene_indices = c()
for (i in 1:length(gad_df$gene_id)){
  index = which(gad_df$gene_id[i] == rownames(sct@assays$RNA))
  if (length(index)>0){
    gene_indices = c(gene_indices,index)
    select = c(select,i)
  }
}
gad_df$gene_indices = gene_indices
gad_score[gad_df$gene_indices,1] = as.numeric(gad_df$gad_stage)
norm_data = t(as.data.frame(sct@assays$RNA@data))
gad_data = data.frame(norm_data %*% gad_score)
names(gad_data) = c('norm_gad')

gad_data$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
gad_data$label = as.character(sct[['Label']][,1])
gad_data$ident = as.character(sct@active.ident)
gad_data$index = 1:dim(gad_data)[1]


umap_out2 = cbind(umap_out, gad_data$norm_gad)
umap_out2$max_val = pmax(gad_data$norm_gad)

ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `gad_data$norm_gad`, fill = `gad_data$norm_gad`)) + 
  #  geom_point(pch = 21, aes(fill = `gad_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5) +
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #  scale_fill_gradientn(colors =  (brewer.pal(n = 9, name =  "BuGn") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "BuGn")), guide = '' ) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) +   #scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greens")[1]), high = (brewer.pal(n = 9, name = "Greens")[9])) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'chartreuse1') + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'chartreuse1') + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
  #  scale_colour_gradient(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 

ggsave(paste(figs.out,"gadAB_umap_all",'.png'))
#### Exponential ####

keep_rows  = which(data2$Label =="Exponential" | data2$Label =="Ciproa" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
data = data2[keep_rows,]

# Make minimum 5 observations
k = colSums(data!=0) 
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
#counts = data.frame(log(m, base = 10),data$Label)
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
ggplot(data = counts, aes(x = `count`, fill = `annotation`,color = `annotation`)) + 
  geom_histogram(alpha = 0.3 ,position = 'identity') + 
  #  geom_point(aes(color = factor(Label) ))+
  #  labs(x = "Total mRNA", y = "rRNA", color = "Label") +
  #  ggtitle(sprintf("Size Scaling of all cells")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 


# Filter out the operons with less than 10 obs total
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
to_filter = which(k > 10)
new_names = names(to_filter)

new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
var_x = apply(sc_precursor,2,var)
mean_x = apply(sc_precursor,2,mean)
var_x_norm = apply(as.matrix(sct@assays$RNA@data),2,var)
mean_x_norm = apply(sct@assays$RNA@data,2,mean)
#sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sct), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sct)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)




sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
#FindMarkers(sct, ident.1 = 'Exponential')

mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
var_total = data.frame(cbind(1:length(varExplained), varExplained, rep('Exponential',length(varExplained))))

#kurtosi
pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total2= data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Exponential',dim(pca_out)[2])))
plot(1:50, kurtosi(pca_out))
ggplot(data = data.frame(cbind(1:length(varExplained), varExplained)), aes (x = V1, y = varExplained)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_MG1655_var_explained','.png'))

plot(1:50, kurtosi(pca_out))
#ggsave(paste(figs.out,'pca_MG1655_var_explained','.png'))
rm(mat)
rm(sc_precursor)

names(data)[str_detect(names(data),'xhl')]

DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
#sct <- FindClusters(sct, verbose = TRUE, algorithm=1,resolution=0.9)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
g = Idents(sct)
g[which(g == 2)]
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.1,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.1,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1,logfc.threshold = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1,logfc.threshold = 0.1)
DimPlot(sct, label = TRUE) + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_MG1655_scifi3','.png'))
#### UMAP Variability ###

#### Variablee Genes

library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)
umap_out$label = gsub('Post',' ',umap_out$label)
unique(umap_out$label)
#umap_out$label = gsub('3hr Cipro Rep1','90 min Cipro',umap_out$label)


p1 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`)) + 
  geom_point(pch = 21, aes(fill = umi_counts),size = 1) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMAP for UMI counts")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
#p2 =
library(rcartocolor)
p2 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 1,alpha = 1.5) +  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  ggtitle(sprintf("E. coli")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
  #  scale_color_manual(values=as.vector(palette.colors()[c(2,3,4,5,6,8)]))
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,10,5)])
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)
clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)

clusters_total = data.frame(cbind(clusters, rep('Exponential',dim(clusters)[1])))
names(clusters_total) = c('ident','counts','condition')

umap_out = umap_out[umap_out$label != 'Overnight',]
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette1_30pcs.png'))

umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential'))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.4,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(12)])

ggsave(paste(figs.out, 'umap_MG1655_clean_palette2_40_pcs_exponential.png'))


library(scales)
show_col(as.vector(rcartocolor::carto_pal(name='Safe')))
colors = hue_pal()(length(unique(umap_out$ident)))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.05,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_exponential_phase_scifi9.png'))

#### Gent ####
keep_rows  = which(data2$Label =="Gent" | data2$Label =="Exponentiala" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
data = data2[keep_rows,]
k = colSums(data!=0) 
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
to_filter = which(k > 10)
new_names = names(to_filter)

new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 10000)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
var_total_gent = data.frame(cbind(1:length(varExplained), varExplained, rep('Gent',length(varExplained))))
var_total = data.frame(rbind(var_total, var_total_gent))
rm(mat)
rm(sc_precursor)


pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total_gent = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Gent',length(varExplained))))
kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_gent))

sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)
umap_out$label = gsub('Post',' ',umap_out$label)
umap_out$label = gsub('Gent','Gentamycin',umap_out$label)
unique(umap_out$label)
library(rcartocolor)
p3 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  labs(x = "Dim 1", y = "Dim 2", color = "Cluster") +
  ggtitle(sprintf("UMAP for Cluster")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),panel.grid = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Paired")+ guides(colour = guide_legend(override.aes = list(size=8)))
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)
clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)
clusters_gent = data.frame(cbind(clusters, rep('Gent',dim(clusters)[1])))
names(clusters_gent) = c('ident','counts','condition')
clusters_total = rbind(clusters_total, clusters_gent)

library(scales)
colors = hue_pal()(length(unique(umap_out$ident)))
umap_out$ident <- factor(umap_out$ident, levels=as.character(sort(as.numeric(unique(umap_out$ident)))), labels=as.character(sort(as.numeric(unique(umap_out$ident)))))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.6,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,12)])
ggsave(paste(figs.out, 'umap_MG1655_gent_phase_scifi9.png'))
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2.png'))
#scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,10,5)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette1.png'))

levs = as.character(0:18)
umap_out$ident = factor(umap_out$ident, levels = levs)
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.4,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_gent_phase_scifi9.png'))


#### Chlor ####
keep_rows  = which(data2$Label =="Chlor" | data2$Label =="Teet" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
data = data2[keep_rows,]
library(dplyr)
k = colSums(data!=0) 
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
to_filter = which(k > 10)
new_names = names(to_filter)

new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 10000)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)

rm(mat)
rm(sc_precursor)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
pca <- sct[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
var_total_chlor = data.frame(cbind(1:length(varExplained), varExplained, rep('Chlor',length(varExplained))))
var_total = data.frame(rbind(var_total, var_total_chlor))

pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total_chlor = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Chlor',length(varExplained))))
kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_chlor))


ggplot(data = data.frame(cbind(1:length(varExplained), varExplained)), aes (x = V1, y = varExplained)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
rm(mat)
rm(sc_precursor)
sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.05,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,logfc.threshold = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6,logfc.threshold = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7,logfc.threshold = 0.1)
cluster8.markers <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1,logfc.threshold = 0.1)
cluster9.markers <- FindMarkers(sct, ident.1 = 9, min.pct = 0.1,logfc.threshold = 0.1)
cluster10.markers <- FindMarkers(sct, ident.1 = 10, min.pct = 0.1, logfc.threshold = 0.1)
cluster11.markers <- FindMarkers(sct, ident.1 = 11, min.pct = 0.1, logfc.threshold = 0.1)
cluster14.markers <- FindMarkers(sct, ident.1 = 14, min.pct = 0.1, logfc.threshold = 0.1)
cluster15.markers <- FindMarkers(sct, ident.1 = 15, min.pct = 0.1, logfc.threshold = 0.1)
cluster16.markers <- FindMarkers(sct, ident.1 = 16, min.pct = 0.1, logfc.threshold = 0.1)
cluster13.markers <- FindMarkers(sct, ident.1 = 13, min.pct = 0.1, logfc.threshold = 0.1)
cluster12.markers <- FindMarkers(sct, ident.1 = 12, min.pct = 0.1, logfc.threshold = 0.1)
cluster17.markers <- FindMarkers(sct, ident.1 = 17, min.pct = 0.1, logfc.threshold = 0.1)
#### PCA ####

#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)
unique(umap_out$label)
library(rcartocolor)
p3 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  labs(x = "Dim 1", y = "Dim 2", color = "Cluster") +
  ggtitle(sprintf("UMAP for Cluster")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),panel.grid = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Paired")+ guides(colour = guide_legend(override.aes = list(size=8)))
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)
clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)
clusters_chlor = data.frame(cbind(clusters, rep('Chlor',dim(clusters)[1])))
names(clusters_chlor) = c('ident','counts','condition')
clusters_total = rbind(clusters_total, clusters_chlor)


library(scales)
colors = hue_pal()(length(unique(umap_out$ident)))
umap_out$ident <- factor(umap_out$ident, levels=as.character(sort(as.numeric(unique(umap_out$ident)))), labels=as.character(sort(as.numeric(unique(umap_out$ident)))))

levs = as.character(0:19)
umap_out$ident = factor(umap_out$ident, levels = levs)
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.1,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(2,3,8,12)])
ggsave(paste(figs.out, 'umap_MG1655_chlor_scifi9.png'))



levs = as.character(0:18)
umap_out$ident = factor(umap_out$ident, levels = levs)
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.05,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_chlor_scifi9.png'))
#### Tet ####
keep_rows  = which(data2$Label =="Tet" | data2$Label =="Teet" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
data = data2[keep_rows,]
library(dplyr)
#data = sample_n(data, 300)
k = colSums(data!=0) 
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
to_filter = which(k > 10)
new_names = names(to_filter)

new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)

rm(mat)
rm(sc_precursor)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
pca <- sct[["pca"]]

# Get the total variance:

pca_out = data.frame(sct[['pca']]@cell.embeddings)
#plot(1:50,sort(abs(kurtosi(pca_out)),decreasing = TRUE))
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
var_total_tet = data.frame(cbind(1:length(varExplained), varExplained, rep('Tet',length(varExplained))))
kurtosi_total_tet = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Tet',length(varExplained))))
var_total = data.frame(rbind(var_total, var_total_tet))
kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_tet))

#sct = JackStraw(sct)

print(sct[["pca"]], dims = 1:25, nfeatures = 5)

DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(3,9)) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
ggsave(paste(figs.out, 'pca_pc39_tet.png'))
kurtosi(pca_out$PC_3)
rm(mat)
rm(sc_precursor)
sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.05,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,logfc.threshold = 0.1)
#cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6,logfc.threshold = 0.1)
#cluster6.markers <- FindMarkers(sct, ident.1 = 6,ident.2=4, min.pct = 0.1)
#cluster7.markers <- FindMarkers(sct, ident.1 = 7, ident.2 = 1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7,logfc.threshold = 0.1)
cluster8.markers <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1,logfc.threshold = 0.1)
cluster9.markers <- FindMarkers(sct, ident.1 = 9, min.pct = 0.1,logfc.threshold = 0.1)
cluster10.markers <- FindMarkers(sct, ident.1 = 10, min.pct = 0.1, logfc.threshold = 0.1)
cluster11.markers <- FindMarkers(sct, ident.1 = 11, min.pct = 0.1, logfc.threshold = 0.1,ident.2 = 3)
cluster14.markers <- FindMarkers(sct, ident.1 = 14, min.pct = 0.1, logfc.threshold = 0.1)
cluster15.markers <- FindMarkers(sct, ident.1 = 15, min.pct = 0.1, logfc.threshold = 0.1)
cluster16.markers <- FindMarkers(sct, ident.1 = 16, min.pct = 0.1, logfc.threshold = 0.1)
cluster13.markers <- FindMarkers(sct, ident.1 = 13, min.pct = 0.1, logfc.threshold = 0.1)
cluster12.markers <- FindMarkers(sct, ident.1 = 12, min.pct = 0.1, logfc.threshold = 0.1)
cluster17.markers <- FindMarkers(sct, ident.1 = 17, min.pct = 0.1, logfc.threshold = 0.1)
#my_levels <- c(0,2,3,1)
#levels(sct) = my_levels
#markers %>%
#  group_by(cluster) %>%  top_n(n = 2, wt = avg_log2FC)
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 5, wt = avg_log2FC) -> top10
#DoHeatmap(sct, features = top10$gene,draw.lines = TRUE)  + scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ 
#  theme(text = element_text(size = 18))
#ggsave(paste(figs.out,'MG1655_heatmap_exp_3hr','.png'))
#### PCA ####
#### get the PC's ###

pca_out = data.frame(sct[['pca']]@cell.embeddings)
plot(1:50,kurtosi(pca_out),decreasing = TRUE)
pca_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
pca_out$label = as.character(sct[['Label']][,1])
pca_out$ident = as.character(sct@active.ident)

ggplot(pca_out,aes(x = `PC_3`, y = `PC_9`)) + 
  geom_point( size = 0.05,alpha = 0.8) +  
  geom_density2d()+
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +theme_bw()+
  ggtitle(sprintf("")) + theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),
                               #panel.grid=element_blank(),panel.border = element_blank(),
#                              panel.background = element_blank(), 
#                              axis.text = element_blank(),
                              #        axis. = element_blank(),
                              #        axis.line = element_line(colour = "black"),
                              #        axis.text=element_text(size=20),axis.text.x=element_blank(),
  #                            axis.ticks.x=element_blank(),
  #                            axis.text.y=element_blank(),
  #                            axis.ticks.y=element_blank() ,
                              axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))
  #  scale_color_simpsons()s
#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)
umap_out$label = gsub('Tetracyclineracycline','Tetracycline',umap_out$label)
unique(umap_out$label)
library(rcartocolor)
#ggsave(paste(figs.out, 'umap_cluster_MG1655_total.png'))
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)

clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)
clusters_tet = data.frame(cbind(clusters, rep('Tet',dim(clusters)[1])))
names(clusters_tet) = c('ident','counts','condition')
clusters_total = rbind(clusters_total, clusters_tet)


library(scales)
colors = hue_pal()(length(unique(umap_out$ident)))
umap_out$ident <- factor(umap_out$ident, levels=as.character(sort(as.numeric(unique(umap_out$ident)))), labels=as.character(sort(as.numeric(unique(umap_out$ident)))))

levs = as.character(0:19)
umap_out$ident = factor(umap_out$ident, levels = levs)


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.05,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_tet_scifi9.png'))

#umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.05,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(7,1,11,9,2,3,8,12)])
ggsave(paste(figs.out, 'umap_MG1655_tet_scifi9.png'))



#### Nal ####

keep_rows  = which(data2$Label =="Exponentiaal" | data2$Label =="Nal" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
data = data2[keep_rows,]
#data = sample_n(data, 10000, replace = TRUE)
k = colSums(data!=0) 
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
to_filter = which(k > 10)
new_names = names(to_filter)

new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
var_total_nal = data.frame(cbind(1:length(varExplained), varExplained, rep('Nal',length(varExplained))))
var_total = data.frame(rbind(var_total, var_total_nal))
var_total$varExplained = as.numeric(as.character(var_total$varExplained))
var_total$V1 = as.numeric(as.character(var_total$V1))
ggplot(data = data.frame(cbind(1:length(varExplained), varExplained)), aes (x = V1, y = varExplained)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_MG1655_var_explained','.png'))

ggplot(data = var_total, aes (x = V1, y = varExplained, col = V3)) +
  geom_line() + 
  geom_point(size = 1) +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_MG1655_var_explained','.png'))

rm(mat)
rm(sc_precursor)


pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total_nal = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Nal',length(varExplained))))

kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_nal))

rm(mat)
rm(sc_precursor)

sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.05,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,logfc.threshold = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.05,logfc.threshold = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.05,logfc.threshold = 0.1)
#### PCA ####

#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)
umap_out$label = gsub('Nal','Nalidixic acid',umap_out$label)
unique(umap_out$label)
library(rcartocolor)
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)

clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)
clusters_nal = data.frame(cbind(clusters, rep('Nal',dim(clusters)[1])))
names(clusters_nal) = c('ident','counts','condition')
clusters_total = rbind(clusters_total, clusters_nal)


levs = as.character(0:18)
umap_out$ident = factor(umap_out$ident, levels = levs)
library(scales)
colors = hue_pal()(length(unique(umap_out$ident)))
umap_out$ident <- factor(umap_out$ident, levels=as.character(sort(as.numeric(unique(umap_out$ident)))), labels=as.character(sort(as.numeric(unique(umap_out$ident)))))


umap_out$ident = factor(umap_out$ident, levels = levs)
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.4,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_nal_scifi9.png'))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.4,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(8,12)])
ggsave(paste(figs.out, 'umap_MG1655_nalidixic_acid.png'))




####


#### Ciprofloxacin ####


keep_rows  = which(data2$Label =="Cipro" | data2$Label =="Naal" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
data = data2[keep_rows,]
k = colSums(data!=0) 
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
to_filter = which(k > 10)
new_names = names(to_filter)

new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")

pca <- sct[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
var_total_cipro = data.frame(cbind(1:length(varExplained), varExplained, rep('Cipro',length(varExplained))))
var_total = data.frame(rbind(var_total, var_total_cipro))
var_total$varExplained = as.numeric(as.character(var_total$varExplained))
var_total$V1 = as.numeric(as.character(var_total$V1))

pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total_cipro = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Cipro',length(varExplained))))

kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_cipro))


ggplot(data = var_total, aes (x = V1, y = varExplained, col = V3)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.05,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,logfc.threshold = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.05,logfc.threshold = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.05,logfc.threshold = 0.1)
#my_levels <- c(0,2,3,1)
#levels(sct) = my_levels
#markers %>%
#  group_by(cluster) %>%  top_n(n = 2, wt = avg_log2FC)
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 5, wt = avg_log2FC) -> top10
#DoHeatmap(sct, features = top10$gene,draw.lines = TRUE)  + scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ 
#  theme(text = element_text(size = 18))
#ggsave(paste(figs.out,'MG1655_heatmap_exp_3hr','.png'))
#### PCA ####

#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)
umap_out$label = gsub('Cipro','Ciprofloxacin',umap_out$label)
unique(umap_out$label)
library(rcartocolor)
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)

clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)
clusters_cipro = data.frame(cbind(clusters, rep('Cipro',dim(clusters)[1])))
names(clusters_cipro) = c('ident','counts','condition')
clusters_total = rbind(clusters_total, clusters_cipro)


levs = as.character(0:3)
umap_out$ident = factor(umap_out$ident, levels = levs)
library(scales)
colors = hue_pal()(length(unique(umap_out$ident)))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.2,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_cipro_scifi9.png'))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.5,alpha = 1.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(3,8,12)])
ggsave(paste(figs.out, 'umap_MG1655_cipro.png'))


#### Cefazolin  ####

data2 = data_total
unique(data2$Label)
#keep_rows  = which(data2$Label =="4HR_COCULTURE" | data2$Label =="STATIONARY" | data2$Label == "Coculture_EXP"| data2$Label == "EXP")
#keep_rows  = which( data2$Label == "Cef"| data2$Label == "Cyclo" | data2$Label =="EXP" | data2$Label =="STATIONARY")
keep_rows  = which( data2$Label == "Cef"| data2$Label == "Cyclao" | data2$Label =="EXP" | data2$Label =="STATIONARY")
#keep_rows  = which(data2$Identity == 'MG1655')
data = data2[keep_rows,]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity'))])
counts = data.frame(log(m, base = 10),data$Label)
names(counts) = c('count','annotation')
ggplot(data = counts, aes(x = `count`, fill = `annotation`,color = `annotation`)) + 
  geom_histogram(alpha = 0.3 ,position = 'identity') + 
  #  geom_point(aes(color = factor(Label) ))+
  #  labs(x = "Total mRNA", y = "rRNA", color = "Label") +
  #  ggtitle(sprintf("Size Scaling of all cells")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 


# Filter out the operons with less than 10 obs total
k = colSums(data[-which(names(data) %in% c('Label','Identity'))])
to_filter = which(k > 10)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]

rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sct), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sct)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)

sct <- RunPCA(sct, verbose = FALSE,features = all.genes)

mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
var_total_cef = data.frame(cbind(1:length(varExplained), varExplained, rep('Cef',length(varExplained))))
var_total = data.frame(rbind(var_total, var_total_cef))


pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total_cef = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Cef',length(varExplained))))

kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_cef))



var_total$varExplained = as.numeric(as.character(var_total$varExplained))
var_total$V1 = as.numeric(as.character(var_total$V1))
ggplot(data = data.frame(cbind(1:length(varExplained), varExplained)), aes (x = V1, y = varExplained)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_MG1655_var_explained','.png'))

ggplot(data = var_total, aes (x = V1, y = varExplained, col = V3)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
print(sct[["pca"]], dims = 1:20, nfeatures = 5)
VizDimLoadings(sct, dims = c(8,15), reduction = "pca",nfeatures = 20, balanced = TRUE)
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(3,7)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/PCA_Relative_Count.png')

sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=1,resolution = 0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_MG1655_scifi7','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scBactoSeq/Figures/2020-05-12/UMAP_Relative_Count.png')
#DimPlot(sct, label = TRUE, reduction = 'tsne')
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_MG1655_scifi3','.png'))
DimPlot(sct, label = TRUE ,  reduction = 'pca') + ggtitle('Clusters (PCA)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
g = Idents(sct)
g[which(g == 2)]
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1)
#cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.1)
#cluster3.markers <- FindMarkers(sct, ident.1 = 3,ident.2=4, min.pct = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1)
#cluster5.markers <- FindMarkers(sct, ident.1 = 5, ident.2 = 8,min.pct = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1)
#cluster6.markers <- FindMarkers(sct, ident.1 = 6,ident.2=2, min.pct = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.1)
cluster8.markers <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1)
cluster9.markers <- FindMarkers(sct, ident.1 = 9, min.pct = 0.1)
cluster10.markers <- FindMarkers(sct, ident.1 = 10, min.pct = 0.1)
cluster11.markers <- FindMarkers(sct, ident.1 = 11, min.pct = 0.1)

#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
umap_out$label = gsub('_',' ',umap_out$label)
umap_out$label = gsub('Cef','Cefazolin',umap_out$label)
unique(umap_out$label)
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)

clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)
clusters_cef = data.frame(cbind(clusters, rep('Cef',dim(clusters)[1])))
names(clusters_cef) = c('ident','counts','condition')
clusters_total = rbind(clusters_total, clusters_cef)


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.7,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_cefazin.png'))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.7,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,11,9,2,3,8,12)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2.png'))

ggsave(paste(figs.out, 'umap_MG1655_cef_scifi9.png'))
#### Cyclo ###

data2 = data_total
unique(data2$Label)
#keep_rows  = which(data2$Label =="4HR_COCULTURE" | data2$Label =="STATIONARY" | data2$Label == "Coculture_EXP"| data2$Label == "EXP")
#keep_rows  = which( data2$Label == "Cyclo"| data2$Label == "Cyclo" | data2$Label =="EXP" | data2$Label =="STATIONARY")
keep_rows  = which( data2$Label == "Cyclo"| data2$Label == "Cyclao" | data2$Label =="EXP" | data2$Label =="STATIONARY")
#keep_rows  = which(data2$Identity == 'MG1655')
data = data2[keep_rows,]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity'))])
counts = data.frame(log(m, base = 10),data$Label)
names(counts) = c('count','annotation')
ggplot(data = counts, aes(x = `count`, fill = `annotation`,color = `annotation`)) + 
  geom_histogram(alpha = 0.3 ,position = 'identity') + 
  #  geom_point(aes(color = factor(Label) ))+
  #  labs(x = "Total mRNA", y = "rRNA", color = "Label") +
  #  ggtitle(sprintf("Size Scaling of all cells")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 


# Filter out the operons with less than 10 obs total
k = colSums(data[-which(names(data) %in% c('Label','Identity'))])
to_filter = which(k > 10)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]

rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
#sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 10000)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sct), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sct)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)

sct <- RunPCA(sct, verbose = FALSE,features = all.genes)

mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
var_total_cyclo = data.frame(cbind(1:length(varExplained), varExplained, rep('Cyclo',length(varExplained))))
var_total = data.frame(rbind(var_total, var_total_cyclo))
var_total$varExplained = as.numeric(as.character(var_total$varExplained))

pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total_cyclo = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Cyclo',length(varExplained))))

kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_cyclo))


var_total$V1 = as.numeric(as.character(var_total$V1))
ggplot(data = var_total, aes (x = V1, y = varExplained, col = V3)) +
  geom_line() + 
  geom_point(size = 0.5) +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
rm(mat)
rm(sct2)
rm(sc_precursor)
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=1,resolution = 0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_MG1655_scifi7','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scBactoSeq/Figures/2020-05-12/UMAP_Relative_Count.png')
#DimPlot(sct, label = TRUE, reduction = 'tsne')
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_MG1655_scifi3','.png'))
DimPlot(sct, label = TRUE ,  reduction = 'pca') + ggtitle('Clusters (PCA)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
g = Idents(sct)
g[which(g == 2)]
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1)
#cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.1)
#cluster3.markers <- FindMarkers(sct, ident.1 = 3,ident.2=4, min.pct = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1)
#cluster5.markers <- FindMarkers(sct, ident.1 = 5, ident.2 = 8,min.pct = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1)
#cluster6.markers <- FindMarkers(sct, ident.1 = 6,ident.2=2, min.pct = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.1)
cluster8.markers <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1)
cluster9.markers <- FindMarkers(sct, ident.1 = 9, min.pct = 0.1)
cluster10.markers <- FindMarkers(sct, ident.1 = 10, min.pct = 0.1)
cluster11.markers <- FindMarkers(sct, ident.1 = 11, min.pct = 0.1)
cluster14.markers <- FindMarkers(sct, ident.1 = 14, min.pct = 0.1)
cluster13.markers <- FindMarkers(sct, ident.1 = 13, min.pct = 0.1)
cluster12.markers <- FindMarkers(sct, ident.1 = 12, min.pct = 0.1)

#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('Cyclo','Cycloserine',umap_out$label)
library(rcartocolor)
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)

clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)
clusters_cyclo = data.frame(cbind(clusters, rep('Cyclo',dim(clusters)[1])))
names(clusters_cyclo) = c('ident','counts','condition')
clusters_total = rbind(clusters_total, clusters_cyclo)


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.8,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_cycloserine_scifi9.png'))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.7,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(11,9,2,3,8,12)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2.png'))


ggsave(paste(figs.out, 'umap_MG1655_cyclo_scifi9.png'))


#### Erythro ####


keep_rows  = which(data2$Label =="Erythro" | data2$Label =="Naal" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
data = data2[keep_rows,]
k = colSums(data!=0) 
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
to_filter = which(k > 10)
new_names = names(to_filter)

new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)

mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
var_total_erythro = data.frame(cbind(1:length(varExplained), varExplained, rep('Erythro',length(varExplained))))
var_total = data.frame(rbind(var_total, var_total_erythro))
var_total$varExplained = as.numeric(as.character(var_total$varExplained))
var_total$V1 = as.numeric(as.character(var_total$V1))


pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total_erythro = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Erythro',length(varExplained))))

kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_erythro))



rm(mat)
rm(sc_precursor)
sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.05,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,logfc.threshold = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.05,logfc.threshold = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.05,logfc.threshold = 0.1)
#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
umap_out$label = gsub('Erythro','Erythromycin',umap_out$label)
unique(umap_out$label)
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)

clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)
clusters_erythro = data.frame(cbind(clusters, rep('Erythro',dim(clusters)[1])))
names(clusters_erythro) = c('ident','counts','condition')
clusters_total = rbind(clusters_total, clusters_erythro)



ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.3,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_erythro_scifi9.png'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.4,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(9,2,3,8,12)])
ggsave(paste(figs.out, 'umap_MG1655_erythro.png'))


#### Lambda30 ####
keep_rows  = which(data2$Label =="Lambda30" | data2$Label =="Lambda30s" |  data2$Label == "Exponentiala"|data2$Label=="3hr_Cef_Rep1a")
data = data2[keep_rows,]
k = colSums(data!=0) 
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
to_filter = which(k > 10)
new_names = names(to_filter)

new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)

pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total_lambda30 = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Lambda30',length(varExplained))))
kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_lambda30))


rm(mat)
rm(sc_precursor)
sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.05,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,logfc.threshold = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.05,logfc.threshold = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.05,logfc.threshold = 0.1)
RidgePlot(sct, features = c('Lambda-A','Lambda-B'), ncol = 2)

library(dplyr)

library(dplyr)

av.exp.sct <- AverageExpression(sct, return.seurat = TRUE)
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
#sct@active.ident = factor(sct@meta.data,labels = sct@meta.data$Label, levels = sct@meta.data$Label)
markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.5)
markers = markers[markers$avg_log2FC > 0.5,]
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 3, wt = avg_log2FC) -> top10

markers = markers %>%
  add_count(cluster) %>%
  group_by(cluster) %>%
  mutate(n = 8) %>%
  #  group_by(rating, .add = TRUE) %>%
  #In old dplyr add = TRUE
  #group_by(rating, add = TRUE) %>%
  sample_n(n, replace = TRUE) %>%
  select(-n)

markers = unique(markers)
#markers = markers[1:8,]
genes = markers$gene
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
g = merge(g, markers)
#g = g[which(g$variable%in% genes),]
g = g[order(g$cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'EC-',replacement = '', x = g$gene)
g$gene = gsub(pattern = 'Lambda-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = g$gene,labels = g$gene)
#g$Cluster = factor(g$Cluster, levels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'),labels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'))
#g$Cluster = factor(g$Cluster, levels = c(3,2,1))
g = g[order(g$cluster),]

library(myriad)
ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.0), oob = scales::squish)+
  labs(fill='', vjust = 1) +
  
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=11,color = 'black'),
        axis.title=element_text(size=13), legend.text = element_text(size = 11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(size=11, vjust = 3)) + scale_colour_hue(l=40)


#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
unique(umap_out$label)
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)

clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)
clusters_lambda30 = data.frame(cbind(clusters, rep('Lambda30',dim(clusters)[1])))
names(clusters_lambda30) = c('ident','counts','condition')

#clusters_total = clusters_total[clusters_total$condition != 'Lambda30',]
clusters_total = rbind(clusters_total, clusters_lambda30)


g = aggregate(umi_counts~ident,umap_out,FUN = length)
freqs = c(0.5, 0.5)
freqs = g$umi_counts/dim(umap_out)[1]
-sum(freqs * log2(freqs))


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.3,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_lambda_30_scifi9.png'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.4,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
    scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(5,6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(8,12)])
ggsave(paste(figs.out, 'umap_MG1655_lambda30.png'))

### Lambda 90 ####
keep_rows  = which(data2$Label =="Lambda90" | data2$Label =="Lambda90s" |  data2$Label == "Exponentiala"|data2$Label=="3hr_Cef_Rep1a")
data = data2[keep_rows,]
k = colSums(data!=0) 
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
to_filter = which(k > 10)
new_names = names(to_filter)

new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)

pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total_lambda90 = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Lambda90',length(varExplained))))
kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_lambda90))


rm(mat)
rm(sc_precursor)
sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.05,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,logfc.threshold = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.05,logfc.threshold = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.05,logfc.threshold = 0.1)
RidgePlot(sct, features = c('Lambda-A','Lambda-B'), ncol = 2)

library(dplyr)

library(dplyr)

av.exp.sct <- AverageExpression(sct, return.seurat = TRUE)
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
#sct@active.ident = factor(sct@meta.data,labels = sct@meta.data$Label, levels = sct@meta.data$Label)
markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.5)
markers = markers[markers$avg_log2FC > 0.5,]
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 3, wt = avg_log2FC) -> top10

markers = markers %>%
  add_count(cluster) %>%
  group_by(cluster) %>%
  mutate(n = 8) %>%
  #  group_by(rating, .add = TRUE) %>%
  #In old dplyr add = TRUE
  #group_by(rating, add = TRUE) %>%
  sample_n(n, replace = TRUE) %>%
  select(-n)

markers = unique(markers)
#markers = markers[1:8,]
genes = markers$gene
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
g = merge(g, markers)
#g = g[which(g$variable%in% genes),]
g = g[order(g$cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'EC-',replacement = '', x = g$gene)
g$gene = gsub(pattern = 'Lambda-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = g$gene,labels = g$gene)
#g$Cluster = factor(g$Cluster, levels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'),labels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'))
#g$Cluster = factor(g$Cluster, levels = c(3,2,1))
g = g[order(g$cluster),]

library(myriad)
ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.0), oob = scales::squish)+
  labs(fill='', vjust = 1) +
  
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=11,color = 'black'),
        axis.title=element_text(size=13), legend.text = element_text(size = 11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(size=11, vjust = 3)) + scale_colour_hue(l=40)


#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
unique(umap_out$label)
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)

clusters = aggregate(umi_counts~ident,umap_out,FUN = length)
clusters$umi_counts = clusters$umi_counts/sum(clusters$umi_counts)
clusters_lambda90 = data.frame(cbind(clusters, rep('Lambda90',dim(clusters)[1])))
names(clusters_lambda90) = c('ident','counts','condition')

#clusters_total = clusters_total[clusters_total$condition != 'Lambda90',]
clusters_total = rbind(clusters_total, clusters_lambda90)


g = aggregate(umi_counts~ident,umap_out,FUN = length)
freqs = c(0.5, 0.5)
freqs = g$umi_counts/dim(umap_out)[1]
-sum(freqs * log2(freqs))


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.7,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 11), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_cluster_lambda_90_scifi9.png'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.7,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(8,12)])
ggsave(paste(figs.out, 'umap_MG1655_lambda90.png'))
#### Variability analysis ###


#### PCA ####
write.csv(var_total,file = paste(figs.out,'var_total.csv', sep = ''))
var_total = read.csv(paste(figs.out,'var_total.csv', sep = ''))
var_total$varExplained = as.numeric(as.character(var_total$varExplained))
#var_total$label = var_total$V3

var_total$label = gsub('Gent','Gentamycin',var_total$label)
var_total$label = gsub('Cipro','Ciprofloxacin',var_total$label)
var_total$label = gsub('Cef','Cefazolin',var_total$label)
var_total$label = gsub('Nal','Nalidixic acid',var_total$label)
#var_total$label = gsub('Cycloserineserine','Cycloserine',var_total$label)
var_total$label = gsub('Cyclo','Cycloserine',var_total$label)
#var_total$label = gsub('Erythromycinmycin','Erythromycin',var_total$label)
var_total$label = gsub('Erythro','Erythromycin',var_total$label)
#var_total$label = gsub('Tetracyclineracycline','Tetracycline',var_total$label)
var_total$label = gsub('Tet','Tetracycline',var_total$label)
var_total$label <- factor(var_total$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential'))




ggplot(data = var_total, aes (x = V1, y = varExplained, col = label)) +
  geom_line(size = 1, alpha = 0.8) + 
  geom_point(size = 0.8, alpha = 0.8) +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("Variance Explained")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),
        #panel.grid=element_blank(),panel.border = element_blank(),
        #        panel.background = element_blank(), 
        #        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        #        axis.ticks.x=element_blank(),
        #        axis.text.y=element_blank(),
        #        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
ggsave(paste(figs.out, 'pca_mg1655_varexplaine.png'))

var_total2 = var_total[var_total$label!= "Cefazolin",]
var_total2 = var_total2[var_total2$label!= "Cycloserine",]
ggplot(data = var_total2, aes (x = V1, y = varExplained, col = label)) +
  geom_line(size = 1, alpha = 0.8) + 
  geom_point(size = 0.8, alpha = 0.8) +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("Variance Explained")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),
        #panel.grid=element_blank(),panel.border = element_blank(),
        #        panel.background = element_blank(), 
        #        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        #        axis.ticks.x=element_blank(),
        #        axis.text.y=element_blank(),
        #        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,3,8,5,6,12)])
ggsave(paste(figs.out, 'pca_mg1655_varexplaine_no_wall.png'))

#### kurtosi ####

#write.csv(x = kurtosi_total, file = paste(figs.out,'kurtosi_total.csv', sep = ''))
kurtosi_total = read.csv(paste(figs.out,'kurtosi_total.csv', sep = ''))
kurtosi_total$label = kurtosi_total$X3

kurtosi_total$label = gsub('Gent','Gentamycin',kurtosi_total$label)
kurtosi_total$label = gsub('Cipro','Ciprofloxacin',kurtosi_total$label)
kurtosi_total$label = gsub('Cef','Cefazolin',kurtosi_total$label)
kurtosi_total$label = gsub('Nal','Nalidixic acid',kurtosi_total$label)
#kurtosi_total$label = gsub('Cycloserineserine','Cycloserine',kurtosi_total$label)
kurtosi_total$label = gsub('Cyclo','Cycloserine',kurtosi_total$label)
#kurtosi_total$label = gsub('Erythromycinmycin','Erythromycin',kurtosi_total$label)
kurtosi_total$label = gsub('Erythro','Erythromycin',kurtosi_total$label)
#kurtosi_total$label = gsub('Tetracyclineracycline','Tetracycline',kurtosi_total$label)
kurtosi_total$label = gsub('Tet','Tetracycline',kurtosi_total$label)
kurtosi_total$label <- factor(kurtosi_total$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'))



kurtosi_total$kurtosiExplained = as.numeric(as.character(kurtosi_total$X2))
kurtosi_total$X1 = as.numeric(as.character(kurtosi_total$X1))

#kurtosi_total_max = kurtosi_total %>%group_by()



kurtosi_max = kurtosi_total[,c('label','kurtosiExplained')] %>% group_by(label) %>% summarise(max = max(kurtosiExplained, na.rm=TRUE))
kurtosi_sum = kurtosi_total[,c('label','kurtosiExplained')] %>% group_by(label) %>% summarise(sum = sum(kurtosiExplained, na.rm=TRUE))
#kurtosi_total = merge(kurtosi_total,kurtosi_max)
kurtosi_total = merge(kurtosi_total,kurtosi_sum)
kurtosi_total$kurtosi_norm = kurtosi_total$kurtosiExplained/kurtosi_total$max


#kurtosi_total$PC = as.vector(sapply(kurtosi_total$X, substr, start = 1, stop = 5))
#kurtosi_total$X = kurtosi_total$X
kurtosi_top5 = kurtosi_total[which(kurtosi_total$X1 %in% 1:5),]
kurtosi_top5 = kurtosi_top5[,c('label','kurtosiExplained', 'sum')] %>% group_by(label) %>% summarise(mean_kurt = mean(kurtosiExplained), mean_sum = mean(sum))
#kurtosi_top5 = kurtosi_top5[which(kurtosi_total$label == 'Tetracycline' | kurtosi_total$label == 'Erythromycin' | kurtosi_total$label == 'Gentamycin' | kurtosi_total$label == 'Chlor' | kurtosi_total$label == 'Exponential' ),]
kurtosi_rows = which(kurtosi_total$label == 'Tetracycline' )#| kurtosi_total$label == 'Erythromycin' | kurtosi_total$label == 'Gentamycin' | kurtosi_total$label == 'Chlor' | kurtosi_total$label == 'Exponential' )
kurtosi_rows = intersect(which(kurtosi_total$X1 %in% 1:5), which(kurtosi_total$label == 'Erythromycin'))
kurtosi_pick = kurtosi_total[kurtosi_rows,]
#mechanisms

kurtosi_top5['mechanism'] = c('30S Ribosome', '30S Ribosome',"Cell Wall", "Cell Wall",'50S Ribosome', '50S Ribosome', 'DNA Damage','DNA Damage','None','Phage','Phage')
kurtosi_top5['cidality'] = c('Cidal', 'Static',"Cidal", "Cidal",'Cidal', 'Static', 'Cidal','Cidal','None','Phage','Phage')
kurtosi_top5$mean_kurt = as.numeric(kurtosi_top5$mean_kurt)
kurtosi_top5$mechanism = factor(kurtosi_top5$mechanism, levels = c('30S Ribosome', "Cell Wall", '50S Ribosome', 'DNA Damage','None','Phage'))
drug_mapper = kurtosi_top5[c('label','mechanism','cidality')]
kurtosi_top5v2 = kurtosi_total[which(kurtosi_total$X1 %in% 1:3),]
kurtosi_top5v2 = merge(kurtosi_top5v2,drug_mapper)

kurtosi_total = merge(kurtosi_total, drug_mapper)

ggplot(data = kurtosi_total[which(kurtosi_total$mechanism %in% c('30S Ribosome','50S Ribosome')),], aes (x = X1, y = kurtosiExplained, col = label)) +
  geom_line(size = 2, alpha = 0.8, aes(linetype = cidality)) + 
  # geom_point(size = 0.8, alpha = 0.8) +
  labs(x = "Rank", y = "kurtosi") +
  ggtitle(sprintf("kurtosi for ribosome antibiotics")) +
  scale_y_log10()+
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=0.01,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12)])
ggsave(paste(figs.out, 'kurtosi_MG1655_conditions.png'))

#library(ggpattern)
ggplot(data = kurtosi_top5, aes (x = mechanism, y = mean_sum, fill = label)) +
  geom_bar(stat='identity',position='dodge') + 
#  geom_line(size = 1, alpha = 0.8) + 
#  geom_bar_pattern(aes(fill = label),stat='identity', position=position_dodge()) + 
#  geom_bar_pattern(aes(fill = label),stat='identity', position=position_dodge()) + 
  #scale_pattern_manual(values = c(Nerd = "stripe", NotNerd = "none")) +
#   geom_point(size = 0.8, alpha = 0.8) +
  labs(x = "Rank", y = "Kurtosis") +
#  scale_y_log10() + 
  ggtitle(sprintf("kurtosi")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,12,5,6)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,5,6,12)])
ggsave(paste(figs.out,'total_kurtosis','.png'))

kurtosi_top5v2$label = factor(gsub(x = kurtosi_top5v2$label, pattern = 'Chlor','Chloramphenicol'))
kurtosi_top5v2$label = factor(kurtosi_top5v2$label,levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'))
#kurtosi_top5v2$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential'))



ggplot(data = kurtosi_top5v2[which(kurtosi_top5v2$mechanism %in% c('30S Ribosome','50S Ribosome')),], aes (x = mechanism, y = kurtosiExplained, fill = label)) +
  #  geom_line(size = 1, alpha = 0.8) + 
  geom_boxplot(aes(color = cidality), lwd = 1)+
#  geom_bar(aes(fill = label),stat='identity', position=position_dodge()) + 
  #   geom_point(size = 0.8, alpha = 0.8) +
  labs(x = "Mechanism", y = "Kurtosis") +
    scale_y_log10() + 
#  ggtitle(sprintf("Kurtosis of Top3 Ranked PCs")) +
  ggtitle(sprintf("")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold.italic"),axis.text = element_text(size = 0.1),
        axis.title=element_text(size=0.1,face="bold"), legend.text = element_text(size = 12), legend.title = element_text(size=0.1,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=2)))+
#  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,12,5,6)])
  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12,5,6)])+
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(12,5,12,5,6)])
ggsave(paste(figs.out,'kurtosis_ribosomes','.png'))
ggplot(data = kurtosi_top5v2[which(kurtosi_top5v2$mechanism %in% c('30S Ribosome','50S Ribosome')),], aes (x = mechanism, y = kurtosiExplained, fill = label)) +
  #  geom_line(size = 1, alpha = 0.8) + 
  geom_boxplot(aes(color = cidality), lwd = 1)+
#  geom_bar(aes(fill = label),stat='identity', position=position_dodge()) + 
  #   geom_point(size = 0.8, alpha = 0.8) +
  labs(x = "Mechanism", y = "Kurtosis") +
    scale_y_log10() + 
#  ggtitle(sprintf("Kurtosis of Top3 Ranked PCs")) +
  ggtitle(sprintf("")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 10, size = 10, face="bold.italic"),axis.text = element_text(size = 10),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 12), legend.title = element_text(size=0.1,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=2)))+
#  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,12,5,6)])
  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12,5,6)])+
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(12,5,12,5,6)])
ggsave(paste(figs.out,'kurtosis_ribosomes_with_words','.png'))

ggplot(data = kurtosi_pick, aes (x = X1, y = kurtosiExplained, col = label)) +
  geom_line(size = 1, alpha = 0.8) + 
  # geom_point(size = 0.8, alpha = 0.8) +
  labs(x = "Rank", y = "Kurtosis") +
  scale_y_log10() + 
  ggtitle(sprintf("kurtosi")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,5,6,12)])
#ggsave(paste(figs.out,'pca_MG1655_var_explained','.png'))

ggplot(data = kurtosi_pick, aes (x = X1, y = kurtosiExplained, col = label)) +
  geom_line(size = 1, alpha = 0.8) + 
  # geom_point(size = 0.8, alpha = 0.8) +
  labs(x = "Rank", y = "Kurtosis") +
  scale_y_log10() + 
  ggtitle(sprintf("kurtosi")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,5,6,12)])
#ggsave(paste(figs.out,'pca_MG1655_var_explained','.png'))

# plot the `` 
#### clusters total ###

write.csv(x = clusters_total, file = paste(figs.out,'clusters_total.csv', sep = ''))
#### Clusters ###

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  # ggtitle(sprintf("E. coli")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
#    scale_color_simpsons()
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe')))

ggsave(paste(figs.out, 'umap_MG1655_clean_palette2_30_pcs_clusters.png'))


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  #  ggtitle(sprintf("E. coli")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
  #  scale_color_simpsons()
  #  scale_color_manual(values=as.vector(palette.colors()[c(2,9,4,5,6,8)]))
  scale_color_manual(values=as.vector(palette.colors()))
ggsave(paste(figs.out, 'umap_MG1655_clean_palette1_30_pcs_clusters.png'))


#Histogram of gadA
gads = data.frame('scRNA-Seq',data$EC_gadA)
write.csv(x=gads,file = '/Volumes/AdamsonLab/brucewang/data/counts/gada_count.csv')
#### all lambdas ####
unique(data_total$Label)
data2 = data_total
#keep_rows  = which( data2$Label == "Lambda90"|data2$Label == "Lambda30" | data2$Label == 'Easxponential')
keep_rows  = which( data2$Label == "Lambda90"|data2$Label == "Lambda30" | data2$Label == 'Exponential')
data = data2[keep_rows,]
k = colSums(data[-which(names(data) %in% c('Label','Identity'))])
#k = colSums(data[which(names(data) %in% ec_names)])
to_filter = which(k > 10)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
data = data[new_names]
#k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
k2 = rowSums(data[,which(names(data) %in% ec_names)])
to_filter = which(k2 >= 10)
data = data[to_filter,]
#sc_precursor = data[,-which(names(data) %in% c('Label','Identity','EC_ECU_04345'))]
#sc_precursor = data[,-which(names(data) %in% c('Label','Identity',lambda_names))]
sc_precursor = data[,-which(names(data) %in% c('Label','Identity','EC_ylcI','EC_rzoR','EC_ynfO','EC_nohD'))]

rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
#sct = NormalizeData(sc, scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sct), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sct)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)

sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
sct <- RunUMAP(sct, dims = 1:40, verbose = TRUE)

sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution = 0.5)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'ident') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#umap
umap_out_all = data.frame(sct[['umap']]@cell.embeddings)
umap_out_all$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out_all$label = as.character(sct[['Label']][,1])
umap_out_all$ident = as.character(sct@active.ident)
umap_out_all$label = gsub('_',' ',umap_out$label)
library(ggsci)
aggregate(umi_counts~ident,umap_out,mean)
sizes = aggregate(umi_counts~ident,umap_out,length)
total = sum(aggregate(umi_counts~ident,umap_out,length)['umi_counts'])
sizes['umi_counts']/total
ggplot(umap_out_all,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.2,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(5,6)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(8,12)])
#ggsave(paste(figs.out, 'umap_MG1655_lambda.png'))
#ggsave(file=paste(figs.out, 'umap_MG1655_lambda_lognorm.svg'), dpi=300)
#ggsave(paste(figs.out, 'umap_MG1655_lambda_lognorm.png'))


#levs = as.character(0:18)
umap_out_all$ident = factor(umap_out_all$ident)

#umap_out$ident = factor(as.numeric(umap_out$ident)+1, levels =levs)

ggplot(umap_out_all,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.05,alpha = 1.5) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  # ggtitle(sprintf("E. coli")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
  scale_color_manual(values=as.vector(palette.colors())[c(9,3,7,5,6,8)])#+ guides(colour = guide_legend(override.aes = list(size=1)))

#  scale_color_manual(values=as.vector(palette.colors())[c(1,2,4,5,6,8)])#+ guides(colour = guide_legend(override.aes = list(size=1)))


#  scale_color_simpsons()s
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,5)])
#ggsave(paste(figs.out, 'umap_MG1655_lambda_clusters.png'))
#ggsave(paste(figs.out, 'umap_MG1655_lambda_clusters_lognorm.png'))
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2_30_pcs.png'))

idents = umap_out_all$ident
#### Kbet on both ###


library(kBET)
lytic = as.vector(idents)
lytic[lytic != "3"] = "0"
library(kBET)
pca_out_all = data.frame(sct[['pca']]@cell.embeddings)
pca.data = list()
pca.data$x = pca_out
batch.silhouette_both <- batch_sil(pca.data, lytic,nPCs = 50)
#batch.estimate = batch.estimate = kBET(df = pca_out,lytic,do.pca = FALSE)
#plot.data = data.frame(class=rep(c('observed', 'expected'), 
#                                 each=length(batch.estimate$stats$kBET.observed)), 
#                       data =  c(batch.estimate$stats$kBET.observed,
#                                 batch.estimate$stats$kBET.expected))
#g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
#  labs(x='Test', y='Rejection rate',title='kBET test results') +
#  theme_bw() +  
#  scale_y_continuous(limits=c(0,1))
#
#write.csv(x = plot.data, file = paste(figs.out, 'both_kbet.csv',sep = ''))
#### Color by Lambda ####


# Create a weighted plot.. 
sct2 = sct
phage_score = matrix(data = 0,ncol=2,nrow = dim(sct2@assays$RNA)[1])
phage_score[,1] = as.numeric(str_detect(rownames(sct2@assays$RNA),'Lambda'))
#phage_score[,1] = as.numeric(str_detect(rownames(sct2@assays$RNA),'Nin'))
phage_score[,2] = as.numeric(str_detect(rownames(sct2@assays$RNA),'EC-'))
#phage_score = as.numeric(phage_score)
#norm_data = t(as.data.frame(sct2@assays$RNA@data))
all_data = t(as.data.frame(GetAssayData(sct2, 'counts')))
norm_data = all_data/rowSums(all_data)
phage_data = data.frame(t(as.matrix(sct2@assays$RNA@counts)) %*% phage_score)
phage_data3 = data.frame(norm_data %*% phage_score)*100
phage_data3$ident = idents

quantile(phage_data3[phage_data3$ident == 2,]$X1,probs = seq(0,0.1,0.01))
quantile(phage_data3[phage_data3$ident == 1,]$X1,probs = seq(0,0.1,0.01))
#phage_data3 = data.frame(norm_data %*% phage_score)




library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
#umap_out$label = gsub('3hr Cipro Rep1','90 min Cipro',umap_out$label)

umap_out2 = cbind(umap_out, phage_data3$X1,phage_data3$X2)
umap_out2$max_val = pmax(phage_data3$X1,phage_data3$X2)


library(dplyr)
ggplot(umap_out2%>%arrange(`phage_data3$X1`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `phage_data3$X1`, fill = `phage_data3$X1`)) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point( size = 0.4,alpha = 0.8) +  
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Blues"))[9]) + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Blues"))[9]) + 
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) #+guides(colour = guide_legend(override.aes = list(size=4)))
#ggsave(paste(figs.out,"Lambda_percentage_umap",'.png'))
#ggsave(paste(figs.out,"Lambda_percentage_umap_lognorm",'.png'))


### Heatmap ###

library(dplyr)

av.exp.sct <- AverageExpression(sct, return.seurat = TRUE)
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
#sct@active.ident = factor(sct@meta.data$Label,labels = sct@meta.data$Label, levels = sct@meta.data$Label)
genes = rownames(sct@assays$RNA)[str_detect(rownames(sct@assays$RNA),'Lambda-')]
genes = sort(genes)

#genes = markers$gene
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
#g = merge(g)
g = g[which(g$gene%in% genes),]
g = g[order(g$Cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'EC-',replacement = '', x = g$gene)
g$gene = gsub(pattern = 'Lambda-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = sort(g$gene),labels = sort(g$gene))

#g$Cluster = factor(g$Cluster, levels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'),labels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'))
g$Cluster = factor(as.numeric(g$Cluster)+1, levels = c(5,4,3,2,1))
#g$Cluster = factor(g$Cluster, levels = c(5,4,3,2,1))

g = g[order(g$Cluster),]

library(myriad)
ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.0), oob = scales::squish)+
  labs(fill='', vjust = 1) +
  
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=11,color = 'black'),
        axis.title=element_text(size=13), legend.text = element_text(size = 11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(size=11, vjust = 3)) + scale_colour_hue(l=40)

#ggsave(paste(figs.out,"Lambda_heatmap",'.png'))
ggsave(paste(figs.out,"Lambda_heatmap_test",'.svg'),  device = "svg")

ggsave(paste(figs.out,"Lambda_heatmap",'.png'))



ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.0), oob = scales::squish)+
  labs(fill='', vjust = 1) +
  
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=0.05,color = 'black'),
        axis.title=element_text(size=0.05), legend.text = element_text(size = 0.05), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.05),
        legend.title = element_text(size=0.05, vjust = 3)) + scale_colour_hue(l=40)


#ggsave(paste(figs.out,"Lambda_heatmap_nowords",'.png'))

#ggsave(paste(figs.out,"Lambda_heatmap",'.png'))

#### Coli genes only ####
unique(data_total$Label)
data2 = data_total
keep_rows  = which( data2$Label == "Lambda90"|data2$Label == "Lambda30" | data2$Label == 'Exponential')
data = data2[keep_rows,]
m = rowSums(data2[-which(names(data2) %in% c('Label','Identity'))])
counts = data.frame(log(m, base = 10),data2$Label)
names(counts) = c('count','annotation')
ggplot(data = counts, aes(x = `count`, fill = `annotation`,color = `annotation`)) + 
  geom_histogram(alpha = 0.3 ,position = 'identity') + 
  #  geom_point(aes(color = factor(Label) ))+
  #  labs(x = "Total mRNA", y = "rRNA", color = "Label") +
  #  ggtitle(sprintf("Size Scaling of all cells")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 


# Filter out the operons with less than 10 obs total
k = colSums(data[-which(names(data) %in% c('Label','Identity'))])
to_filter = which(k > 10)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
k2 = rowSums(data[,which(names(data) %in% ec_names)])

to_filter = which(k2 >= 10)
data = data[to_filter,]
#data = data[c(intersect(new_names,lambda_names),'Label','Identity')]
## Make Seurat
#aggregate(EC_lpp~Label,data2,sum)


#sc_precursor = data[,-which(names(data) %in% c('Label','Identity','EC_ECU_04345'))]
sc_precursor = data[,-which(names(data) %in% c('Label','Identity',lambda_names,'EC_ylcI','EC_rzoR','EC_ynfO','EC_nohD'))]
#sc_precursor = data[,-which(names(data) %in% c('Label','Identity'))]

rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
#sct = NormalizeData(sc, scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sct), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sct)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)

sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
#FindMarkers(sct, ident.1 = 'Exponential')

print(sct[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_MG1655_scifi7','.png'))
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/PCA_Relative_Count.png')

sct <- RunUMAP(sct, dims = 1:40, verbose = TRUE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution = 0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#### Project the new idents ####
levels(sct@active.ident) = c(levels(sct@active.ident),"4",'5','6')
idents3 = rep('0', length(idents))
idents3[which(phage_data$X1 > 10)] = '1'
#idents3[which(phage_data$X1 > 20)] = '2'

names(idents3) = 1:length(idents3)
names(idents) = 1:length(idents)
sct@active.ident = factor(idents)
#sct@active.ident = factor(idents3)
#cluster0.markers <- FindMarkers(sct, ident.1 = "0", min.pct = 0.1,logfc.threshold = 0.1)


markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.05)
#markers = markers[markers$avg_log2FC > 0.5,]
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 3, wt = avg_log2FC) -> top10

markers = markers %>%
  add_count(cluster) %>%
  group_by(cluster) %>%
  mutate(n = 10) %>%
  #  group_by(rating, .add = TRUE) %>%
  #In old dplyr add = TRUE
  #group_by(rating, add = TRUE) %>%
  sample_n(n, replace = TRUE) %>%
  select(-n)

markers = unique(markers)
#markers = markers[1:8,]
genes = markers$gene

#### Heatmap ####
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 



library(dplyr)

av.exp.sct <- AverageExpression(sct, return.seurat = TRUE)
#av.exp <- AverageExpression(sct, return.seurat = FALSE)

#av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
#sct@active.ident = factor(sct@meta.data$Label,labels = sct@meta.data$Label, levels = sct@meta.data$Label)
#genes = rownames(sct@assays$RNA)[str_detect(rownames(sct@assays$RNA),'Lambda-')]
#genes = sort(genes)

#genes = markers$gene
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
#g = merge(g)
g = g[which(g$gene%in% genes),]
g = g[order(g$Cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'EC-',replacement = '', x = g$gene)
g$gene = gsub(pattern = 'Lambda-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = sort(g$gene),labels = sort(g$gene))

#g$Cluster = factor(g$Cluster, levels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'),labels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'))
g$Cluster = factor(as.numeric(g$Cluster)+1, levels = c(6,5,4,3,2,1))
#g$Cluster = factor(g$Cluster, levels = c(5,4,3,2,1))

g = g[order(g$Cluster),]

library(myriad)
ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.0), oob = scales::squish)+
  labs(fill='', vjust = 1) +
  
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=11,color = 'black'),
        axis.title=element_text(size=13), legend.text = element_text(size = 11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(size=11, vjust = 3)) + scale_colour_hue(l=40)

#ggsave(paste(figs.out,"coli_clusters_heatmap",'.png'))


cluster0.markers <- FindMarkers(sct, ident.1 = 0,min.pct = 0.1,logfc.threshold = 0.1)
#cluster1.markers <- FindMarkers(sct, ident.1 = 1, ident.2 = 0,min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, ident.2 = 0,min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2,min.pct = 0.05,logfc.threshold = 0.05)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.05,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,logfc.threshold = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.05,logfc.threshold = 0.1)


ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.0), oob = scales::squish)+
  labs(fill='', vjust = 1) +
  
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=0.05,color = 'black'),
        axis.title=element_text(size=0.05), legend.text = element_text(size = 0.05), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.05),
        legend.title = element_text(size=0.05, vjust = 3)) + scale_colour_hue(l=40)

#ggsave(paste(figs.out,"coli_clusters_heatmap_no_words",'.png'))

#### Go Term ####

library(topGO)
ontology = 'BP'
# install.packages("org.Mm.eg.db") 
# Take from the markers
#cluster2.marker_subset = cluster2.markers[which(cluster2.markers$avg_log2FC > 0),]
#rownames(cluster2.marker_subset) = gsub(pattern = 'EC-',replacement = '',x=rownames(cluster2.marker_subset))
use_markers = cluster3.markers
expressed.genes <- rownames(use_markers[intersect(which(use_markers$avg_log2FC > 0), 
                                                  which(use_markers$p_val_adj< 0.05)),])
down.genes <- rownames(use_markers[intersect(which(use_markers$avg_log2FC < 0), 
                                             which(use_markers$p_val_adj< 0.05)),])
expressed.genes = gsub(pattern = 'EC-',replacement = '',x = expressed.genes)
#expressed.genes = expressed.genes[which(expressed.genes != "ompC")]
#expressed.genes = expressed.genes[which(expressed.genes != "ompA")]
down.genes = gsub(pattern = 'EC-',replacement = '',x = down.genes)
all.genes = gsub(pattern = 'EC-',replacement = '', x = rownames(sct@assays$RNA))

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
geneList2 <- ifelse(all.genes %in% down.genes, 1, 0)
names(geneList) <- all.genes
names(geneList2) <- all.genes
GOdata <- new("topGOdata",
              ontology = ontology, # use biological process ontology
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.EcK12.eg.db", ID = "symbol")
#results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher", cutOff = 0.05)
goEnrichment = GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
#goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment$Fisher <- as.numeric(goEnrichment$Fisher)
goEnrichment <- goEnrichment[goEnrichment$Fisher<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","Fisher")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))


GOdata <- new("topGOdata",
              ontology = ontology, # use biological process ontology
              allGenes = geneList2,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.EcK12.eg.db", ID = "symbol")
#results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
resultFisher2 <- runTest(GOdata, algorithm = "elim", statistic = "fisher",cutOff = 0.05)
goEnrichment2 = GenTable(GOdata, Fisher = resultFisher2, topNodes = 20, numChar = 60)
#goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment2$Fisher <- as.numeric(goEnrichment2$Fisher)
goEnrichment2 <- goEnrichment2[goEnrichment2$Fisher<0.05,]
goEnrichment2 <- goEnrichment2[,c("GO.ID","Term","Fisher")]
goEnrichment2$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment2$Term)
goEnrichment2$Term <- gsub("\\.\\.\\.$", "", goEnrichment2$Term)
goEnrichment2$Term <- paste(goEnrichment2$GO.ID, goEnrichment2$Term, sep=", ")
goEnrichment2$Term <- factor(goEnrichment2$Term, levels=rev(goEnrichment2$Term))


goEnrichment$Fish = -log10(goEnrichment$Fisher)
goEnrichment2$Fish = log10(goEnrichment2$Fisher)
goEnrichment = goEnrichment[1:8,]
goEnrichment2 = goEnrichment2[1:8,]
goEnrichment = rbind(goEnrichment, goEnrichment2)

goEnrichment$Term = gsub(pattern = 'GO:[0-9]+, ', replacement = '',x = goEnrichment$Term)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
library(stringr)
goEnrichment$Term = firstup(goEnrichment$Term)
goEnrichment = goEnrichment[order(goEnrichment$Fish),]
#goEnrichment$Term <- factor(goEnrichment$Term, levels=(goEnrichment$Term))
#goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
#goEnrichment$Fish <- factor(goEnrichment$Fish, levels=unique(sort(goEnrichment$Fish)))
ggplot(goEnrichment, aes(x=reorder(Term,Fish), y=Fish, fill = `Fish`)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p value)") +
  ggtitle("") +  
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
#  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Fisher)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    ##    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=8, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10,  hjust=1.10, color = 'black'),
    axis.text.y=element_text(angle=0, size=10,  vjust=0.5, color = 'black'),
    axis.title=element_text(size=0.01, face="bold"),
    ##    legend.key=element_blank(),     #removes the border
    ##    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    #    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

ggsave(paste(figs.out,"go_term_lytic_coli",'.svg'),  device = "svg")
write.csv(x = cluster1.markers, file = paste(figs.out,"genes_lytic_genes",'.csv',sep = ''))



#### Kbet ####

library(kBET)
lytic = as.vector(idents)
lytic[lytic != "3"] = "0"
library(kBET)
pca_out_coli = data.frame(sct[['pca']]@cell.embeddings)
#batch.estimate = batch.estimate = kBET(df = pca_out,lytic,do.pca = FALSE)
# get a list
pca.data = list()
pca.data$x = pca_out
batch.silhouette_coli <- batch_sil(pca.data, lytic,nPCs = 50)
plot.data = data.frame(class=rep(c('observed', 'expected'), 
                                 each=length(batch.estimate$stats$kBET.observed)), 
                       data =  c(batch.estimate$stats$kBET.observed,
                                 batch.estimate$stats$kBET.expected))
g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='Test', y='Rejection rate',title='kBET test results') +
  theme_bw() +  
  scale_y_continuous(limits=c(0,1))

write.csv(x = plot.data, file = paste(figs.out, 'coli_kbet.csv',sep = ''))
#### PCA ####
umap_out_coli = data.frame(sct[['umap']]@cell.embeddings)
umap_out_coli$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out_coli$label = as.character(sct[['Label']][,1])
umap_out_coli$ident = as.character(sct@active.ident)
umap_out_coli$ident = idents
umap_out_coli$label = gsub('_',' ',umap_out$label)
library(ggsci)
aggregate(umi_counts~ident,umap_out,mean)
sizes = aggregate(umi_counts~ident,umap_out,length)
total = sum(aggregate(umi_counts~ident,umap_out,length)['umi_counts'])
sizes['umi_counts']/total

#umap_out = umap_out[umap_out$label != 'Overnight',]
#umap_out = umap_out[umap_out$label != 'Lambda90',]
#umap_out = umap_out[umap_out$label != 'Lambda30',]


ggplot(umap_out_coli,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.2,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,5,6)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(8,12)])
ggsave(paste(figs.out, 'umap_MG1655_lambda_no_lambda.png'))
#ggsave(paste(figs.out, 'umap_MG1655_lambda_no_lambda_lognormalize.png'))


umap_out_coli$ident = factor(umap_out_coli$ident, levels = levs)

umap_out_coli$ident = factor(as.numeric(umap_out_coli$ident)+1, levels =levs)

ggplot(umap_out_coli,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.2,alpha = 0.7) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  # ggtitle(sprintf("E. coli")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
  scale_color_manual(values=as.vector(palette.colors())[c(9,3,7,5,6,8)])#+ guides(colour = guide_legend(override.aes = list(size=1)))


#  scale_color_simpsons()s
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,5)])
ggsave(paste(figs.out, 'umap_MG1655_no_lambda.png'))
#ggsave(paste(figs.out, 'umap_MG1655_lambda_clusters_lognorm_coli_only.png'))
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2_30_pcs.png'))

# Create a weighted plot.. 
umap_out2 = cbind(umap_out_coli, phage_data3$X1,phage_data3$X2)
umap_out2$max_val = pmax(phage_data3$X1,phage_data3$X2)


ggplot(umap_out2%>%arrange(`phage_data3$X1`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `phage_data3$X1`, fill = `phage_data3$X1`)) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point( size = 0.4,alpha = 0.7) +  
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Blues"))[9]) + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Blues"))[9]) + 
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) #+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(paste(figs.out,"Lambda_percentage_umap_no_lambda_lognormalize",'.png'))

#### Lambda genes only ####
unique(data_total$Label)
data2 = data_total
keep_rows  = which( data2$Label == "Lambda90"|data2$Label == "Lambda30" | data2$Label == 'Expaonential')
data = data2[keep_rows,]
m = rowSums(data2[-which(names(data2) %in% c('Label','Identity'))])
counts = data.frame(log(m, base = 10),data2$Label)
names(counts) = c('count','annotation')
ggplot(data = counts, aes(x = `count`, fill = `annotation`,color = `annotation`)) + 
  geom_histogram(alpha = 0.3 ,position = 'identity') + 
  #  geom_point(aes(color = factor(Label) ))+
  #  labs(x = "Total mRNA", y = "rRNA", color = "Label") +
  #  ggtitle(sprintf("Size Scaling of all cells")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 


# Filter out the operons with less than 10 obs total
k = colSums(data[-which(names(data) %in% c('Label','Identity'))])
to_filter = which(k > 10)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
k2 = rowSums(data[,which(names(data) %in% ec_names)])

to_filter = which(k2 >= 10)
data = data[to_filter,]
#data = data[c(intersect(new_names,lambda_names),'Label','Identity')]
## Make Seurat
#aggregate(EC_lpp~Label,data2,sum)


#sc_precursor = data[,-which(names(data) %in% c('Label','Identity','EC_ECU_04345'))]
sc_precursor = data[,-which(names(data) %in% c('Label','Identity',ec_names))]
#sc_precursor = data[,-which(names(data) %in% c('Label','Identity'))]

rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
#sct = NormalizeData(sc, scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sct), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sct)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)

sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
#FindMarkers(sct, ident.1 = 'Exponential')

print(sct[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_MG1655_scifi7','.png'))
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/PCA_Relative_Count.png')

sct <- RunUMAP(sct, dims = 1:5, verbose = TRUE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:5, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution = 0.2)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#### PCA ####
umap_out_lambda = data.frame(sct[['umap']]@cell.embeddings)
umap_out_lambda$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out_lambda$label = as.character(sct[['Label']][,1])
umap_out_lambda$ident = as.character(sct@active.ident)
umap_out_lambda$ident = idents
umap_out_lambda$label = gsub('_',' ',umap_out$label)
library(ggsci)
aggregate(umi_counts~ident,umap_out,mean)
sizes = aggregate(umi_counts~ident,umap_out,length)
total = sum(aggregate(umi_counts~ident,umap_out,length)['umi_counts'])
sizes['umi_counts']/total

#umap_out = umap_out[umap_out$label != 'Overnight',]
#umap_out = umap_out[umap_out$label != 'Lambda90',]
#umap_out = umap_out[umap_out$label != 'Lambda30',]


ggplot(umap_out_lambda,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.2,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(5,6)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(8,12)])
#ggsave(paste(figs.out, 'umap_MG1655_lambda_no_coli.png'))
ggsave(paste(figs.out, 'umap_MG1655_lambda_no_coli_lognormalize.png'))


levs = as.character(0:18)
#umap_out$ident = factor(umap_out$ident, levels = levs)

umap_out$ident = factor(idents, levels =levs)

ggplot(umap_out_lambda,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.05,alpha = 1.5) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  # ggtitle(sprintf("E. coli")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
  scale_color_manual(values=as.vector(palette.colors())[c(9,3,7,5,6,8)])#+ guides(colour = guide_legend(override.aes = list(size=1)))

  #scale_color_manual(values=as.vector(palette.colors())[c(1,2,4,5,6,8,9,10,11,3,7)])#+ guides(colour = guide_legend(override.aes = list(size=1)))


#  scale_color_simpsons()s
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,5)])
ggsave(paste(figs.out, 'umap_MG1655_lambda_clusters_no_coli_lambda.png'))
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2_30_pcs.png'))
#### Silhouette lambda ####

library(kBET)
lytic = as.vector(idents)
lytic[lytic != "3"] = "0"
random_size = length(which(lytic == 3))
random_labels = rep("0", length(lytic))
random_indices = sample(x = 1:length(lytic), size = random_size, replace = FALSE)
random_labels[random_indices] = "3"
#random_labels[which(idents == "4")] = "3"
library(kBET)
pca_out_lambda = data.frame(sct[['pca']]@cell.embeddings)
#batch.estimate = batch.estimate = kBET(df = pca_out,lytic,do.pca = FALSE)
# get a list
ggplots = list()
idents3 = rep('0', length(idents))
idents3[which(phage_data$X1 > 20.9)] = '1'

for (cluster in unique(idents)){
#for (cluster in c('2')){
  print(as.numeric(cluster))
  reps = 20
  batch.silhouette_lambda_null = c()
  batch.silhouette_coli_null = c()
  batch.silhouette_both_null = c()
  batch.silhouette_lambda_lytic = c()
  batch.silhouette_coli_lytic = c()
  batch.silhouette_both_lytic = c()
  data_size = dim(pca_out_all)[1]/2
  idents2 = idents
  lytic = as.vector(idents2)
  random_size = length(which(lytic == cluster))
  random_labels = rep("9", length(lytic))
  random_indices = sample(x = 1:length(lytic), size = random_size, replace = FALSE)
  random_labels[random_indices] = as.character(cluster)

  for (i in 1:reps){
    rand_sample = sample(x = 1:dim(pca_out_all)[1], size = data_size, replace = FALSE)
    print(rand_sample[1:5])
    idents2 = idents[rand_sample]
    lytic = as.vector(idents2)
    #  lytic[lytic != "3"] = "0"
    lytic[lytic != as.character(cluster)] = "9"
    #  random_size = length(which(lytic == 3))
#    random_size = length(which(lytic == cluster))
#    random_labels = rep("9", length(lytic))
#    random_indices = sample(x = 1:length(lytic), size = random_size, replace = FALSE)
#    #  random_labels[random_indices] = as.characte
#    random_labels[random_indices] = as.character(cluster)
    random_labels2 = random_labels[rand_sample]
    print(random_indices[1:5])
    pca.data_lambda = list()
    pca.data_lambda$x = pca_out_lambda[rand_sample,]
    pca.data_coli = list()
    pca.data_coli$x = pca_out_coli[rand_sample,]
    pca.data_both = list()
    pca.data_both$x = pca_out_all[rand_sample,]
    #  random_labels = rep("9", length(lytic))
    #  random_indices = sample(x = 1:length(lytic), size = random_size, replace = FALSE)
    random_labels[random_indices] = as.character(cluster)
    batch.silhouette_lambda_null = c(batch.silhouette_lambda_null,batch_sil(pca.data_lambda, random_labels2,nPCs = 5) )
    batch.silhouette_lambda_lytic = c(batch.silhouette_lambda_lytic,batch_sil(pca.data_lambda, lytic,nPCs = 5) )
    batch.silhouette_coli_null = c(batch.silhouette_coli_null,batch_sil(pca.data_coli, random_labels2,nPCs = 40) )
    batch.silhouette_coli_lytic = c(batch.silhouette_coli_lytic,batch_sil(pca.data_coli, lytic,nPCs = 40) )
    batch.silhouette_both_null = c(batch.silhouette_both_null,batch_sil(pca.data_both, random_labels2,nPCs = 40) )
    batch.silhouette_both_lytic = c(batch.silhouette_both_lytic,batch_sil(pca.data_both, lytic,nPCs = 40) )
  }
  
  
  #batch.silhouette_lambda_lytic <- batch_sil(pca.data_lambda, lytic,nPCs = 5)
  ##batch.silhouette_lambda_null <- batch_sil(pca.data_lambda, random_labels,nPCs = 5)
  #batch.silhouette_coli_lytic <- batch_sil(pca.data_coli, lytic,nPCs = 40)
  ##batch.silhouette_coli_null <- batch_sil(pca.data_coli, random_labels,nPCs = 40)
  #batch.silhouette_both_lytic <- batch_sil(pca.data_both, lytic,nPCs = 40)
  #batch.silhouette_both_null <- batch_sil(pca.data_both, random_labels,nPCs = 40)
  
  
  silhouettes = c(batch.silhouette_lambda_lytic,batch.silhouette_lambda_null,
                  batch.silhouette_coli_lytic,batch.silhouette_coli_null,
                  batch.silhouette_both_lytic,batch.silhouette_both_null
  )
  species = c(rep('Lambda',length(batch.silhouette_lambda_lytic)),rep('Lambda',length(batch.silhouette_lambda_null)),
              rep('E. coli',length(batch.silhouette_coli_lytic)), rep('E. coli',length(batch.silhouette_coli_null)),
              rep('Both',length(batch.silhouette_both_lytic)), rep('Both',length(batch.silhouette_both_null))
  )
  cluster_label = c(rep(paste('Cluster ', cluster, sep = ''),length(batch.silhouette_lambda_lytic)),rep('Random',length(batch.silhouette_lambda_null)),
                    rep(paste('Cluster ', cluster, sep = ''),length(batch.silhouette_coli_lytic)),rep('Random',length(batch.silhouette_coli_null)),
                    rep(paste('Cluster ', cluster, sep = ''),length(batch.silhouette_both_null)),rep('Random',length(batch.silhouette_both_null))
  )
  
  to_plot = data.frame(cbind(silhouettes,species,cluster_label))
  to_plot$silhouettes = as.numeric(as.character(silhouettes))
  to_plot$species2 = factor(to_plot$species,levels = c('Both','Lambda','E. coli'), labels = c('Both','Lambda','E. coli'))
  
  t.test(x = batch.silhouette_both_lytic, y = batch.silhouette_coli_lytic, alternative = 'g')
  
  ggplot() + 
    geom_bar(data = to_plot, aes(x = `species2`, y= `silhouettes`, fill = `cluster_label`), position = "dodge", stat="summary", color = 'black')+
    geom_point(data =to_plot, aes(x = `species2`, y= `silhouettes`, fill = `cluster_label`), position =  position_dodge(width = .9))+
    #  geom_errorbar(data = g, aes(x = `Fixation`, y= `CT`, fill = `Sample`,ymin = 0, ymax = 30), position =  position_dodge(width = .9))+
    #  scale_fill_brewer(palette = "Set2") + 
    scale_fill_manual(values = as.vector(palette.colors()[c(5,9,10,11,3,7)]))+
    scale_color_manual(values = 'black')+
    #+ guides(colour = guide_legend(override.aes = list(size=1)))
    #  geom_line(data = subset(roc_curve_frame, X3 %like% "Max"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
    # geom_line(stat = "identity",aes(color = factor(X3)), linetype = 'dotted') +
    labs(x = "Genome", y = "Silhouette", color = "Method") +
    #ggtitle(sprintf("5s rRNA Depletion")) +
    #  coord_cartesian(xlim= c(0,1.0))+
    #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
    #  coord_equal(ratio=1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=12),
          axis.title=element_text(size=12,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=12,face="bold")) + scale_colour_hue(l=40)
  ggplots[[cluster]] = to_plot
  ggsave(paste(figs.out, 'silhouette_score_with_words_cluster_',cluster, '.png', sep = ''))
}
m = bind_rows(ggplots, .id = "column_label")
ggplot() + 
  geom_bar(data = m, aes(x = `species2`, y= `silhouettes`, fill = `cluster_label`), position = "dodge", stat="summary", color = 'black')+
  geom_point(data =m, aes(x = `species2`, y= `silhouettes`, fill = `cluster_label`), position =  position_dodge(width = .9))+
  #  geom_errorbar(data = g, aes(x = `Fixation`, y= `CT`, fill = `Sample`,ymin = 0, ymax = 30), position =  position_dodge(width = .9))+
  #  scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values=as.vector(palette.colors())[c(9,3,7,5,6,8)])+#+ guides(colour = guide_legend(override.aes = list(size=1)))

#  scale_fill_manual(values = as.vector(palette.colors()[c(5,3,7,2,9)]))+
  scale_color_manual(values = 'black')+
  #+ guides(colour = guide_legend(override.aes = list(size=1)))
  #  geom_line(data = subset(roc_curve_frame, X3 %like% "Max"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  # geom_line(stat = "identity",aes(color = factor(X3)), linetype = 'dotted') +
  labs(x = "Genome", y = "Silhouette", color = "Method") +
  #ggtitle(sprintf("5s rRNA Depletion")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=12,face="bold")) + scale_colour_hue(l=40)
ggsave(paste(figs.out, 'silhouette_score_with_words_cluster_all', '.png', sep = ''))

cluster2 = m[m$column_label == 2,]
ggplot() + 
  geom_bar(data = cluster2, aes(x = `species2`, y= `silhouettes`, fill = `cluster_label`), position = "dodge", stat="summary", color = 'black')+
  geom_point(data =cluster2, aes(x = `species2`, y= `silhouettes`, fill = `cluster_label`), position =  position_dodge(width = .9))+
  #  geom_errorbar(data = g, aes(x = `Fixation`, y= `CT`, fill = `Sample`,ymin = 0, ymax = 30), position =  position_dodge(width = .9))+
  #  scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values=as.vector(palette.colors())[c(7,6,8)])+#+ guides(colour = guide_legend(override.aes = list(size=1)))
  scale_color_manual(values = 'black')+
  #+ guides(colour = guide_legend(override.aes = list(size=1)))
  #  geom_line(data = subset(roc_curve_frame, X3 %like% "Max"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  # geom_line(stat = "identity",aes(color = factor(X3)), linetype = 'dotted') +
  labs(x = "Genome", y = "Silhouette", color = "Method") +
  #ggtitle(sprintf("5s rRNA Depletion")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=12,face="bold")) + scale_colour_hue(l=40)
ggsave(paste(figs.out, 'silhouette_score_with_words_cluster_cluster2', '.png', sep = ''))
t.test(x = subset(subset(cluster2, cluster_label == "Cluster 2"),species == "E. coli")['silhouettes'],
       y = subset(subset(cluster2, cluster_label == "Random"),species == "E. coli")['silhouettes'], alternative = 'g')

reps = 10
batch.silhouette_lambda_null = c()
batch.silhouette_coli_null = c()
batch.silhouette_both_null = c()
batch.silhouette_lambda_lytic = c()
batch.silhouette_coli_lytic = c()
batch.silhouette_both_lytic = c()
data_size = dim(pca_out_all)[1]/2
cluster = 2

#kBET(df = pca_out_lambda[rand_sample,],random_labels,do.pca = FALSE)
#coli_lytic = kBET(df = pca_out_coli[rand_sample,],lytic,do.pca = FALSE)
#coli_random = kBET(df = pca_out_coli[rand_sample,],random_labels,do.pca = FALSE)
#lambda_lytic = kBET(df = pca_out_lambda[rand_sample,],lytic,do.pca = FALSE)
#lambda_random = kBET(df = pca_out_lambda[rand_sample,],random_labels,do.pca = FALSE)
#all_lytic = kBET(df = pca_out_all[rand_sample,],lytic,do.pca = FALSE)
#all_random = kBET(df = pca_out_all[rand_sample,],random_labels,do.pca = FALSE)
for (i in 1:reps){
  rand_sample = sample(x = 1:dim(pca_out_all)[1], size = data_size, replace = FALSE)
  print(rand_sample[1:5])
  idents2 = idents[rand_sample]
  lytic = as.vector(idents2)
#  lytic[lytic != "3"] = "0"
  lytic[lytic != as.character(cluster)] = "9"
#  random_size = length(which(lytic == 3))
  random_size = length(which(lytic == cluster))
  random_labels = rep("9", length(lytic))
  random_indices = sample(x = 1:length(lytic), size = random_size, replace = FALSE)
#  random_labels[random_indices] = as.characte
  random_labels[random_indices] = as.character(cluster)
  print(random_indices[1:5])
  pca.data_lambda = list()
  pca.data_lambda$x = pca_out_lambda[rand_sample,]
  pca.data_coli = list()
  pca.data_coli$x = pca_out_coli[rand_sample,]
  pca.data_both = list()
  pca.data_both$x = pca_out_all[rand_sample,]
#  random_labels = rep("9", length(lytic))
#  random_indices = sample(x = 1:length(lytic), size = random_size, replace = FALSE)
  random_labels[random_indices] = as.character(cluster)
  batch.silhouette_lambda_null = c(batch.silhouette_lambda_null,batch_sil(pca.data_lambda, random_labels,nPCs = 5) )
  batch.silhouette_lambda_lytic = c(batch.silhouette_lambda_lytic,batch_sil(pca.data_lambda, lytic,nPCs = 5) )
  batch.silhouette_coli_null = c(batch.silhouette_coli_null,batch_sil(pca.data_coli, random_labels,nPCs = 40) )
  batch.silhouette_coli_lytic = c(batch.silhouette_coli_lytic,batch_sil(pca.data_coli, lytic,nPCs = 40) )
  batch.silhouette_both_null = c(batch.silhouette_both_null,batch_sil(pca.data_both, random_labels,nPCs = 40) )
  batch.silhouette_both_lytic = c(batch.silhouette_both_lytic,batch_sil(pca.data_both, lytic,nPCs = 40) )
}


#batch.silhouette_lambda_lytic <- batch_sil(pca.data_lambda, lytic,nPCs = 5)
##batch.silhouette_lambda_null <- batch_sil(pca.data_lambda, random_labels,nPCs = 5)
#batch.silhouette_coli_lytic <- batch_sil(pca.data_coli, lytic,nPCs = 40)
##batch.silhouette_coli_null <- batch_sil(pca.data_coli, random_labels,nPCs = 40)
#batch.silhouette_both_lytic <- batch_sil(pca.data_both, lytic,nPCs = 40)
#batch.silhouette_both_null <- batch_sil(pca.data_both, random_labels,nPCs = 40)


silhouettes = c(batch.silhouette_lambda_lytic,batch.silhouette_lambda_null,
             batch.silhouette_coli_lytic,batch.silhouette_coli_null,
             batch.silhouette_both_lytic,batch.silhouette_both_null
                )
species = c(rep('Lambda',length(batch.silhouette_lambda_lytic)),rep('Lambda',length(batch.silhouette_lambda_null)),
            rep('E. coli',length(batch.silhouette_coli_lytic)), rep('E. coli',length(batch.silhouette_coli_null)),
            rep('Both',length(batch.silhouette_both_lytic)), rep('Both',length(batch.silhouette_both_null))
                )
cluster_label = c(rep(paste('Cluster ', cluster, sep = ''),length(batch.silhouette_lambda_lytic)),rep('Random',length(batch.silhouette_lambda_null)),
                  rep(paste('Cluster ', cluster, sep = ''),length(batch.silhouette_coli_lytic)),rep('Random',length(batch.silhouette_coli_null)),
                  rep(paste('Cluster ', cluster, sep = ''),length(batch.silhouette_both_null)),rep('Random',length(batch.silhouette_both_null))
                  )

to_plot = data.frame(cbind(silhouettes,species,cluster_label))
to_plot$silhouettes = as.numeric(as.character(silhouettes))
to_plot$species2 = factor(to_plot$species,levels = c('Both','Lambda','E. coli'), labels = c('Both','Lambda','E. coli'))

t.test(x = batch.silhouette_both_lytic, y = batch.silhouette_coli_lytic, alternative = 'g')

ggplot() + 
  geom_bar(data = to_plot, aes(x = `species2`, y= `silhouettes`, fill = `cluster_label`), position = "dodge", stat="summary", color = 'black')+
  geom_point(data =to_plot, aes(x = `species2`, y= `silhouettes`, fill = `cluster_label`), position =  position_dodge(width = .9))+
  #  geom_errorbar(data = g, aes(x = `Fixation`, y= `CT`, fill = `Sample`,ymin = 0, ymax = 30), position =  position_dodge(width = .9))+
#  scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values = as.vector(palette.colors()[c(5,9,10,11,3,7)]))+
  scale_color_manual(values = 'black')+
  #+ guides(colour = guide_legend(override.aes = list(size=1)))
  #  geom_line(data = subset(roc_curve_frame, X3 %like% "Max"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  # geom_line(stat = "identity",aes(color = factor(X3)), linetype = 'dotted') +
  labs(x = "Genome", y = "Silhouette", color = "Method") +
  #ggtitle(sprintf("5s rRNA Depletion")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=12,face="bold")) + scale_colour_hue(l=40)
ggsave(paste(figs.out, 'silhouette_score_with_words_cluster_',cluster, '.png', sep = ''))

ggplot() + 
  geom_bar(data = to_plot, aes(x = `species2`, y= `silhouettes`, fill = `cluster_label`), position = "dodge", stat="summary", fun.y = 'mean', color = 'black')+
  #  geom_point(data = g, aes(x = `Sample`, y= `Delta_Ct`, fill = `Sample`), position =  position_dodge(width = .9))+
  #  geom_errorbar(data = g, aes(x = `Fixation`, y= `CT`, fill = `Sample`,ymin = 0, ymax = 30), position =  position_dodge(width = .9))+
  #  scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values = as.vector(palette.colors()[c(5,9,10,11,3,7)]))+
  scale_color_manual(values = 'black')+
  #+ guides(colour = guide_legend(override.aes = list(size=1)))
  #  geom_line(data = subset(roc_curve_frame, X3 %like% "Max"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  # geom_line(stat = "identity",aes(color = factor(X3)), linetype = 'dotted') +
  labs(x = "", y = " ", color = "") +
  #ggtitle(sprintf("5s rRNA Depletion")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 0.01, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=0.01),
        axis.title=element_text(size=0.01,face="bold"), legend.text = element_text(size = 0.01), legend.title = element_text(size=0.01,face="bold")) + scale_colour_hue(l=40)
ggsave(paste(figs.out, 'silhouette_score_no_words.png'))





batch.estimate_all = batch.estimate = kBET(df = pca_out_all,lytic,do.pca = FALSE)
batch.estimate_coli = batch.estimate = kBET(df = pca_out_coli,lytic,do.pca = FALSE)
batch.estimate_lambda = batch.estimate = kBET(df = pca_out_lambda,lytic,do.pca = FALSE)
#batch.estimate_all = batch.estimate = kBET(df = pca_out_all,random_labels,do.pca = FALSE)

#### Color by Lambda ####


# Create a weighted plot.. 
umap_out2 = cbind(umap_out_lambda, phage_data3$X1,phage_data3$X2)
umap_out2$max_val = pmax(phage_data3$X1,phage_data3$X2)


ggplot(umap_out2%>%arrange(`phage_data3$X1`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `phage_data3$X1`, fill = `phage_data3$X1`)) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point( size = 0.4,alpha = 0.8) +  
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Blues"))[9]) + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Blues"))[9]) + 
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) #+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(paste(figs.out,"Lambda_percentage_umap_no_coli_lognormalize",'.png'))


#### Lambda only ####

### Heatmap ###
library(dplyr)

av.exp.sct <- AverageExpression(sct, return.seurat = TRUE)
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
#sct@active.ident = factor(sct@meta.data$Label,labels = sct@meta.data$Label, levels = sct@meta.data$Label)
markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.5)
markers = markers[markers$avg_log2FC > 0.5,]
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 3, wt = avg_log2FC) -> top10

markers = markers %>%
  add_count(cluster) %>%
  group_by(cluster) %>%
  mutate(n = 8) %>%
  #  group_by(rating, .add = TRUE) %>%
  #In old dplyr add = TRUE
  #group_by(rating, add = TRUE) %>%
  sample_n(n, replace = TRUE) %>%
  select(-n)

markers = unique(markers)
#markers = markers[1:8,]
#genes = markers$gene
genes = lambda_names
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
#g = merge(g, markers)
g = g[which(g$variable%in% genes),]
g = g[order(g$cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'EC-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = g$gene,labels = g$gene)
#g$Cluster = factor(g$Cluster, levels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'),labels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'))
g$Cluster = factor(as.numeric(g$Cluster)+1, levels = c(5,4,3,2,1))
#g$Cluster = factor(g$Cluster, levels = c(5,4,3,2,1))

g = g[order(g$cluster),]

library(myriad)
ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.0), oob = scales::squish)+
  labs(fill='', vjust = 1) +
  
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=11,color = 'black'),
        axis.title=element_text(size=13), legend.text = element_text(size = 11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(size=11, vjust = 3)) + scale_colour_hue(l=40)

ggsave(paste(figs.out,"Lambda_heatmap",'.png'))
#### Change the labels ####

lytic = union(which(umap_out$ident == "4"),which(umap_out$ident == "3"))
lysogenic = union(which(umap_out$ident == "0"),union(which(umap_out$ident == "1"),which(umap_out$ident == "2")))
levels(sct@active.ident) = c(levels(sct@active.ident),"Lytic","Lysogenic")
sct@active.ident[lytic] = "1"
sct@active.ident[lysogenic] = "2"
#sct@active.ident[look_into1] = "6"
#sct@active.ident[look_into2] = "7"



#### All lambda genes ####

library(dplyr)

av.exp.sct <- AverageExpression(sct, return.seurat = TRUE)
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
#sct@active.ident = factor(sct@meta.data$Label,labels = sct@meta.data$Label, levels = sct@meta.data$Label)
genes = rownames(sct@assays$RNA)[str_detect(rownames(sct@assays$RNA),'Lambda-')]
genes = sort(genes)

#genes = markers$gene
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
#g = merge(g)
g = g[which(g$gene%in% genes),]
g = g[order(g$Cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'EC-',replacement = '', x = g$gene)
g$gene = gsub(pattern = 'Lambda-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = sort(g$gene),labels = sort(g$gene))

#g$Cluster = factor(g$Cluster, levels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'),labels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'))
g$Cluster = factor(as.numeric(g$Cluster)+1, levels = c(5,4,3,2,1))
#g$Cluster = factor(g$Cluster, levels = c(5,4,3,2,1))

g = g[order(g$Cluster),]

library(myriad)
ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.0), oob = scales::squish)+
  labs(fill='', vjust = 1) +
  
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=11,color = 'black'),
        axis.title=element_text(size=13), legend.text = element_text(size = 11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(size=11, vjust = 3)) + scale_colour_hue(l=40)

ggsave(paste(figs.out,"Lambda_heatmap",'.png'))



ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.0), oob = scales::squish)+
  labs(fill='', vjust = 1) +
  
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=0.05,color = 'black'),
        axis.title=element_text(size=0.05), legend.text = element_text(size = 0.05), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.05),
        legend.title = element_text(size=0.05, vjust = 3)) + scale_colour_hue(l=40)


ggsave(paste(figs.out,"Lambda_heatmap_nowords",'.png'))


umap_out3 = data.frame(sct[['umap']]@cell.embeddings)
umap_out3$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out3$label = as.character(sct[['Label']][,1])
umap_out3$ident = as.character(sct@active.ident)

umap_out3 = cbind(umap_out3, phage_data3$X1,phage_data3$X2)
umap_out2$max_val = pmax(phage_data3$X1,phage_data3$X2)


ggplot(umap_out3%>%arrange(`phage_data3$X1`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `phage_data3$X1`, fill = `phage_data3$X1`)) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5) +
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #  scale_fill_gradientn(colors =  (brewer.pal(n = 9, name =  "BuGn") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "BuGn")), guide = '' ) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) +   #scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greens")[1]), high = (brewer.pal(n = 9, name = "Greens")[9])) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Blues"))[9]) + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Blues"))[9]) + 
  #  scale_colour_gradient(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 
ggsave(paste(figs.out,"removing_lambda_genes",'.png'))


#ggplot(umap_out2%>%arrange(`phage_data3$X2`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `phage_data3$X2`, fill = `phage_data3$X2`)) + 
#  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
#  geom_point(size =1.5) +
#  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
#  labs(y = "", x = "") +
#  ggtitle(sprintf("")) +
#  #  scale_fill_gradientn(colors =  (brewer.pal(n = 9, name =  "BuGn") ), guide = "colourbar") + 
#  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "BuGn")), guide = '' ) + 
#  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) +   #scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greens")[1]), high = (brewer.pal(n = 9, name = "Greens")[9])) + 
#  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Blues"))[9]) + 
#  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Blues"))[9]) + 
#  #  scale_colour_gradient(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
#  #  coord_cartesian(xlim= c(0,1.0))+
#  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
#  #  coord_equal(ratio=1) +
#  theme_bw() +
#  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank() ,
#        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
#  labs(fill = "", colour = "") 



#### Variance vector ####
var_matrix = matrix(ncol = 2)
mean_matrix = matrix(ncol = 2)
var_norm_matrix = matrix(ncol = 2)
mean_norm_matrix = matrix(ncol = 2)
var_pca_matrix = matrix(ncol = 2)
mean_pca_matrix = matrix(ncol = 2)
for (label in unique(data_total$Label)){
  print(label)
    unique(data_total$Label)
    keep_rows  = which(data2$Label ==label)
    data = data2[keep_rows,]
    k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
    to_filter = which(k2 >= 10)
    data = data[to_filter,]
    sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
    rownames(sc_precursor) = 1:dim(sc_precursor)[1]
    sc_precursor = t(sc_precursor)
    sc_meta = data[which(names(data) %in% c('Label'))]
    sc = CreateSeuratObject(
      #  sc_precursor,
      sc_precursor,
      project = "SeuratProject",
      assay = "RNA",
      min.cells = 0,
      min.features = 0,
      names.field = 1,
      names.delim = "_",
      meta.data = sc_meta
    )
    sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
    var_x = apply(sc_precursor,1,kurtosi)
    mean_x = apply(sc_precursor,1,mean)
    var_x_norm = apply(as.matrix(sct@assays$RNA@data),1,kurtosi)
    mean_x_norm = apply(sct@assays$RNA@data,1,kurtosi)
    var_matrix = rbind(var_matrix, cbind(var_x, rep(label,length(var_x))))
    mean_matrix = rbind(mean_matrix, cbind(mean_x, rep(label,length(mean_x))))
    var_norm_matrix = rbind(var_norm_matrix, cbind(var_x_norm, rep(label,length(var_x_norm))))
    mean_norm_matrix = rbind(mean_norm_matrix, cbind(mean_x_norm, rep(label,length(mean_x_norm))))
    
    sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
    all.genes <- rownames(sct)
    sct <- ScaleData(sct, features = all.genes)
    sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
    pca <- sct[["pca"]]
    
    pca_out = data.frame(sct[['pca']]@cell.embeddings)
    kurtosi_tota_x = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep(label,dim(pca_out)[1])))
    var_total = data.frame(rbind(var_total, var_total_tet))
    kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_tet))
}


genes = names(data)[-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
var_frame = as.data.frame(var_matrix)
var_norm_frame = as.data.frame(var_norm_matrix)
var_norm_frame[,1] = as.numeric(as.character(var_norm_frame[,1]))
var_frame[,1] = as.numeric(as.character(var_frame[,1]))
names(var_frame) = c('Variance','Condition')
names(var_norm_frame) = c('Variance','Condition')


mean_frame = as.data.frame(mean_matrix)
mean_norm_frame = as.data.frame(mean_norm_matrix)
mean_norm_frame[,1] = as.numeric(as.character(mean_norm_frame[,1]))
mean_frame[,1] = as.numeric(as.character(mean_frame[,1]))
names(mean_frame) = c('Mean','Condition')
names(mean_norm_frame) = c('Mean','Condition')

fano = var_frame
fano[,1] = fano[,1]/mean_frame[,1]

var_frame$Mean = mean_frame$Mean
var_norm_frame$Mean = mean_norm_frame$Mean

ggplot( data = fano[keep_rows,],aes(x=Variance, fill=Condition, color = Condition)) +
#  geom_point()+
  scale_y_log10()+
  geom_freqpoly(binwidth = 0.1)+
#  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") +  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 
var_frame$var_sqrt = sqrt(var_frame$Variance)
var_frame$fano = var_frame$Variance/var_frame$Mean
var_frame$log_fano = log(var_frame$fano)
var_frame$log_mean = log(var_frame$Mean )
var_frame$norm_var = var_norm_frame$Variance
var_frame$log_norm_var = log(var_norm_frame$Variance)
var_frame$norm_mean = mean_norm_frame$Mean
var_frame$log_norm_mean = log(mean_norm_frame$Mean)
var_frame$norm_fano = var_frame$norm_var/var_frame$norm_mean
var_frame$log_norm_fano = log(var_frame$norm_fano)
var_frame$sqrt_var = sqrt(var_frame$norm_var)
var_frame = var_frame[!is.na(var_frame$Condition),]
var_frame = var_frame[which(var_frame$Condition != 'Overnight'),]
var_frame$gene = rep(genes,11)
keep_rows  = which(var_frame$Condition =='Cef' | var_frame$Condition == "Lambda30" |  var_frame$Condition == "Exponential" | var_frame$Condition == 'Chlora')
#test = var_frame[which(var_frame$norm_fano > 2),]
#test = test[which(test$norm_mean > 0.1),]
#test = test[which(test$Condition == 'Lambda30'),]
#tail(test,n = 15)
#test = test[which(test$gene == 'Lambda_A'),]
ggplot( data = var_frame[keep_rows,],aes(x=Mean, y = Variance,fill=Condition, color = Condition)) +
 geom_point()+
#  xlim(0.1, 3)+
  #geom_line()+
#  scale_y_log10() + 
#  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") +  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") #+xlim(0.1,2.5)

ggplot( data = var_frame[keep_rows,],aes(x=Mean, y = fano,fill=Condition, color = Condition)) +
 geom_point()+
#  geom_line()+
#  scale_y_log10() + 
#  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
#  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") +  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 


keep_rows  = which(var_frame$Condition =='Teta' | var_frame$Condition == "Lambda30" |  var_frame$Condition == "Exponential" | var_frame$Condition == 'Chlora')
ggplot( data = var_frame[keep_rows,],aes(x=Variance, fill=Condition, color = Condition)) +
#  geom_point()+
  #  xlim(0.1, 3)+
  #geom_line()+scale_y_log10()+
  geom_freqpoly(binwidth = 0.1)+
    scale_y_log10() + 
  #  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  #  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") +  theme_bw() +
  ylab('Count') + xlab('Variance') + 
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") #+xlim(0.1,2.5)
ggsave(paste(figs.out,'lambda_exponential_variance','.png'))
ggplot( data = var_frame[keep_rows,],aes(x=Mean, y = Variance,fill=Condition, color = Condition)) +
  geom_point()+
  #  xlim(0.1, 3)+
  #geom_line()+
  #  scale_y_log10() + 
  #  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  #  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  ggtitle('Variance/Mean') + 
  labs(fill="") +  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") #+xlim(0.1,2.5)

ggsave(paste(figs.out,'lambda_exponential_mean_variance','.png'))
ggplot( data = var_frame[keep_rows,],aes(x=Mean, y = fano,fill=Condition, color = Condition)) +
  geom_point(alpha = 0.5)+
  #  xlim(0.1, 3)+
  #geom_line()+
  #  scale_y_log10() + 
  #  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  #  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  ggtitle('Var/Mean to Mean') + 
  labs(fill="") +  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") #+xlim(0.1,2.5)

ggsave(paste(figs.out,'lambda_exponential_mean_fano','.png'))


ggplot( data = var_frame[keep_rows,],aes(x=norm_mean, y = norm_var,fill=Condition, color = Condition)) +
  geom_point(alpha = 0.6)+
  #  xlim(0.1, 3)+
  #geom_line()+
  #  scale_y_log10() + 
  #  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  #  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  ggtitle('Normalized Variance/Mean') + 
  labs(fill="") +  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") #+xlim(0.1,2.5)

ggsave(paste(figs.out,'lambda_exponential_mean_norm_var','.png'))


keep_rows  = which(var_frame$Condition =='Tet' | var_frame$Condition == "Exponential" |  var_frame$Condition == "Exponential" | var_frame$Condition == 'Chlora')
ggplot( data = var_frame[keep_rows,],aes(x=Variance, fill=Condition, color = Condition)) +
#  geom_point()+
  #  xlim(0.1, 3)+
  #geom_line()+scale_y_log10()+
  geom_freqpoly(binwidth = 0.1)+
    scale_y_log10() + 
  #  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  #  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") +  theme_bw() +
  ylab('Count') + xlab('Variance') + 
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") #+xlim(0.1,2.5)
ggsave(paste(figs.out,'lambda_exponential_variance','.png'))



# Matrix plot
library(dplyr)

av.exp.sct <- AverageExpression(sct, return.seurat = TRUE)
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
#sct@active.ident = factor(sct@meta.data$Label,labels = sct@meta.data$Label, levels = sct@meta.data$Label)
markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.5)
markers = markers[markers$avg_log2FC > 0.5,]
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 3, wt = avg_log2FC) -> top10

markers = markers %>%
  add_count(cluster) %>%
  group_by(cluster) %>%
  mutate(n = 8) %>%
  #  group_by(rating, .add = TRUE) %>%
  #In old dplyr add = TRUE
  #group_by(rating, add = TRUE) %>%
  sample_n(n, replace = TRUE) %>%
  select(-n)

markers = unique(markers)
#markers = markers[1:8,]
genes = markers$gene
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
g = merge(g, markers)
#g = g[which(g$variable%in% genes),]
g = g[order(g$cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'EC-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = g$gene,labels = g$gene)
#g$Cluster = factor(g$Cluster, levels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'),labels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'))
g$Cluster = factor(as.numeric(g$Cluster)+1, levels = c(5,4,3,2,1))
#g$Cluster = factor(g$Cluster, levels = c(5,4,3,2,1))

g = g[order(g$cluster),]

library(myriad)
ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.0), oob = scales::squish)+
  labs(fill='', vjust = 1) +
  
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=11,color = 'black'),
        axis.title=element_text(size=13), legend.text = element_text(size = 11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(size=11, vjust = 3)) + scale_colour_hue(l=40)

#### Cidal vs static ###
#### Tet ####
keep_rows  = which(data2$Label =="Tet" | data2$Label =="Chlor" |  data2$Label == "Erythro"|data2$Label=="Gent")
data = data2[keep_rows,]
library(dplyr)
#data = sample_n(data, 300)
k = colSums(data!=0) 
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]
#### Counts/gene
#### Counts/gene
m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
counts = data.frame(m,data$Label)
names(counts) = c('count','annotation')
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
to_filter = which(k > 10)
new_names = names(to_filter)

new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data[which(names(data) %in% c('Label'))]
sc = CreateSeuratObject(
  #  sc_precursor,
  sc_precursor,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = sc_meta
)
sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
#sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 10000)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)

rm(mat)
rm(sc_precursor)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
pca <- sct[["pca"]]

# Get the total variance:

pca_out = data.frame(sct[['pca']]@cell.embeddings)
#plot(1:50,sort(abs(kurtosi(pca_out)),decreasing = TRUE))
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
var_total_tet = data.frame(cbind(1:length(varExplained), varExplained, rep('Tet',length(varExplained))))
kurtosi_total_tet = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Tet',length(varExplained))))
var_total = data.frame(rbind(var_total, var_total_tet))
kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_tet))

#sct = JackStraw(sct)

print(sct[["pca"]], dims = 1:25, nfeatures = 5)

DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(1,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
kurtosi(pca_out$PC_3)
rm(mat)
rm(sc_precursor)
sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.05,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,logfc.threshold = 0.1)
#cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6,logfc.threshold = 0.1)
#cluster6.markers <- FindMarkers(sct, ident.1 = 6,ident.2=4, min.pct = 0.1)
#cluster7.markers <- FindMarkers(sct, ident.1 = 7, ident.2 = 1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7,logfc.threshold = 0.1)
cluster8.markers <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1,logfc.threshold = 0.1)
cluster9.markers <- FindMarkers(sct, ident.1 = 9, min.pct = 0.1,logfc.threshold = 0.1)
cluster10.markers <- FindMarkers(sct, ident.1 = 10, min.pct = 0.1, logfc.threshold = 0.1)
cluster11.markers <- FindMarkers(sct, ident.1 = 11, min.pct = 0.1, logfc.threshold = 0.1,ident.2 = 3)
cluster14.markers <- FindMarkers(sct, ident.1 = 14, min.pct = 0.1, logfc.threshold = 0.1)
cluster15.markers <- FindMarkers(sct, ident.1 = 15, min.pct = 0.1, logfc.threshold = 0.1)
cluster16.markers <- FindMarkers(sct, ident.1 = 16, min.pct = 0.1, logfc.threshold = 0.1)
cluster13.markers <- FindMarkers(sct, ident.1 = 13, min.pct = 0.1, logfc.threshold = 0.1)
cluster12.markers <- FindMarkers(sct, ident.1 = 12, min.pct = 0.1, logfc.threshold = 0.1)
cluster17.markers <- FindMarkers(sct, ident.1 = 17, min.pct = 0.1, logfc.threshold = 0.1)
#my_levels <- c(0,2,3,1)
#levels(sct) = my_levels
#markers %>%
#  group_by(cluster) %>%  top_n(n = 2, wt = avg_log2FC)
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 5, wt = avg_log2FC) -> top10
#DoHeatmap(sct, features = top10$gene,draw.lines = TRUE)  + scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ 
#  theme(text = element_text(size = 18))
#ggsave(paste(figs.out,'MG1655_heatmap_exp_3hr','.png'))
#### PCA ####
#### get the PC's ###

pca_out = data.frame(sct[['pca']]@cell.embeddings)
plot(1:50,sort(abs(kurtosi(pca_out)),decreasing = TRUE))
pca_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
pca_out$label = as.character(sct[['Label']][,1])
pca_out$ident = as.character(sct@active.ident)

#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)


umap_out$label = gsub('Gent','Gentamycin',umap_out$label)
umap_out$label = gsub('Cipro','Ciprofloxacin',umap_out$label)
umap_out$label = gsub('Cef','Cefazolin',umap_out$label)
umap_out$label = gsub('Nal','Nalidixic acid',umap_out$label)
#umap_out$label = gsub('Cycloserineserine','Cycloserine',umap_out$label)
umap_out$label = gsub('Cyclo','Cycloserine',umap_out$label)
#umap_out$label = gsub('Erythromycinmycin','Erythromycin',umap_out$label)
umap_out$label = gsub('Erythro','Erythromycin',umap_out$label)
#umap_out$label = gsub('Tetracyclineracycline','Tetracycline',umap_out$label)
umap_out$label = gsub('Tet','Tetracycline',umap_out$label)
umap_out$label = gsub('Chlor','Chloramphenicol',umap_out$label)
unique(umap_out$label)
#umap_out$label = gsub('8hr','T360',umap_out$label)
#colors = data.frame(read.csv(paste(figs.out,'color_map.csv'), sep = ','))[,-1]
#umap_out = umap_out[umap_out$label=='STATIONARY',]
#umap_out = merge(umap_out, colors,by.x = 'label',by.y = "umap_out.label")
library(ggsci)

library(scales)
colors = hue_pal()(length(unique(umap_out$ident)))
umap_out$ident <- factor(umap_out$ident, levels=as.character(sort(as.numeric(unique(umap_out$ident)))), labels=as.character(sort(as.numeric(unique(umap_out$ident)))))

aggregate(umi_counts~ident,umap_out,mean)
sizes = aggregate(umi_counts~ident,umap_out,length)
total = sum(aggregate(umi_counts~ident,umap_out,length)['umi_counts'])
sizes['umi_counts']/total
umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Erythromycin","Chloramphenicol"), labels=c("Gentamycin", "Tetracycline", "Erythromycin","Chloramphenicol"))

#umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.05,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,3,8,12)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2.png'))
#scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,10,5)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette1.png'))
ggsave(paste(figs.out, 'umap_MG1655_clean_palette_40pcs_statics_all.png'))


levs = as.character(0:18)
umap_out$ident = factor(umap_out$ident, levels = levs)
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.05,alpha = 0.8) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 8), legend.title = element_text(size=15,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_MG1655_statics_clusters.png'))
#scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,10,5)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette1.png'))
#ggsave(paste(figs.out, 'umap_MG1655_clusters_40pcs.png'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.05,alpha = 1.5) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
  # ggtitle(sprintf("E. coli")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))#+
#  scale_color_simpsons()s
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,5)])
ggsave(paste(figs.out, 'umap_MG1655_tet_chlor_clusters.png'))
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2_30_pcs.png'))

#### Tet downsamples ####
library(grDevices)
col1 = as.vector(rcartocolor::carto_pal(name='Safe'))[c(7)]
#### Tet ####
kurtosi_total = matrix(ncol = 3)
nrows =c(500,1000,2500,5000,7500,10000,15000,20000,27374)
nsamples = length(nrows)
colors = colorRampPalette(c(col1, 'red4'))(nsamples)
for (i in 1:length(nrows)){
  keep_rows  = which(data2$Label =="Tet" | data2$Label =="Teet" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
  n_rows = nrows[i]
  data = data2[keep_rows,]
  data = data[sample(dim(data)[1], n_rows,FALSE),]
  library(dplyr)
  #data = sample_n(data, 300)
  k = colSums(data!=0) 
  to_filter = which(k > 5)
  new_names = names(to_filter)
  new_names = c(new_names)
  data = data[new_names]
  #### Counts/gene
  #### Counts/gene
  m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
  counts = data.frame(m,data$Label)
  names(counts) = c('count','annotation')
  k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
  to_filter = which(k > 10)
  new_names = names(to_filter)
  new_names = c(new_names, "Label",'Identity')
  new_names = c(new_names)
  data = data[new_names]
  k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
  to_filter = which(k2 >= 10)
  data = data[to_filter,]
  ## Make Seurat
  
  sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
  rownames(sc_precursor) = 1:dim(sc_precursor)[1]
  sc_precursor = t(sc_precursor)
  sc_meta = data[which(names(data) %in% c('Label'))]
  sc = CreateSeuratObject(
    #  sc_precursor,
    sc_precursor,
    project = "SeuratProject",
    assay = "RNA",
    min.cells = 0,
    min.features = 0,
    names.field = 1,
    names.delim = "_",
    meta.data = sc_meta
  )
  sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
  #sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 10000)
  sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
  top10 <- head(VariableFeatures(sct), 10)
  all.genes <- rownames(sct)
  
  rm(mat)
  rm(sc_precursor)
  sct <- ScaleData(sct, features = all.genes)
  sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
  pca_out = data.frame(sct[['pca']]@cell.embeddings)
  kurtosi_total_tet = cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep(n_rows,dim(pca_out)[1]))
  kurtosi_total =rbind(kurtosi_total, kurtosi_total_tet)
  
  #sct = JackStraw(sct)
  
  print(sct[["pca"]], dims = 1:25, nfeatures = 5)
  
  DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(3,9)) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
    theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#  ggsave(paste(figs.out, 'pca_pc39_tet.png'))
  kurtosi(pca_out$PC_3)
  rm(mat)
  rm(sc_precursor)
  sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
  sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
  sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
  #### UMAP Variability ###
  library(ggsci)
  umap_out = data.frame(sct[['umap']]@cell.embeddings)
  umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
  umap_out$label = as.character(sct[['Label']][,1])
  umap_out$ident = as.character(sct@active.ident)
  #umap_out = umap_out[umap_out$label=='STATIONARY',]
  umap_out$label = gsub('_',' ',umap_out$label)
  umap_out$label = gsub('Tet',n_rows,umap_out$label)
  unique(umap_out$label)
  library(rcartocolor)
  #ggsave(paste(figs.out, 'umap_cluster_MG1655_total.png'))
  
  #umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'))
  ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
    geom_point( size =(1/i) * 1,alpha = 0.8) +  
    #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
    labs(x = "", y = "", color = "Condition") +
    ggtitle(sprintf("")) +
    #  coord_cartesian(xlim= c(0,1.0))+
    #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
    #  coord_equal(ratio=1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text = element_blank(),
          #        axis. = element_blank(),
          #        axis.line = element_line(colour = "black"),
          #        axis.text=element_text(size=20),axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank() ,
          axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
    #  scale_color_simpsons()s
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
    scale_color_manual(values=colors[i])
  ggsave(paste(figs.out, 'umap_MG1655_tet_sampled_',n_rows,'.png'))
  
  
  
}


kurtosi_total = kurtosi_total[-1,]
kurtosi_total = data.frame(kurtosi_total)
names(kurtosi_total) = c('X1','kurtosi','label')

#kurtosi_total$label <- factor(kurtosi_total$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'))



kurtosi_total$kurtosi = as.numeric(as.character(kurtosi_total$kurtosi))
kurtosi_total$X1 = as.numeric(as.character(kurtosi_total$X1))

#kurtosi_total_max = kurtosi_total %>%group_by()



kurtosi_max = kurtosi_total[,c('label','kurtosi')] %>% group_by(label) %>% summarise(max = max(kurtosiExplained, na.rm=TRUE))
kurtosi_sum = kurtosi_total[,c('label','kurtosi')] %>% group_by(label) %>% summarise(sum = sum(kurtosiExplained, na.rm=TRUE))
#kurtosi_total = merge(kurtosi_total,kurtosi_max)
kurtosi_total = merge(kurtosi_total,kurtosi_sum)
kurtosi_total$kurtosi_norm = kurtosi_total$kurtosiExplained/kurtosi_total$max


#kurtosi_total$PC = as.vector(sapply(kurtosi_total$X, substr, start = 1, stop = 5))
#kurtosi_total$X = kurtosi_total$X
kurtosi_top5 = kurtosi_total[which(kurtosi_total$X1 %in% 1:5),]
kurtosi_top5 = kurtosi_top5[,c('label','kurtosiExplained', 'sum')] %>% group_by(label) %>% summarise(mean_kurt = mean(kurtosiExplained), mean_sum = mean(sum))
#kurtosi_top5 = kurtosi_top5[which(kurtosi_total$label == 'Tetracycline' | kurtosi_total$label == 'Erythromycin' | kurtosi_total$label == 'Gentamycin' | kurtosi_total$label == 'Chlor' | kurtosi_total$label == 'Exponential' ),]
kurtosi_rows = which(kurtosi_total$label == 'Tetracycline' )#| kurtosi_total$label == 'Erythromycin' | kurtosi_total$label == 'Gentamycin' | kurtosi_total$label == 'Chlor' | kurtosi_total$label == 'Exponential' )
kurtosi_rows = intersect(which(kurtosi_total$X1 %in% 1:5), which(kurtosi_total$label == 'Erythromycin'))
kurtosi_pick = kurtosi_total[kurtosi_rows,]
#mechanisms

kurtosi_top5['mechanism'] = c('30S Ribosome', '30S Ribosome',"Cell Wall", "Cell Wall",'50S Ribosome', '50S Ribosome', 'DNA Damage','DNA Damage','None','Phage','Phage')
kurtosi_top5['cidality'] = c('Cidal', 'Static',"Cidal", "Cidal",'Cidal', 'Static', 'Cidal','Cidal','None','Phage','Phage')
kurtosi_top5$mean_kurt = as.numeric(kurtosi_top5$mean_kurt)
kurtosi_top5$mechanism = factor(kurtosi_top5$mechanism, levels = c('30S Ribosome', "Cell Wall", '50S Ribosome', 'DNA Damage','None','Phage'))
drug_mapper = kurtosi_top5[c('label','mechanism','cidality')]
kurtosi_top5v2 = kurtosi_total[which(kurtosi_total$X1 %in% 1:3),]
kurtosi_top5v2 = merge(kurtosi_top5v2,drug_mapper)

kurtosi_total = merge(kurtosi_total, drug_mapper)

kurtosi_total$label2 = as.character(kurtosi_total$label)
ggplot(data = kurtosi_total, aes (x = X1, y = kurtosi, group = label, color = label)) +
   geom_point(size = 0.8, alpha = 0.8) +
  geom_line(position = 'identity', size = 1, alpha = 0.8 ) + 
  labs(x = "Rank", y = "kurtosi") +
  ggtitle(sprintf("kurtosi for ribosome antibiotics")) +
  scale_y_log10()+
  theme_bw() +
  scale_color_gradient2(
    low = "blue", 
    mid = "gray", 
    high = "brown", 
    midpoint = 2500
  )+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=0.01,face="bold"))# + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12)])
ggsave(paste(figs.out, 'kurtosi_MG1655_conditions.png'))



#### Erythro subsample ####
kurtosi_total = matrix(ncol = 3)
nrows =c(1000,2500,5000,7500,10000,15000,20000,30000)
#nrows =c(500,1000)#,2500,5000,7500,10000,15000,20000)
nsamples = length(nrows)
#col1 = 
col1 = as.vector(rcartocolor::carto_pal(name='Safe'))[c(7)]

#colors = colorRampPalette(c("skyblue1", col1))(nsamples)
colors = colorRampPalette(c("skyblue1", 'gray','brown'))(nsamples)
for (i in 1:length(nrows)){
  keep_rows  = which(data2$Label =="Tet" | data2$Label =="Teet" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
  n_rows = nrows[i]
  data = data2[keep_rows,]
  data = data[sample(dim(data)[1], n_rows,FALSE),]
  library(dplyr)
  #data = sample_n(data, 300)
  k = colSums(data!=0) 
  to_filter = which(k > 5)
  new_names = names(to_filter)
  new_names = c(new_names)
  data = data[new_names]
  #### Counts/gene
  #### Counts/gene
  m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
  counts = data.frame(m,data$Label)
  names(counts) = c('count','annotation')
  k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
  to_filter = which(k > 10)
  new_names = names(to_filter)
  new_names = c(new_names, "Label",'Identity')
  new_names = c(new_names)
  data = data[new_names]
  k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
  to_filter = which(k2 >= 10)
  data = data[to_filter,]
  ## Make Seurat
  
  sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
  rownames(sc_precursor) = 1:dim(sc_precursor)[1]
  sc_precursor = t(sc_precursor)
  sc_meta = data[which(names(data) %in% c('Label'))]
  sc = CreateSeuratObject(
    #  sc_precursor,
    sc_precursor,
    project = "SeuratProject",
    assay = "RNA",
    min.cells = 0,
    min.features = 0,
    names.field = 1,
    names.delim = "_",
    meta.data = sc_meta
  )
#  sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
  sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
  sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
  top10 <- head(VariableFeatures(sct), 10)
  all.genes <- rownames(sct)
  
  rm(mat)
  rm(sc_precursor)
  sct <- ScaleData(sct, features = all.genes)
  sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
  pca_out = data.frame(sct[['pca']]@cell.embeddings)
  kurtosi_total_tet = cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep(n_rows,dim(pca_out)[1]))
  kurtosi_total =rbind(kurtosi_total, kurtosi_total_tet)
  
  pcs = names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]
  
  #sct = JackStraw(sct)
  
  print(sct[["pca"]], dims = 1:25, nfeatures = 5)
  
  DimPlot(sct, label = TRUE,group.by = 'Label',label.size = 0.1,reduction = 'pca', cols =  colors[i],dims =as.numeric(gsub('PC_','',names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]))) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
    theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  #  ggsave(paste(figs.out, 'pca_pc39_tet.png'))
  ggsave(paste(figs.out, 'pca_MG1655_tet_sampled_scifi1039',n_rows,'.png'))
  kurtosi(pca_out$PC_3)
  rm(mat)
  rm(sc_precursor)
  sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
  sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
  sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
  #### UMAP Variability ###
  library(ggsci)
  umap_out = data.frame(sct[['umap']]@cell.embeddings)
  umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
  umap_out$label = as.character(sct[['Label']][,1])
  umap_out$ident = as.character(sct@active.ident)
  #umap_out = umap_out[umap_out$label=='STATIONARY',]
  umap_out$label = gsub('_',' ',umap_out$label)
  umap_out$label = gsub('Tet',n_rows,umap_out$label)
  unique(umap_out$label)
  library(rcartocolor)
  #ggsave(paste(figs.out, 'umap_cluster_MG1655_total.png'))
  
  #umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'))
  ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
    geom_point( size =(1/i) * 1,alpha = 0.8) +  
    #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
    labs(x = "", y = "", color = "Condition") +
    ggtitle(sprintf("")) +
    #  coord_cartesian(xlim= c(0,1.0))+
    #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
    #  coord_equal(ratio=1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text = element_blank(),
          #        axis. = element_blank(),
          #        axis.line = element_line(colour = "black"),
          #        axis.text=element_text(size=20),axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank() ,
          axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
    #  scale_color_simpsons()s
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
    scale_color_manual(values=colors[i])
  ggsave(paste(figs.out, 'umap_MG1655_tet_sampled_scifi10_lognorm',n_rows,'.png'))
  
  
  
}

kurtosi_total = kurtosi_total[-1,]
kurtosi_total = data.frame(kurtosi_total)
names(kurtosi_total) = c('X1','kurtosi','label')

#kurtosi_total$label <- factor(kurtosi_total$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'))



kurtosi_total$kurtosi = as.numeric(as.character(kurtosi_total$kurtosi))
kurtosi_total$X1 = as.numeric(as.character(kurtosi_total$X1))
kurtosi_total2 = subset(kurtosi_total, X1 < 15)
#ggplot(data = subset(kurtosi_total, label > 1000), aes (x = X1, y = kurtosi, group = label, color = label)) +
ggplot(data =kurtosi_total2, aes (x = X1, y = kurtosi, group = label, color = label)) +
   geom_point(size = 0.8, alpha = 0.8) +
  geom_line(position = 'identity', size = 0.5, alpha = 1.8 ) + 
  labs(x = "Rank", y = "kurtosi") +
  ggtitle(sprintf("kurtosi for ribosome antibiotics")) +
#  scale_y_log10()+
  theme_bw() +
  scale_color_gradient2(
    low = "skyblue1", 
    mid = "gray", 
    high = 'brown', 
    midpoint = (max(nrows) + min(nrows))/2
  )+
  theme(plot.title = element_text(hjust = 0.5, size = 0.01, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 9), legend.title = element_text(size=0.01,face="bold"))# + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12)])
#ggsave(paste(figs.out, 'umap_MG1655_tet_scifi10_ramp_with_words','.png'))
ggsave(paste(figs.out, 'umap_MG1655_tet_scifi10_ramp_with_words_lognorm','.png'))

###

#### Gad on the whole thing 

#gad_vec = c('EC_gadA','EC_gadB','EC_gadC')
gad_vec = c('EC_tfaQ','EC_tfaR','EC_pinQ','EC_pinR','EC_ynaE','EC_ydfK','EC_insI2','EC_insI3')
gad_df = data.frame(gad_vec, rep(1,length(gad_vec)))
names(gad_df) = c('gene_id','gad_stage')
gad_df$gene_id = gsub('_','-',gad_df$gene_id)
gad_score = matrix(data = 0,ncol=1,nrow = dim(sct@assays$RNA)[1])
gene_indices = c()
for (i in 1:length(gad_df$gene_id)){
  index = which(gad_df$gene_id[i] == rownames(sct@assays$RNA))
  if (length(index)>0){
    gene_indices = c(gene_indices,index)
    select = c(select,i)
  }
}
gad_df$gene_indices = gene_indices
gad_score[gad_df$gene_indices,1] = as.numeric(gad_df$gad_stage)
norm_data = t(as.data.frame(sct@assays$RNA@data))
gad_data = data.frame(norm_data %*% gad_score)
names(gad_data) = c('norm_gad')

gad_data$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
gad_data$label = as.character(sct[['Label']][,1])
gad_data$ident = as.character(sct@active.ident)
gad_data$index = 1:dim(gad_data)[1]


umap_out2 = cbind(umap_out, gad_data$norm_gad)
umap_out2$max_val = pmax(gad_data$norm_gad)

ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `gad_data$norm_gad`, fill = `gad_data$norm_gad`)) + 
  #  geom_point(pch = 21, aes(fill = `gad_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5) +
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #  scale_fill_gradientn(colors =  (brewer.pal(n = 9, name =  "BuGn") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "BuGn")), guide = '' ) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) +   #scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greens")[1]), high = (brewer.pal(n = 9, name = "Greens")[9])) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'chartreuse1') + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'chartreuse1') + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
  #  scale_colour_gradient(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 

ggsave(paste(figs.out,"gadAB_umap_all",'.png'))
ggplot(data =kurtosi_total2, aes (x = X1, y = kurtosi, group = label, color = label)) +
  geom_point(size = 1.3, alpha = 0.8) +
  geom_line(position = 'identity', size = 1.2, alpha = 1.8 ) + 
  labs(x = "Rank", y = "kurtosi") +
  ggtitle(sprintf("kurtosi for ribosome antibiotics")) +
  #  scale_y_log10()+
  theme_bw() +
  scale_color_gradient2(
    low = "skyblue1", 
    mid = "gray", 
    high = 'brown', 
    midpoint = (max(nrows) + min(nrows))/2
  )+
  theme(plot.title = element_text(hjust = 0.5, size = 0.01, face="bold.italic"),
        axis.title=element_text(size=0.05,face="bold"), legend.text = element_text(size = 0.5), legend.title = element_text(size=0.01,face="bold"),
        axis.text = element_blank(),
        panel.grid =element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
        axis.text.y=element_blank())
#        axis.ticks.y=element_blank() ,)# + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12)])
#ggsave(paste(figs.out, 'umap_MG1655_tet_scifi10_ramp_no_words','.png'))
ggsave(paste(figs.out, 'umap_MG1655_tet_scifi10_ramp_no_words_lognorm','.png'))
#library(ggpattern)
ggplot(data = kurtosi_top5, aes (x = mechanism, y = mean_sum, fill = label)) +
  geom_bar(stat='identity',position='dodge') + 
  #  geom_line(size = 1, alpha = 0.8) + 
  #  geom_bar_pattern(aes(fill = label),stat='identity', position=position_dodge()) + 
  #  geom_bar_pattern(aes(fill = label),stat='identity', position=position_dodge()) + 
  #scale_pattern_manual(values = c(Nerd = "stripe", NotNerd = "none")) +
  #   geom_point(size = 0.8, alpha = 0.8) +
  labs(x = "Rank", y = "Kurtosis") +
  #  scale_y_log10() + 
  ggtitle(sprintf("kurtosi")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,5,6,12)])
ggsave(paste(figs.out,'total_kurtosis','.png'))

kurtosi_top5v2$label = factor(gsub(x = kurtosi_top5v2$label, pattern = 'Chlor','Chloramphenicol'))
kurtosi_top5v2$label = factor(kurtosi_top5v2$label,levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'))
#kurtosi_top5v2$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential'))



ggplot(data = kurtosi_top5v2[which(kurtosi_top5v2$mechanism %in% c('30S Ribosome','50S Ribosome')),], aes (x = mechanism, y = kurtosiExplained, fill = label)) +
  #  geom_line(size = 1, alpha = 0.8) + 
  geom_boxplot(aes(color = cidality), lwd = 1)+
  #  geom_bar(aes(fill = label),stat='identity', position=position_dodge()) + 
  #   geom_point(size = 0.8, alpha = 0.8) +
  labs(x = "Mechanism", y = "Kurtosis") +
  scale_y_log10() + 
  #  ggtitle(sprintf("Kurtosis of Top3 Ranked PCs")) +
  ggtitle(sprintf("")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold.italic"),axis.text = element_text(size = 0.1),
        axis.title=element_text(size=0.1,face="bold"), legend.text = element_text(size = 12), legend.title = element_text(size=0.1,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=2)))+
  #  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,12,5,6)])
  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12,5,6)])+
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(12,5,12,5,6)])
ggsave(paste(figs.out,'kurtosis_ribosomes','.png'))
ggplot(data = kurtosi_top5v2[which(kurtosi_top5v2$mechanism %in% c('30S Ribosome','50S Ribosome')),], aes (x = mechanism, y = kurtosiExplained, fill = label)) +
  #  geom_line(size = 1, alpha = 0.8) + 
  geom_boxplot(aes(color = cidality), lwd = 1)+
  #  geom_bar(aes(fill = label),stat='identity', position=position_dodge()) + 
  #   geom_point(size = 0.8, alpha = 0.8) +
  labs(x = "Mechanism", y = "Kurtosis") +
  scale_y_log10() + 
  #  ggtitle(sprintf("Kurtosis of Top3 Ranked PCs")) +
  ggtitle(sprintf("")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 10, size = 10, face="bold.italic"),axis.text = element_text(size = 10),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 12), legend.title = element_text(size=0.1,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=2)))+
  #  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,12,5,6)])
  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12,5,6)])+
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(12,5,12,5,6)])
ggsave(paste(figs.out,'kurtosis_ribosomes_with_words','.png'))

ggplot(data = kurtosi_pick, aes (x = X1, y = kurtosiExplained, col = label)) +
  geom_line(size = 1, alpha = 0.8) + 
  # geom_point(size = 0.8, alpha = 0.8) +
  labs(x = "Rank", y = "Kurtosis") +
  scale_y_log10() + 
  ggtitle(sprintf("kurtosi")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,5,6,12)])
#ggsave(paste(figs.out,'pca_MG1655_var_explained','.png'))


#### static vs cidal ####

#### Static subsample ####
data_total = 2
kurtosi_total = matrix(ncol = 3)
#nrows =c(1000,2500,5000,7500,8027,10000,15000,20000,30000)
nrows =c(1000,2500,5000,7500,10000,15000,25000,35000,45000,55000,65000,75000)
#nrows =c(500,1000)#,2500,5000,7500,10000,15000,20000)
nsamples = length(nrows)
#col1 = 
col1 = as.vector(rcartocolor::carto_pal(name='Safe'))[c(7)]

#colors = colorRampPalette(c("skyblue1", col1))(nsamples)
colors = colorRampPalette(c("skyblue1", 'gray','brown'))(nsamples)
feature_colors = colorRampPalette(c("gray", 'gray','yellow','red'))(5)
library(scales)
show_col(feature_colors)

for (i in 1:length(nrows)){
  keep_rows  = which(data2$Label =="Tet" | data2$Label =="Chlor" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
  n_rows = nrows[i]
  data = data2[keep_rows,]
  data = data[sample(dim(data)[1], n_rows,FALSE),]
  library(dplyr)
  #data = sample_n(data, 300)
  k = colSums(data!=0) 
  to_filter = which(k > 5)
  new_names = names(to_filter)
  new_names = unique(c(new_names,'EC_insI-2'))
  data = data[new_names]
  #### Counts/gene
  #### Counts/gene
  m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
  counts = data.frame(m,data$Label)
  names(counts) = c('count','annotation')
  k = colSums(data[-which(names(data) %in% c('Label','Identity'))])
  to_filter = which(k > 10)
  new_names = names(to_filter)
  new_names = c(new_names, "Label",'Identity','EC_insI-2')
  new_names = unique(c(new_names))
  data = data[new_names]
  k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
  to_filter = which(k2 >= 10)
  data = data[to_filter,]
  ## Make Seurat
  
  sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
  rownames(sc_precursor) = 1:dim(sc_precursor)[1]
  sc_precursor = t(sc_precursor)
  sc_meta = data[which(names(data) %in% c('Label'))]
  sc = CreateSeuratObject(
    #  sc_precursor,
    sc_precursor,
    project = "SeuratProject",
    assay = "RNA",
    min.cells = 0,
    min.features = 0,
    names.field = 1,
    names.delim = "_",
    meta.data = sc_meta
  )
  #sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
  sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
  sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
  top10 <- head(VariableFeatures(sct), 10)
  all.genes <- rownames(sct)
  
  rm(mat)
  rm(sc_precursor)
  sct <- ScaleData(sct, features = all.genes)
#  sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
#  sct <- RunPCA(sct, verbose = FALSE,features = all.genes, npcs = 50)
  sct <- RunPCA(sct, verbose = FALSE,features = all.genes, npcs = 100)
  pca_out = data.frame(sct[['pca']]@cell.embeddings)
  kurtosi_total_tet = cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep(n_rows,dim(pca_out)[1]))
  kurtosi_total =rbind(kurtosi_total, kurtosi_total_tet)
  
  pcs = names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]
  
  #sct = JackStraw(sct)
  
  print(sct[["pca"]], dims = 1:25, nfeatures = 5)
  
  DimPlot(sct, label = TRUE,label.size = 0.1,reduction = 'pca', cols =  colors[i],dims =as.numeric(gsub('PC_','',names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]))) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
    theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  #  ggsave(paste(figs.out, 'pca_pc39_tet.png'))
  pca_out = data.frame(sct[['pca']]@cell.embeddings)
  
  pca_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
  pca_out$label = as.character(sct[['Label']][,1])
  pca_out$label = gsub('Tet',n_rows,pca_out$label)
  pca_out$label = gsub('Chlor',n_rows,pca_out$label)
  ggplot(pca_out,aes_string(x = pcs[1], y = pcs[2],col = "label")) + 
    geom_point( size =(1/i) * 1,alpha = 0.8) +  
    #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
    labs(x = "", y = "", color = "Condition") +
    ggtitle(sprintf("")) +
    #  coord_cartesian(xlim= c(0,1.0))+
    #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
    #  coord_equal(ratio=1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text = element_blank(),
          #        axis. = element_blank(),
          #        axis.line = element_line(colour = "black"),
          #        axis.text=element_text(size=20),axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank() ,
          axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
    #  scale_color_simpsons()s
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
    scale_color_manual(values=colors[i])
  ggsave(paste(figs.out, 'pca_MG1655_tet_chlor_sampled_scifi_lognorm',n_rows,pcs[1],pcs[2],'.png'))
  kurtosi(pca_out$PC_3)
  rm(mat)
  rm(sc_precursor)
  sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
  sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
  
  
  #### UMAP Variability ###
  library(ggsci)
  umap_out = data.frame(sct[['umap']]@cell.embeddings)
  umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
  umap_out$label = as.character(sct[['Label']][,1])
  umap_out$ident = as.character(sct@active.ident)
  #umap_out = umap_out[umap_out$label=='STATIONARY',]
  umap_out$label = gsub('_',' ',umap_out$label)
  unique(umap_out$label)
  library(rcartocolor)
  #ggsave(paste(figs.out, 'umap_cluster_MG1655_total.png'))
  #umap_out$label = gsub('Tet',n_rows,umap_out$label)
  umap_out$label = gsub('Tet',n_rows,umap_out$label)
  umap_out$label = gsub('Chlor',n_rows,umap_out$label)
  
  unique(umap_out$label)
  library(rcartocolor)
  #ggsave(paste(figs.out, 'umap_cluster_MG1655_total.png'))
  
  #umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'))
  ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
    geom_point( size =(1/i) * 1,alpha = 0.8) +  
    #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
    labs(x = "", y = "", color = "Condition") +
    ggtitle(sprintf("")) +
    #  coord_cartesian(xlim= c(0,1.0))+
    #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
    #  coord_equal(ratio=1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text = element_blank(),
          #        axis. = element_blank(),
          #        axis.line = element_line(colour = "black"),
          #        axis.text=element_text(size=20),axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank() ,
          axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
    #  scale_color_simpsons()s
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
    scale_color_manual(values=colors[i])
  ggsave(paste(figs.out, 'umap_MG1655_tet_chlor_sampled_scifi10_lognorm',n_rows,'.png'))
  gad_vec = c('EC_insI-2')
  gad_df = data.frame(gad_vec, rep(1,length(gad_vec)))
  names(gad_df) = c('gene_id','gad_stage')
  gad_df$gene_id = gsub('_','-',gad_df$gene_id)
  gad_score = matrix(data = 0,ncol=1,nrow = dim(sct@assays$RNA)[1])
  gene_indices = c()
  for (j in 1:length(gad_df$gene_id)){
    index = which(gad_df$gene_id[j] == rownames(sct@assays$RNA))
    if (length(index)>0){
      gene_indices = c(gene_indices,index)
      select = c(select,j)
    }
  }
  gad_df$gene_indices = gene_indices
  gad_score[gad_df$gene_indices,1] = as.numeric(gad_df$gad_stage)
  norm_data = t(as.data.frame(sct@assays$RNA@data))
  gad_data = data.frame(norm_data %*% gad_score)
  names(gad_data) = c('norm_gad')
  
  gad_data$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
  gad_data$label = as.character(sct[['Label']][,1])
  gad_data$ident = as.character(sct@active.ident)
  gad_data$index = 1:dim(gad_data)[1]
  umap_out2 = cbind(umap_out, gad_data$norm_gad)
  umap_out2$max_val = pmax(gad_data$norm_gad)
  
  ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `gad_data$norm_gad`, fill = `gad_data$norm_gad`)) + 
    #  geom_point(pch = 21, aes(fill = `gad_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
    geom_point( size =(1/i) * 1,alpha = 0.8) +  
    
    #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
    labs(y = "", x = "") +
    ggtitle(sprintf("")) +
#    scale_fill_gradient2(  low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) + 
#    scale_colour_gradient2(low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) +     theme_bw() +
    scale_fill_gradientn(  colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) + 
    scale_colour_gradientn(colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) +     theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text = element_blank(),
          #        axis. = element_blank(),
          #        axis.line = element_line(colour = "black"),
          #        axis.text=element_text(size=20),axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank() ,
          axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=4)))+
    labs(fill = "", colour = "") 
#  FeaturePlot(sct,features = c( 'EC-insI-2'), reduction = 'umap',order = TRUE, pt.size =(1/(2*i)) * 1 )& scale_colour_gradientn(colours = feature_colors)
  ggsave(paste(figs.out, 'featureplot_MG1655_tet_chlor_sampled_scifi_lognorm',n_rows,'.png'))
  
  #umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'))
#  ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
#    geom_point( size =(1/i) * 1,alpha = 0.8) +  
#    #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
#    labs(x = "", y = "", color = "Condition") +
#    ggtitle(sprintf("")) +
#    #  coord_cartesian(xlim= c(0,1.0))+
#    #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
#    #  coord_equal(ratio=1) +
#    theme_bw() +
#    theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
#          panel.background = element_blank(), 
#          axis.text = element_blank(),
#          #        axis. = element_blank(),
#          #        axis.line = element_line(colour = "black"),
#          #        axis.text=element_text(size=20),axis.text.x=element_blank(),
#          axis.ticks.x=element_blank(),
#          axis.text.y=element_blank(),
#          axis.ticks.y=element_blank() ,
#          axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
#    scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(2,7,9,2,3,8,12)])
#  
#    #  scale_color_simpsons()s
#    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
#    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
#    #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
#    scale_color_manual(values=colors[i])
  
  
  
}

gad_vec = c('EC_tfaR')
gad_df = data.frame(gad_vec, rep(1,length(gad_vec)))
names(gad_df) = c('gene_id','gad_stage')
gad_df$gene_id = gsub('_','-',gad_df$gene_id)
gad_score = matrix(data = 0,ncol=1,nrow = dim(sct@assays$RNA)[1])
gene_indices = c()
for (j in 1:length(gad_df$gene_id)){
  index = which(gad_df$gene_id[j] == rownames(sct@assays$RNA))
  if (length(index)>0){
    gene_indices = c(gene_indices,index)
    select = c(select,j)
  }
}
gad_df$gene_indices = gene_indices
gad_score[gad_df$gene_indices,1] = as.numeric(gad_df$gad_stage)
norm_data = t(as.data.frame(sct@assays$RNA@data))
gad_data = data.frame(norm_data %*% gad_score)
names(gad_data) = c('norm_gad')

gad_data$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
gad_data$label = as.character(sct[['Label']][,1])
gad_data$ident = as.character(sct@active.ident)
gad_data$index = 1:dim(gad_data)[1]
umap_out2 = cbind(umap_out, gad_data$norm_gad)
pca_out2 = cbind(pca_out, gad_data$norm_gad)
umap_out2$max_val = pmax(gad_data$norm_gad)
pca_out2$max_val = pmax(gad_data$norm_gad)

ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `gad_data$norm_gad`, fill = `gad_data$norm_gad`)) + 
  #  geom_point(pch = 21, aes(fill = `gad_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point( size =(1/i) * 1,alpha = 0.8) +  
  
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #    scale_fill_gradient2(  low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) + 
  #    scale_colour_gradient2(low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) +     theme_bw() +
  scale_fill_gradientn(  colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) + 
  scale_colour_gradientn(colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) +     theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=4)))+
  labs(fill = "", colour = "") 
#  FeaturePlot(sct,features = c( 'EC-insI-2'), reduction = 'umap',order = TRUE, pt.size =(1/(2*i)) * 1 )& scale_colour_gradientn(colours = feature_colors)
ggsave(paste(figs.out, 'featureplot_MG1655_tet_chlor_tfaR',n_rows,'.png'))

ggplot(pca_out2%>%arrange(`max_val`), aes_string(x = pcs[1], y=pcs[2])) + 
  #  geom_point(pch = 21, aes(fill = `gad_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point( size =(1/i) * 1,alpha = 0.8,aes(colour = `gad_data$norm_gad`, fill = `gad_data$norm_gad`)) +  
  
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #    scale_fill_gradient2(  low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) + 
  #    scale_colour_gradient2(low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) +     theme_bw() +
  scale_fill_gradientn(  colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) + 
  scale_colour_gradientn(colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) +     theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=4)))+
  labs(fill = "", colour = "") 
#  FeaturePlot(sct,features = c( 'EC-insI-2'), reduction = 'umap',order = TRUE, pt.size =(1/(2*i)) * 1 )& scale_colour_gradientn(colours = feature_colors)
ggsave(paste(figs.out, 'pca_featureplot_MG1655_tet_chlor_tfaR',n_rows,'.png'))


gad_vec = c('EC_pinR')
gad_df = data.frame(gad_vec, rep(1,length(gad_vec)))
names(gad_df) = c('gene_id','gad_stage')
gad_df$gene_id = gsub('_','-',gad_df$gene_id)
gad_score = matrix(data = 0,ncol=1,nrow = dim(sct@assays$RNA)[1])
gene_indices = c()
for (j in 1:length(gad_df$gene_id)){
  index = which(gad_df$gene_id[j] == rownames(sct@assays$RNA))
  if (length(index)>0){
    gene_indices = c(gene_indices,index)
    select = c(select,j)
  }
}
gad_df$gene_indices = gene_indices
gad_score[gad_df$gene_indices,1] = as.numeric(gad_df$gad_stage)
gad_data = data.frame(norm_data %*% gad_score)
names(gad_data) = c('norm_gad')

gad_data$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
gad_data$label = as.character(sct[['Label']][,1])
gad_data$ident = as.character(sct@active.ident)
gad_data$index = 1:dim(gad_data)[1]
umap_out2 = cbind(umap_out, gad_data$norm_gad)
pca_out2 = cbind(pca_out, gad_data$norm_gad)
umap_out2$max_val = pmax(gad_data$norm_gad)
pca_out2$max_val = pmax(gad_data$norm_gad)

ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `gad_data$norm_gad`, fill = `gad_data$norm_gad`)) + 
  #  geom_point(pch = 21, aes(fill = `gad_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point( size =(1/i) * 1,alpha = 0.8) +  
  
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #    scale_fill_gradient2(  low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) + 
  #    scale_colour_gradient2(low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) +     theme_bw() +
  scale_fill_gradientn(  colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) + 
  scale_colour_gradientn(colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) +     theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=4)))+
  labs(fill = "", colour = "") 
#  FeaturePlot(sct,features = c( 'EC-insI-2'), reduction = 'umap',order = TRUE, pt.size =(1/(2*i)) * 1 )& scale_colour_gradientn(colours = feature_colors)
ggsave(paste(figs.out, 'featureplot_MG1655_tet_chlor_pinR',n_rows,'.png'))

ggplot(pca_out2%>%arrange(`max_val`), aes_string(x = pcs[1], y=pcs[2])) + 
  #  geom_point(pch = 21, aes(fill = `gad_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point( size =(1/i) * 1,alpha = 0.8,aes(colour = `gad_data$norm_gad`, fill = `gad_data$norm_gad`)) +  
  
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #    scale_fill_gradient2(  low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) + 
  #    scale_colour_gradient2(low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) +     theme_bw() +
  scale_fill_gradientn(  colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) + 
  scale_colour_gradientn(colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) +     theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=4)))+
  labs(fill = "", colour = "") 
#  FeaturePlot(sct,features = c( 'EC-insI-2'), reduction = 'umap',order = TRUE, pt.size =(1/(2*i)) * 1 )& scale_colour_gradientn(colours = feature_colors)
ggsave(paste(figs.out, 'pca_featureplot_MG1655_tet_chlor_pinR',n_rows,'.png'))

kurts = kurtosi(data.frame(sct[['pca']]@cell.embeddings))
kurt_frame = data.frame(cbind(1:100, kurts))
names(kurt_frame) = c('PC','kurtosis')
ggplot(kurt_frame, aes(x = PC, y = kurtosis)) + 
  #  geom_point(pch = 21, aes(fill = `gad_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point( size = 2,alpha = 0.8) +  
  
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "Kurtosis", x = "Principal Component") +
  ggtitle(sprintf("")) +
  #    scale_fill_gradient2(  low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) + 
  #    scale_colour_gradient2(low = (brewer.pal(n = 9, name = "Greys"))[2], mid = 'yellow', high = 'red',limits = c(0, 5), oob = scales::squish) +     theme_bw() +
  scale_fill_gradientn(  colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) + 
  scale_colour_gradientn(colors = feature_colors,limits = c(0, 2.0), oob = scales::squish) +     theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(), 
#        axis.text = element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=4)))+
  labs(fill = "", colour = "") 
  #FeaturePlot(sct,features = c( 'EC-insI-2'), reduction = 'umap',order = TRUE, pt.size =(1/(2*i)) * 1 )& scale_colour_gradientn(colours = feature_colors)
ggsave(paste(figs.out, 'kurtosis_unordered',n_rows,'.svg'), device = 'svg')



#kurtosi(sct@assays$RNA@data[which(rownames(sct@assays$RNA@data) == 'EC-insI-2'),])
#x = kurtosi(t(as.matrix(sct@assays$RNA@data)))
#hist.data = hist(sct@assays$RNA@data[which(rownames(sct@assays$RNA@data) == 'EC-ydfK'),], plot = F)
#hist.data$counts = log10(hist.data$counts + 1)
#plot(hist.data, ylab='log2(Frequency)')

kurtosi_total = kurtosi_total[-1,]
kurtosi_total = data.frame(kurtosi_total)
names(kurtosi_total) = c('X1','kurtosi','label')

#kurtosi_total$label <- factor(kurtosi_total$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'))



kurtosi_total$kurtosi = as.numeric(as.character(kurtosi_total$kurtosi))
kurtosi_total$X1 = as.numeric(as.character(kurtosi_total$X1))
kurtosi_total2 = subset(kurtosi_total, X1 < 15)
#ggplot(data = subset(kurtosi_total, label > 1000), aes (x = X1, y = kurtosi, group = label, color = label)) +
kurtosi_total2 = kurtosi_total2[kurtosi_total2$label %in% c(1000,2500,5000,10000,15000,35000,55000,75000),]
ggplot(data =kurtosi_total2, aes (x = X1, y = kurtosi, group = label, color = label)) +
  geom_point(size = 0.8, alpha = 0.8) +
  geom_line(position = 'identity', size = 0.5, alpha = 1.8 ) + 
  labs(x = "Rank", y = "kurtosi") +
  ggtitle(sprintf("kurtosi for ribosome antibiotics")) +
  #  scale_y_log10()+
  theme_bw() +
  scale_color_gradient2(
    low = "skyblue1", 
    mid = "gray", 
    high = 'brown', 
    midpoint = (75000 + 1000)/2
  )+
  theme(plot.title = element_text(hjust = 0.5, size = 0.01, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 9), legend.title = element_text(size=0.01,face="bold"))# + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12)])
#ggsave(paste(figs.out, 'umap_MG1655_tet_scifi10_ramp_with_words','.png'))
ggsave(paste(figs.out, 'umap_MG1655_static_scifi10_ramp_with_words_lognormv2','.svg'), device = 'svg', dpi = 300)


ggplot(data =kurtosi_total2, aes (x = X1, y = kurtosi, group = label, color = label)) +
  geom_point(size = 1.3, alpha = 0.8) +
  geom_line(position = 'identity', size = 1.2, alpha = 1.8 ) + 
  labs(x = "Rank", y = "kurtosi") +
  ggtitle(sprintf("kurtosi for ribosome antibiotics")) +
  #  scale_y_log10()+
  theme_bw() +
  scale_color_gradient2(
    low = "skyblue1", 
    mid = "gray", 
    high = 'brown', 
    midpoint = (75000 + 1000)/2
  )+
  theme(plot.title = element_text(hjust = 0.5, size = 0.01, face="bold.italic"),
#        axis.title=element_text(size=0.05,face="bold"), 
        legend.text = element_text(size = 10), legend.title = element_text(size=0.01,face="bold"),
#        axis.text = element_blank(),
        panel.grid =element_blank())
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        #        axis.ticks.x=element_blank(),
#        axis.text.y=element_blank())
#        axis.ticks.y=element_blank() ,)# + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12)])
#ggsave(paste(figs.out, 'umap_MG1655_tet_scifi10_ramp_no_words','.png'))
ggsave(paste(figs.out, 'umap_MG1655_static_scifi10_ramp_with_words_lognormv3','.svg'), device = 'svg', dpi = 300)
ggsave(paste(figs.out, 'umap_MG1655_static_scifi10_ramp_with_words_lognormv3','.png'),  dpi = 300)

#ggsave(paste(figs.out, 'umap_MG1655_static_scifi10_ramp_no_words_lognormv2','.png'))
#library(ggpattern)
#### Erythro subsample ####
kurtosi_total2 = matrix(ncol = 3)
nrows =c(1000,2500,5000,7500,8027)
#nrows =c(500,1000)#,2500,5000,7500,10000,15000,20000)
nsamples = length(nrows)
#col1 = 
col1 = as.vector(rcartocolor::carto_pal(name='Safe'))[c(7)]

#colors = colorRampPalette(c("skyblue1", col1))(nsamples)
colors = colorRampPalette(c("skyblue1", 'gray','brown'))(nsamples)
for (i in 1:length(nrows)){
  keep_rows  = which(data2$Label =="Gent" | data2$Label =="Erythro" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
  n_rows = nrows[i]
  data = data2[keep_rows,]
  data = data[sample(dim(data)[1], n_rows,FALSE),]
  library(dplyr)
  #data = sample_n(data, 300)
  k = colSums(data!=0) 
  to_filter = which(k > 5)
  new_names = names(to_filter)
  new_names = c(new_names)
  data = data[new_names]
  #### Counts/gene
  #### Counts/gene
  m = rowSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
  counts = data.frame(m,data$Label)
  names(counts) = c('count','annotation')
  k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345','BS-trnO-Ala','BS-trnA-Ala'))])
  to_filter = which(k > 10)
  new_names = names(to_filter)
  new_names = c(new_names, "Label",'Identity')
  new_names = c(new_names)
  data = data[new_names]
  k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
  to_filter = which(k2 >= 10)
  data = data[to_filter,]
  ## Make Seurat
  
  sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala'))]
  rownames(sc_precursor) = 1:dim(sc_precursor)[1]
  sc_precursor = t(sc_precursor)
  sc_meta = data[which(names(data) %in% c('Label'))]
  sc = CreateSeuratObject(
    #  sc_precursor,
    sc_precursor,
    project = "SeuratProject",
    assay = "RNA",
    min.cells = 0,
    min.features = 0,
    names.field = 1,
    names.delim = "_",
    meta.data = sc_meta
  )
#  sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
  sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
  sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
  top10 <- head(VariableFeatures(sct), 10)
  all.genes <- rownames(sct)
  
  rm(mat)
  rm(sc_precursor)
  sct <- ScaleData(sct, features = all.genes)
  sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
  pca_out = data.frame(sct[['pca']]@cell.embeddings)
  kurtosi_total_tet = cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep(n_rows,dim(pca_out)[1]))
  kurtosi_total2 =rbind(kurtosi_total2, kurtosi_total_tet)
  
  pcs = names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]
  
  #sct = JackStraw(sct)
  
  print(sct[["pca"]], dims = 1:25, nfeatures = 5)
  
  DimPlot(sct, label = TRUE,label.size = 0.1,reduction = 'pca', cols =  colors[i],dims =as.numeric(gsub('PC_','',names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]))) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
    theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  #  ggsave(paste(figs.out, 'pca_pc39_tet.png'))
  ggsave(paste(figs.out, 'pca_MG1655_gent_erythro_sampled_scifi1039_lognorm',n_rows,'.png'))
  rm(mat)
  rm(sc_precursor)
  sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
  sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
  #### UMAP Variability ###
  library(ggsci)
  umap_out = data.frame(sct[['umap']]@cell.embeddings)
  umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
  umap_out$label = as.character(sct[['Label']][,1])
  umap_out$ident = as.character(sct@active.ident)
  #umap_out = umap_out[umap_out$label=='STATIONARY',]
  umap_out$label = gsub('_',' ',umap_out$label)
  unique(umap_out$label)
  library(rcartocolor)
  #ggsave(paste(figs.out, 'umap_cluster_MG1655_total.png'))
  
  #umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'))
  ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
    geom_point( size =(1/i) * 1,alpha = 0.8) +  
    #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
    labs(x = "", y = "", color = "Condition") +
    ggtitle(sprintf("")) +
    #  coord_cartesian(xlim= c(0,1.0))+
    #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
    #  coord_equal(ratio=1) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text = element_blank(),
          #        axis. = element_blank(),
          #        axis.line = element_line(colour = "black"),
          #        axis.text=element_text(size=20),axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank() ,
          axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 15), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))+
    scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(9,4,9,2,3,8,12)])
  
  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,5,6,8,7)])
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
  #    scale_color_manual(values=colors[i])
  ggsave(paste(figs.out, 'umap_MG1655_gent_erythro_sampled_scifi10_lognorm',n_rows,'.png'))
  
  
  
}


kurtosi_total = kurtosi_total[-1,]
kurtosi_total2 = kurtosi_total2[-1,]
kurtosi_total = data.frame(kurtosi_total)
kurtosi_total2 = data.frame(kurtosi_total2)
kurtosi_total$moa = 'Static'
kurtosi_total2$moa = 'Cidal'
names(kurtosi_total) = c('X1','kurtosi','label','moa')
names(kurtosi_total2) = c('X1','kurtosi','label','moa')

kurtosi_total3 = rbind(kurtosi_total, kurtosi_total2)

#kurtosi_total$label <- factor(kurtosi_total$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential','Lambda30','Lambda90'))



kurtosi_total3$kurtosi = as.numeric(as.character(kurtosi_total3$kurtosi))
kurtosi_total3$X1 = as.numeric(as.character(kurtosi_total3$X1))
kurtosi_total4 = subset(kurtosi_total3, X1 < 15)
kurtosi_total5 = subset(kurtosi_total4, label < 10000)
kurtosi_total5 = subset(kurtosi_total5, label > 8000)
kurtosi_total6 = subset(kurtosi_total4,   moa == 'Static')



ggplot(data =kurtosi_total5, aes (x = X1, y = kurtosi, group = moa, color = moa)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_line(position = 'identity', size = 1.5, alpha = 1.8 ) + 
  labs(x = "Rank", y = "kurtosi") +
  ggtitle(sprintf("kurtosi for ribosome antibiotics")) +
  #  scale_y_log10()+
  theme_bw() +
#  scale_color_gradient2(
#    low = "skyblue1", 
#    mid = "gray", 
#    high = 'brown', 
#    midpoint = (max(nrows) + min(nrows))/2
#  )+
  theme(plot.title = element_text(hjust = 0.5, size = 0.01, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 9), legend.title = element_text(size=0.01,face="bold"))# + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::car
ggsave(paste(figs.out,'kurtosis_static_cidal_lognorm','.png'))



ggplot(data =kurtosi_total5, aes (x = X1, y = kurtosi, group = moa, color = moa)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_line(position = 'identity', size = 1.2, alpha = 1.8 ) + 
  labs(x = "Rank", y = "kurtosi") +
  ggtitle(sprintf("kurtosi for ribosome antibiotics")) +
  #  scale_y_log10()+
  theme_bw() +
  #  scale_color_gradient2(
  #    low = "skyblue1", 
  #    mid = "gray", 
  #    high = 'brown', 
  #    midpoint = (max(nrows) + min(nrows))/2
  #  )+
  theme(plot.title = element_text(hjust = 0.5, size = 0.01, face="bold.italic"),
        axis.title=element_text(size=0.05,face="bold"), legend.text = element_text(size = 0.5), legend.title = element_text(size=0.01,face="bold"),
        axis.text = element_blank(),
        panel.grid =element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        #        axis.ticks.x=element_blank(),
        axis.text.y=element_blank())+scale_color_manual(values = c('black','gray'))
ggsave(paste(figs.out,'kurtosis_static_cidal_lognorm_nowords','.png'))


DimPlot(sct, label = TRUE,label.size = 0.1,reduction = 'pca', group.by = 'Label', dims =as.numeric(gsub('PC_','',names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]))) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

ggplot(data =kurtosi_total6, aes (x = X1, y = kurtosi, group = label, color = label)) +
  geom_point(size = 0.8, alpha = 0.8) +
  geom_line(position = 'identity', size = 0.5, alpha = 1.8 ) + 
  labs(x = "Rank", y = "kurtosi") +
  ggtitle(sprintf("kurtosi for ribosome antibiotics")) +
  #  scale_y_log10()+
  theme_bw() +
  scale_color_gradient2(
    low = "skyblue1", 
    mid = "gray", 
    high = 'brown', 
#    midpoint = (max(nrows) + min(nrows))/2
    midpoint = (30000 + min(nrows))/2
  )+
  theme(plot.title = element_text(hjust = 0.5, size = 0.01, face="bold.italic"),
        axis.title=element_text(size=10,face="bold"), legend.text = element_text(size = 9), legend.title = element_text(size=0.01,face="bold"))# + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12)])
ggsave(paste(figs.out,'kurtosis_ribosomes_with_words_lognorm','.png'))


ggplot(data =kurtosi_total6, aes (x = X1, y = kurtosi, group = label, color = label)) +
  geom_point(size = 0.8, alpha = 0.8) +
  geom_line(position = 'identity', size = 0.5, alpha = 1.8 ) + 
  labs(x = "Rank", y = "kurtosi") +
  ggtitle(sprintf("kurtosi for ribosome antibiotics")) +
  #  scale_y_log10()+
  theme_bw() +
  scale_color_gradient2(
    low = "skyblue1", 
    mid = "gray", 
    high = 'brown', 
    #    midpoint = (max(nrows) + min(nrows))/2
    midpoint = (30000 + min(nrows))/2
  )+
  theme(plot.title = element_text(hjust = 0.5, size = 0.01, face="bold.italic"),
        axis.title=element_text(size=0.05,face="bold"), legend.text = element_text(size = 0.5), legend.title = element_text(size=0.01,face="bold"),
        axis.text = element_blank(),
        panel.grid =element_blank(),
        #        axis. = element_blank(),
        #        axis.line = element_line(colour = "black"),
        #        axis.text=element_text(size=20),axis.text.x=element_blank(),
        #        axis.ticks.x=element_blank(),
        axis.text.y=element_blank())
#        axis.ticks.y=element_blank() ,)# + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=4)))
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,5,6,12)])
#  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,9,2,12)])

ggsave(paste(figs.out,'kurtosis_ribosomes_no_words_lognorm','.png'))



FeaturePlot(sct,features = c('Lambda-NinB', 'Lambda-A'), reduction = 'umap',order = TRUE)

