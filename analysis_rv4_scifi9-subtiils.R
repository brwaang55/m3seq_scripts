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
library(countreg)
figs.out = '~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/scifiseq9_figs//'
dir.create(figs.out)
#Colors
Breaks=seq(0,60,1)
#dir.create('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/slides/Figures/2020-07-20/')
Colors=rev(brewer.pal(11,"Spectral"))
colors=colorRampPalette(Colors)(120)
#Exponential + Stationary
#data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/petri_seq_with_anti.csv')
data = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/expression_filtered_scifi9_BS168_USA300_25.csv')
cols = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/gene_index_scifi9_BS168_USA300_25.csv')
annotation = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/cell_index_scifi9_BS168_USA300_25.csv')
#colnames(data) = cols[,1]annotation = merge(annotation,supp_annotation,by.x = 'r1',by.y = 'BARCODE_1')
annotation$treatment = annotation$treatment
annotation = annotation[order(annotation$cell_index),]
#colnames(data) = cols[,1]
annotation$treatment = gsub('Lib1_HB_Tn5_BS168_STAPH_','',annotation$treatment)

data['Label'] = annotation$treatment
data['Identity'] = annotation$identity
#data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/Data/BS168_genes_no_ribo_stationary.csv')
"Label" %in% names(data)
#data = data[,-1]
#data = data[,-(dim(data)[2])]
cols = names(data)

#### Past ####
data_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/expression_filtered_scifi10_BS168_USA300_25.csv')
cols_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/gene_index_scifi10_BS168_USA300_25.csv')
annotation_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/cell_index_scifi10_BS168_USA300_25.csv')
annotation_past$treatment = gsub('Lib2_HB_Tn5_BS168_STAPH_','',annotation_past$treatment)
data_past['Label'] = annotation_past$treatment
data_past['Identity'] = annotation_past$identity
#annotation_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/cell_index_scifi7_len40.csv')
#annotation_past$treatment = str_replace(annotation_past$treatment,'2HR_CEFOZOLIN','2hr_Cipro')
#annotation_past$treatment = str_replace(annotation_past$treatment,'2HR_CIPRO','2hr_Cef')
#annotation_past$treatment = str_replace(annotation_past$treatment,'4HR_CEFOZOLIN','4hr_Cef')
#colnames(data) = cols[,1]
#annotation_past$treatment = gsub('Lib1_HB_Tn5_BS168_STAPH_','Past_',annotation_past$treatment)
data_past['Label'] = annotation_past$treatment
data_past['Identity'] = annotation_past$identity
#data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/Data/subtilis_genes_no_ribo_stationary.csv')
cols_past = names(data_past)
# Make minimum 5 observations
### Make this only E. BS168 ####
new_names = names(data)
ec_names = new_names[str_detect(new_names,'BS_')]
data = data[ec_names]
data_past = data_past[ec_names]

#data_past = data_past[ec_names]
data['Label'] = annotation$treatment
data['Identity'] = annotation$identity
data_past['Label'] = annotation_past$treatment
data_past['Identity'] = annotation_past$identity
#keep_rows  = which(data_past$Identity == 'coli')
#keep_rows  = which(data$Identity == 'coli' & data$Label !="CRISPRI" & data$Label != "M9"& data$Label != "30MIN_FIX")
#data_past = data_past[keep_rows,]
#data_total = rbind(data,data_past)
data_total = rbind(data,data_past)
data_past = 9
#data_total = rbind(data)
#data_total = rbind(data_past)
data2 = data_total
data_total = 9
#### Make this only  B. subtilis ####
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
kurtosi_total= data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Exponential',dim(pca_out)[2])))
plot(1:50, kurtosi(pca_out))
ggplot(data = data.frame(cbind(1:length(kurtosi_total), kurtosi_total)), aes (x = V1, y = kurtosi_total)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_BS168_combined_var_explained','.png'))

plot(1:50, kurtosi(pca_out))
#ggsave(paste(figs.out,'pca_BS168_combined_var_explained','.png'))
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
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1,logfc.threshold = 0.1)
DimPlot(sct, label = TRUE) + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,label.size = 0.1,reduction = 'pca', dims =as.numeric(gsub('PC_','',names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]))) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#  ggsave(paste(figs.out, 'pca_pc39_tet.png'))
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_BS168_combined_scifi3','.png'))
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
#ggsave(paste(figs.out, 'umap_BS168_combined_clean_palette1_30pcs.png'))

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

ggsave(paste(figs.out, 'umap_BS168_combined_clean_palette2_40_pcs_exponential.png'))


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
ggsave(paste(figs.out, 'umap_BS168_combined_cluster_exponential_phase_scifi9.png'))

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
DimPlot(sct, label = TRUE,label.size = 0.1,reduction = 'pca', ,dims =as.numeric(gsub('PC_','',names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]))) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#  ggsave(paste(figs.out, 'pca_pc39_tet.png'))
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.1,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.1,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1,logfc.threshold = 0.1)
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
ggsave(paste(figs.out, 'umap_BS168_combined_gent_phase_scifi9.png'))
#ggsave(paste(figs.out, 'umap_BS168_combined_clean_palette2.png'))
#scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,10,5)])
#ggsave(paste(figs.out, 'umap_BS168_combined_clean_palette1.png'))

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
ggsave(paste(figs.out, 'umap_BS168_combined_cluster_gent_phase_scifi9.png'))


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
DimPlot(sct, label = TRUE,label.size = 0.1,reduction = 'pca', dims =as.numeric(gsub('PC_','',names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]))) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
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
cluster11.markers <- FindMarkers(sct, ident.1 = 11, min.pct = 0.1, logfc.threshold = 0.1,ident.2 = 3)
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
ggsave(paste(figs.out, 'umap_BS168_combined_chlor_scifi9.png'))



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
ggsave(paste(figs.out, 'umap_BS168_combined_cluster_chlor_scifi9.png'))
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
#ggsave(paste(figs.out,'BS168_combined_heatmap_exp_3hr','.png'))
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
#ggsave(paste(figs.out, 'umap_cluster_BS168_combined_total.png'))
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
ggsave(paste(figs.out, 'umap_BS168_combined_cluster_tet_scifi9.png'))

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
ggsave(paste(figs.out, 'umap_BS168_combined_tet_scifi9.png'))



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
#ggsave(paste(figs.out,'pca_BS168_combined_var_explained','.png'))

ggplot(data = var_total, aes (x = V1, y = varExplained, col = V3)) +
  geom_line() + 
  geom_point(size = 1) +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_BS168_combined_var_explained','.png'))

rm(mat)
rm(sc_precursor)


pca_out = data.frame(sct[['pca']]@cell.embeddings)
kurtosi_total_nal = data.frame(cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep('Nal',length(varExplained))))

kurtosi_total = data.frame(rbind(kurtosi_total, kurtosi_total_nal))
kurtosi_total$kurtosis = as.numeric(as.character(kurtosi_total$X2))
kurtosi_total$X1 = as.numeric(as.character(kurtosi_total$X1))
ggplot(data = kurtosi_total, aes (x = X1, y = kurtosis, col = X3)) +
  geom_line() + 
  geom_point(size = 1) +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  #  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_BS168_combined_var_explained','.png'))
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
ggsave(paste(figs.out, 'umap_BS168_combined_cluster_nal_scifi9.png'))
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
ggsave(paste(figs.out, 'umap_BS168_combined_nalidixic_acid.png'))




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


DimPlot(sct, label = TRUE,label.size = 0.1,reduction = 'pca', dims =as.numeric(gsub('PC_','',names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]))) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#  ggsave(paste(figs.out, 'pca_pc39_tet.png'))

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
#ggsave(paste(figs.out,'BS168_combined_heatmap_exp_3hr','.png'))
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
ggsave(paste(figs.out, 'umap_BS168_combined_cluster_cipro_scifi9.png'))
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
ggsave(paste(figs.out, 'umap_BS168_combined_cipro.png'))


#### Cefazolin  ####

unique(data2$Label)
#keep_rows  = which(data2$Label =="4HR_COCULTURE" | data2$Label =="STATIONARY" | data2$Label == "Coculture_EXP"| data2$Label == "EXP")
#keep_rows  = which( data2$Label == "Cef"| data2$Label == "Cyclo" | data2$Label =="EXP" | data2$Label =="STATIONARY")
keep_rows  = which( data2$Label == "Cef"| data2$Label == "Cyclao" | data2$Label =="EXP" | data2$Label =="STATIONARY")
#keep_rows  = which(data2$Identity == 'BS168_combined')
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
#ggsave(paste(figs.out,'pca_BS168_combined_var_explained','.png'))

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
#ggsave(paste(figs.out,'umap_BS168_combined_scifi7','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scBactoSeq/Figures/2020-05-12/UMAP_Relative_Count.png')
#DimPlot(sct, label = TRUE, reduction = 'tsne')
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_BS168_combined_scifi3','.png'))
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
ggsave(paste(figs.out, 'umap_BS168_combined_cluster_cefazin.png'))
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
#ggsave(paste(figs.out, 'umap_BS168_combined_clean_palette2.png'))

ggsave(paste(figs.out, 'umap_BS168_combined_cef_scifi9.png'))
#### Cyclo ###

unique(data2$Label)
#keep_rows  = which(data2$Label =="4HR_COCULTURE" | data2$Label =="STATIONARY" | data2$Label == "Coculture_EXP"| data2$Label == "EXP")
#keep_rows  = which( data2$Label == "Cyclo"| data2$Label == "Cyclo" | data2$Label =="EXP" | data2$Label =="STATIONARY")
keep_rows  = which( data2$Label == "Cyclo"| data2$Label == "Cyclao" | data2$Label =="EXP" | data2$Label =="STATIONARY")
#keep_rows  = which(data2$Identity == 'BS168_combined')
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
#ggsave(paste(figs.out,'umap_BS168_combined_scifi7','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scBactoSeq/Figures/2020-05-12/UMAP_Relative_Count.png')
#DimPlot(sct, label = TRUE, reduction = 'tsne')
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_BS168_combined_scifi3','.png'))
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
ggsave(paste(figs.out, 'umap_BS168_combined_cluster_cycloserine_scifi9.png'))
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
#ggsave(paste(figs.out, 'umap_BS168_combined_clean_palette2.png'))


ggsave(paste(figs.out, 'umap_BS168_combined_cyclo_scifi9.png'))


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
ggsave(paste(figs.out, 'umap_BS168_combined_cluster_erythro_scifi9.png'))

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
ggsave(paste(figs.out, 'umap_BS168_combined_erythro.png'))

library(dplyr)

av.exp.sct <- AverageExpression(sct, group.by = 'ident',return.seurat = TRUE)
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top10
#markers = markers[1:8,]
genes = markers$gene
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
g = merge(g, markers)
#g = g[which(g$variable%in% genes),]
g = g[order(g$cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'BS-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = g$gene,labels = g$gene)
#g = g[order(g$cluster),]

library(myriad)
ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile() +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.5), oob = scales::squish)+
  xlab('Gene')+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11), panel.grid =  element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=11,color = 'black'),
        axis.title=element_text(size=13), legend.text = element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(size=11)) + scale_colour_hue(l=40)

####
# Create a weighted plot.. 
sct2 = sct
phage_anno = data.frame(read.csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/scBactoSeq/data/spbeta_genes.csv'))
phage_anno$gene_id = gsub('_','-',phage_anno$gene_id)
phage_anno$sp_st1 = as.numeric(phage_anno$spbeta_stage==1)
phage_anno$sp_st2 = as.numeric(phage_anno$spbeta_stage==2)
phage_anno$sp_st3 = as.numeric(phage_anno$spbeta_stage==3)
phage_anno$pb_st1 = as.numeric(phage_anno$pbs_stage==1)
phage_anno$pb_st2 = as.numeric(phage_anno$pbs_stage==2)
phage_anno$pb_st3 = as.numeric(phage_anno$pbs_stage==3)
gene_indices = c()
select = c()
for (i in 1:length(phage_anno$gene_id)){
  index = which(phage_anno$gene_id[i] == rownames(sct2@assays$RNA))
  if (length(index)>0){
    gene_indices = c(gene_indices,index)
    select = c(select,i)
  }
}
phage_anno2 = phage_anno[select,]
phage_anno2$gene_indices = gene_indices
phage_score = matrix(data = 0,ncol=2,nrow = dim(sct2@assays$RNA)[1])
phage_score[phage_anno2$gene_indices,1] = as.numeric(phage_anno2$spbeta_score)
phage_score[phage_anno2$gene_indices,2] = as.numeric(phage_anno2$pbsx_score)
phage_score2 = matrix(data = 0,ncol=8,nrow = dim(sct2@assays$RNA)[1])
phage_score2[phage_anno2$gene_indices,1] = phage_anno2$sp_st1
phage_score2[phage_anno2$gene_indices,2] = phage_anno2$sp_st2
phage_score2[phage_anno2$gene_indices,3] = phage_anno2$sp_st3
phage_score2[phage_anno2$gene_indices,4] = phage_anno2$pb_st1
phage_score2[phage_anno2$gene_indices,5] = phage_anno2$pb_st2
phage_score2[phage_anno2$gene_indices,6] = phage_anno2$pb_st3
phage_score2[phage_anno2$gene_indices,7] = phage_anno2$spbeta_stage
phage_score2[phage_anno2$gene_indices,8] = phage_anno2$pbs_stage
norm_data = t(as.data.frame(sct2@assays$RNA@data))
#phage_data = data.frame(norm_data %*% phage_score)
#phage_data2 = data.frame(norm_data %*% phage_score2)
phage_data = data.frame(t(as.matrix(sct2@assays$RNA@counts)) %*% phage_score)
phage_data2 = data.frame(t(as.matrix(sct2@assays$RNA@counts)) %*% phage_score2)
phage_data3 = data.frame(norm_data %*% phage_score)
divide1 = phage_data[,1]
divide2 = phage_data[,2]
divide1[which(divide1 < 2)] = 0 
divide2[which(divide2 < 2)] = 0 
phage_data2[,7] = phage_data2[,7]/divide1
phage_data2[,8] = phage_data2[,8]/divide2
is.na(phage_data2)<-sapply(phage_data2, is.infinite)
phage_data2[is.na(phage_data2)] = 0
#### Cell Stage ####
names(phage_data2) = c('sp_st1','sp_st2','sp_st3','pb_st1','pb_st2','pb_st3','sp_stage','pb_stage')
names(phage_data3) = c('norm_sp','norm_pbs')
phage_data = cbind(phage_data, phage_data2)
phage_data = cbind(phage_data, phage_data3)

phage_data$umi_counts = as.numeric(as.character(sct2[['nCount_RNA']][,1]))
phage_data$label = as.character(sct2[['Label']][,1])
phage_data$ident = as.character(sct2@active.ident)
phage_data$index = 1:dim(phage_data)[1]
phage_data$norm_pbs2 = phage_data$norm_pbs/max(phage_data$norm_pbs)
phage_data$norm_sp2 = phage_data$norm_sp/max(phage_data$norm_sp)
#phage_data$x1_stage = phage_data2$X1
#phage_data$x2_stage = phage_data2$X2

library(dplyr)
#ggplot(phage_data%>%arrange(sp_stage),aes(x = `X1`, y = `X2`)) + 
p1 = ggplot(phage_data%>%arrange(pb_stage),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = pb_stage),size =3 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +#+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "PbsX stage") 


p2 = ggplot(phage_data%>%arrange(sp_stage),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = sp_stage),size =3 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-0,0.20) + ylim(-0,0.20)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Sp\U03B2 stage") 

p2
aggregate(norm_sp~ident,phage_data,min)
aggregate(norm_pbs~ident,phage_data,min)
pbsx=5.88
spb=18.86
#pbsx=0
#spb=0
length(which(phage_data$norm_pbs > pbsx))/dim(phage_data)[1]
length(which(phage_data$norm_sp > spb))/dim(phage_data)[1]
length(intersect(which(phage_data$norm_pbs > pbsx),which(phage_data$norm_sp > spb)))/dim(phage_data)[1]
length(which(phage_data$norm_pbs > pbsx))/dim(phage_data)[1]*length(which(phage_data$norm_sp > spb))/dim(phage_data)[1]
p1 + p2
#ggsave(paste(figs.out,'stage_plots','.png'))
phage_data
p3 = ggplot(phage_data,aes(x = `norm_sp`, y = `norm_pbs`,fill = `ident`)) + 
  geom_point(pch = 21, size =2.5, alpha = 0.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
#  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_color_manual(values=as.vector(palette.colors()))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Cluster") 

p3
ggsave(paste(figs.out,'phages_overlaid_cluster_nal','.png'))
#Use 5th percentile
quantile(phage_data[phage_data$ident == 2,]$norm_pbs,probs = seq(0,0.1,0.01))
quantile(phage_data[phage_data$ident == 1,]$norm_sp,probs = seq(0,0.1,0.01))
is_pbsx = which(phage_data$norm_pbs > pbsx)
is_spb = which(phage_data$norm_sp > spb)
is_both = intersect(is_pbsx, is_spb)
label = rep('Not Induced',dim(phage_data)[1])
label[is_pbsx] = 'PbsX'
label[is_spb] = 'Sp\U03B2'
label[is_both] = 'Both'
phage_data$phage = label


ggplot(phage_data,aes(x = `norm_sp`, y = `norm_pbs`,fill = `phage`, color = `phage`)) + 
#  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "", y = "") +
#  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_fill_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  scale_color_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 

ggsave(paste(figs.out,'phages_overlaid_coloring_nal','.png'))


p4 = ggplot(phage_data,aes(x = `norm_pbs2`, y = `norm_sp2`)) + 
  geom_point(pch = 21, aes(fill = pca_out$label),size =2 ) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
#  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-0,0.40) + ylim(-0,0.40)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(fill = "Louvain Cluster") 


p4
#ggsave(paste(figs.out,'normed_phages_overlaid_cluster','.png'))

p3 = ggplot(phage_data,aes(x = `norm_pbs`, y = `norm_sp`,fill = umi_counts)) + 
  geom_point(pch = 21, size =2 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +#+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Cluster") 


ggplot(phage_data%>%arrange(sp_stage),aes(x = `X1`, y = `X2`)) + 
  geom_point(pch = 21, aes(fill = sp_stage),size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
#  geom_density_2d()+
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw()+ #   xlim(-1,250) + ylim(-1,30)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
      axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  
  

ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_stage`, y = `pb_stage`)) + 
  geom_point(pch = 21,size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
#  geom_density_2d()+
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw()+   #+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  



#ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_st2`, y = `pb_st2`)) + 
#ggplot(phage_data%>%arrange(pb_st2),aes(x = `pb_st2`, y = `pb_st3`)) + 
ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_st2`, y = `sp_st3`)) + 
  geom_point(pch = 21, aes(fill = pb_st2),size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-1,10) + ylim(-1,10)+
theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
      axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  



ggplot(phage_data,aes(x = `X1`, y = `X2`,col = `ident`)) + 
  geom_point( size = 1.5) +  labs(x = "Sp\U03B2", y = "PbsX", color = "Class", parse = TRUE) +
  ggtitle(sprintf("Phage-scored cells")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))


ggplot(phage_data,aes(x = `X1`, y = `X2`,col = `ident`)) + 
  geom_point( size = 1) +  labs(x = "Sp\U03B2", y = "PbsX", color = "Class", parse = TRUE) +
  ggtitle(sprintf("Phage-scored cells")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-10,4) + ylim(-2,4)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
####

library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out2 = cbind(umap_out, phage_data$norm_pbs,phage_data$norm_sp)
umap_out2$max_val = pmax(phage_data$norm_pbs,phage_data$norm_sp)
umap_out2$phage = phage_data$phage
temp_phage = matrix(phage_data[c('norm_pbs','norm_sp')],ncol=2)
colnames(temp_phage) = paste0('PHAGE_',1:2)
sct[["phage"]] <- CreateDimReducObject(embeddings = temp_phage, key = "PHAGE_", assay = DefaultAssay(sct))

look_into = intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_pbs2 > 0.02))
look_into = union(which(phage_data$norm_sp2 > 0.05),which(phage_data$norm_pbs2 > 0.05))
look_into1 = which(phage_data$norm_sp2 > 0.5)
look_into2 = intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_sp2 < 0.5))
look_into3 = intersect(intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_pbs2 > 0.02)),intersect(which(phage_data$norm_sp2 > 0.05),which(phage_data$norm_sp2 < 0.5)))
look_into4 =  which(phage_data$norm_pbs2 > 0.1)
levels(sct@active.ident) = c(levels(sct@active.ident),"6","7",'8','9','30')
#sct@active.ident[look_into1] = "6"
#sct@active.ident[look_into2] = "7"
#sct@active.ident[look_into3] = "8"
#sct@active.ident[look_into4] = "9"
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1,ident.2 = 9)
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1,ident.2 = 8)
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1,ident.2 = 7)
#cluster4.markers2 <- FindMarkers(sct, ident.1 =6, min.pct = 0.05)
#up = intersect(rownames(cluster4.markers1[cluster4.markers1$avg_log2FC > 0,]),rownames(cluster4.markers2[cluster4.markers2$avg_log2FC > 0,]))
#down = intersect(rownames(cluster4.markers1[cluster4.markers1$avg_log2FC < 0,]),rownames(cluster4.markers2[cluster4.markers2$avg_log2FC < 0,]))
#sct@assays$RNA
#cluster4.markers1[up,]
#DimPlot(sct,reduction='pca')
phage_data5 = cbind(phage_data, norm_data)
ggplot(phage_data5%>%arrange(`BS-xtmA`),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = `BS-xtmA`),size =3 ) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
    scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-0,0.40) + ylim(-0,0.40)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(fill = "rnjB expression") 
ggsave(paste(figs.out,'rnjB_phageplot','.png'))

library(dplyr)
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `phage_data$norm_pbs`, fill = `phage_data$norm_pbs`)) + 
#  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5) +
#  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
#  scale_fill_gradientn(colors =  (brewer.pal(n = 9, name =  "BuGn") ), guide = "colourbar") + 
#  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "BuGn")), guide = '' ) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) +   #scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greens")[1]), high = (brewer.pal(n = 9, name = "Greens")[9])) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
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
ggsave(paste(figs.out,"pbsx_percentage_umap_nal",'.png'))
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `phage_data$norm_pbs`, fill = `phage_data$norm_pbs`)) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5) +
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #  scale_fill_gradientn(colors =  (brewer.pal(n = 9, name =  "BuGn") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "BuGn")), guide = '' ) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) +   #scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greens")[1]), high = (brewer.pal(n = 9, name = "Greens")[9])) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
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
#ggsave(paste(figs.out,"pbsx_percentage_umap_40pcs",'.png'))
#ggsave(paste(figs.out,"pbsx_percentage_umap_40pcs_cipro_only",'.png'))
library(prismatic)
new_col = prismatic::clr_lighten('red2',shift=0.1)
#new_col = prismatic::clr_saturate('red3',shift=0.1)
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2` )) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5, aes( colour = `phage_data$norm_sp`, fill = `phage_data$norm_sp`) )+
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #scale_colour_gradientn(low = (brewer.pal(n = 9, name = "Reds")[1]), high = (brewer.pal(n = 9, name = "Reds")[9]), guide = 'colourbar') + 
  #scale_fill_gradientn(colors =  (brewer.pal(n = 9, name = "Reds") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "Reds")), guide = '' ) + 
#  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
#  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red') + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red') + 
#  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
#  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
#    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
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
  labs(fill = "", colour = '') 


ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2` )) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5, aes( colour = `phage_data$norm_sp`, fill = `phage_data$norm_sp`) )+
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #scale_colour_gradientn(low = (brewer.pal(n = 9, name = "Reds")[1]), high = (brewer.pal(n = 9, name = "Reds")[9]), guide = 'colourbar') + 
  #scale_fill_gradientn(colors =  (brewer.pal(n = 9, name = "Reds") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "Reds")), guide = '' ) + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red') + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red') + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
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
  labs(fill = "", colour = '') 

#ggsave(paste(figs.out,"sp_percentage_umap_v2",'.png'))
ggsave(paste(figs.out,"sp_percentage_umap_nal",'.png'))
#ggsave(paste(figs.out,"sp_percentage_umap_30pcs",'.png'))

ggplot(umap_out2,aes(x = `UMAP_1`, y = `UMAP_2`,fill = `phage`, color = `phage`)) + 
  #  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "", y = "") +
  #  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
#  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_fill_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  scale_color_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        #        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        #        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 
ggsave(paste(figs.out,"nal_umap_phage_plot",'.png'))

####
# Create a weighted plot.. 
sct2 = sct
phage_anno = data.frame(read.csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/scBactoSeq/data/spbeta_genes.csv'))
phage_anno$gene_id = gsub('_','-',phage_anno$gene_id)
phage_anno$sp_st1 = as.numeric(phage_anno$spbeta_stage==1)
phage_anno$sp_st2 = as.numeric(phage_anno$spbeta_stage==2)
phage_anno$sp_st3 = as.numeric(phage_anno$spbeta_stage==3)
phage_anno$pb_st1 = as.numeric(phage_anno$pbs_stage==1)
phage_anno$pb_st2 = as.numeric(phage_anno$pbs_stage==2)
phage_anno$pb_st3 = as.numeric(phage_anno$pbs_stage==3)
gene_indices = c()
select = c()
for (i in 1:length(phage_anno$gene_id)){
  index = which(phage_anno$gene_id[i] == rownames(sct2@assays$RNA))
  if (length(index)>0){
    gene_indices = c(gene_indices,index)
    select = c(select,i)
  }
}
phage_anno2 = phage_anno[select,]
phage_anno2$gene_indices = gene_indices
phage_score = matrix(data = 0,ncol=2,nrow = dim(sct2@assays$RNA)[1])
phage_score[phage_anno2$gene_indices,1] = as.numeric(phage_anno2$spbeta_score)
phage_score[phage_anno2$gene_indices,2] = as.numeric(phage_anno2$pbsx_score)
phage_score2 = matrix(data = 0,ncol=8,nrow = dim(sct2@assays$RNA)[1])
phage_score2[phage_anno2$gene_indices,1] = phage_anno2$sp_st1
phage_score2[phage_anno2$gene_indices,2] = phage_anno2$sp_st2
phage_score2[phage_anno2$gene_indices,3] = phage_anno2$sp_st3
phage_score2[phage_anno2$gene_indices,4] = phage_anno2$pb_st1
phage_score2[phage_anno2$gene_indices,5] = phage_anno2$pb_st2
phage_score2[phage_anno2$gene_indices,6] = phage_anno2$pb_st3
phage_score2[phage_anno2$gene_indices,7] = phage_anno2$spbeta_stage
phage_score2[phage_anno2$gene_indices,8] = phage_anno2$pbs_stage
norm_data = t(as.data.frame(sct2@assays$RNA@data))
#phage_data = data.frame(norm_data %*% phage_score)
#phage_data2 = data.frame(norm_data %*% phage_score2)
phage_data = data.frame(t(as.matrix(sct2@assays$RNA@counts)) %*% phage_score)
phage_data2 = data.frame(t(as.matrix(sct2@assays$RNA@counts)) %*% phage_score2)
phage_data3 = data.frame(norm_data %*% phage_score)
divide1 = phage_data[,1]
divide2 = phage_data[,2]
divide1[which(divide1 < 2)] = 0 
divide2[which(divide2 < 2)] = 0 
phage_data2[,7] = phage_data2[,7]/divide1
phage_data2[,8] = phage_data2[,8]/divide2
is.na(phage_data2)<-sapply(phage_data2, is.infinite)
phage_data2[is.na(phage_data2)] = 0
#### Cell Stage ####
names(phage_data2) = c('sp_st1','sp_st2','sp_st3','pb_st1','pb_st2','pb_st3','sp_stage','pb_stage')
names(phage_data3) = c('norm_sp','norm_pbs')
phage_data = cbind(phage_data, phage_data2)
phage_data = cbind(phage_data, phage_data3)

phage_data$umi_counts = as.numeric(as.character(sct2[['nCount_RNA']][,1]))
phage_data$label = as.character(sct2[['Label']][,1])
phage_data$ident = as.character(sct2@active.ident)
phage_data$index = 1:dim(phage_data)[1]
phage_data$norm_pbs2 = phage_data$norm_pbs/max(phage_data$norm_pbs)
phage_data$norm_sp2 = phage_data$norm_sp/max(phage_data$norm_sp)
#phage_data$x1_stage = phage_data2$X1
#phage_data$x2_stage = phage_data2$X2

library(dplyr)
#ggplot(phage_data%>%arrange(sp_stage),aes(x = `X1`, y = `X2`)) + 
p1 = ggplot(phage_data%>%arrange(pb_stage),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = pb_stage),size =3 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +#+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "PbsX stage") 


p2 = ggplot(phage_data%>%arrange(sp_stage),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = sp_stage),size =3 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-0,0.20) + ylim(-0,0.20)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Sp\U03B2 stage") 

p2
length(which(phage_data$norm_pbs > 3.7))/dim(phage_data)[1]
length(which(phage_data$norm_sp > 8.13))/dim(phage_data)[1]
aggregate(norm_sp~ident,phage_data,min)
aggregate(norm_pbs~ident,phage_data,min)
pbsx = 19.609
spb = 11.634
length(intersect(which(phage_data$norm_pbs > pbsx),which(phage_data$norm_sp > spb)))/dim(phage_data)[1]
length(which(phage_data$norm_sp > spb))/dim(phage_data)[1] * length(which(phage_data$norm_pbs > pbsx))/dim(phage_data)[1]
length(which(phage_data$norm_pbs > pbsx))/dim(phage_data)[1]
length(which(phage_data$norm_sp > spb))/dim(phage_data)[1]
length(intersect(which(phage_data$norm_pbs > pbsx),which(phage_data$norm_sp > spb)))/dim(phage_data)[1]
length(which(phage_data$norm_sp > spb))/dim(phage_data)[1] * length(which(phage_data$norm_pbs > pbsx))/dim(phage_data)[1]
p1 + p2
#ggsave(paste(figs.out,'stage_plots','.png'))

p3 = ggplot(phage_data,aes(x = `norm_sp`, y = `norm_pbs`,fill = `ident`)) + 
  geom_point(pch = 21, size =2.5 , alpha = 0.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_color_manual(values=as.vector(palette.colors()))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Cluster") 

p3
ggsave(paste(figs.out,'phages_overlaid_cluster_cipro','.png'))
#Use 5th percentile
quantile(phage_data[phage_data$ident == 3,]$norm_pbs,probs = seq(0,0.1,0.01))
#quantile(phage_data[phage_data$norm_pbs > 1,]$norm_pbs,probs = seq(0,0.5,0.05))

quantile(phage_data[phage_data$ident == 4,]$norm_sp,probs = seq(0,0.1,0.01))
is_pbsx = which(phage_data$norm_pbs > pbsx)
is_spb = which(phage_data$norm_sp > spb)
is_both = intersect(is_pbsx, is_spb)
label = rep('Not Induced',dim(phage_data)[1])
label[is_pbsx] = 'PbsX'
label[is_spb] = 'Sp\U03B2'
label[is_both] = 'Both'
phage_data$phage = label


ggplot(phage_data,aes(x = `norm_sp`, y = `norm_pbs`,fill = `phage`, color = `phage`)) + 
  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  #  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_fill_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  scale_color_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 

ggsave(paste(figs.out,'phages_overlaid_coloring_cipro','.png'))




p4 = ggplot(phage_data,aes(x = `norm_pbs2`, y = `norm_sp2`)) + 
  geom_point(pch = 21, aes(fill = pca_out$label),size =2 ) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-0,0.40) + ylim(-0,0.40)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(fill = "Louvain Cluster") 


p4
#ggsave(paste(figs.out,'normed_phages_overlaid_cluster','.png'))

p3 = ggplot(phage_data,aes(x = `norm_pbs`, y = `norm_sp`,fill = umi_counts)) + 
  geom_point(pch = 21, size =2 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +#+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Cluster") 


ggplot(phage_data%>%arrange(sp_stage),aes(x = `X1`, y = `X2`)) + 
  geom_point(pch = 21, aes(fill = sp_stage),size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
  #  geom_density_2d()+
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw()+ #   xlim(-1,250) + ylim(-1,30)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  


ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_stage`, y = `pb_stage`)) + 
  geom_point(pch = 21,size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
  #  geom_density_2d()+
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw()+   #+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  



#ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_st2`, y = `pb_st2`)) + 
#ggplot(phage_data%>%arrange(pb_st2),aes(x = `pb_st2`, y = `pb_st3`)) + 
ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_st2`, y = `sp_st3`)) + 
  geom_point(pch = 21, aes(fill = pb_st2),size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  



ggplot(phage_data,aes(x = `X1`, y = `X2`,col = `ident`)) + 
  geom_point( size = 1.5) +  labs(x = "Sp\U03B2", y = "PbsX", color = "Class", parse = TRUE) +
  ggtitle(sprintf("Phage-scored cells")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))


ggplot(phage_data,aes(x = `X1`, y = `X2`,col = `ident`)) + 
  geom_point( size = 1) +  labs(x = "Sp\U03B2", y = "PbsX", color = "Class", parse = TRUE) +
  ggtitle(sprintf("Phage-scored cells")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-10,4) + ylim(-2,4)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
####

library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out2 = cbind(umap_out, phage_data$norm_pbs,phage_data$norm_sp)
umap_out2$max_val = pmax(phage_data$norm_pbs,phage_data$norm_sp)
umap_out2$phage = phage_data$phage
temp_phage = matrix(phage_data[c('norm_pbs','norm_sp')],ncol=2)
colnames(temp_phage) = paste0('PHAGE_',1:2)
sct[["phage"]] <- CreateDimReducObject(embeddings = temp_phage, key = "PHAGE_", assay = DefaultAssay(sct))

look_into = intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_pbs2 > 0.02))
look_into = union(which(phage_data$norm_sp2 > 0.05),which(phage_data$norm_pbs2 > 0.05))
look_into1 = which(phage_data$norm_sp2 > 0.5)
look_into2 = intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_sp2 < 0.5))
look_into3 = intersect(intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_pbs2 > 0.02)),intersect(which(phage_data$norm_sp2 > 0.05),which(phage_data$norm_sp2 < 0.5)))
look_into4 =  which(phage_data$norm_pbs2 > 0.1)
levels(sct@active.ident) = c(levels(sct@active.ident),"6","7",'8','9','30')
#sct@active.ident[look_into1] = "6"
#sct@active.ident[look_into2] = "7"
#sct@active.ident[look_into3] = "8"
#sct@active.ident[look_into4] = "9"
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1,ident.2 = 9)
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1,ident.2 = 8)
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1,ident.2 = 7)
#cluster4.markers2 <- FindMarkers(sct, ident.1 =6, min.pct = 0.05)
#up = intersect(rownames(cluster4.markers1[cluster4.markers1$avg_log2FC > 0,]),rownames(cluster4.markers2[cluster4.markers2$avg_log2FC > 0,]))
#down = intersect(rownames(cluster4.markers1[cluster4.markers1$avg_log2FC < 0,]),rownames(cluster4.markers2[cluster4.markers2$avg_log2FC < 0,]))
#sct@assays$RNA
#cluster4.markers1[up,]
#DimPlot(sct,reduction='pca')
phage_data5 = cbind(phage_data, norm_data)
ggplot(phage_data5%>%arrange(`BS-xtmA`),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = `BS-xtmA`),size =3 ) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-0,0.40) + ylim(-0,0.40)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(fill = "rnjB expression") 
ggsave(paste(figs.out,'rnjB_phageplot','.png'))

library(dplyr)
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `phage_data$norm_pbs`, fill = `phage_data$norm_pbs`)) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5) +
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #  scale_fill_gradientn(colors =  (brewer.pal(n = 9, name =  "BuGn") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "BuGn")), guide = '' ) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) +   #scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greens")[1]), high = (brewer.pal(n = 9, name = "Greens")[9])) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
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
ggsave(paste(figs.out,"pbsx_percentage_umap_cipro",'.png'))
#ggsave(paste(figs.out,"pbsx_percentage_umap_40pcs",'.png'))
#ggsave(paste(figs.out,"pbsx_percentage_umap_40pcs_cipro_only",'.png'))
library(prismatic)
new_col = prismatic::clr_lighten('red2',shift=0.1)
#new_col = prismatic::clr_saturate('red3',shift=0.1)
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2` )) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5, aes( colour = `phage_data$norm_sp`, fill = `phage_data$norm_sp`) )+
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #scale_colour_gradientn(low = (brewer.pal(n = 9, name = "Reds")[1]), high = (brewer.pal(n = 9, name = "Reds")[9]), guide = 'colourbar') + 
  #scale_fill_gradientn(colors =  (brewer.pal(n = 9, name = "Reds") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "Reds")), guide = '' ) + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red') + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red') + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
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
  labs(fill = "", colour = '') 

#ggsave(paste(figs.out,"sp_percentage_umap_v2",'.png'))
ggsave(paste(figs.out,"sp_percentage_umap_cipro",'.png'))
#ggsave(paste(figs.out,"sp_percentage_umap_30pcs",'.png'))

ggplot(umap_out2,aes(x = `UMAP_1`, y = `UMAP_2`,fill = `phage`, color = `phage`)) + 
  #  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  geom_point(pch = 21, size =0.7, alpha = 0.8) +  labs(x = "", y = "") +
  #  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  #  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_fill_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  scale_color_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        #        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        #        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 
ggsave(paste(figs.out,"cipro_umap_phage_plot",'.png'))



#### All ####


#keep_rows  = which(data2$Label =="Tet" | data2$Label =="Naal" |  data2$Label == "Ciproa"|data2$Label=="3hr_Cef_Rep1a")
data = data2
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
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
rm(mat)
rm(sc_precursor)
rm(sc)
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
cluster8.markers <- FindMarkers(sct, ident.1 = 8, min.pct = 0.05,logfc.threshold = 0.1)
cluster9.markers <- FindMarkers(sct, ident.1 = 9, min.pct = 0.05,logfc.threshold = 0.1)
cluster10.markers <- FindMarkers(sct, ident.1 = 10, min.pct = 0.05,logfc.threshold = 0.1)
cluster11.markers <- FindMarkers(sct, ident.1 = 11, min.pct = 0.05,logfc.threshold = 0.1)
cluster12.markers <- FindMarkers(sct, ident.1 = 12, min.pct = 0.05,logfc.threshold = 0.1)
cluster13.markers <- FindMarkers(sct, ident.1 = 13, min.pct = 0.05,logfc.threshold = 0.1)


#my_levels <- c(0,2,3,1)
#levels(sct) = my_levels
#markers %>%
#  group_by(cluster) %>%  top_n(n = 2, wt = avg_log2FC)
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 5, wt = avg_log2FC) -> top10
#DoHeatmap(sct, features = top10$gene,draw.lines = TRUE)  + scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ 
#  theme(text = element_text(size = 18))
#ggsave(paste(figs.out,'BS168_heatmap_exp_3hr','.png'))
##### Matrix Plot ####

#### Matrixplot
library(dplyr)

#av.exp.sct <- AverageExpression(sct, group.by = 'Label',return.seurat = TRUE)
av.exp.sct <- AverageExpression(sct, group.by = 'ident',return.seurat = TRUE)
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
sct@active.ident = factor(sct@meta.data$Label,labels = sct@meta.data$Label, levels = sct@meta.data$Label)
markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.5,group.by = 'Label')
markers2 <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.05,group.by = 'Label')
markers2 %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top10

markers = markers %>%
  add_count(cluster) %>%
  group_by(cluster) %>%
  mutate(n = 8) %>%
  #  group_by(rating, .add = TRUE) %>%
  #In old dplyr add = TRUE
  #group_by(rating, add = TRUE) %>%
  sample_n(n, replace = TRUE) #%>%
#  select(-n)

markers2[markers2$gene == "BS-recA",]
markers = rbind(markers, markers2[markers2$gene == "BS-recA",])
markers = unique(markers)
#markers = markers[markers$avg_log2FC > 0.5,]
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 3, wt = avg_log2FC) -> top10

#markers = markers[1:8,]
genes = markers$gene
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
g = merge(g, markers)
#g = g[which(g$variable%in% genes),]
g = g[order(g$cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'BS-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = g$gene,labels = g$gene)
g$Cluster = factor(g$Cluster, levels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'),labels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'))
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
unique(umap_out$label)
library(rcartocolor)
#  scale_color_simpsons()#scale_color_brewer(palette="Set2")
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
#p1 + p2
#ggsave(paste(figs.out, 'umap_subtilis_all_label.png',sep=''))
#ggsave(paste(figs.out, 'umap_cluster_BS168_total.png'))
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)

umap_out = umap_out[umap_out$label != 'Overnight',]

umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chlor",'Ciprofloxacin','Nalidixic acid','Exponential'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.2,alpha = 0.8, stroke = 0.2) +  
  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("B. subtilis")) +
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
  #  scale_color_simpsons()s  #  scale_color_simpsons()s
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(4,7,1,11,9,2,3,8,12)])
#ggsave(paste(figs.out, 'umap_BS168_clean_palette2.png'))
ggsave(paste(figs.out, 'umap_BS168_clean_palette_all.png'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.05,alpha = 0.6) +   labs(x = "", y = "", color = "label") +
  ggtitle(sprintf("B. subtilis All")) +
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
  #  scale_color_simpsons()s
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,11,9,12,4,8,7)])
ggsave(paste(figs.out, 'umap_BS168_drugs_all_scifi9.png'))

library(scales)
colors = hue_pal()(length(unique(umap_out$ident)))
umap_out$ident <- factor(umap_out$ident, levels=as.character(sort(as.numeric(unique(umap_out$ident)))), labels=as.character(sort(as.numeric(unique(umap_out$ident)))))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.05,alpha = 0.6) +   labs(x = "", y = "", color = "label") +
  ggtitle(sprintf("B. subtilis All")) +
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
#  scale_color_brewer(palette="scales")
  #  scale_color_simpsons()s
  scale_color_manual(values=colors)
ggsave(paste(figs.out, 'umap_BS168_clusters_all_scifi9.png'))

#### Add phages ####

# Create a weighted plot.. 
sct2 = sct
phage_anno = data.frame(read.csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/scBactoSeq/data/spbeta_genes.csv'))
phage_anno$gene_id = gsub('_','-',phage_anno$gene_id)
phage_anno$sp_st1 = as.numeric(phage_anno$spbeta_stage==1)
phage_anno$sp_st2 = as.numeric(phage_anno$spbeta_stage==2)
phage_anno$sp_st3 = as.numeric(phage_anno$spbeta_stage==3)
phage_anno$pb_st1 = as.numeric(phage_anno$pbs_stage==1)
phage_anno$pb_st2 = as.numeric(phage_anno$pbs_stage==2)
phage_anno$pb_st3 = as.numeric(phage_anno$pbs_stage==3)
gene_indices = c()
select = c()
for (i in 1:length(phage_anno$gene_id)){
  index = which(phage_anno$gene_id[i] == rownames(sct2@assays$RNA))
  if (length(index)>0){
    gene_indices = c(gene_indices,index)
    select = c(select,i)
  }
}
phage_anno2 = phage_anno[select,]
phage_anno2$gene_indices = gene_indices
phage_score = matrix(data = 0,ncol=2,nrow = dim(sct2@assays$RNA)[1])
phage_score[phage_anno2$gene_indices,1] = as.numeric(phage_anno2$spbeta_score)
phage_score[phage_anno2$gene_indices,2] = as.numeric(phage_anno2$pbsx_score)
phage_score2 = matrix(data = 0,ncol=8,nrow = dim(sct2@assays$RNA)[1])
phage_score2[phage_anno2$gene_indices,1] = phage_anno2$sp_st1
phage_score2[phage_anno2$gene_indices,2] = phage_anno2$sp_st2
phage_score2[phage_anno2$gene_indices,3] = phage_anno2$sp_st3
phage_score2[phage_anno2$gene_indices,4] = phage_anno2$pb_st1
phage_score2[phage_anno2$gene_indices,5] = phage_anno2$pb_st2
phage_score2[phage_anno2$gene_indices,6] = phage_anno2$pb_st3
phage_score2[phage_anno2$gene_indices,7] = phage_anno2$spbeta_stage
phage_score2[phage_anno2$gene_indices,8] = phage_anno2$pbs_stage
norm_data = t(as.data.frame(sct2@assays$RNA@data))
#phage_data = data.frame(norm_data %*% phage_score)
#phage_data2 = data.frame(norm_data %*% phage_score2)
phage_data = data.frame(t(as.matrix(sct2@assays$RNA@counts)) %*% phage_score)
phage_data2 = data.frame(t(as.matrix(sct2@assays$RNA@counts)) %*% phage_score2)
phage_data3 = data.frame(norm_data %*% phage_score)
divide1 = phage_data[,1]
divide2 = phage_data[,2]
divide1[which(divide1 < 2)] = 0 
divide2[which(divide2 < 2)] = 0 
phage_data2[,7] = phage_data2[,7]/divide1
phage_data2[,8] = phage_data2[,8]/divide2
is.na(phage_data2)<-sapply(phage_data2, is.infinite)
phage_data2[is.na(phage_data2)] = 0
#### Cell Stage ####
names(phage_data2) = c('sp_st1','sp_st2','sp_st3','pb_st1','pb_st2','pb_st3','sp_stage','pb_stage')
names(phage_data3) = c('norm_sp','norm_pbs')
phage_data = cbind(phage_data, phage_data2)
phage_data = cbind(phage_data, phage_data3)

phage_data$umi_counts = as.numeric(as.character(sct2[['nCount_RNA']][,1]))
phage_data$label = as.character(sct2[['Label']][,1])
phage_data$ident = as.character(sct2@active.ident)
phage_data$index = 1:dim(phage_data)[1]
phage_data$norm_pbs2 = phage_data$norm_pbs/max(phage_data$norm_pbs)
phage_data$norm_sp2 = phage_data$norm_sp/max(phage_data$norm_sp)
#phage_data$x1_stage = phage_data2$X1
#phage_data$x2_stage = phage_data2$X2

library(dplyr)
#ggplot(phage_data%>%arrange(sp_stage),aes(x = `X1`, y = `X2`)) + 
p1 = ggplot(phage_data%>%arrange(pb_stage),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = pb_stage),size =3 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +#+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "PbsX stage") 


p2 = ggplot(phage_data%>%arrange(sp_stage),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = sp_stage),size =3 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-0,0.20) + ylim(-0,0.20)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Sp\U03B2 stage") 

p2
aggregate(norm_sp~ident,phage_data,min)
aggregate(norm_pbs~ident,phage_data,min)
pbsx=2.94
spb=16.67
#pbsx=0
#spb=0
length(which(phage_data$norm_pbs > pbsx))/dim(phage_data)[1]
length(which(phage_data$norm_sp > spb))/dim(phage_data)[1]
length(intersect(which(phage_data$norm_pbs > pbsx),which(phage_data$norm_sp > spb)))/dim(phage_data)[1]
length(which(phage_data$norm_pbs > pbsx))/dim(phage_data)[1]*length(which(phage_data$norm_sp > spb))/dim(phage_data)[1]
p1 + p2
#ggsave(paste(figs.out,'stage_plots','.png'))
phage_data
p3 = ggplot(phage_data,aes(x = `norm_sp`, y = `norm_pbs`,fill = `ident`)) + 
  geom_point(pch = 21, size =2.5, alpha = 0.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_color_manual(values=as.vector(palette.colors()))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Cluster") 

p3
quantile(phage_data[phage_data$ident == 2,]$norm_pbs,probs = seq(0,0.1,0.01))
quantile(phage_data[phage_data$ident == 1,]$norm_sp,probs = seq(0,0.1,0.01))
is_pbsx = which(phage_data$norm_pbs > pbsx)
is_spb = which(phage_data$norm_sp > spb)
is_both = intersect(is_pbsx, is_spb)
label = rep('Not Induced',dim(phage_data)[1])
label[is_pbsx] = 'PbsX'
label[is_spb] = 'Sp\U03B2'
label[is_both] = 'Both'
phage_data$phage = label


ggplot(phage_data,aes(x = `norm_sp`, y = `norm_pbs`,fill = `phage`, color = `phage`)) + 
  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  #  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_fill_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  scale_color_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 



p4 = ggplot(phage_data,aes(x = `norm_pbs2`, y = `norm_sp2`)) + 
  geom_point(pch = 21, aes(fill = pca_out$label),size =2 ) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-0,0.40) + ylim(-0,0.40)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(fill = "Louvain Cluster") 


p4
#ggsave(paste(figs.out,'normed_phages_overlaid_cluster','.png'))

p3 = ggplot(phage_data,aes(x = `norm_pbs`, y = `norm_sp`,fill = umi_counts)) + 
  geom_point(pch = 21, size =2 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +#+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Cluster") 


ggplot(phage_data%>%arrange(sp_stage),aes(x = `X1`, y = `X2`)) + 
  geom_point(pch = 21, aes(fill = sp_stage),size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
  #  geom_density_2d()+
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw()+ #   xlim(-1,250) + ylim(-1,30)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  


ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_stage`, y = `pb_stage`)) + 
  geom_point(pch = 21,size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
  #  geom_density_2d()+
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw()+   #+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  



#ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_st2`, y = `pb_st2`)) + 
#ggplot(phage_data%>%arrange(pb_st2),aes(x = `pb_st2`, y = `pb_st3`)) + 
ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_st2`, y = `sp_st3`)) + 
  geom_point(pch = 21, aes(fill = pb_st2),size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  



ggplot(phage_data,aes(x = `X1`, y = `X2`,col = `ident`)) + 
  geom_point( size = 1.5) +  labs(x = "Sp\U03B2", y = "PbsX", color = "Class", parse = TRUE) +
  ggtitle(sprintf("Phage-scored cells")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))


ggplot(phage_data,aes(x = `X1`, y = `X2`,col = `ident`)) + 
  geom_point( size = 1) +  labs(x = "Sp\U03B2", y = "PbsX", color = "Class", parse = TRUE) +
  ggtitle(sprintf("Phage-scored cells")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-10,4) + ylim(-2,4)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
####

library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out2 = cbind(umap_out, phage_data$norm_pbs,phage_data$norm_sp)
umap_out2$max_val = pmax(phage_data$norm_pbs,phage_data$norm_sp)
temp_phage = matrix(phage_data[c('norm_pbs','norm_sp')],ncol=2)
colnames(temp_phage) = paste0('PHAGE_',1:2)
sct[["phage"]] <- CreateDimReducObject(embeddings = temp_phage, key = "PHAGE_", assay = DefaultAssay(sct))

look_into = intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_pbs2 > 0.02))
look_into = union(which(phage_data$norm_sp2 > 0.05),which(phage_data$norm_pbs2 > 0.05))
look_into1 = which(phage_data$norm_sp2 > 0.5)
look_into2 = intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_sp2 < 0.5))
look_into3 = intersect(intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_pbs2 > 0.02)),intersect(which(phage_data$norm_sp2 > 0.05),which(phage_data$norm_sp2 < 0.5)))
look_into4 =  which(phage_data$norm_pbs2 > 0.1)
levels(sct@active.ident) = c(levels(sct@active.ident),"6","7",'8','9','30')
#sct@active.ident[look_into1] = "6"
#sct@active.ident[look_into2] = "7"
#sct@active.ident[look_into3] = "8"
#sct@active.ident[look_into4] = "9"
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1,ident.2 = 9)
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1,ident.2 = 8)
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1,ident.2 = 7)
#cluster4.markers2 <- FindMarkers(sct, ident.1 =6, min.pct = 0.05)
#up = intersect(rownames(cluster4.markers1[cluster4.markers1$avg_log2FC > 0,]),rownames(cluster4.markers2[cluster4.markers2$avg_log2FC > 0,]))
#down = intersect(rownames(cluster4.markers1[cluster4.markers1$avg_log2FC < 0,]),rownames(cluster4.markers2[cluster4.markers2$avg_log2FC < 0,]))
#sct@assays$RNA
#cluster4.markers1[up,]
#DimPlot(sct,reduction='pca')
phage_data5 = cbind(phage_data, norm_data)
ggplot(phage_data5%>%arrange(`BS-xtmA`),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = `BS-xtmA`),size =3 ) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-0,0.40) + ylim(-0,0.40)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(fill = "rnjB expression") 
ggsave(paste(figs.out,'rnjB_phageplot','.png'))

library(dplyr)
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `phage_data$norm_pbs`, fill = `phage_data$norm_pbs`)) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =0.5,alpha = 0.8) +
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #  scale_fill_gradientn(colors =  (brewer.pal(n = 9, name =  "BuGn") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "BuGn")), guide = '' ) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) +   #scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greens")[1]), high = (brewer.pal(n = 9, name = "Greens")[9])) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
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
ggsave(paste(figs.out,"pbsx_percentage_umap_all",'.png'))
#ggsave(paste(figs.out,"pbsx_percentage_umap_40pcs",'.png'))
#ggsave(paste(figs.out,"pbsx_percentage_umap_40pcs_cipro_only",'.png'))
library(prismatic)
new_col = prismatic::clr_lighten('red2',shift=0.1)
#new_col = prismatic::clr_saturate('red3',shift=0.1)
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2` )) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =0.5, alpha = 0.8, aes( colour = `phage_data$norm_sp`, fill = `phage_data$norm_sp`) )+
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
#  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
#  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  #scale_colour_gradientn(low = (brewer.pal(n = 9, name = "Reds")[1]), high = (brewer.pal(n = 9, name = "Reds")[9]), guide = 'colourbar') + 
  #scale_fill_gradientn(colors =  (brewer.pal(n = 9, name = "Reds") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "Reds")), guide = '' ) + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red2') + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red2') + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
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
  labs(fill = "", colour = '') 

#ggsave(paste(figs.out,"sp_percentage_umap_v2",'.png'))
ggsave(paste(figs.out,"sp_percentage_umap_all",'.png'))
#ggsave(paste(figs.out,"sp_percentage_umap_30pcs",'.png'))




exp_diff = FindMarkers(sct, ident.1 = 0,ident.2 = 4)
stat_diff = FindMarkers(sct, ident.1 = 0,ident.2 = 11)
sch_diff = FindMarkers(sct, ident.1 = 1,ident.2 = 10)
#FeaturePlot(sct,features = c('1','guaBA','gadBC','hipBA'), reduction = 'umap')
#FeaturePlot(sct,features = c('BS-ykuO','BS-sigH','BS-suhB','BS-hemE'), reduction = 'umap')
#FeaturePlot(sct,features = c('BS-tufA','BS-putA','BS-sdiA','BS-bcsA'), reduction = 'umap')
#FeaturePlot(sct,features = c('BS-rplP','BS-dmlA','BS-sdiA','BS-tap'), reduction = 'umap',size  = 4)
#ggsave(paste(figs.out,'BS168_geneplot_scifi3','.png'))

DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_BS168_scifi3','.png'))

#FeaturePlot(sct,features =  c('acoA'), reduction = 'umap') & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-srfAB','BS-sigF','BS-cwlC','BS-BSU-04345'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-rpsQ','BS-trnS-Leu1','BS-floT','BS-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-gapB','BS-mntA','BS-floT','BS-pftA'), reduction = 'pca')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-opuCB','BS-ydjP','BS-sigB','BS-cwlC'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-rpsQ','BS-manA','BS-ykzE','BS-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-rppA','BS-groEL','BS-xtmA','BS-yonO'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_oosmootic1','.png'))
FeaturePlot(sct,features =  c('BS-dcuA','BS-tnaA','BS-katG','BS-malK'), reduction = 'pca')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_oosmootic1_pca','.png'))
#FeaturePlot(sct,features =  c('BS-ydbS','BS-katE','BS-clpE','BS-cwlC'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-spxA','BS-yhgE','BS-nfeDB','BS-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-ydbS','BS-katE','BS-clpE','BS-cwlC'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_oosmootic2','.png'))
FeaturePlot(sct,features =  c('BS-ydbS','BS-katE','BS-clpE','BS-cwlC'), reduction = 'pca')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-spxA','BS-yhgE','BS-nfeDB','BS-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_oosmootic2_pca.png'))
FeaturePlot(sct,features =  c('BS-yuaI','BS-gapB','BS-groES','BS-cwlC'), reduction = 'pca')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-spxA','BS-yhgE','BS-nfeDB','BS-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#ggsave(paste(figs.out,'BS168_geneplot_exp_heat_cef','.png'))
levels(sct)  = c(0,1,2,3,4,6,5,7)
markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
library(dplyr)
markers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC) -> top10
new_genes = top10$gene
#new_genes = insert(top10$gene,ats = 18,values = c('BS-spo0M','BS-cwlC'))

DoHeatmap(sct, new_genes,draw.lines = TRUE)  + scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_heatmap_osmotic','.png'))


#### PCA ####


pca_out = data.frame(sct[['pca']]@cell.embeddings)
pca_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
pca_out$label = as.character(sct[['Label']][,1])
pca_out$ident = as.character(sct@active.ident)
library(scales)
colors = hue_pal()(length(unique(pca_out$label)))
col_key = data.frame(cbind(colors, unique(pca_out$label)))
colnames(col_key) = c('color','label')
#pca_out = merge(pca_out,col_key,by.x = 'label', by.y = 'label')
#plot3d(pca_out[,2],pca_out[,3],pca_out[,4], col = as.vector(pca_out$color),xlab = 'PC1',ylab = 'PC2',
#       zlab = 'PC3',bty = "b2", main ="BS168 PCA",size = 3,alpha = 0.8,type = 'p', specular="black" )
#legend3d("topright", legend = paste(col_key$label), pch = 16, col = col_key$color, cex=2, inset=c(0.08))


p1 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`)) + 
  geom_point(pch = 21, aes(fill = umi_counts),size = 1.0) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
p2 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`,col = `label`)) + 
  geom_point( size = 1.0) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("PCA - BS168 - 8hr Cipro")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")
p1 + p2
p2
ggsave(paste(figs.out, 'pca_8hr_BS168.png'))

#### UMAP Variability ###

#### Variablee Genes


umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)


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
p2 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 1,alpha = 1.5) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMAP - BS168 - 8hr Cipro")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")
p3 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMAP for Cluster")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Spectral")
p1 + p2
ggsave(paste(figs.out, 'umap_umi_counts_BS168_osmotic.png'))
p2
ggsave(paste(figs.out, 'umap_samples_BS168_8hr_cipro.png'))










#######
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.5,alpha = 0.8) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMAP for Class - E. BS168")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Spectral")
ggsave(paste(figs.out, 'umap_umi_counts_BS168.png'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.1,alpha = 0.1) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  geom_point(data = subset(umap_out, label == "2HR_OSMOTIC"),aes(x = `UMAP_1`, y = `UMAP_2`), size = 1,alpha = 0.8) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("BS168 Osmotic Shock (0.3M NaCl)")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Spectral")
ggsave(paste(figs.out, 'BS168_umap_umi_counts_osmotic.png'))


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.1,alpha = 0.1) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  geom_point(data = subset(umap_out, label == "2hr_Cef"),aes(x = `UMAP_1`, y = `UMAP_2`), size = 1,alpha = 0.8) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("E. BS168 Osmotic Shock")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Spectral")

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.1,alpha = 0.1) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  geom_point(data = subset(umap_out, label == "2HR_CIPRO"),aes(x = `UMAP_1`, y = `UMAP_2`), size = 1,alpha = 0.8) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("E. BS168 Cefozolin")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(),  panel.grid = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Spectral")
ggsave(paste(figs.out, 'BS168_umap_umi_counts_osmotic.png'))



ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.1,alpha = 0.1) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  geom_point(data = subset(umap_out, label == "4hr_Cef"),aes(x = `UMAP_1`, y = `UMAP_2`), size = 1,alpha = 0.8) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("E. BS168 Cefozolin")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(),  panel.grid = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Spectral")
#ggsave(paste(figs.out, 'BS168_umap_umi_counts_cefozolin.png'))


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.1,alpha = 0.1) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  geom_point(data = subset(umap_out, label == "2hr_Cef"),aes(x = `UMAP_1`, y = `UMAP_2`), size = 1,alpha = 0.8) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("E. BS168 Cefozolin")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(),  panel.grid = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Spectral")
#ggsave(paste(figs.out, 'BS168_umap_umi_counts_cefozolin.png'))

#### Cipro & Nal ####

keep_rows  = which(data2$Label =="Cipro" | data2$Label =="Nal" |  data2$Label == "Exponaential"|data2$Label=="3hr_Cef_Rep1a")
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

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_trnO-Ala','BS_trnA-Ala','BS_trnE-Gly','BS_trnSL-Gly1'))]
var_nal = apply(sc_precursor,2,var)
mean_nal = apply(sc_precursor,2,mean)
fano_nal = var_nal/mean_nal^2
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
#sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 10000)
sct = FindVariableFeatures(sct, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(sct), 10)
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)
sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
rm(mat)
rm(sc_precursor)
sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.3)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,label.size = 0.1,reduction = 'pca', dims =as.numeric(gsub('PC_','',names(sort(abs(kurtosi(pca_out)),decreasing = TRUE))[1:2]))) + ggtitle('PCA - Tetracycline "High Kurtosis"') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
kurtosi_total_all = cbind(1:length(sort(abs(kurtosi(pca_out)),decreasing = TRUE)), sort(abs(kurtosi(pca_out)),decreasing = TRUE), rep("All",dim(pca_out)[1]))

cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, ident.2 = 3,min.pct = 0.1,logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05,logfc.threshold = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.05,logfc.threshold = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,logfc.threshold = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.05,logfc.threshold = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.05,logfc.threshold = 0.1)
cluster8.markers <- FindMarkers(sct, ident.1 = 8, min.pct = 0.05,logfc.threshold = 0.1)


zeroes = which(sct@active.ident == 0)
ones = which(sct@active.ident == 1)
twos = which(sct@active.ident == 2)
threes = which(sct@active.ident == 3)
fours = which(sct@active.ident == 4)
fives = which(sct@active.ident == 5)
sixes = which(sct@active.ident == 6)
sevens = which(sct@active.ident == 7)

sct@active.ident[threes] = 5
sct@active.ident[fives] = 3
sct@active.ident[twos] = 4
sct@active.ident[fours] = 2

#my_levels <- c(0,2,3,1)
#levels(sct) = my_levels
#markers %>%
#  group_by(cluster) %>%  top_n(n = 2, wt = avg_log2FC)
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 5, wt = avg_log2FC) -> top10
#DoHeatmap(sct, features = top10$gene,draw.lines = TRUE)  + scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ 
#  theme(text = element_text(size = 18))
#ggsave(paste(figs.out,'BS168_heatmap_exp_3hr','.png'))
#### PCA ####

library(dplyr)

av.exp.sct <- AverageExpression(sct, group.by = 'ident',return.seurat = TRUE)
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
markers2 <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers2 %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top10


markers = markers2 %>%
  add_count(cluster) %>%
  group_by(cluster) %>%
  mutate(n = 6) %>%
  #  group_by(rating, .add = TRUE) %>%
  #In old dplyr add = TRUE
  #group_by(rating, add = TRUE) %>%
  sample_n(n, replace = TRUE) #%>%
#  select(-n)

markers2[markers2$gene == "BS-yonO",]
markers = rbind(markers, markers2[markers2$gene == "BS-yonO",])
markers = unique(markers)
markers = markers[markers$gene != 'BS-mtbP',]
#markers = markers[1:8,]
genes = markers$gene
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
g = merge(g, markers)
#g = g[which(g$variable%in% genes),]
g = g[order(g$cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'BS-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = g$gene,labels = g$gene)
#g$Cluster = factor(as.numeric(g$Cluster)+1, levels = c(10,9,3,2,1,8,7,6,5,4))
g$Cluster = factor(as.numeric(g$Cluster)+1, levels = c(8,7,6,5,4,3,2,1))

#g = g[order(g$cluster),]

library(myriad)
ggplot(g, aes(x = gene, Cluster, fill = value)) +
  geom_tile(lwd = 0.5, color = 'black') +
  coord_equal()+
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")),
                       limits = c(-1, 1.5), oob = scales::squish)+
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

#ggsave(paste(figs.out,'cipro_nal_matrixplot_with_words_lognorm','.png'))
ggsave(paste(figs.out,'cipro_nal_matrixplot_with_words_lognorm','.svg'), device = 'svg')

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

ggsave(paste(figs.out,'subtilis_prophages_matrixplot_with_no_words_lognorm','.png'))
##

#### UMAP Variability ###
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
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
umap_out$label = gsub('Cycloserineserine','Cycloserine',umap_out$label)
umap_out$label = gsub('Erythromycinmycin','Erythromycin',umap_out$label)
umap_out$label = gsub('Tetracyclineracycline','Tetracycline',umap_out$label)
unique(umap_out$label)
library(rcartocolor)
#  scale_color_simpsons()#scale_color_brewer(palette="Set2")
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
#p1 + p2
#ggsave(paste(figs.out, 'umap_subtilis_all_label.png',sep=''))
#ggsave(paste(figs.out, 'umap_cluster_BS168_total.png'))
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.5,alpha = 1.0) +   labs(x = "", y = "", color = "label") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
   labs(color='Condition') +
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
  #  scale_color_simpsons()s
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(3,8,12)])
ggsave(paste(figs.out, 'umap_BS168_cipro_nal_scifi9_lognorm.png'))

umap_out$ident = as.character(as.numeric(umap_out$ident)+1)
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.5,alpha = 1.0) +   labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +  
  #  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,11,2,9,7,4,3,8,12)]) + 
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
  scale_color_manual(values=as.vector(palette.colors())[c(1,2,4,5,6,8,9)])#+ guides(colour = guide_legend(override.aes = list(size=1)))
ggsave(paste(figs.out, 'umap_BS168_cluster_nal_cipro_lognorm.png'))
####
# Create a weighted plot.. 
sct2 = sct
phage_anno = data.frame(read.csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/scBactoSeq/data/spbeta_genes.csv'))
phage_anno$gene_id = gsub('_','-',phage_anno$gene_id)
phage_anno$sp_st1 = as.numeric(phage_anno$spbeta_stage==1)
phage_anno$sp_st2 = as.numeric(phage_anno$spbeta_stage==2)
phage_anno$sp_st3 = as.numeric(phage_anno$spbeta_stage==3)
phage_anno$pb_st1 = as.numeric(phage_anno$pbs_stage==1)
phage_anno$pb_st2 = as.numeric(phage_anno$pbs_stage==2)
phage_anno$pb_st3 = as.numeric(phage_anno$pbs_stage==3)
gene_indices = c()
select = c()
for (i in 1:length(phage_anno$gene_id)){
  index = which(phage_anno$gene_id[i] == rownames(sct2@assays$RNA))
  if (length(index)>0){
    gene_indices = c(gene_indices,index)
    select = c(select,i)
  }
}
phage_anno2 = phage_anno[select,]
phage_anno2$gene_indices = gene_indices
phage_score = matrix(data = 0,ncol=2,nrow = dim(sct2@assays$RNA)[1])
phage_score[phage_anno2$gene_indices,1] = as.numeric(phage_anno2$spbeta_score)
phage_score[phage_anno2$gene_indices,2] = as.numeric(phage_anno2$pbsx_score)
phage_score2 = matrix(data = 0,ncol=8,nrow = dim(sct2@assays$RNA)[1])
phage_score2[phage_anno2$gene_indices,1] = phage_anno2$sp_st1
phage_score2[phage_anno2$gene_indices,2] = phage_anno2$sp_st2
phage_score2[phage_anno2$gene_indices,3] = phage_anno2$sp_st3
phage_score2[phage_anno2$gene_indices,4] = phage_anno2$pb_st1
phage_score2[phage_anno2$gene_indices,5] = phage_anno2$pb_st2
phage_score2[phage_anno2$gene_indices,6] = phage_anno2$pb_st3
phage_score2[phage_anno2$gene_indices,7] = phage_anno2$spbeta_stage
phage_score2[phage_anno2$gene_indices,8] = phage_anno2$pbs_stage
all_data = t(as.data.frame(GetAssayData(sct2, 'counts')))
norm_data = all_data/rowSums(all_data)*100

#phage_data = data.frame(norm_data %*% phage_score)
#phage_data2 = data.frame(norm_data %*% phage_score2)
phage_data = data.frame(norm_data %*% phage_score)
phage_data2 = data.frame(norm_data %*% phage_score2)
phage_data3 = data.frame(norm_data %*% phage_score)
divide1 = phage_data[,1]
divide2 = phage_data[,2]
divide1[which(divide1 < 2)] = 0 
divide2[which(divide2 < 2)] = 0 
phage_data2[,7] = phage_data2[,7]/divide1
phage_data2[,8] = phage_data2[,8]/divide2
is.na(phage_data2)<-sapply(phage_data2, is.infinite)
phage_data2[is.na(phage_data2)] = 0
#### Cell Stage ####
names(phage_data2) = c('sp_st1','sp_st2','sp_st3','pb_st1','pb_st2','pb_st3','sp_stage','pb_stage')
names(phage_data3) = c('norm_sp','norm_pbs')
phage_data = cbind(phage_data, phage_data2)
phage_data = cbind(phage_data, phage_data3)

phage_data$umi_counts = as.numeric(as.character(sct2[['nCount_RNA']][,1]))
phage_data$label = as.character(sct2[['Label']][,1])
phage_data$ident = as.character(sct2@active.ident)
phage_data$index = 1:dim(phage_data)[1]
phage_data$norm_pbs2 = phage_data$norm_pbs/max(phage_data$norm_pbs)
phage_data$norm_sp2 = phage_data$norm_sp/max(phage_data$norm_sp)
#phage_data$x1_stage = phage_data2$X1
#phage_data$x2_stage = phage_data2$X2

library(dplyr)
#ggplot(phage_data%>%arrange(sp_stage),aes(x = `X1`, y = `X2`)) + 
p1 = ggplot(phage_data%>%arrange(pb_stage),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = pb_stage),size =3 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +#+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "PbsX stage") 


p2 = ggplot(phage_data%>%arrange(sp_stage),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = sp_stage),size =3 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-0,0.20) + ylim(-0,0.20)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Sp\U03B2 stage") 

p2
aggregate(norm_sp~ident,phage_data,min)
aggregate(norm_pbs~ident,phage_data,min)
#pbsx=7.39
# Old
pbsx=11.02
spb=16.67
# Combined
pbsx=8.09
spb=14.70
# Log norm
pbsx=8.40
spb=15.0

#2.57, vs 2.52
#pbsx=4.87
#spb=13.25
#pbsx=0
#spb=0
length(which(phage_data$norm_pbs > pbsx))/dim(phage_data)[1]
#17.09
length(which(phage_data$norm_sp > spb))/dim(phage_data)[1]
#14.43
length(intersect(which(phage_data$norm_pbs > pbsx),which(phage_data$norm_sp > spb)))/dim(phage_data)[1]
#2.4357
length(which(phage_data$norm_pbs > pbsx))/dim(phage_data)[1]*length(which(phage_data$norm_sp > spb))/dim(phage_data)[1]
#2.467
length(which(phage_data$norm_pbs > pbsx))/dim(phage_data)[1] - length(intersect(which(phage_data$norm_pbs > pbsx),which(phage_data$norm_sp > spb)))/dim(phage_data)[1]
#14.659
length(which(phage_data$norm_sp > spb))/dim(phage_data)[1]- length(intersect(which(phage_data$norm_pbs > pbsx),which(phage_data$norm_sp > spb)))/dim(phage_data)[1]
#11.990
100-11.99-14.659 -2.437



p1 + p2
#ggsave(paste(figs.out,'stage_plots','.png'))
phage_data
p3 = ggplot(phage_data,aes(x = `norm_sp`, y = `norm_pbs`,fill = `ident`)) + 
  geom_point(pch = 21, size =2.5, alpha = 0.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_color_manual(values=as.vector(palette.colors()))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Cluster") 

p3
ggsave(paste(figs.out,'phages_overlaid_cluster_cipro_nal','.png'))
#Use 5th percentile
#quantile(phage_data[phage_data$ident == 3,]$norm_pbs,probs = seq(0,0.1,0.01))
#quantile(phage_data[phage_data$ident == 4,]$norm_sp,probs = seq(0,0.1,0.01))
quantile(phage_data[phage_data$ident == 4,]$norm_pbs,probs = seq(0,0.1,0.01))
quantile(phage_data[phage_data$ident == 5,]$norm_sp,probs = seq(0,0.1,0.01))
is_pbsx = which(phage_data$norm_pbs > pbsx)
is_spb = which(phage_data$norm_sp > spb)
is_both = intersect(is_pbsx, is_spb)
label = rep('Not Induced',dim(phage_data)[1])
label[is_pbsx] = 'PbsX'
label[is_spb] = 'Sp\U03B2'
label[is_both] = 'Both'
phage_data$phage = label


ggplot(phage_data,aes(x = `norm_sp`, y = `norm_pbs`,fill = `phage`, color = `phage`)) + 
  #  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "", y = "") +
  #  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_fill_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  scale_color_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        #        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        #        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 

ggsave(paste(figs.out,'phages_overlaid_coloring_nal_cipro_combined_lognorm','.png'))
#ggsave(paste(figs.out,'phages_overlaid_coloring_nal_cipro_combined','.png'))


p4 = ggplot(phage_data,aes(x = `norm_pbs2`, y = `norm_sp2`)) + 
  geom_point(pch = 21, aes(fill = pca_out$label),size =2 ) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-0,0.40) + ylim(-0,0.40)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(fill = "Louvain Cluster") 


p4
#ggsave(paste(figs.out,'normed_phages_overlaid_cluster','.png'))

p3 = ggplot(phage_data,aes(x = `norm_pbs`, y = `norm_sp`,fill = umi_counts)) + 
  geom_point(pch = 21, size =2 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +#+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Cluster") 


ggplot(phage_data%>%arrange(sp_stage),aes(x = `X1`, y = `X2`)) + 
  geom_point(pch = 21, aes(fill = sp_stage),size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
  #  geom_density_2d()+
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw()+ #   xlim(-1,250) + ylim(-1,30)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  


ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_stage`, y = `pb_stage`)) + 
  geom_point(pch = 21,size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
  #  geom_density_2d()+
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw()+   #+ xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  



#ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_st2`, y = `pb_st2`)) + 
#ggplot(phage_data%>%arrange(pb_st2),aes(x = `pb_st2`, y = `pb_st3`)) + 
ggplot(phage_data%>%arrange(pb_st2),aes(x = `sp_st2`, y = `sp_st3`)) + 
  geom_point(pch = 21, aes(fill = pb_st2),size =2 ) +  labs(x = "Sp\U03B2", y = "PbsX") +
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))  



ggplot(phage_data,aes(x = `X1`, y = `X2`,col = `ident`)) + 
  geom_point( size = 1.5) +  labs(x = "Sp\U03B2", y = "PbsX", color = "Class", parse = TRUE) +
  ggtitle(sprintf("Phage-scored cells")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-1,10) + ylim(-1,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))


ggplot(phage_data,aes(x = `X1`, y = `X2`,col = `ident`)) + 
  geom_point( size = 1) +  labs(x = "Sp\U03B2", y = "PbsX", color = "Class", parse = TRUE) +
  ggtitle(sprintf("Phage-scored cells")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-10,4) + ylim(-2,4)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
####

library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out2 = cbind(umap_out, phage_data$norm_pbs,phage_data$norm_sp)
umap_out2$max_val = pmax(phage_data$norm_pbs,phage_data$norm_sp)
umap_out2$phage = phage_data$phage
temp_phage = matrix(phage_data[c('norm_pbs','norm_sp')],ncol=2)
colnames(temp_phage) = paste0('PHAGE_',1:2)
sct[["phage"]] <- CreateDimReducObject(embeddings = temp_phage, key = "PHAGE_", assay = DefaultAssay(sct))

look_into = intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_pbs2 > 0.02))
look_into = union(which(phage_data$norm_sp2 > 0.05),which(phage_data$norm_pbs2 > 0.05))
look_into1 = which(phage_data$norm_sp2 > 0.5)
look_into2 = intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_sp2 < 0.5))
look_into3 = intersect(intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_pbs2 > 0.02)),intersect(which(phage_data$norm_sp2 > 0.05),which(phage_data$norm_sp2 < 0.5)))
look_into4 =  which(phage_data$norm_pbs2 > 0.1)
levels(sct@active.ident) = c(levels(sct@active.ident),"6","7",'8','9','30')
#sct@active.ident[look_into1] = "6"
#sct@active.ident[look_into2] = "7"
#sct@active.ident[look_into3] = "8"
#sct@active.ident[look_into4] = "9"
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1,ident.2 = 9)
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1,ident.2 = 8)
#cluster4.markers1 <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1,ident.2 = 7)
#cluster4.markers2 <- FindMarkers(sct, ident.1 =6, min.pct = 0.05)
#up = intersect(rownames(cluster4.markers1[cluster4.markers1$avg_log2FC > 0,]),rownames(cluster4.markers2[cluster4.markers2$avg_log2FC > 0,]))
#down = intersect(rownames(cluster4.markers1[cluster4.markers1$avg_log2FC < 0,]),rownames(cluster4.markers2[cluster4.markers2$avg_log2FC < 0,]))
#sct@assays$RNA
#cluster4.markers1[up,]
#DimPlot(sct,reduction='pca')
phage_data5 = cbind(phage_data, norm_data)
ggplot(phage_data5%>%arrange(`BS-xtmA`),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = `BS-xtmA`),size =3 ) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-0,0.40) + ylim(-0,0.40)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+
  labs(fill = "rnjB expression") 
ggsave(paste(figs.out,'rnjB_phageplot','.png'))

library(dplyr)

library(dplyr)
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `phage_data$norm_pbs`, fill = `phage_data$norm_pbs`)) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =0.6,alpha = 1.2) +
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #  scale_fill_gradientn(colors =  (brewer.pal(n = 9, name =  "BuGn") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "BuGn")), guide = '' ) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) +   #scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greens")[1]), high = (brewer.pal(n = 9, name = "Greens")[9])) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Greens"))[9]) + 
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
ggsave(paste(figs.out,"pbsx_percentage_umap_nal_cipro_combined",'.png'))
#ggsave(paste(figs.out,"pbsx_percentage_umap_40pcs",'.png'))
#ggsave(paste(figs.out,"pbsx_percentage_umap_40pcs_cipro_only",'.png'))
library(prismatic)
new_col = prismatic::clr_lighten('red2',shift=0.1)
#new_col = prismatic::clr_saturate('red3',shift=0.1)
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2` )) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =0.6, alpha = 1.2, aes( colour = `phage_data$norm_sp`, fill = `phage_data$norm_sp`) )+
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  #scale_colour_gradientn(low = (brewer.pal(n = 9, name = "Reds")[1]), high = (brewer.pal(n = 9, name = "Reds")[9]), guide = 'colourbar') + 
  #scale_fill_gradientn(colors =  (brewer.pal(n = 9, name = "Reds") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "Reds")), guide = '' ) + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red2') + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red2') + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
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
  labs(fill = "", colour = '') 

#ggsave(paste(figs.out,"sp_percentage_umap_v2",'.png'))
ggsave(paste(figs.out,"sp_percentage_umap_nal_cipro_combined",'.png'))
#ggsave(paste(figs.out,"sp_percentage_umap_30pcs",'.png'))


#ggsave(paste(figs.out,"pbsx_percentage_umap_40pcs",'.png'))
#ggsave(paste(figs.out,"pbsx_percentage_umap_40pcs_cipro_only",'.png'))
library(prismatic)
new_col = prismatic::clr_lighten('red2',shift=0.1)
#new_col = prismatic::clr_saturate('red3',shift=0.1)
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2` )) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5, aes( colour = `phage_data$norm_sp`, fill = `phage_data$norm_sp`) )+
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #scale_colour_gradientn(low = (brewer.pal(n = 9, name = "Reds")[1]), high = (brewer.pal(n = 9, name = "Reds")[9]), guide = 'colourbar') + 
  #scale_fill_gradientn(colors =  (brewer.pal(n = 9, name = "Reds") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "Reds")), guide = '' ) + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red') + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red') + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
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
  labs(fill = "", colour = '') 


ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2` )) + 
  #  geom_point(pch = 21, aes(fill = `phage_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =1.5, aes( colour = `phage_data$norm_sp`, fill = `phage_data$norm_sp`) )+
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #scale_colour_gradientn(low = (brewer.pal(n = 9, name = "Reds")[1]), high = (brewer.pal(n = 9, name = "Reds")[9]), guide = 'colourbar') + 
  #scale_fill_gradientn(colors =  (brewer.pal(n = 9, name = "Reds") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "Reds")), guide = '' ) + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = (brewer.pal(n = 9, name = "Reds"))[9]) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red') + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'red') + 
  #  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = new_col) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) + 
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
  labs(fill = "", colour = '') 

#ggsave(paste(figs.out,"sp_percentage_umap_v2",'.png'))
#ggsave(paste(figs.out,"sp_percentage_umap_nal_cipro",'.png'))
#ggsave(paste(figs.out,"sp_percentage_umap_30pcs",'.png'))

ggplot(umap_out2,aes(x = `UMAP_1`, y = `UMAP_2`,fill = `phage`, color = `phage`)) + 
  #  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "Sp\U03B2 (%)", y = "PbsX (%)") +
  geom_point(pch = 21, size =1.0, alpha = 1.8) +  labs(x = "", y = "") +
  #  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
  #  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  #  theme_bw() + xlim(-2,100) + ylim(-2,100)+
  scale_fill_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+ guides(colour = guide_legend(override.aes = list(size=1)))+
  scale_color_manual(values=c('burlywood4','Gray','#238B45','#EF3B2C'))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        #        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        #        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 
ggsave(paste(figs.out,"nal_cipro_umap_phage_plot",'.png'))

