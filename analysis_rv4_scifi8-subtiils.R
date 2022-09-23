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
figs.out = '~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/scifiseq_figs//'
dir.create(figs.out)
#Colors
Breaks=seq(0,60,1)
#dir.create('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/slides/Figures/2020-07-20/')
Colors=rev(brewer.pal(11,"Spectral"))
colors=colorRampPalette(Colors)(120)
#Exponential + Stationary
#data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/petri_seq_with_anti.csv')
data = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/expression_filtered_scifi8_NISSLE_BS168_25.csv')
cols = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/gene_index_scifi8_NISSLE_BS168_25.csv')
annotation = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/cell_index_scifi8_NISSLE_BS168_25.csv')
supp_annotation = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/scBactoSeq/data/Barcodes/bw_scifi8_2samples_supplemental.csv')[c('LIBRARY_NAME','BARCODE_1')]
annotation = merge(annotation,supp_annotation,by.x = 'r1',by.y = 'BARCODE_1')
annotation$treatment = annotation$LIBRARY_NAME
annotation = annotation[order(annotation$cell_index),]
#colnames(data) = cols[,1]
annotation$treatment = gsub('KLENOW_Pre_Lib1_Tn5_BS168_Nissle_','Post',annotation$treatment)
data['Label'] = annotation$treatment
data['Identity'] = annotation$identity
#data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/Data/BS168_genes_no_ribo_stationary.csv')
"Label" %in% names(data)
#data = data[,-1]
#data = data[,-(dim(data)[2])]
cols = names(data)

### Make this only E. BS168 ####
new_names = names(data)
ec_names = new_names[str_detect(new_names,'BS_')]
data = data[ec_names]

data['Label'] = annotation$treatment
data['Identity'] = annotation$identity
data2 = data
keep_rows  = which(data2$Identity == 'BS168')
#keep_rows  = which(data$Identity == 'BS168' & data$Label !="CRISPRI" & data$Label != "M9"& data$Label != "30MIN_FIX")
data2 = data2[keep_rows,]
#### Second data import ####
#data_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/expression_filtered_scifi7_len40.csv')
#cols_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/gene_index_scifi7_len40.csv')
#annotation_past = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/cell_index_scifi7_len40.csv')
#annotation_past$treatment = str_replace(annotation_past$treatment,'2HR_CEFOZOLIN','2hr_Cipro')
#annotation_past$treatment = str_replace(annotation_past$treatment,'2HR_CIPRO','2hr_Cef')
#annotation_past$treatment = str_replace(annotation_past$treatment,'4HR_CEFOZOLIN','4hr_Cef')
##colnames(data) = cols[,1]
#annotation_past$treatment = gsub('KLENOW_Lib1_Tn5_','Past_',annotation_past$treatment)
#data_past['Label'] = annotation_past$treatment
#data_past['Identity'] = annotation_past$identity
##data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/Data/subtilis_genes_no_ribo_stationary.csv')
#cols_past = names(data_past)
## Make minimum 5 observations
#
#### Make this only  B. subtilis ####
##ec_names = new_names[str_detect(new_names,'BS_')]
#data_past = data_past[ec_names]
#data_past['Label'] = annotation_past$treatment
#data_past['Identity'] = annotation_past$identity
#keep_rows  = which(data_past$Identity == 'subtilis')
##keep_rows  = which(data$Identity == 'coli' & data$Label !="CRISPRI" & data$Label != "M9"& data$Label != "30MIN_FIX")
#data_past = data_past[keep_rows,]
##keep_rows  = which(data$Identity == 'BS168'  & data$Label != "PENTA_EXP"& data$Label != "30MIN_FIX")
#unique(data_past$Label)
#data_total = rbind(data2,data_past)
data_total = rbind(data2)
#write.csv(data_total,file='~/Documents/Data/subtilis_scifi7_8.csv')

unique(data_total$Label)

#keep_rows  = which(data2$Label =="Exponential" | data2$Label =="Stationary" |  data2$Label == "3hr_Cipro_Rep1"|data2$Label=="3hr_Cef_Rep1")
#keep_rows  = which(data_total$Label =="PostExponential" | data_total$Label =="PostStationary" |  data_total$Label == "Post3hr_Cipro_Rep1"|data_total$Label=="Past_EXP"|data_total$Label=="Past_STATIONARY")
#keep_rows  = which(data_total$Label =="PostExponential" | data_total$Label =="PostStationary" |data_total$Label=="Past_EaXP"|data_total$Label=="Past_aSTATIONARY")
#keep_rows  = which(data_total$Label =="PostExponential" )
keep_rows  = which(data_total$Label =="Post8hr_Cipro_Rep1" | data_total$Label == "Post8hr_Cipro_Rep2" )
#keep_rows  = which(data2$Label =="STATIONARY" | data2$Label == "EXP")
#keep_rows  = which(data2$Identity == 'BS168')
#data = data2[keep_rows,]
data = data_total
data = data[keep_rows,]

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
k = colSums(data[-which(names(data) %in% c('Label','Identity','BS-BSU-04345'))])
to_filter = which(k > 10)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data = data[to_filter,]
#data = data[keep_rows,]
## Make Seurat

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
sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
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

print(sct[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sct, dims = 1:2, reduction = "pca",nfeatures = 20, balanced = TRUE)
#ggsave(paste(figs.out,'dim_loadings_BS168','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/pca_loadings.png')
#DimHeatmap(sct, dims = 1:3, cells = 500, balanced = TRUE)
ElbowPlot(sct,ndims = 15)


#mat <- Seurat::GetAssayData(sct, assay = "SCT", slot = "scale.data")
mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
ggplot(data = data.frame(cbind(1:length(varExplained), varExplained)), aes (x = V1, y = varExplained)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_BS168_var_explained','.png'))

#sct = JackStraw(sct)

#data_to_write_out <- t(as.data.frame(as.matrix(sct@assays$SCT@scale.data)))
#fwrite(x = data_to_write_out, row.names = FALSE, file = "~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scale_transformed.csv")

names(data)[str_detect(names(data),'xhl')]

DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_BS168_scifi7','.png'))
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/PCA_Relative_Count.png')

sct <- RunUMAP(sct, dims = 1:10, verbose = FALSE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:10, verbose = TRUE)
#sct <- FindClusters(sct, verbose = TRUE, algorithm=1,resolution=0.9)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.6)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_BS168_scifi7','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scBactoSeq/Figures/2020-05-12/UMAP_Relative_Count.png')
#DimPlot(sct, label = TRUE, reduction = 'tsne')
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_BS168_scifi3','.png'))
DimPlot(sct, label = TRUE ,  reduction = 'pca') + ggtitle('Clusters (PCA)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
g = Idents(sct)
g[which(g == 2)]
#cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1)
#cluster1.markers <- FindMarkers(sct, ident.1 = 1, ident.2=0,min.pct = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.1)
#cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5,min.pct = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1)
#cluster6.markers <- FindMarkers(sct, ident.1 = 6,ident.2=4, min.pct = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.1)
cluster8.markers <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1)
cluster9.markers <- FindMarkers(sct, ident.1 = 9, min.pct = 0.1)
cluster10.markers <- FindMarkers(sct, ident.1 = 10, min.pct = 0.1)
cluster11.markers <- FindMarkers(sct, ident.1 = 11, min.pct = 0.1)

cluster14.markers <- FindMarkers(sct, ident.1 = 14, min.pct = 0.1)
cluster13.markers <- FindMarkers(sct, ident.1 = 13, min.pct = 0.1)
cluster12.markers <- FindMarkers(sct, ident.1 = 12, min.pct = 0.1)
exp_diff = FindMarkers(sct, ident.1 = 0,ident.2 = 4)
stat_diff = FindMarkers(sct, ident.1 = 0,ident.2 = 11)
sch_diff = FindMarkers(sct, ident.1 = 1,ident.2 = 10)
#FeaturePlot(sct,features = c('1','guaBA','gadBC','hipBA'), reduction = 'umap')
#FeaturePlot(sct,features = c('BS-ykuO','BS-sigH','BS-suhB','BS-hemE'), reduction = 'umap')
#FeaturePlot(sct,features = c('BS-tufA','BS-putA','BS-sdiA','BS-bcsA'), reduction = 'umap')
#FeaturePlot(sct,features = c('BS-rplP','BS-dmlA','BS-sdiA','BS-tap'), reduction = 'umap',size  = 4)
#ggsave(paste(figs.out,'BS168_geneplot_scifi3','.png'))
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
#ggsave(paste(figs.out,'umap_clustered_BS168_scifi3','.png'))

#FeaturePlot(sct,features =  c('acoA'), reduction = 'umap') & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-gadB','BS-aceA','BS-ssrA','BS-sulA'), reduction = 'pca')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-gapB','BS-mntA','BS-floT','BS-pftA'), reduction = 'pca')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-yddK','BS-recA','BS-spo0M','BS-xkdG'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
greens = (brewer.pal(n = 11, name = "Greens"))
reds = (brewer.pal(n = 11, name = "Reds"))
p = FeaturePlot(object = sct,
            features = c("BS-xtmA", "BS-yonO"),
#           cols = c("white", "#00441B", "#67000D", "pink"),  
           cols = c(gray(level = 0.90), greens[7], reds[7], "pink"),  
          min.cutoff = "0",max.cutoff = "q90",pt.size = 1.6,
            reduction = "umap", order = TRUE, blend.threshold = 0.01,combine = FALSE,
            blend = TRUE)
p[[3]] & NoAxes()
#&DarkTheme()
ggsave(paste(figs.out,'BS168_geneplot_phages','.png'))
#
p[[4]] & labs(x = '',y='') &theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) 
ggsave(paste(figs.out,'BS168_geneplot_color_map','.png'))
#&DarkTheme() 

FeaturePlot(sct,features =  c('BS-dcuA','BS-tnaA','BS-manA','BS-malK'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-spxA','BS-yhgE','BS-nfeDB','BS-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-xtmA','BS-yonO','BS-cwlC','BS-yosX'), reduction = 'umap',order=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(0,5))
#FeaturePlot(sct,features =  c('BS-xtmA','BS-yonO'), reduction = 'umap',order=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-xtmA','BS-yonO'), reduction = 'umap',order=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Greens")))
m1 = FeaturePlot(sct,features =  c('BS-xtmA'), reduction = 'umap',order=TRUE)& scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "Greens")))
m2 = FeaturePlot(sct,features =  c('BS-yonO'), reduction = 'umap',order=TRUE)& scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "Reds")))
#ggsave(paste(figs.out,'BS168_geneplot_exp_stationary_with_batches_umap','.png'))
ggsave(paste(figs.out,'BS168_geneplot_phages','.png'))
FeaturePlot(sct,features =  c('BS-gapB','BS-spo0M','BS-cwlC','BS-fusA'), reduction = 'pca',order=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_exp_stationary_pca_with_baatches','.png'))
FeaturePlot(sct,features =  c('BS-dcuA','BS-tnaA','BS-sppA','BS-malK'), reduction = 'pca')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-spxA','BS-yhgE','BS-nfeDB','BS-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_exp_stationary_pca','.png'))
markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
library(dplyr)

my_levels <- c(0,2,3,1,4,5)
levels(sct) = my_levels
markers %>%
  group_by(cluster) %>%  top_n(n = 2, wt = avg_log2FC)


markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(sct, features = top10$gene,draw.lines = TRUE)  + scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ 
  theme(text = element_text(size = 18))
ggsave(paste(figs.out,'BS168_heatmap_exp_3hr','.png'))
#### PCA ####


pca_out = data.frame(sct[['pca']]@cell.embeddings)
pca_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
pca_out$label = as.character(sct[['Label']][,1])
pca_out$ident = as.character(sct@active.ident)
library(scales)
colors = hue_pal()(length(unique(pca_out$label)))
col_key = data.frame(cbind(colors, unique(pca_out$label)))
colnames(col_key) = c('color','label')
pca_out$label = gsub('Past_EXP','Batch1_OD_0.7',pca_out$label)
pca_out$label = gsub('Past_STATIONARY','Batch1_OD_1.65',pca_out$label)
pca_out$label = gsub('PostExponential','Batch2_OD_0.6',pca_out$label)
pca_out$label = gsub('PostStationary','Batch2_OD_2.4',pca_out$label)
#pca_out = merge(pca_out,col_key,by.x = 'label', by.y = 'label')
#plot3d(pca_out[,2],pca_out[,3],pca_out[,4], col = as.vector(pca_out$color),xlab = 'PC1',ylab = 'PC2',
#       zlab = 'PC3',bty = "b2", main ="BS168 PCA",size = 3,alpha = 0.8,type = 'p', specular="black" )
#legend3d("topright", legend = paste(col_key$label), pch = 16, col = col_key$color, cex=2, inset=c(0.08))


p1 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`)) + 
  geom_point(pch = 21, aes(fill = umi_counts),size = 1.2) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMI counts - B. subtilis")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
p2 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`,col = `label`)) + 
  geom_point( size = 1.2) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("Exp/Stationary- 2 Batches - B. subtilis")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
p1 + p2
p2
ggsave(paste(figs.out, 'pca_umi_counts_BS168_EXP_with_bathces.png'))
#pca_out2 = pca_out[pca_out$label=='STATIONARY',]

p1 = ggplot(pca_out2,aes(x = `PC_1`, y = `PC_2`)) + 
  geom_point(pch = 21, aes(fill = umi_counts),size = 1.2) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMI counts - E. BS168")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
p1
ggsave(paste(figs.out, 'pca_BS168_EXP.png'))
p3 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("PCA for Cluster")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Paired")+ guides(colour = guide_legend(override.aes = list(size=10)))
p3
ggsave(paste(figs.out, 'pca_cluster_BS168_EXP.png'))

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
umap_out$label = gsub('3hr Cipro Rep1','90 min Cipro',umap_out$label)
umap_out$label = gsub('3hr Cef Rep1','90 min Cef',umap_out$label)
umap_out$label = gsub('90 min','T90',umap_out$label)
umap_out$label = gsub('8hr','T360',umap_out$label)
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
  ggtitle(sprintf("B. subtilis")) +
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
#  scale_color_manual(values=as.vector(palette.colors()[c(8,3,4,5,6,7)]))
  
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
p2
ggsave(paste(figs.out, 'umap_subtilis_all_label.png',sep=''))
p3
ggsave(paste(figs.out, 'umap_cluster_BS168_total.png'))
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)
g = ggplot_build(p2)
colors = g$data[[1]]['colour']
color_matrix = cbind(colors, umap_out$label)
color_matrix = unique(color_matrix)
write.csv(color_matrix,paste(figs.out,'color_map.csv'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 1,alpha = 1.5) +  
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
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
#  scale_color_simpsons()
    scale_color_manual(values=as.vector(palette.colors()[c(2,9,4,5,6,8)]))
  
  #scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,10,5)])
ggsave(paste(figs.out, 'umap_BS168_clean_palette1.png'))
#ggsave(paste(figs.out, 'umap_BS168_clean_palette1_30pcs.png'))


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 1,alpha = 1.5) +  
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
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
  #  scale_color_simpsons()s
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,10,5)])
#ggsave(paste(figs.out, 'umap_BS168_clean_palette2.png'))
ggsave(paste(figs.out, 'umap_BS168_clean_palette2_30_pcs.png'))


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +   labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("B. subtilis")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + 
  scale_color_manual(values=as.vector(palette.colors()))+ guides(colour = guide_legend(override.aes = list(size=8)))
ggsave(paste(figs.out, 'umap_BS168_cluster_clean.png'))

#### Cefazolin ####
data2 = data_total 
unique(data2$Label)
phage_anno = data.frame(read.csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/scBactoSeq/data/spbeta_genes.csv'))
#phage_anno$gene_id = gsub('_','-',phage_anno$gene_id)a
#keep_rows  = which(data2$Label =="4HR_COCULTURE" | data2$Label =="STATIONARY" | data2$Label == "Coculture_EXP"| data2$Label == "EXP")
#keep_rows  = which(data2$Label =="Post3hr_Cipro_Rep1" | data2$Label == "Post8hr_Cipro_Rep2"| data2$Label == "Post8hr_Cipro_Rep1" | data2$Label =="PostExponential"|data2$Label =="PostStationary")
#keep_rows  = which(data2$Label =="Post3hr_Cipro_Rep1" | data2$Label == "Post8hr_Cipro_Rep1"| data2$Label == "Post8hr_Cipro_Rep2")
##keep_rows  = which(data2$Label =="Post3hr_Cipro_Rep1"  | data2$Label == "Pastd_EXP"| data2$Label == "Past_4HdR_CIPRO")
#keep_rows  = which(data2$Label =="Post3hr_Cipro_Rep1"  | data2$Label == "Pastd_EXP"| data2$Label == "Past_4HdR_CIPRO")
#keep_rows  = which(data2$Identity == 'BS168')
#data = data2[keep_rows,]
data = data2
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


# Filter out the operons with less than 5 obs total (smaller sample sizie)
k = colSums(data[-which(names(data) %in% c('Label','Identity'))])
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
#phage_anno3 = phage_anno[phage_anno$pbsx_score == 1,]
phage_anno3 = phage_anno
#new_names = setdiff(new_names,phage_anno3$gene_id)
data = data[new_names]
k2 = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
#to_filter = which(k2 >= 10)
to_filter = which(k2 >=15 )
data = data[to_filter,]
## Make Seurat

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_BSU_04345'))]
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
#sct = NormalizeData(sc, scale.factor = 100)
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

sct <- RunPCA(sct, verbose = FALSE,features = all.genes)
#FindMarkers(sct, ident.1 = 'Exponential')

print(sct[["pca"]], dims = 1:10, nfeatures = 5)
VizDimLoadings(sct, dims = c(1,3), reduction = "pca",nfeatures = 22, balanced = FALSE)
#ggsave(paste(figs.out,'dim_loadings_BS168','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/pca_loadings.png')
#DimHeatmap(sct, dims = 1:3, cells = 500, balanced = TRUE)
ElbowPlot(sct,ndims = 15)


#mat <- Seurat::GetAssayData(sct, assay = "SCT", slot = "scale.data")
mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
ggplot(data = data.frame(cbind(1:length(varExplained), varExplained)), aes (x = V1, y = varExplained)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_BS168_var_explained','.png'))

#sct = JackStraw(sct)

#data_to_write_out <- t(as.data.frame(as.matrix(sct@assays$SCT@scale.data)))
#fwrite(x = data_to_write_out, row.names = FALSE, file = "~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scale_transformed.csv")

names(data)[str_detect(names(data),'spo')]

DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_BS168_scifi7','.png'))
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(1,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/PCA_Relative_Count.png')

sct <- RunUMAP(sct, dims = 1:40, verbose = FALSE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:40, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=4,resolution = 0.5)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_BS168_scifi7','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scBactoSeq/Figures/2020-05-12/UMAP_Relative_Count.png')
#DimPlot(sct, label = TRUE, reduction = 'tsne')
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_BS168_scifi3','.png'))
DimPlot(sct, label = TRUE ,  reduction = 'pca',dims=c(1,3)) + ggtitle('Clusters (PCA)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
g = Idents(sct)
g[which(g == 2)]
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, ident.2=0,min.pct = 0.1)
#cluster1.markers <- FindMarkers(sct, ident.1 = 1,ident.2=0 ,min.pct = 0.1)
#cluster2.markers <- FindMarkers(sct, ident.1 = 2, ident.2 = 3,min.pct = 0.05)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.05)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05)
cluster3.markers <- FindMarkers(sct, ident.1 = 3,ident.2=1, min.pct = 0.1)
#pstS <- FindMarkers(sct, ident.1 = 4, min.pct = 0.1,ident.2 = 1)
#cluster4.markers <- FindMarkers(sct, ident.1 = 4,ident.2=0 ,min.pct = 0.05)
cluster4.markers <- FindMarkers(sct, ident.1 = 4,min.pct = 0.05)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.05,ident.2=4)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, ident.2 = 0,min.pct = 0.05)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6,ident.2=2, min.pct = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.1)
cluster8.markers <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1)
cluster9.markers <- FindMarkers(sct, ident.1 = 9, min.pct = 0.1)
cluster10.markers <- FindMarkers(sct, ident.1 = 10, min.pct = 0.1)
cluster11.markers <- FindMarkers(sct, ident.1 = 11, min.pct = 0.1)
cluster14.markers <- FindMarkers(sct, ident.1 = 14, min.pct = 0.1)
cluster13.markers <- FindMarkers(sct, ident.1 = 13, min.pct = 0.1)
cluster12.markers <- FindMarkers(sct, ident.1 = 12, min.pct = 0.1)
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
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca',dims=c(1,5)) + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,reduction = 'pca',dims=c(1,4)) + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_BS168_scifi3','.png'))
names(data)[str_detect(names(data),'yku')]
#FeaturePlot(sct,features =  c('acoA'), reduction = 'umap') & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-recX','BS-recA','BS-yonO','BS-xkdG'), reduction = 'umap',order=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-xtmA','BS-rnjB','BS-dutB','BS-conQ'), reduction = 'umap',order=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-xtmA','BS-xkdB','BS-xkdE','BS-xpf'), reduction = 'pca')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_cipro_pbsX','.png'))
FeaturePlot(sct,features =  c('BS-ykuG','BS-ybfG','BS-xkdD'), reduction = 'pca',order=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-yonF','BS-xkdN','BS-yonO','BS-xpf'), reduction = 'pca',order = TRUE, dims=c(1,2))& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#yobI -- pifA abortive infection system. 
FeaturePlot(sct,features =  c('BS-rpoE','BS-yopP'), reduction = 'pca',order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-uvrX','BS-yoqL'), reduction = 'pca',order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-pghB','BS-ydcL'), reduction = 'pca',order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
# Doesn't seem llike induction here requires aimX..? Don't see much of it
#FeaturePlot(sct,features =  c('BS-yoqL','BS-yorS','BS-yonO','BS-yomX'), reduction = 'pca',dims=c(1,2),order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#xre binds to malP region... why..? jk that's comopletely wrong
#xkdA is an immA analogue, xre is the repressor
FeaturePlot(sct,features =  c('BS-yoqL','BS-yonO','BS-xkdA','BS-xtmA'), reduction = 'pca',dims=c(1,4),order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-yomX','BS-yonO','BS-yoyJ','BS-ykuG'), reduction = 'pca',dims=c(2,4),order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-yoqL','BS-yoqH','BS-yorF','BS-yosL'), reduction = 'pca',dims=c(1,4),order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-ywqAA','BS-yosL','BS-yorF','BS-yonF'), reduction = 'pca',dims=c(1,4),order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_cipro','.png'))
FeaturePlot(sct,features =  c('BS-nrdEB','BS-ybfG','BS-yorF','BS-yonF'), reduction = 'pca',dims=c(1,4),order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-xtmA','BS-ybfG','BS-yorF','BS-yoqH'), reduction = 'pca',dims=c(1,2),order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_cipro_xtmA','.png'))
FeaturePlot(sct,features =  c('BS-yoqL','BS-yorC','BS-yorF','BS-yomX'), reduction = 'pca',dims=c(4,5),order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-yoqL','BS-xtmA','BS-yorF','BS-yomX'), reduction = 'pca',dims=c(1,4),order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_cipro_xtmA','.png'))
VizDimLoadings(sct, dims = c(1,4), reduction = "pca",nfeatures = 30, balanced = FALSE)
FeaturePlot(sct,features =  c('BS-xis','BS-yoqC','BS-xkdF','BS-ybfG'), dims=c(1,3),reduction = 'pca',order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-ydiF','BS-BSU-21058'), reduction = 'pca',order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-yuaI','BS-gapB','BS-groES','BS-cwlC'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('BS-spxA','BS-yhgE','BS-nfeDB','BS-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'BS168_geneplot_cipro','.png'))
#mylevels  =  c(0,2,1,3,7,4,6,5)
#levels(sct) = mylevels
markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
library(dplyr)

markers %>%
  group_by(cluster) %>%  top_n(n = 2, wt = avg_log2FC)

markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
library(R.utils)
new_genes = insert(top10$gene,ats = 5,values = c('BS-nicK','BS-BSU-21058'))
new_genes =top10$gene
DoHeatmap(sct, features = new_genes,draw.lines = TRUE)  + scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+theme(text = element_text(size = 20))
ggsave(paste(figs.out,'BS168_heatmap_exp_heat_cef','.png'))

#### PCA ####
#sct@reductions$pca@cell.embeddings[,2] = -sct@reductions$pca@cell.embeddings[,2]

pca_out = data.frame(sct[['pca']]@cell.embeddings)
pca_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
pca_out$label = as.character(sct[['Label']][,1])
pca_out$ident = as.character(sct@active.ident)
pca_out$index = 1:dim(pca_out)[1]
pca_1 = pca_out[pca_out$ident == 3,]
pca_2 = pca_out[pca_out$ident == 2,]
pca_1 = pca_out[which(pca_out$PC_2 >1),]
pca_2 = pca_out[which(pca_out$PC_1 < -1),]
pca_0 = pca_out[pca_out$ident == 0,]
hist(pca_2$PC_2)
hist(pca_1$PC_1)
look_into = pca_2[which(pca_2$PC_2 > 1),]$index
look_into0 = pca_0[which(pca_0$PC_2 > 0.5),]$index
levels(sct@active.ident) = c(levels(sct@active.ident),"4")
levels(sct@active.ident) = c(levels(sct@active.ident),"5")
sct@active.ident[look_into] = "4"
sct@active.ident[look_into0] = "5"
norm_data = as.data.frame(t(as.data.frame(sct@assays$RNA@data)))
further = norm_data[look_into,]
see = colSums(further)
library(scales)
#colors = hue_pal()(length(unique(pca_out$label)))
#col_key = data.frame(cbind(colors, unique(pca_out$label)))
colors = hue_pal()(length(unique(pca_out$ident)))
colors = brewer.pal(n = length(unique(pca_out$ident)), name = "Set2")
col_key = data.frame(cbind(colors, unique(pca_out$ident)))
colnames(col_key) = c('color','label')
#pca_out = merge(pca_out,col_key,by.x = 'label', by.y = 'label')
pca_out = merge(pca_out,col_key,by.x = 'ident', by.y = 'label')
library(rgl)
plot3d(pca_out[,2],pca_out[,3],pca_out[,5], col = as.vector(pca_out$color),xlab = 'PC1',ylab = 'PC2',
       zlab = 'PC4',bty = "b2", main ="BS168 PCA",size = 10,alpha = 0.8,type = 'p', specular="black" )
legend3d("topright", legend = paste(col_key$label), pch = 16, col = col_key$color, cex=2, inset=c(0.08))


p1 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_4`)) + 
  geom_point(pch = 21, aes(fill = umi_counts),size = 1) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
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
  geom_point( size = 1) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("Class - E. BS168")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
p1 + p2
p2
ggsave(paste(figs.out, 'pca_umi_counts_BS168_cipro_pca.png'))
p3 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`,col = `ident`)) + 
  geom_point( size = 2) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("Clusters - Subtilis with Cipro")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-10,4) + ylim(-2,4)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="scales")+ guides(colour = guide_legend(override.aes = list(size=10)))
p1
p3
ggsave(paste(figs.out, 'pca_BS168_ciproo3hr.png'))

#### UMAP Variability ###

#### Variablee Genes


umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)

p1 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`)) + 
  geom_point(pch = 21, aes(fill = umi_counts),size = 1) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMAP for UMI counts")) +
  #scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) )+
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
p2 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 1,alpha = 1.0) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMAP -- 3hr Cipro")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 1.0) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2") + guides(colour = guide_legend(override.aes = list(size=10)))
p3 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMAP for Cluster")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 1.0) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
p1 + p2
p2
ggsave(paste(figs.out, 'umap_umi_counts_BS168_drug.png'))
p3
ggsave(paste(figs.out, 'umap_cluster_BS168_ciproo3hr.png'))
sizes = aggregate(umi_counts~ident,umap_out,length)
total = sum(aggregate(umi_counts~ident,umap_out,length)['umi_counts'])
sizes['umi_counts']/total

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
phage_data = data.frame(t(sc_precursor) %*% phage_score)
phage_data2 = data.frame(t(sc_precursor) %*% phage_score2)
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
p1 + p2
ggsave(paste(figs.out,'stage_plots','.png'))

p3 = ggplot(phage_data,aes(x = `norm_pbs`, y = `norm_sp`,fill = `label`)) + 
  geom_point(pch = 21, size =2 ) +  labs(y = "Sp\U03B2 (%)", x = "PbsX (%)") +
  ggtitle(sprintf("Sp\U03B2 vs PbsX percentages")) +
#  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + xlim(-2,90) + ylim(-2,100)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=10)))+ labs(fill = "Cluster") 

p3
ggsave(paste(figs.out,'phages_overlaid_cluster','.png'))



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
ggsave(paste(figs.out,'normed_phages_overlaid_cluster','.png'))

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


ggplot(phage_data,aes(x = `x1_stage`, y = `x2_stage`,col = `ident`)) + 
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

umap_out2 = cbind(umap_out, phage_data$norm_pbs,phage_data$norm_sp)
umap_out2$max_val = pmax(phage_data$norm_pbs,phage_data$norm_sp)
temp_phage = matrix(phage_data[c('norm_pbs','norm_sp')],ncol=2)
colnames(temp_phage) = paste0('PHAGE_',1:2)
sct[["phage"]] <- CreateDimReducObject(embeddings = temp_phage, key = "PHAGE_", assay = DefaultAssay(sct))

#look_into = intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_pbs2 > 0.02))
#look_into = union(which(phage_data$norm_sp2 > 0.05),which(phage_data$norm_pbs2 > 0.05))
#look_into1 = which(phage_data$norm_sp2 > 0.5)
#look_into2 = intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_sp2 < 0.5))
#look_into3 = intersect(intersect(which(phage_data$norm_sp2 > 0.1),which(phage_data$norm_pbs2 > 0.02)),intersect(which(phage_data$norm_sp2 > 0.05),which(phage_data$norm_sp2 < 0.5)))
#look_into4 =  which(phage_data$norm_pbs2 > 0.1)
#levels(sct@active.ident) = c(levels(sct@active.ident),"6","7",'8','9','10')
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
ggplot(phage_data5%>%arrange(`BS-yoqH`),aes(x = `norm_pbs`, y = `norm_sp`)) + 
  geom_point(pch = 21, aes(fill = `BS-yonO`),size =3 ) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
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
#ggsave(paste(figs.out,"pbsx_percentage_umap",'.png'))
ggsave(paste(figs.out,"pbsx_percentage_umap_30pcs",'.png'))
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

ggsave(paste(figs.out,"sp_percentage_umap_v2",'.png'))
#ggsave(paste(figs.out,"sp_percentage_umap_30pcs",'.png'))

FeaturePlot(sct,features =  c('BS-yotJ','BS-sspC'), reduction = 'pca',dims=c(1,2),order = TRUE,pt.size=2)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('BS-rnjB'), reduction = 'pca',order = TRUE,pt.size=4)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'rnjB_featureplot','.png'))
names(data)[str_detect(names(data),'ykz')]

pca_out$ident = as.character(sct@active.ident)

p3 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`,col = `ident`)) + 
  geom_point( size = 2) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("New Clusters - Subtilis with Cipro")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() + #xlim(-10,4) + ylim(-2,4)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
p3
ggsave(paste(figs.out, 'pca_BS168_ciproo3hr_new.png'))

#### Osmootic  Shoock and stuff now ####

data2 = data_total
unique(data2$Label)
#keep_rows  = which(data2$Label =="4HR_COCULTURE" | data2$Label =="STATIONARY" | data2$Label == "Coculture_EXP"| data2$Label == "EXP")
keep_rows  = which( data2$Label == "Post8hr_Cipro_Rep1"| data2$Label == "Past_2hr_Cipro" | data2$Label =="EXP" | data2$Label =="STATIONARY")
#keep_rows  = which(data2$Identity == 'BS168')
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

sc_precursor = data[,-which(names(data) %in% c('Label','Identity','BS_BSU_04345'))]
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
VizDimLoadings(sct, dims = 1:2, reduction = "pca",nfeatures = 20, balanced = TRUE)
#ggsave(paste(figs.out,'dim_loadings_BS168','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/pca_loadings.png')
#DimHeatmap(sct, dims = 1:3, cells = 500, balanced = TRUE)
ElbowPlot(sct,ndims = 15)


#mat <- Seurat::GetAssayData(sct, assay = "SCT", slot = "scale.data")
mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
ggplot(data = data.frame(cbind(1:length(varExplained), varExplained)), aes (x = V1, y = varExplained)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_BS168_var_explained','.png'))

#sct = JackStraw(sct)

#data_to_write_out <- t(as.data.frame(as.matrix(sct@assays$SCT@scale.data)))
#fwrite(x = data_to_write_out, row.names = FALSE, file = "~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scale_transformed.csv")

names(data)[str_detect(names(data),'yok')]

DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_BS168_scifi7','.png'))
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/PCA_Relative_Count.png')

sct <- RunUMAP(sct, dims = 1:10, verbose = FALSE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:10, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=1,resolution = 0.5)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_BS168_scifi7','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scBactoSeq/Figures/2020-05-12/UMAP_Relative_Count.png')
#DimPlot(sct, label = TRUE, reduction = 'tsne')
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_BS168_scifi3','.png'))
DimPlot(sct, label = TRUE ,  reduction = 'pca') + ggtitle('Clusters (PCA)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
g = Idents(sct)
g[which(g == 2)]
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, ident.2=0,min.pct = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, ident.2=0,min.pct = 0.1)
#cluster3.markers <- FindMarkers(sct, ident.1 = 3,ident.2=4, min.pct = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4,ident.2=0, min.pct = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, ident.2 = 8,min.pct = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6,ident.2=2, min.pct = 0.1)
cluster7.markers <- FindMarkers(sct, ident.1 = 7, min.pct = 0.1)
cluster8.markers <- FindMarkers(sct, ident.1 = 8, min.pct = 0.1)
cluster9.markers <- FindMarkers(sct, ident.1 = 9, min.pct = 0.1)
cluster10.markers <- FindMarkers(sct, ident.1 = 10, min.pct = 0.1)
cluster11.markers <- FindMarkers(sct, ident.1 = 11, min.pct = 0.1)
cluster14.markers <- FindMarkers(sct, ident.1 = 14, min.pct = 0.1)
cluster13.markers <- FindMarkers(sct, ident.1 = 13, min.pct = 0.1)
cluster12.markers <- FindMarkers(sct, ident.1 = 12, min.pct = 0.1)
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
FeaturePlot(sct,features =  c('BS-dcuA','BS-tnaA','BS-katG','BS-malK'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
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




#### 
library(topGO)
library(gage)
bs = kegg.gsets(species = 'bsu')
#ec = kegg.gsets(species = 'eco')
cluster9.markers




