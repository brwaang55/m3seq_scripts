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
library(reticulate)
figs.out = '~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/scifiseq_figs/'
dir.create(figs.out)
#Colors
Breaks=seq(0,60,1)
#dir.create('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/slides/Figures/2020-07-20/')
Colors=rev(brewer.pal(11,"Spectral"))
colors=colorRampPalette(Colors)(120)
#Exponential + Stationary
#data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/petri_seq_with_anti.csv')
data = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/expression_filtered_scifi8_MG1655_25.csv')
cols = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/gene_index_scifi8_MG1655_25.csv')
annotation = read_csv('/Volumes/AdamsonLab/brucewang/data/scifi3_saved_data/cell_index_scifi8_MG1655_25.csv')
#colnames(data) = cols[,1]
#colnames(data) = cols[,1]
annotation$treatment = gsub('KLENOW_Pre_Lib1_Tn5_MG1655_','Post_',annotation$treatment)
data['Label'] = annotation$treatment
data['Identity'] = annotation$identity
#data = read_csv('~/Dropbox (Princeton)/Notes/research/adamson_gitai_wingreen/Data/MG1655_genes_no_ribo_stationary.csv')
"Label" %in% names(data)
#data = data[,-1]
#### data_past ###
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
# Make minimum 5 observations

### Make this only  B. subtilis ####
ec_names = names(data)[str_detect(names(data),'EC_')]
data = data[ec_names]
data_past = data_past[ec_names]
data['Label'] = annotation$treatment
data['Identity'] = annotation$identity
#data_past['Label'] = annotation_past$treatment
#data_past['Identity'] = annotation_past$identity
#keep_rows  = which(data_past$Identity == 'coli')
#keep_rows  = which(data$Identity == 'coli' & data$Label !="CRISPRI" & data$Label != "M9"& data$Label != "30MIN_FIX")
#data_past = data_past[keep_rows,]
data2 = data
keep_rows  = which(data2$Identity == 'MG1655')
data2 = data[keep_rows,]
#keep_rows  = which(data$Identity == 'BS168'  & data$Label != "PENTA_EXP"& data$Label != "30MIN_FIX")
#unique(data_past$Label)
#data_total = rbind(data2,data_past)
data_total = rbind(data2)
#data_total = rbind(data2,data_past)
#write.csv(data_total,file='~/Documents/Data/subtilis_scifi7_8.csv')
#write.csv(data_total,file='~/Documents/Data/subtilis_coli7_8.csv')
unique(data_total$Label)
### Make this only E. MG1655 ####
data = data_total
#data = data[,-(dim(data)[2])]
cols = names(data)

#### Make counts of data ####
umis_per_cell = rowSums(data[,-which(names(data) %in% c('Label','Identity','anti-rRNA'))])
labels = data$Label
umis_per_cell_frame = data.frame(cbind(labels, umis_per_cell))
table(umis_per_cell_frame$labels)
library(dplyr)
umis_per_cell_frame %>% group_by(labels) %>% summarise(median = median(umis_per_cell))
# Make minimum 5 observations
k = colSums(data!=0)
to_filter = which(k > 5)
new_names = names(to_filter)
new_names = c(new_names)
data = data[new_names]

unique(data$Label)
lambda_names = c('EC_gad','EC_gad','EC_gadC')
#keep_rows  = which(data$Identity == 'MG1655' & data$Label !="CRISPRI" & data$Label != "M9"& data$Label != "30MIN_FIX")
#keep_rows  = which(data$Label =="Post_Exponential" | data$Label =="Post_Stationary"  | data$Label =="Post_3hr_Cipro_Rep1" )#| data$Label == "Past_STATIONARY" |data$Label == "Past_EXP")
keep_rows  = which( data$Label =="Post_Stationary" )#| data$Label == "Past_STATIONARY" |data$Label == "Past_EXP")
#keep_rows  = which(data$Label =="Post_Exponential" | data$Label =="Post_Stationary"  | data$Label =="Past_STATIONARY" )#| data$Label == "Past_STATIONARY" |data$Label == "Past_EXP")
#keep_rows  = which(data$Label =="Post_Expaonential" | data$Label =="Post_Stationary"  | data$Label =="Past_STATIaONARY" )#| data$Label == "Past_STATIONARY" |data$Label == "Past_EXP")
#keep_rows  = which(data$Label =="Post_8hr_Cipro_Rep1" )#| data$Label == "Past_STATIONARY" |data$Label == "Past_EXP")
#keep_rows  = which(data2$Label =="Stationary" | data2$Label == "3hr_Cipro_Rep1")
#keep_rows  = which( data2$Label == "3hr_Cipro_Rep1"| data2$Label == "8hr_Cipro_Rep1"| data2$Label == "8hr_Cipro_Rep2")
data = data[keep_rows,]
#keep_rows  = which(data$Identity == 'MG1655'  & data$Label != "PENTA_EXP"& data$Label != "30MIN_FIX")
unique(data2$Label)
#keep_rows  = which(data2$Label =="STATIONARY" | data2$Label == "EXP")
#keep_rows  = which(data2$Identity == 'MG1655')
#data = data2[keep_rows,]
to_filter = which(k > 10)
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
#sct = NormalizeData(sc,normalization.method = "RC", scale.factor = 100)
sct = NormalizeData(sc,normalization.method = "LogNormalize", scale.factor = 100)
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

print(sct[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sct, dims = 1:2, reduction = "pca",nfeatures = 20, balanced = TRUE)
#ggsave(paste(figs.out,'dim_loadings_MG1655','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/pca_loadings.png')
#DimHeatmap(sct, dims = 1:3, cells = 500, balanced = TRUE)
ElbowPlot(sct,ndims = 15)


#mat <- Seurat::GetAssayData(sct, assay = "SCT", slot = "scale.data")
mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance
varExplained = eigValues / sum(eigValues)
ggplot(data = data.frame(cbind(1:length(varExplained), varExplained)), aes (x = V1, y = varExplained)) +
  geom_line() + 
  geom_point() +  labs(x = "PC", y = "Variance Explained") +
  ggtitle(sprintf("PCA Variance Explained - Scaling by UMI Count")) +
  theme_bw() +
#  ylim(0,0.0225)+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_MG1655_var_explained','.png'))

#sct = JackStraw(sct)

#data_to_write_out <- t(as.data.frame(as.matrix(sct@assays$SCT@scale.data)))
#fwrite(x = data_to_write_out, row.names = FALSE, file = "~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scale_transformed.csv")

names(data)[str_detect(names(data),'adi')]

DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_MG1655_scifi7','.png'))
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/PCA_Relative_Count.png')

#1:10 normally
sct <- RunUMAP(sct, dims = 1:10, verbose = FALSE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:10, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution=0.2)
DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_MG1655_scifi7','.png'))
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scBactoSeq/Figures/2020-05-12/UMAP_Relative_Count.png')
#DimPlot(sct, label = TRUE, reduction = 'tsne')

#ggsave(paste(figs.out,'umap_clustered_MG1655_scifi3','.png'))
#DimPlot(sct, label = TRUE ,  reduction = 'pca') + ggtitle('Clusters (PCA)- Relative Counts') + theme_bw() +
#  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
g = Idents(sct)
#g[which(g == 2)]
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1, logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1, logfc.threshold = 0.1)
#cluster1.markers <- FindMarkers(sct, ident.1 = 1, ident.2=0,min.pct = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1, logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.05, logfc.threshold = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.1, logfc.threshold = 0.1)


RidgePlot(sct, features = c('EC-gadA','EC-gadB'), ncol = 2)
DotPlot(sct, features = c('EC-gadA','EC-gadB'))+  scale_colour_gradient2(low = "#FF00FF",  high = 'chartreuse1')+ylab('') + xlab('')


#### Matrixplot
library(dplyr)

av.exp.sct <- AverageExpression(sct, return.seurat = TRUE)
av.exp = Seurat::GetAssayData(av.exp.sct, assay = "RNA", slot = "scale.data")
#av.exp <- AverageExpression(sct)
av.df <- as.data.frame(t(av.exp))
#sct@active.ident = factor(sct@meta.data$Label,labels = sct@meta.data$Label, levels = sct@meta.data$Label)
#markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.5,group.by = 'Label')
markers <- FindAllMarkers(sct, min.pct = 0.05, logfc.threshold = 0.1)
markers = markers[abs(markers$avg_log2FC) > 0.1,]
#markers %>%
#  group_by(cluster) %>%
#  top_n(n = 3, wt = avg_log2FC) -> top10

library(dplyr)
#markers = markers %>%
#  add_count(cluster) %>%
#  group_by(cluster) %>%
#  mutate(n = 2) %>%
#  #  group_by(rating, .add = TRUE) %>%
#  #In old dplyr add = TRUE
#  #group_by(rating, add = TRUE) %>%
#  sample_n(n, replace = TRUE) %>%
#  select(-n)
markers =   Reduce(rbind,                                 # Top N highest values by group
       by(markers,
          markers["cluster"],
          head,
          n = 6))

markers = unique(markers)
#markers = markers[1:8,]
genes = markers$gene
av.df$cluster = rownames(av.df)
g = reshape2::melt(av.df, 'cluster')
names(g) = c("Cluster",'gene','value')
cluster0.markers
markers = markers[markers$cluster!=0, ]
g = merge(g, markers)
#g = g[which(g$variable%in% genes),]
g = g[order(g$cluster),]
#cor.df <- tidyr::gather(data = av.df , 'gene', 'expression')
g$gene = gsub(pattern = 'EC-',replacement = '', x = g$gene)
g$gene = factor(g$gene, levels = g$gene,labels = g$gene)
#g$Cluster = factor(g$Cluster, levels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'),labels = c("Gent", "Tet", "Cef","Cyclo","Erythro","Chlor",'Cipro','Nal','Exponential','Overnight'))
g$Cluster = as.numeric(as.character(g$Cluster))+1
g$Cluster = factor(g$Cluster, levels = c(3,2,1))
#g$cluster = factor(g$cluster + 1)
g = g[order(g$cluster),]

library(myriad)
ggplot(g, aes(x = gene, y= Cluster, fill = value)) +
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


#ggsave(paste(figs.out,'stationary_matrixplot_with_words','.png'))
ggsave(paste(figs.out,'stationary_matrixplot_with_words_lognorm','.png'))
ggsave(paste(figs.out,"go_term_stationary_coli",'.svg'),  device = "svg")




#FeaturePlot(sct, features = c('EC-rplF','EC-rpoC', 'rplG'), order = TRUE)
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

#ggsave(paste(figs.out,'stationary_matrixplot_with_no_words','.png'))
ggsave(paste(figs.out,'stationary_matrixplot_with_no_words_lognorm','.png'))
##

#g[which(g == 2)]
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1, logfc.threshold = 0.05)
#cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.05)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1, logfc.threshold = 0.1)
#cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1, logfc.threshold = 0.1)
#cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1, logfc.threshold = 0.1)
##cluster1.markers <- FindMarkers(sct, ident.1 = 1, ident.2=0,min.pct = 0.1)
#cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1, logfc.threshold = 0.1)
library(topGO)
ontology = 'BP'
# install.packages("org.Mm.eg.db") 
# Take from the markers
#cluster2.marker_subset = cluster2.markers[which(cluster2.markers$avg_log2FC > 0),]
#rownames(cluster2.marker_subset) = gsub(pattern = 'EC-',replacement = '',x=rownames(cluster2.marker_subset))
expressed.genes <- rownames(cluster1.markers[intersect(which(cluster1.markers$avg_log2FC > 0), 
                            which(cluster1.markers$p_val< 0.05)),])
down.genes <- rownames(cluster1.markers[intersect(which(cluster1.markers$avg_log2FC < 0), 
                                                       which(cluster1.markers$p_val< 0.05)),])
expressed.genes = gsub(pattern = 'EC-',replacement = '',x = expressed.genes)
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
#goEnrichment$Term <- factor(goEnrichment$Term, levels=(goEnrichment$Term))
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
ggplot(goEnrichment[1:9,], aes(x=Term, y=Fish, fill = `Fish`)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p value)") +
  ggtitle("") +  
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Fisher)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
        panel.grid = element_blank(),
    legend.position='none',
##    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=8, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=0.01, face="bold", hjust=1.10, color = 'black'),
    axis.text.y=element_text(angle=0, size=0.01, face="bold", vjust=0.5, color = 'black'),
    axis.title=element_text(size=0.01, face="bold"),
##    legend.key=element_blank(),     #removes the border
##    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
#    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
#ggsave(paste(figs.out,'biological_process_gad_lognorm_v2','.png'), dpi = 300)
ggsave(paste(figs.out,'biological_process_gad_lognorm_v2','.svg'), dpi = 300, device = 'svg')
ggplot(goEnrichment[1:9,], aes(x=Term, y=Fish, fill = `Fish`)) +
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
ggsave(paste(figs.out,'biological_process_gad_with_text_lognorm_v2','.svg'), dpi = 300, device = 'svg')

#GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
#### Cluster 3


library(topGO)
ontology = 'BP'
# install.packages("org.Mm.eg.db") 
# Take from the markers
#cluster2.marker_subset = cluster2.markers[which(cluster2.markers$avg_log2FC > 0),]
#rownames(cluster2.marker_subset) = gsub(pattern = 'EC-',replacement = '',x=rownames(cluster2.marker_subset))
expressed.genes <- rownames(cluster3.markers[intersect(which(cluster2.markers$avg_log2FC > 0), 
                                                       which(cluster2.markers$p_val< 0.05)),])
down.genes <- rownames(cluster3.markers[intersect(which(cluster2.markers$avg_log2FC < 0), 
                                                  which(cluster2.markers$p_val< 0.05)),])
expressed.genes = gsub(pattern = 'EC-',replacement = '',x = expressed.genes)
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
#goEnrichment$Term <- factor(goEnrichment$Term, levels=(goEnrichment$Term))
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
ggplot(goEnrichment, aes(x=Term, y=Fish, fill = `Fish`)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p value)") +
  ggtitle("") +  
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Fisher)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=8, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=0, face="bold", hjust=1.10, color = 'black'),
    axis.text.y=element_text(angle=0, size=10, face="bold", vjust=0.5, color = 'black'),
    axis.title=element_text(size=0, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
ggsave(paste(figs.out,'biological_process_cluster3_lognorm','.png'))
#GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
#####


#### Gad on the whole thing 
library(dplyr)
umap_out = data.frame(sct[['umap']]@cell.embeddings)

rps = names(data)[str_detect(names(data),'rps')]
rpl = names(data)[str_detect(names(data),'rpl')]
#rp_vec = c(rps, rpl)
#gad_vec = rp_vec
gad_vec = c('EC_gadA','EC_gadB','EC_gadC')
#gad_vec = c('EC_gadA','EC_gadB')
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
#umap_out2$UMAP_1 = -umap_out2$UMAP_1
#umap_out2$UMAP_2 = -umap_out2$UMAP_2
ggplot(umap_out2%>%arrange(`max_val`),aes(x = `UMAP_1`, y = `UMAP_2`, colour = `gad_data$norm_gad`, fill = `gad_data$norm_gad`)) + 
  #  geom_point(pch = 21, aes(fill = `gad_data$norm_pbs`),size =1.5) +  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  geom_point(size =0.8, alpha = 1) +
  #  labs(y = "Relative Sp\U03B2 (%)", x = "Relative PbsX (%)") +
  labs(y = "", x = "") +
  ggtitle(sprintf("")) +
  #  scale_fill_gradientn(colors =  (brewer.pal(n = 9, name =  "BuGn") ), guide = "colourbar") + 
  #  scale_colour_gradientn(colours =(brewer.pal(n = 9, name = "BuGn")), guide = '' ) + 
  #    scale_fill_gradientn(colours =(brewer.pal(n = 9, name = "Greens")) ) +   #scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greens")[1]), high = (brewer.pal(n = 9, name = "Greens")[9])) + 
  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'chartreuse1',limits = c(0, 5), oob = scales::squish) + 
  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'chartreuse1',limits = c(0, 5), oob = scales::squish) + 
#  scale_fill_gradient(  low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'chartreuse1') + 
#  scale_colour_gradient(low = (brewer.pal(n = 9, name = "Greys"))[2], high = 'chartreuse1') + 
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

ggsave(paste(figs.out,"gadAB_umap_stationary_mg1655",'.png'))

#org.EcK12.eg.db
#### PCA ####


pca_out = data.frame(sct[['pca']]@cell.embeddings)
pca_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
pca_out$label = as.character(sct[['Label']][,1])
pca_out$label = gsub('Past_EXP','Batch1_OD_0.6',pca_out$label)
pca_out$label = gsub('Past_STATIONARY','Batch1_OD_1.5',pca_out$label)
pca_out$label = gsub('Post_Exponential','Batch2_OD_0.6',pca_out$label)
pca_out$label = gsub('Post_Stationary','Batch2_OD_2.8',pca_out$label)
pca_out$ident = as.character(sct@active.ident)
library(scales)
colors = hue_pal()(length(unique(pca_out$label)))
col_key = data.frame(cbind(colors, unique(pca_out$label)))
colnames(col_key) = c('color','label')
#pca_out = merge(pca_out,col_key,by.x = 'label', by.y = 'label')
#plot3d(pca_out[,2],pca_out[,3],pca_out[,4], col = as.vector(pca_out$color),xlab = 'PC1',ylab = 'PC2',
#       zlab = 'PC3',bty = "b2", main ="MG1655 PCA",size = 3,alpha = 0.8,type = 'p', specular="black" )
#legend3d("topright", legend = paste(col_key$label), pch = 16, col = col_key$color, cex=2, inset=c(0.08))


p1 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`)) + 
  geom_point(pch = 21, aes(fill = umi_counts),size = 0.8) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMI counts - E. MG1655")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
p2 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`,col = `label`)) + 
  geom_point( size = 0.5,alpha=0.8) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("MG1655 PCA-Batch1 & 2")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
p1 + p2
p2
ggsave(paste(figs.out, 'pca_MG1655_batchesexp_stat.png'))
#pca_out2 = pca_out[pca_out$label=='STATIONARY',]

p1 = ggplot(pca_out2,aes(x = `PC_1`, y = `PC_2`)) + 
  geom_point(pch = 21, aes(fill = umi_counts),size = 1.2) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMI counts - E. MG1655")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
p1
ggsave(paste(figs.out, 'pca_umi_counts_MG1655_stat_exp.png'))
p3 = ggplot(pca_out,aes(x = `PC_3`, y = `PC_5`,col = `ident`)) + 
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
ggsave(paste(figs.out, 'pca_cluster_MG1655_exp_stat.png'))

#### UMAP Variability ###

#### Variablee Genes


umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out$label = gsub('Post_Stationary','Batch2_OD_2.8',umap_out$label)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)
umap_out$label = gsub('Post',' ',umap_out$label)
umap_out$label = gsub('3hr Cipro Rep1','90 min Cipro',umap_out$label)
umap_out$label = gsub('3hr Cef Rep1','90 min Cef',umap_out$label)
umap_out$label = gsub('90 min','T90',umap_out$label)
umap_out$label = gsub('8hr','T360',umap_out$label)

#umap_out$label = gsub('8hr','T360',umap_out$label)
#colors = data.frame(read.csv(paste(figs.out,'color_map.csv'), sep = ','))[,-1]
#umap_out = umap_out[umap_out$label=='STATIONARY',]
colors
#umap_out = merge(umap_out, colors,by.x = 'label',by.y = "umap_out.label")
library(ggsci)
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
  geom_point( size = 1,alpha = 1.5) +  
  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
#  labs(x = "", y = "", color = "Condition") +
  ggtitle(sprintf("E. coli")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),panel.grid = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
#        axis. = element_blank(),
#        axis.line = element_line(colour = "black"),
#        axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
  scale_color_manual(values =unique(umap_out$colour))
 # scale_color_simpsons()#scale_color_brewer(palette="Set2")
p3 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.4,alpha = 1.5) +  labs(x = "UMAP1", y = "UMAP2", color = "Cluster") +
  ggtitle(sprintf("UMAP for Cluster")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),panel.grid = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Paired")+ guides(colour = guide_legend(override.aes = list(size=8)))
p2
ggsave(paste(figs.out, 'umap_MG1655.png'))
p3
#ggsave(paste(figs.out, 'umap_cluster_MG1655.png'))
#ggsave(paste(figs.out, 'umap_cluster_3hr_cipro.png'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 0.6,alpha = 1.7) +  labs(x = "", y = "", color = "Cluster") +
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
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=8)))+
  scale_color_manual(values=as.vector(palette.colors()[c(2,4,3,5,8,1,6)]))+
  guides(colour = guide_legend(override.aes = list(size=8)))
ggsave(paste(figs.out, 'umap_cluster_MG1655_figure_log_normalize.png'))
DotPlot(sct, features = c('EC-gadA','EC-gadB'), scale.max = 50, col.min = -1, col.max = 1.5)+  
  scale_colour_gradient2(low = "black",  high = 'chartreuse1')+ylab('') + xlab('') +  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
                                                                                          panel.background = element_blank(), 
                                                                                          axis.text = element_blank(),
                                                                                          #        axis. = element_blank(),
                                                                                          #        axis.line = element_line(colour = "black"),
                                                                                          #        axis.text=element_text(size=20),axis.text.x=element_blank(),
                                                                                          #                                                                                          axis.ticks.x=element_blank(),
                                                                                          #                                                                                          axis.text.y=element_blank(),
                                                                                          #                                                                                          axis.ticks.y=element_blank() ,
                                                                                          axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold"))
ggsave(paste(figs.out,'MG1655_stationary_dotplot','.png'))

aggregate(umi_counts~ident,umap_out,mean)
sizes = aggregate(umi_counts~ident,umap_out,length)
total = sum(aggregate(umi_counts~ident,umap_out,length)['umi_counts'])
sizes['umi_counts']/total
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 1,alpha = 1.5) +  
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
ggsave(paste(figs.out, 'umap_MG1655_clean.png'))

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +   labs(x = "", y = "", color = "Cluster") +
  ggtitle(sprintf("E. coli")) +
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
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Paired")+ guides(colour = guide_legend(override.aes = list(size=8)))
ggsave(paste(figs.out, 'umap_MG1655_cluster_clean.png'))


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
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
  scale_color_manual(values=as.vector(palette.colors()[c(2,9,4,5,8)]))

#scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,10,5)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette1.png'))
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette_30pcs.png'))


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 1,alpha = 1.5) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
 # ggtitle(sprintf("B. subtilis")) +
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
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(1,2,3,4,5)])
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2.png'))
#ggsave(paste(figs.out, 'umap_MG1655_clean_palette2_30_pcs.png'))


#### Clusters ###

ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Condition") +
  # ggtitle(sprintf("B. subtilis")) +
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
#### Exponential & stationary ####

aggregate(EC_lpp~Label,data_total,sum)

unique(data_total$Label)
data2 = data_total
#keep_rows  = which(data2$Label =="4HR_COCULTURE" | data2$Label =="STATIONARY" | data2$Label == "Coculture_EXP"| data2$Label == "EXP")
#keep_rows  = which( data2$Label == "Post_3hr_Cipro_Rep1" | data2$Label == "Post_8hr_Cipro_Rep1" |data2$Label == "Post_8hr_Cipro_Rep2"  )
#keep_rows  = which( data2$Label == "Post_3hr_Cipro_Rep1" | data2$Label == "Post_8hr_Cipro_Rep1" |data2$Label == "Post_8hr_Cipro_Rep2"  )
#keep_rows  = which( data2$Label == "Post_3hr_Cipro_Rep1" | data2$Label == "Post_8hr_Cipro_Rep1" |data2$Label == "Post_8hr_Cipro_Rep2"  )
#keep_rows  = which( data2$Label == "Post_8hr_Cipro_Rep1" | data2$Label == "Post_8hr_Cipro_Rep2" )
keep_rows  = which( data2$Label == "Post_Stationary" | data2$Label == "Post_Exponential")
#keep_rows  = which( )
#keep_rows  = grepl('Post',data2$Label)
#keep_rows  = which(data2$Identity == 'MG1655')
data2 = data2[keep_rows,]
#data2 = data2
#### Counts/gene
#### Counts/gene
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
k = colSums(data2[-which(names(data2) %in% c('Label','Identity'))])
to_filter = which(k > 10)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data2 = data2[new_names]
k2 = rowSums(data2[,-which(names(data2) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data2 = data2[to_filter,]
## Make Seurat
aggregate(EC_lpp~Label,data2,sum)


sc_precursor = data2[,-which(names(data2) %in% c('Label','Identity','EC_ECU_04345'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data2[which(names(data2) %in% c('Label'))]
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
#FindMarkers(sct, ident.1 = 'Exponential')

print(sct[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sct, dims = 1:2, reduction = "pca",nfeatures = 20, balanced = TRUE)
FeaturePlot(sct, features=c('EC-rpsC'), order = TRUE)
ElbowPlot(sct,ndims = 15)
mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]

total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
sct <- RunUMAP(sct, dims = 1:10, verbose = TRUE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:10, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution = 0.3)
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
g = Idents(sct)
g[which(g == 2)]
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1, logfc.threshold = 0.05)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.1)
#### Print dat umap ####
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)
unique(umap_out$label)
library(rcartocolor)
#ggsave(paste(figs.out, 'umap_cluster_BS168_combined_total.png'))
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)


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
ggsave(paste(figs.out, 'umap_exponential_stationary_scifi8.png'))

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
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(7,9,2,3,8,12)])
ggsave(paste(figs.out, 'umap_exponential_stationary_scifi8.png'))
ggsave(paste(figs.out,"umap_exponential_stationary_scifi8.svg"),  device = "svg")
#### Markers ####
library(topGO)
ontology = 'BP'
# install.packages("org.Mm.eg.db") 
# Take from the markers
#cluster2.marker_subset = cluster2.markers[which(cluster2.markers$avg_log2FC > 0),]
#rownames(cluster2.marker_subset) = gsub(pattern = 'EC-',replacement = '',x=rownames(cluster2.marker_subset))
expressed.genes <- rownames(cluster1.markers[intersect(which(cluster1.markers$avg_log2FC > 0.5), 
                                                       which(cluster1.markers$p_val< 0.05)),])
down.genes <- rownames(cluster1.markers[intersect(which(cluster1.markers$avg_log2FC < 0.5), 
                                                  which(cluster1.markers$p_val< 0.05)),])
expressed.genes = gsub(pattern = 'EC-',replacement = '',x = expressed.genes)
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
#goEnrichment$Term <- factor(goEnrichment$Term, levels=(goEnrichment$Term))
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
ggplot(goEnrichment, aes(x=Term, y=Fish, fill = `Fish`)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p value)") +
  ggtitle("") +  
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Fisher)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    panel.grid = element_blank(),
    legend.position='none',
    ##    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=8, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=0.01, face="bold", hjust=1.10, color = 'black'),
    axis.text.y=element_text(angle=0, size=0.01, face="bold", vjust=0.5, color = 'black'),
    axis.title=element_text(size=0.01, face="bold"),
    ##    legend.key=element_blank(),     #removes the border
    ##    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    #    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
ggsave(paste(figs.out,'biological_process_exponential','.png'), dpi = 300)
ggplot(goEnrichment, aes(x=Term, y=Fish, fill = `Fish`)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p value)") +
  ggtitle("") +  
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Fisher)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    panel.grid = element_blank(),
    
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
ggsave(paste(figs.out,'biological_process_exponential_with_text_lognorm_v2','.png'), dpi = 300)
ggsave(paste(figs.out,"biological_process_exponential_with_text_lognorm_v2",'.svg'),  device = "svg")

#### Stationary go term ####
library(topGO)
ontology = 'BP'
# install.packages("org.Mm.eg.db") 
# Take from the markers
#cluster2.marker_subset = cluster2.markers[which(cluster2.markers$avg_log2FC > 0),]
#rownames(cluster2.marker_subset) = gsub(pattern = 'EC-',replacement = '',x=rownames(cluster2.marker_subset))
expressed.genes <- rownames(cluster0.markers[intersect(which(cluster0.markers$avg_log2FC > 0.5), 
                                                       which(cluster0.markers$p_val< 0.05)),])
down.genes <- rownames(cluster0.markers[intersect(which(cluster0.markers$avg_log2FC < 0.5), 
                                                  which(cluster0.markers$p_val< 0.05)),])
expressed.genes = gsub(pattern = 'EC-',replacement = '',x = expressed.genes)
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
#goEnrichment$Term <- factor(goEnrichment$Term, levels=(goEnrichment$Term))
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
ggplot(goEnrichment[1:9,], aes(x=Term, y=Fish, fill = `Fish`)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p value)") +
  ggtitle("") +  
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Fisher)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    panel.grid = element_blank(),
    legend.position='none',
    ##    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=8, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=0.01, face="bold", hjust=1.10, color = 'black'),
    axis.text.y=element_text(angle=0, size=0.01, face="bold", vjust=0.5, color = 'black'),
    axis.title=element_text(size=0.01, face="bold"),
    ##    legend.key=element_blank(),     #removes the border
    ##    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    #    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
ggsave(paste(figs.out,'biological_process_stationary','.png'), dpi = 300)
ggplot(goEnrichment[1:9,], aes(x=Term, y=Fish, fill = `Fish`)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p value)") +
  ggtitle("") +  
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Fisher)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    panel.grid = element_blank(),
    
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
ggsave(paste(figs.out,'biological_process_stationary_with_text_lognorm_v2','.png'), dpi = 300)
ggsave(paste(figs.out,"biological_process_stationary_with_text_lognorm_v2",'.svg'),  device = "svg")


#GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
#### Cluster 3


library(topGO)
ontology = 'BP'
# install.packages("org.Mm.eg.db") 
# Take from the markers
#cluster2.marker_subset = cluster2.markers[which(cluster2.markers$avg_log2FC > 0),]
#rownames(cluster2.marker_subset) = gsub(pattern = 'EC-',replacement = '',x=rownames(cluster2.marker_subset))
expressed.genes <- rownames(cluster3.markers[intersect(which(cluster2.markers$avg_log2FC > 0), 
                                                       which(cluster2.markers$p_val< 0.05)),])
down.genes <- rownames(cluster3.markers[intersect(which(cluster2.markers$avg_log2FC < 0), 
                                                  which(cluster2.markers$p_val< 0.05)),])
expressed.genes = gsub(pattern = 'EC-',replacement = '',x = expressed.genes)
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
#goEnrichment$Term <- factor(goEnrichment$Term, levels=(goEnrichment$Term))
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
ggplot(goEnrichment, aes(x=Term, y=Fish, fill = `Fish`)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("-log10(p value)") +
  ggtitle("") +  
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Fisher)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=8, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=0, face="bold", hjust=1.10, color = 'black'),
    axis.text.y=element_text(angle=0, size=10, face="bold", vjust=0.5, color = 'black'),
    axis.title=element_text(size=0, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
ggsave(paste(figs.out,'biological_process_cluster3_lognorm','.png'))
#### Biological replicates ####
#### Exponential & stationary ####

aggregate(EC_lpp~Label,data_total,sum)

unique(data_total$Label)
data2 = data_total
#keep_rows  = which(data2$Label =="4HR_COCULTURE" | data2$Label =="STATIONARY" | data2$Label == "Coculture_EXP"| data2$Label == "EXP")
#keep_rows  = which( data2$Label == "Post_3hr_Cipro_Rep1" | data2$Label == "Post_8hr_Cipro_Rep1" |data2$Label == "Post_8hr_Cipro_Rep2"  )
#keep_rows  = which( data2$Label == "Post_3hr_Cipro_Rep1" | data2$Label == "Post_8hr_Cipro_Rep1" |data2$Label == "Post_8hr_Cipro_Rep2"  )
#keep_rows  = which( data2$Label == "Post_3hr_Cipro_Rep1" | data2$Label == "Post_8hr_Cipro_Rep1" |data2$Label == "Post_8hr_Cipro_Rep2"  )
#keep_rows  = which( data2$Label == "Post_8hr_Cipro_Rep1" | data2$Label == "Post_8hr_Cipro_Rep2" )
keep_rows  = which( data2$Label == "Post_8hr_Cipro_Rep1" | data2$Label == "Post_8hr_Cipro_Rep2")
#keep_rows  = which( )
#keep_rows  = grepl('Post',data2$Label)
#keep_rows  = which(data2$Identity == 'MG1655')
data2 = data2[keep_rows,]
#data2 = data2
#### Counts/gene
#### Counts/gene
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
k = colSums(data2[-which(names(data2) %in% c('Label','Identity'))])
to_filter = which(k > 10)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data2 = data2[new_names]
k2 = rowSums(data2[,-which(names(data2) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data2 = data2[to_filter,]
## Make Seurat
aggregate(EC_lpp~Label,data2,sum)


sc_precursor = data2[,-which(names(data2) %in% c('Label','Identity','EC_ECU_04345'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data2[which(names(data2) %in% c('Label'))]
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
#FindMarkers(sct, ident.1 = 'Exponential')

print(sct[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sct, dims = 1:2, reduction = "pca",nfeatures = 20, balanced = TRUE)
FeaturePlot(sct, features=c('EC-rpsC'), order = TRUE)
ElbowPlot(sct,ndims = 15)
mat <- Seurat::GetAssayData(sct, assay = "RNA", slot = "scale.data")
pca <- sct[["pca"]]

total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
sct <- RunUMAP(sct, dims = 1:10, verbose = TRUE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:10, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=3,resolution = 0.3)
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
g = Idents(sct)
g[which(g == 2)]
cluster0.markers <- FindMarkers(sct, ident.1 = 0, min.pct = 0.1, logfc.threshold = 0.05)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1,logfc.threshold = 0.1)
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.1)
#### Print dat umap ####
library(ggsci)
umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out = umap_out[umap_out$label=='STATIONARY',]
umap_out$label = gsub('_',' ',umap_out$label)
unique(umap_out$label)
library(rcartocolor)
#ggsave(paste(figs.out, 'umap_cluster_BS168_combined_total.png'))
aggregate(umi_counts~ident,umap_out,median)
aggregate(umi_counts~ident,umap_out,FUN = length)


library(scales)
colors = hue_pal()(length(unique(umap_out$ident)))
umap_out$ident <- factor(umap_out$ident, levels=as.character(sort(as.numeric(unique(umap_out$ident)))), labels=as.character(sort(as.numeric(unique(umap_out$ident)))))

levs = as.character(0:19)
umap_out$ident = factor(umap_out$ident, levels = levs)

#umap_out$label <- factor(umap_out$label, levels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'), labels=c("Gentamycin", "Tetracycline", "Cefazolin","Cycloserine","Erythromycin","Chloramphenicol",'Ciprofloxacin','Nalidixic acid','Lambda30','Lambda90','Exponential'))
ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `label`)) + 
  geom_point( size = 0.5,alpha = 0.8) +  
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
  scale_color_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(2,3,8,12)])
ggsave(paste(figs.out, 'umap_replicates_scifi8.png'))
ggsave(paste(figs.out,"umap_replicates_scifi8.svg"),  device = "svg")
#### Silhouette score ####


pca_out = data.frame(sct[['pca']]@cell.embeddings)
pca_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
pca_out$label = as.character(sct[['Label']][,1])
pca_out$ident = as.character(sct@active.ident)

library(kBET)
labels = umap_out$label
replicate = as.vector(labels)
replicate[replicate == "Post 8hr Cipro Rep2"] = "0"
replicate[replicate == "Post 8hr Cipro Rep1"] = "1"
random_size = length(which(replicate == 1))
random_labels = rep("0", length(replicate))
random_indices = sample(x = 1:length(replicate), size = random_size, replace = FALSE)
random_labels[random_indices] = "1"
#random_labels[which(idents == "4")] = "3"
library(kBET)
ggplots = list()

cluster = "1"
print(as.numeric(cluster))
reps = 20
batch.silhouette_null = c()
batch.silhouette_replicate = c()
data_size = dim(pca_out)[1]/2
random_size = length(which(replicate == cluster))
random_labels = rep("0", length(replicate))
random_indices = sample(x = 1:length(replicate), size = random_size, replace = FALSE)
random_labels[random_indices] = as.character(cluster)

for (i in 1:reps){
  rand_sample = sample(x = 1:dim(pca_out)[1], size = data_size, replace = FALSE)
  print(rand_sample[1:5])
  idents2 = replicate[rand_sample]
  #  replicate[replicate != "3"] = "0"
  idents2[idents2 != as.character(cluster)] = "0"
  #  random_size = length(which(replicate == 3))
  #    random_size = length(which(replicate == cluster))
  #    random_labels = rep("9", length(replicate))
  #    random_indices = sample(x = 1:length(replicate), size = random_size, replace = FALSE)
  #    #  random_labels[random_indices] = as.characte
  #    random_labels[random_indices] = as.character(cluster)
  random_labels2 = random_labels[rand_sample]
  print(random_indices[1:5])
  pca.data= list()
  pca.data$x = pca_out[rand_sample,]
  random_labels[random_indices] = as.character(cluster)
  batch.silhouette_null = c(batch.silhouette_null,batch_sil(pca.data, random_labels2,nPCs = 10) )
  batch.silhouette_replicate = c(batch.silhouette_replicate,batch_sil(pca.data,idents2,nPCs = 10) )
}


#batch.silhouette_lambda_replicate <- batch_sil(pca.data_lambda, replicate,nPCs = 5)
##batch.silhouette_lambda_null <- batch_sil(pca.data_lambda, random_labels,nPCs = 5)
#batch.silhouette_coli_replicate <- batch_sil(pca.data_coli, replicate,nPCs = 40)
##batch.silhouette_coli_null <- batch_sil(pca.data_coli, random_labels,nPCs = 40)
#batch.silhouette_both_replicate <- batch_sil(pca.data_both, replicate,nPCs = 40)
#batch.silhouette_both_null <- batch_sil(pca.data_both, random_labels,nPCs = 40)


silhouettes = c(batch.silhouette_replicate,batch.silhouette_null
)
cluster_label = c(rep(paste('Replicate1'),length(batch.silhouette_replicate)),rep('Random',length(batch.silhouette_null))
)

to_plot = data.frame(cbind(silhouettes,cluster_label))
to_plot$silhouettes = as.numeric(as.character(silhouettes))

t.test(x = batch.silhouette_replicate, y = batch.silhouette_null, alternative = 'g')

ggplot() + 
  geom_bar(data = to_plot, aes(x = `cluster_label`, y= `silhouettes`, fill = `cluster_label`), position = "dodge", stat="summary", color = 'black')+
  geom_point(data =to_plot, aes(x = `cluster_label`, y= `silhouettes`, fill = `cluster_label`), position =  position_dodge(width = .9))+
  #  geom_errorbar(data = g, aes(x = `Fixation`, y= `CT`, fill = `Sample`,ymin = 0, ymax = 30), position =  position_dodge(width = .9))+
  #  scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values=as.vector(rcartocolor::carto_pal(name='Safe'))[c(40,2,8,12)])+

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
        panel.grid = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=12,face="bold")) + scale_colour_hue(l=40)
#p = 0.8877
ggsave(paste(figs.out, 'silhouette_score_duplicates_ecoli', '.png', sep = ''))
ggsave(paste(figs.out,"silhouette_score_duplicates_ecoli",'.svg'),  device = "svg")

m = bind_rows(ggplots, .id = "column_label")

#### Variablee Genes


umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)
#umap_out$ident = as.character(sct@assays$RNA[[1]])
#norm_lpp_index = which('EC-lpp'== rownames(sct@assays$RNA))
#norm_ompC_index = which('EC-ompC'== rownames(sct@assays$RNA))
#norm_recA_index = which('EC-recA'== rownames(sct@assays$RNA))
#umap_out$EC_lpp_total = as.numeric(as.character(data2$EC_lpp))
#umap_out$EC_lpp_norm = as.numeric(as.character(as.vector(sct@assays$RNA[norm_lpp_index,])))
#umap_out$EC_ompC_total = as.numeric(as.character(data2$EC_ompC))
#umap_out$EC_ompC_norm = as.numeric(as.character(as.vector(sct@assays$RNA[norm_ompC_index,])))
#umap_out$EC_recA_total = as.numeric(as.character(data2$EC_recA))
#umap_out$EC_recA_norm = as.numeric(as.character(as.vector(sct@assays$RNA[norm_recA_index,])))


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
  ggtitle(sprintf("UMAP for 3hr Cipro - MG1655")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
p3 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMAP for 3hr Cipro - MG1655")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
p2
#ggsave(paste(figs.out, 'umap_umi_counts_MG1655_cipro_3hr.png'))
#p3
#ggsave(paste(figs.out, 'umap_cluster_MG1655_cipro_8hr.png'))

#aggregate(EC_lpp_total~ident,umap_out,sum)
#aggregate(EC_lpp_total~ident,umap_out,mean)
#aggregate(EC_ompC_total~ident,umap_out,sum)
#aggregate(EC_ompC_total~ident,umap_out,mean)
#aggregate(EC_recA_total~ident,umap_out,sum)
#aggregate(EC_recA_total~ident,umap_out,mean)
#aggregate(EC_ompC_total~ident,umap_out,sum)
#aggregate(EC_ompC_total~ident,umap_out,mean)
#aggregate(EC_lpp_norm~ident,umap_out,mean)
#aggregate(EC_lpp_norm~ident,umap_out,length)
#aggregate(umi_counts~ident,umap_out,mean)
#aggregate(umi_counts~ident,umap_out,length)

#### gadABC ####
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



umap_out = data.frame(sct[['umap']]@cell.embeddings)
umap_out$umi_counts = as.numeric(as.character(sct[['nCount_RNA']][,1]))
umap_out$label = as.character(sct[['Label']][,1])
umap_out$ident = as.character(sct@active.ident)

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
  theme(plot.title = element_text(hjust = 0.5, size = 20, face="bold.italic"),panel.grid=element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() ,
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 20), legend.title = element_text(size=20,face="bold")) + guides(colour = guide_legend(override.aes = list(size=8)))+
  labs(fill = "", colour = "") 

ggsave(paste(figs.out,"gadA_umap_stationar_only_with_words",'.png'))


ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  
  #  labs(x = "UMAP1", y = "UMAP2", color = "Condition") +
  labs(x = "", y = "", color = "Cluster") +
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
ggsave(paste(figs.out, 'ec_stationary_cluster_10pcs.png'))

#temp_phage = matrix(gad_data[c('norm_pbs','norm_sp')],ncol=2)
#colnames(temp_phage) = paste0('PHAGE_',1:2)
#sct[["phage"]] <- CreateDimReducObject(embeddings = temp_phage, key = "PHAGE_", assay = DefaultAssay(sct))

#look_into = intersect(which(gad_data$norm_sp2 > 0.1),which(gad_data$norm_pbs2 > 0.02))
#look_into = union(which(gad_data$norm_sp2 > 0.05),which(gad_data$norm_pbs2 > 0.05))
#look_into1 = which(gad_data$norm_sp2 > 0.5)
#look_into2 = intersect(which(gad_data$norm_sp2 > 0.1),which(gad_data$norm_sp2 < 0.5))
#look_into3 = intersect(intersect(which(gad_data$norm_sp2 > 0.1),which(gad_data$norm_pbs2 > 0.02)),intersect(which(gad_data$norm_sp2 > 0.05),which(gad_data$norm_sp2 < 0.5)))
#look_into4 =  which(gad_data$norm_pbs2 > 0.1)
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
gad_data5 = cbind(gad_data, norm_data)
#gad_data$x1_stage = gad_data2$X1
#gad_data$x2_stage = gad_data2$X2



#gad_data = data.frame(norm_data %*% phage_score)
#gad_data2 = data.frame(norm_data %*% phage_score2)

#### Keep on exponential/stationary ###



unique(data_total$Label)
data2 = data_total
#keep_rows  = which(data2$Label =="Post_Exponential" | data2$Label =="Post_Stationary" )
#keep_rows  = which( data2$Label == "Post_3hr_Cipro_Rep1" )
#keep_rows  = which(data2$Identity == 'MG1655')
#data2 = data2[keep_rows,]
#### Counts/gene
#### Counts/gene
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
k = colSums(data2[-which(names(data2) %in% c('Label','Identity'))])
to_filter = which(k > 10)
new_names = names(to_filter)
new_names = c(new_names, "Label",'Identity')
new_names = c(new_names)
data2 = data2[new_names]
k2 = rowSums(data2[,-which(names(data2) %in% c('Label','Identity','anti-rRNA'))])
to_filter = which(k2 >= 10)
data2 = data2[to_filter,]
## Make Seurat

sc_precursor = data2[,-which(names(data2) %in% c('Label','Identity','EC_ECU_04345'))]
rownames(sc_precursor) = 1:dim(sc_precursor)[1]
sc_precursor = t(sc_precursor)
sc_meta = data2[which(names(data2) %in% c('Label'))]
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
#ggsave(paste(figs.out,'dim_loadings_MG1655','.png'))
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
#ggsave(paste(figs.out,'pca_MG1655_var_explained','.png'))

#sct = JackStraw(sct)

#data_to_write_out <- t(as.data.frame(as.matrix(sct@assays$SCT@scale.data)))
#fwrite(x = data_to_write_out, row.names = FALSE, file = "~/Dropbox (Princeton)/Notes/research/gitai_wingreen/scale_transformed.csv")

names(data)[str_detect(names(data),'spo')]

DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca') + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'pca_MG1655_scifi7','.png'))
DimPlot(sct, label = TRUE,group.by = 'Label',reduction = 'pca', dims = c(2,3)) + ggtitle('PCA - Normalized Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave('~/Dropbox (Princeton)/Notes/research/gitai_wingreen/slides/Figures/2020-05-12/PCA_Relative_Count.png')

sct <- RunUMAP(sct, dims = 1:10, verbose = TRUE)
#sct <- RunTSNE(sct, dims = 1:10, verbose = FALSE)

sct <- FindNeighbors(sct, dims = 1:10, verbose = TRUE)
sct <- FindClusters(sct, verbose = TRUE, algorithm=1,resolution = 0.35)
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
cluster1.markers <- FindMarkers(sct, ident.1 = 1, min.pct = 0.1)
cluster2.markers <- FindMarkers(sct, ident.1 = 2, min.pct = 0.1)
cluster3.markers <- FindMarkers(sct, ident.1 = 3, min.pct = 0.1,ident.2=0)
#cluster3.markers <- FindMarkers(sct, ident.1 = 3,ident.2=4, min.pct = 0.1)
cluster4.markers <- FindMarkers(sct, ident.1 = 4, min.pct = 0.1,ident.2 = 0)
cluster5.markers <- FindMarkers(sct, ident.1 = 5, min.pct = 0.1)
#cluster5.markers <- FindMarkers(sct, ident.1 = 5, ident.2 = 8,min.pct = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6, min.pct = 0.1)
cluster6.markers <- FindMarkers(sct, ident.1 = 6,ident.2=3, min.pct = 0.1)
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
#FeaturePlot(sct,features = c('EC-ykuO','EC-sigH','EC-suhB','EC-hemE'), reduction = 'umap')
#FeaturePlot(sct,features = c('EC-tufA','EC-putA','EC-sdiA','EC-bcsA'), reduction = 'umap')
#FeaturePlot(sct,features = c('EC-rplP','EC-dmlA','EC-sdiA','EC-tap'), reduction = 'umap',size  = 4)
#ggsave(paste(figs.out,'MG1655_geneplot_scifi3','.png'))

DimPlot(sct, label = TRUE,group.by = 'Label') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE ,  ) + ggtitle('Clusters (UMAP)- Relative Counts') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,reduction = 'pca') + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
DimPlot(sct, label = TRUE,reduction = 'pca',dims = c(1,2)) + ggtitle('UMAP') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.3),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#ggsave(paste(figs.out,'umap_clustered_MG1655_scifi3','.png'))

#FeaturePlot(sct,features =  c('acoA'), reduction = 'umap') & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('EC-srfAB','EC-sigF','EC-cwlC','EC-ECU-04345'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('EC-rpsQ','EC-trnS-Leu1','EC-floT','EC-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('EC-gapB','EC-mntA','EC-floT','EC-pftA'), reduction = 'pca')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('EC-opuCB','EC-ydjP','EC-sigB','EC-cwlC'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('EC-rpsQ','EC-manA','EC-ykzE','EC-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#FeaturePlot(sct,features =  c('EC-spxA','EC-yhgE','EC-nfeDB','EC-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
FeaturePlot(sct,features =  c('EC-lpp','EC-recA','EC-ompC','EC-yhiM'), reduction = 'umap',order=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'MG1655_geneplot_cipro_genemap','.png'))
FeaturePlot(sct,features =  c('EC-gltW','EC-tufA','EC-gadA','EC-gadB'), reduction = 'umap',order=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave(paste(figs.out,'MG1655_geneplot_exp_stationary_post','.png'))
#FeaturePlot(sct,features =  c('EC-spxA','EC-yhgE','EC-nfeDB','EC-pftA'), reduction = 'umap')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
#ggsave(paste(figs.out,'MG1655_geneplot_exp_heat_cef','.png'))
levels(sct)  = c(0,1,2,3)
markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
library(dplyr)
markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top10
new_genes = top10$gene
library(R.utils)
new_genes = insert(top10$gene,ats = 10,values = c('EC-recA'))

DoHeatmap(sct, new_genes,draw.lines = TRUE)  + scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ 
  theme(text = element_text(size = 20))
ggsave(paste(figs.out,'MG1655_heatmap_cipro_3hr','.png'))
means = aggregate(umi_counts~ident,umap_out,mean)
sizes = aggregate(umi_counts~ident,umap_out,length)
sizes$umi_frac = sizes$umi_counts/sum(sizes$umi_counts)
sizes$mean = means$umi_counts
sizes
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
#       zlab = 'PC3',bty = "b2", main ="MG1655 PCA",size = 3,alpha = 0.8,type = 'p', specular="black" )
#legend3d("topright", legend = paste(col_key$label), pch = 16, col = col_key$color, cex=2, inset=c(0.08))


p1 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`)) + 
  geom_point(pch = 21, aes(fill = umi_counts),size = 1.0) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMI counts - E. MG1655")) +
  scale_fill_gradientn(colours =rev(brewer.pal(n = 11, name = "RdYlBu")) ) + 
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+ guides(colour = guide_legend(override.aes = list(size=10)))
p2 = ggplot(pca_out,aes(x = `PC_1`, y = `PC_2`,col = `ident`)) + 
  geom_point( size = 1.0) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("PCA - 3hr Cipro")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
p2
ggsave(paste(figs.out, 'pca_umi_counts_MG1655_3hr.png'))

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
  ggtitle(sprintf("UMAP for Exp/Stationary  - MG1655")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
p3 = ggplot(umap_out,aes(x = `UMAP_1`, y = `UMAP_2`,col = `ident`)) + 
  geom_point( size = 1,alpha = 1.5) +  labs(x = "Dim 1", y = "Dim 2", color = "Class") +
  ggtitle(sprintf("UMAP for Exp/Stationary - MG1655")) +
  #  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)+scale_color_brewer(palette="Set2")+ guides(colour = guide_legend(override.aes = list(size=10)))
p2
ggsave(paste(figs.out, 'umap_umi_counts_MG1655_exp_stationary_3hr.png'))
p3
ggsave(paste(figs.out, 'umap_cluster_MG1655_exp_stationary.png'))
means = aggregate(umi_counts~ident,umap_out,mean)
sizes = aggregate(umi_counts~ident,umap_out,length)
sizes$umi_frac = sizes$umi_counts/sum(sizes$umi_counts)
sizes$mean = means$umi_counts
sizes
