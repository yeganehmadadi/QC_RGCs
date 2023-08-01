
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(DoubletFinder)
#---------------------------------------------- Loading single-cell RNA-seq count data
root.path= ""
file.path= paste0(root.path, "dataset/")
output.path= paste0(root.path, "results/")

dataFileName1 = paste(file.path, "C1.csv", sep="") 
dataFileName2 = paste(file.path, "D1.csv", sep="")  
dataFileName3 = paste(file.path, "E1.csv", sep="") 
dataFileName4 = paste(file.path, "F1.csv", sep="") 
dataFileName5 = paste(file.path, "H1.csv", sep="") 
dataFileName6 = paste(file.path, "I1.csv", sep="")   
dataFileName7 = paste(file.path, "M1.csv", sep="") 
dataFileName8 = paste(file.path, "N1.csv", sep="")

# read data
singleCellDataGenes_C <- read.csv(file= dataFileName1, head=TRUE, row.names = "Ensembl.ID", sep= ",", stringsAsFactors= FALSE)
singleCellDataGenes_D <- read.csv(file= dataFileName2, head=TRUE, row.names = "Ensembl.ID", sep= ",", stringsAsFactors= FALSE) 
singleCellDataGenes_E <- read.csv(file= dataFileName3, head=TRUE, row.names = "Ensembl.ID", sep= ",", stringsAsFactors= FALSE)  
singleCellDataGenes_F <- read.csv(file= dataFileName4, head=TRUE, row.names = "Ensembl.ID", sep= ",", stringsAsFactors= FALSE)  
singleCellDataGenes_H <- read.csv(file= dataFileName5, head=TRUE, row.names = "Ensembl.ID", sep= ",", stringsAsFactors= FALSE)  
singleCellDataGenes_I <- read.csv(file= dataFileName6, head=TRUE, row.names = "Ensembl.ID", sep= ",", stringsAsFactors= FALSE)  
singleCellDataGenes_M <- read.csv(file= dataFileName7, head=TRUE, row.names = "Ensembl.ID", sep= ",", stringsAsFactors= FALSE) 
singleCellDataGenes_N <- read.csv(file= dataFileName8, head=TRUE, row.names = "Ensembl.ID", sep= ",", stringsAsFactors= FALSE) 

#---------------------------------------------- Removing duplicate genes
selectHighestExpressedDuplicatedGenesFun <- function(singleCellDataGenes){
  # ----Smart and quick way to find genes that are more expressed from duplicated genes ---- #
  # creating a temporary label that is a combination of gene and the amount of expression
  # *Note: sep="-" is very important to avoid gene name confusion
  singleCellDataGenes$tmp.label = paste(singleCellDataGenes$Gene.Symbol, rowSums(singleCellDataGenes[, c(2:ncol(singleCellDataGenes))] > 0 )+1000, sep="-")
  # order data based on this new label
  single.tmp = singleCellDataGenes[order(singleCellDataGenes$tmp.label), ]
  # pick the last duplicated label (belongs to highest expressed gene)
  single.tmp.new = single.tmp[!duplicated(single.tmp$Gene.Symbol, fromLast= TRUE),]

  #----------- Test Routine: to make sure cell with greater number of genes selected
  t1= single.tmp[!duplicated(single.tmp$Gene.Symbol, fromLast= TRUE),]
  t2= single.tmp[!duplicated(single.tmp$Gene.Symbol, fromLast= FALSE),]
  print("Two test routines are perfomring:")
  r1= sum(rowSums( t1[, c(2:ncol(singleCellDataGenes))] > 0 ) < rowSums( t2[, c(2:ncol(singleCellDataGenes))] > 0 ) )  # sum should be zero
  r2= sum(duplicated(t1$Gene.Symbol))                                       # final check after removing duplicates
  print(paste("Number of genes with wrong lower expresseion is: ", toString(r1), sep= ""))
  print(paste("Number of duplicated genes is: ", toString(r2), sep=""))
  # ----------
  return ( single.tmp.new[,1:(ncol(singleCellDataGenes)-1)] )
}

singleCellDataGenes_C = selectHighestExpressedDuplicatedGenesFun(singleCellDataGenes_C)
rownames(singleCellDataGenes_C) = singleCellDataGenes_C$Gene.Symbol
singleCellDataGenes_D = selectHighestExpressedDuplicatedGenesFun(singleCellDataGenes_D)
rownames(singleCellDataGenes_D)= singleCellDataGenes_D$Gene.Symbol
singleCellDataGenes_E = selectHighestExpressedDuplicatedGenesFun(singleCellDataGenes_E)
rownames(singleCellDataGenes_E) = singleCellDataGenes_E$Gene.Symbol
singleCellDataGenes_F = selectHighestExpressedDuplicatedGenesFun(singleCellDataGenes_F)
rownames(singleCellDataGenes_F)= singleCellDataGenes_F$Gene.Symbol
singleCellDataGenes_H= selectHighestExpressedDuplicatedGenesFun(singleCellDataGenes_H)
rownames(singleCellDataGenes_H)= singleCellDataGenes_H$Gene.Symbol
singleCellDataGenes_I = selectHighestExpressedDuplicatedGenesFun(singleCellDataGenes_I)
rownames(singleCellDataGenes_I)= singleCellDataGenes_I$Gene.Symbol
singleCellDataGenes_M = selectHighestExpressedDuplicatedGenesFun(singleCellDataGenes_M)
rownames(singleCellDataGenes_M)= singleCellDataGenes_M$Gene.Symbol
singleCellDataGenes_N = selectHighestExpressedDuplicatedGenesFun(singleCellDataGenes_N)
rownames(singleCellDataGenes_N)= singleCellDataGenes_N$Gene.Symbol

Plate_C = singleCellDataGenes_C[,-1] #remove column Gene.Symbol
colnames(Plate_C) = paste("C", colnames(Plate_C), sep="_") #add prefix to column names
Plate_D = singleCellDataGenes_D[,-1] #remove column Gene.Symbol
colnames(Plate_D) = paste("D", colnames(Plate_D), sep="_") #add prefix to column names
Plate_E = singleCellDataGenes_E[,-1] #remove column Gene.Symbol
colnames(Plate_E) = paste("E", colnames(Plate_E), sep="_") #add prefix to column names
Plate_F = singleCellDataGenes_F[,-1] #remove column Gene.Symbol
colnames(Plate_F) = paste("F", colnames(Plate_F), sep="_") #add prefix to column names
Plate_H = singleCellDataGenes_H[,-1] #remove column Gene.Symbol
colnames(Plate_H) = paste("H", colnames(Plate_H), sep="_") #add prefix to column names
Plate_I = singleCellDataGenes_I[,-1] #remove column Gene.Symbol
colnames(Plate_I) = paste("I", colnames(Plate_I), sep="_") #add prefix to column names
Plate_M = singleCellDataGenes_M[,-1] #remove column Gene.Symbol
colnames(Plate_M) = paste("M", colnames(Plate_M), sep="_") #add prefix to column names
Plate_N = singleCellDataGenes_N[,-1] #remove column Gene.Symbol
colnames(Plate_N) = paste("N", colnames(Plate_N), sep="_") #add prefix to column names

#---------------------------------------------- Quality control
#---------------------------------------------- Create Seurat objects
# We can now load the expression matrices into objects and then merge them into a single merged object. 
sdata.C <- CreateSeuratObject(counts= Plate_C, project= "C", min.cells= 5, min.features= 200)
dim(sdata.C)
sdata.D <- CreateSeuratObject(counts= Plate_D, project= "D", min.cells= 5, min.features= 200)
dim(sdata.D)
sdata.E <- CreateSeuratObject(counts= Plate_E, project= "E", min.cells= 5, min.features= 200)
dim(sdata.E)
sdata.F <- CreateSeuratObject(counts= Plate_F, project= "F", min.cells= 5, min.features= 200)
dim(sdata.F)
sdata.H <- CreateSeuratObject(counts= Plate_H, project= "H", min.cells= 5, min.features= 200)
dim(sdata.H)
sdata.I <- CreateSeuratObject(counts= Plate_I, project= "I", min.cells= 5, min.features= 200)
dim(sdata.I)
sdata.M <- CreateSeuratObject(counts= Plate_M, project= "M", min.cells= 5, min.features= 200)
dim(sdata.M)
sdata.N <- CreateSeuratObject(counts= Plate_N, project= "N", min.cells= 5, min.features= 200)
dim(sdata.N)

# Check the metadata in the Seurat objects
head(sdata.C@meta.data)
head(sdata.D@meta.data)

#---------------------------------------------- Create a merged Seurat object
merged_seurat <- merge(sdata.C, c(sdata.D, sdata.E, sdata.F, sdata.H, sdata.I, sdata.M, sdata.N), 
                       add.cell.ids = c("C", "D", "E", "F", "H", "I", "M", "N"))
dim(merged_seurat)
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

# Here it is how the count matrix and the metadata look like for every cell
as.data.frame(merged_seurat@assays$RNA@counts[1:10, 1:2])
head(merged_seurat@meta.data, 10)
write.csv(merged_seurat@assays$RNA@counts, file= "merged_seurat.csv")

#---------------------------------------------- Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

#---------------------------------------------- Compute percent mito ratio
# add information directly to the metadata slot in the Seurat object using the $ operator
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
head(merged_seurat$mitoRatio, 10)

#---------------------------------------------- calculate the proportion gene expression that comes from ribosomal proteins
merged_seurat$riboRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^Rp[sl]")
merged_seurat$riboRatio <- merged_seurat@meta.data$riboRatio / 100
head(merged_seurat$riboRatio, 10)

#---------------------------------------------- Create metadata dataframe
metadata <- merged_seurat@meta.data

#---------------------------------------------- Add cell IDs to metadata
metadata$cells <- rownames(metadata)

#---------------------------------------------- Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

#---------------------------------------------- Create plate column
metadata$plate <- NA
metadata$plate[which(str_detect(metadata$cells, "^C_"))] <- "C"
metadata$plate[which(str_detect(metadata$cells, "^D_"))] <- "D"
metadata$plate[which(str_detect(metadata$cells, "^E_"))] <- "E"
metadata$plate[which(str_detect(metadata$cells, "^F_"))] <- "F"
metadata$plate[which(str_detect(metadata$cells, "^H_"))] <- "H"
metadata$plate[which(str_detect(metadata$cells, "^I_"))] <- "I"
metadata$plate[which(str_detect(metadata$cells, "^M_"))] <- "M"
metadata$plate[which(str_detect(metadata$cells, "^N_"))] <- "N"

#---------------------------------------------- updating metadata to our Seurat object
# assigning the dataframe into the meta.data slot
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

#---------------------------------------------- Assessing the quality metrics
# We will assess various metrics and then decide on which cells are low quality andshould be removed from the analysis
# Cell counts; UMI counts per cell; Genes detected per cell; UMIs vs. genes detected; Mitochondrial counts ratio; Novelty

# Visualize the number of cell counts per plate
metadata %>% 
  ggplot(aes(x=plate, fill=plate)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "white") +
  ggtitle("NCells")

#---------------
# Visualize the distribution of number UMIs (transcripts) per cell via histogram
metadata %>% 
  ggplot(aes(color=plate, x=nUMI, fill= plate)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 250000)

#---------------
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=plate, x=nGene, fill= plate)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 900)

# For high quality data, the proportional histogram should contain a single large peak that represents cells that were encapsulated. 
# If we see a small shoulder to the right of the major peak or a bimodal distribution of the cells, that can indicate a couple of things: 
# It might be that there are a set of cells that failed for some reason. It could also be that there are biologically different types of cells.

# Visualize the distribution of genes detected per cell via box plot
metadata %>% 
  ggplot(aes(x=plate, y=log10(nGene), fill=plate)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

#---------------
# Visualize the distribution of mitochondrial gene expression detected per cell
# This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells
metadata %>% 
  ggplot(aes(color=plate, x=mitoRatio, fill=plate)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.6)

dim(merged_seurat@assays$RNA@counts)
#----------------------------------------------- QC metrics stored in Seurat
# Visualize QC metrics as a violin plot
VlnPlot(merged_seurat, features= c("nGene", "nUMI", "mitoRatio", "riboRatio"), ncol = 2)

#----------------------------------------------- Filtering
#To filter, we will go back to our Seurat object and use the subset() function:

#--------------- Cell-level filtering
# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_merged_seurat <- subset(x = merged_seurat, 
                                 subset= (nUMI >= 250000) & 
                                   (nGene >= 900) & 
                                   (mitoRatio < 0.6)&
                                   (riboRatio >= 0.02))
dim(filtered_merged_seurat@assays$RNA@counts)
write.csv(filtered_merged_seurat@assays$RNA@counts, file= "Cell-level_filtered_merged_seurat.csv")

#--------------- Gene-level filtering
# First we will remove genes that have zero expression in all cells.
# Second we will keep only genes which are expressed in 5 or more cells.

# Extract counts
counts <- GetAssayData(object = filtered_merged_seurat, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 3 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 5
# Only keeping those genes expressed in more than 3 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_merged_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_merged_seurat@meta.data)
dim(filtered_merged_seurat@assays$RNA@counts)
write.csv(filtered_merged_seurat@assays$RNA@counts, file= "Gene-level_filtered_merged_seurat.csv")

#----------------------------------------------- genes contribute the most to such reads
# We can also see which genes contribute the most to such reads. We can for instance plot the percentage of counts per gene.
# # Compute the relative expression of each gene per cell Use sparse matrix.
par(mar = c(4, 8, 2, 1))
C <- filtered_merged_seurat@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

#---------------------------------------------- again plot, after filtering

# Add number of genes per UMI for each cell to metadata
filtered_merged_seurat$log10GenesPerUMI <- log10(filtered_merged_seurat$nFeature_RNA) / log10(filtered_merged_seurat$nCount_RNA)

# Compute percent mito ratio
filtered_merged_seurat$mitoRatio <- PercentageFeatureSet(object = filtered_merged_seurat, pattern = "^mt-")
filtered_merged_seurat$mitoRatio <- filtered_merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
clean_metadata <- filtered_merged_seurat@meta.data

# Add cell IDs to metadata
clean_metadata$cells <- rownames(clean_metadata)

# Rename columns
clean_metadata <- clean_metadata %>%
  dplyr::rename(seq_folder1 = orig.ident,
                nCount = nCount_RNA,
                nFeature = nFeature_RNA)

# Create plate column
clean_metadata$plate <- NA
clean_metadata$plate[which(str_detect(clean_metadata$cells, "^C_"))] <- "C"
clean_metadata$plate[which(str_detect(clean_metadata$cells, "^D_"))] <- "D"
clean_metadata$plate[which(str_detect(clean_metadata$cells, "^E_"))] <- "E"
clean_metadata$plate[which(str_detect(clean_metadata$cells, "^F_"))] <- "F"
clean_metadata$plate[which(str_detect(clean_metadata$cells, "^H_"))] <- "H"
clean_metadata$plate[which(str_detect(clean_metadata$cells, "^I_"))] <- "I"
clean_metadata$plate[which(str_detect(clean_metadata$cells, "^M_"))] <- "M"
clean_metadata$plate[which(str_detect(clean_metadata$cells, "^N_"))] <- "N"

# Add metadata back to Seurat object
filtered_merged_seurat@meta.data <- clean_metadata

# Visualize the number of cell counts per plate
clean_metadata %>% 
  ggplot(aes(x=plate, fill=plate)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "white") +
  ggtitle("NCells")

# Visualize the number UMIs (transcripts) per cell
clean_metadata %>% 
  ggplot(aes(color=plate, x=nCount, fill= plate)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 

# Visualize the distribution of genes detected per cell via histogram
clean_metadata %>% 
  ggplot(aes(color=plate, x=nFeature, fill= plate)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()  

# Visualize the distribution of genes detected per cell via boxplot
clean_metadata %>% 
  ggplot(aes(x=plate, y=log10(nFeature), fill=plate)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the distribution of mitochondrial gene expression detected per cell
clean_metadata %>% 
  ggplot(aes(color=plate, x=mitoRatio, fill= plate)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic()

#---------------------------------------------- Count normalization
seurat_phase <- NormalizeData(filtered_merged_seurat)
dim(seurat_phase)

#---------------------------------------------- identification of the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)
dim(seurat_phase)

#---------------------------------------------- Scale the counts
seurat_phase <- ScaleData(seurat_phase)
dim(seurat_phase)

#---------------------------------------------- Principal Component Analysis (PCA)
# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by plate
DimPlot(seurat_phase,
        reduction = "pca")

# Plot the PCA split by plate
PCAPlot(seurat_phase,
        split.by = "plate") 

#---------------------------------------------- UMAP visualization
# Run UMAP
seurat_phase <- RunUMAP(seurat_phase, 
                        dims = 1:40,
                        reduction = "pca")

# Plot UMAP colored by plate                            
DimPlot(seurat_phase) 

# Plot UMAP split by plates
DimPlot(seurat_phase,
        split.by = "plate")  

#---------------------------------------------- Identify significant PCs 
# Explore heatmap of PCs
DimHeatmap(seurat_phase, 
           dims = 1:18, 
           cells = 500, 
           balanced = TRUE)

# print out the top positive and negative genes by PCA scores driving the PCs.
# Printing out the most variable genes driving PCs
print(x = seurat_phase[["pca"]], 
      dims = 1:18, 
      nfeatures = 5)

# Plot the elbow plot
# The elbow plot is another helpful way to determine how many PCs to use for clustering so that we are capturing majority of the variation in the data.
ElbowPlot(object = seurat_phase, 
          ndims = 40)

#--------------- Doublet Finder
# Extremely high number of detected genes could indicate doublets. 
#However, depending on the cell type composition in your sample, you may have cells with higher number of genes (and also higher counts) from one cell type.
# we run doubletFinder, selecting first 18 PCs and a pK value of 0.9. To optimize the parameters, you can run the paramSweep function in the package.

# define the expected number of doublet cellscells.
nExp <- round(ncol(seurat_phase) * 0.04)  # expect 4% doublets
seurat_phase <- doubletFinder_v3(seurat_phase, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:18)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(seurat_phase@meta.data)[grepl("DF.classification", colnames(seurat_phase@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(seurat_phase, group.by = "plate") + NoAxes(),
                   DimPlot(seurat_phase, group.by = DF.name) + NoAxes())

# We should expect that two cells have more detected genes than a single cell, lets check if our predicted doublets also have more detected genes in general.
VlnPlot(seurat_phase, features = "nGene", group.by = DF.name, pt.size = 0.1)

# lets remove all predicted doublets from our data.
seurat_phase = seurat_phase[, seurat_phase@meta.data[, DF.name] == "Singlet"]
dim(seurat_phase)
write.csv(seurat_phase@assays$RNA@counts, file= "seurat_phase_after_doublet_finder.csv")

#--------------------------
# Create metadata dataframe
cleanFinal_metadata <- seurat_phase@meta.data

# Add cell IDs to metadata
cleanFinal_metadata$cells <- rownames(cleanFinal_metadata)

# Rename columns
cleanFinal_metadata <- cleanFinal_metadata %>%
  dplyr::rename(seq_folder1 = orig.ident,
                nCount = nCount_RNA,
                nFeature = nFeature_RNA)

# Create plate column
cleanFinal_metadata$plate <- NA
cleanFinal_metadata$plate[which(str_detect(cleanFinal_metadata$cells, "^C_"))] <- "C"
cleanFinal_metadata$plate[which(str_detect(cleanFinal_metadata$cells, "^D_"))] <- "D"
cleanFinal_metadata$plate[which(str_detect(cleanFinal_metadata$cells, "^E_"))] <- "E"
cleanFinal_metadata$plate[which(str_detect(cleanFinal_metadata$cells, "^F_"))] <- "F"
cleanFinal_metadata$plate[which(str_detect(cleanFinal_metadata$cells, "^H_"))] <- "H"
cleanFinal_metadata$plate[which(str_detect(cleanFinal_metadata$cells, "^I_"))] <- "I"
cleanFinal_metadata$plate[which(str_detect(cleanFinal_metadata$cells, "^M_"))] <- "M"
cleanFinal_metadata$plate[which(str_detect(cleanFinal_metadata$cells, "^N_"))] <- "N"

# Add metadata back to Seurat object
seurat_phase@meta.data <- cleanFinal_metadata

# Visualize the number of cell counts per plate
cleanFinal_metadata %>% 
  ggplot(aes(x=plate, fill=plate)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "white") +
  ggtitle("NCells")

#----------------------------------------------- genes contribute the most to such reads
# We can also see which genes contribute the most to such reads. We can for instance plot the percentage of counts per gene.
# # Compute the relative expression of each gene per cell Use sparse matrix.
par(mar = c(4, 8, 2, 1))
C <- filtered_merged_seurat@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

#---------------------------------------------- Clustering the cells
# To perform clustering, we determine the genes that are most different in their expression between cells. 
# Then, we use these genes to determine which correlated genes sets are responsible for the largest differences in expression between cells.
# Seurat uses a graph-based clustering approach, which embeds cells in a graph structure, using a K-nearest neighbor (KNN) graph (by default), with edges drawn between cells with similar gene expression patterns. Then, it attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.

# Determine the K-nearest neighbor graph
seurat_phase <- FindNeighbors(object = seurat_phase, dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_phase <- FindClusters(object = seurat_phase, resolution = c(0.4, 0.6, 0.8, 1.0, 1.4, 1.6, 1.8, 2.0))

# Explore resolutions
seurat_phase@meta.data %>% 
  View()

# Assign identity of clusters
Idents(object = seurat_phase) <- "RNA_snn_res.1.8"

# Plot the UMAP
DimPlot(seurat_phase,
        reduction = "umap",
        label = TRUE,
        label.size = 3)

#---------------------------------------------- Segregation of clusters by plate
# Extract identity and plate information from Seurat object to determine the number of cells per cluster per plate
n_cells <- FetchData(seurat_phase, 
                     vars = c("ident", "plate")) %>%
  dplyr::count(ident, plate) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

# UMAP of cells in each cluster by plates
DimPlot(seurat_phase, 
        label = TRUE, 
        split.by = "plate", ncol=3, label.size = 3)

# Determine metrics to plot present in seurat_phase@meta.data
metrics <-  c("nUMI", "nGene", "mitoRatio", "riboRatio")
FeaturePlot(seurat_phase, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.3, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

#---------------------------------------------- Exploration of the PCs driving the different clusters
# Defining the information in the Seurat object of interest
columns <- c(paste0("PC_", 1:20),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the Seurat object
pc_data <- FetchData(seurat_phase, 
                     vars = columns)

# Extract the UMAP coordinates for the first 10 cells
seurat_phase@reductions$umap@cell.embeddings[1:10, 1:2]

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_phase, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
# we can see that the genes driving in each PC exhibit higher expression in which cluster
map(paste0("PC_", 1:20), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# For instance, the genes driving PC_2 exhibit higher expression in clusters 6, 11, and 17 (maybe a bit higher in 15, too). 
# We could look back at our genes driving this PC to get an idea of what the cell types might be.
# Examine PCA results 
print(seurat_phase[["pca"]], dims = 1:20, nfeatures = 5)

#---------------------------------------------- Exploring known cell type markers (base on Tran et. al.)
DimPlot(object = seurat_phase, 
        reduction = "umap", 
        label = TRUE)

# To access the expression levels of all genes, rather than just the 3000 most highly variable genes, 
# we can use the normalized count data stored in the RNA assay slot.

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_phase) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_phase <- NormalizeData(seurat_phase, verbose = FALSE)

FeaturePlot(seurat_phase, 
            reduction = "umap", 
            features = c("Rbpms", "Thy1", "Slc17a6", "Pou4f1", "Pou4f2", "Pou4f3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# Vln plot 
VlnPlot(object = seurat_phase,
        features = c("Rbpms", "Thy1", "Slc17a6", "Pou4f1", "Pou4f2", "Pou4f3"))

# Dot Plot
DotPlot(object = seurat_phase, cols = c("lightgrey", "red"),
        features = c("Rbpms", "Thy1", "Slc17a6", "Pou4f1", "Pou4f2", "Pou4f3", "Tfap2a", "Gad1", "Slc6a9", "Lhx1", "Onecut1", "Vsx2", "Otx2", "Arr3", "Rho", "Rlbp1", "Aqp4", "Pecam1", "Kcnj8", "Fcrls", "P2ry12"))

# --------------------------------- Adding Gene Annotations

annotations <- read.csv("MetaData/00_Master_RGC_Y.csv")

# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat_phase.markers <- FindAllMarkers(seurat_phase)
seurat_phase.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(seurat_phase.markers, "All markers for each clusters.csv", row.names=FALSE)

# --------------------------------- Evaluating marker genes
# Extract top 30 markers per cluster
seurat_phase.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top30

# Visualize top 30 markers per cluster
View(top30)
write.csv(top30, "top30 markers for each clusters.csv", row.names=FALSE)

# --------------------------------- Visualizing marker genes
# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 30 markers for each cluster.
DoHeatmap(object = seurat_phase, features= top30$gene)

DotPlot(object   = seurat_phase,
        features = as.character(unique(top30$gene)),
        #group.by = "cluster",
        assay    = "RNA") +
  coord_flip()

# --------------------------------- Assigning cell type identity to clusters based cluster markers on Tran et al. (Fig. 1E)
# We can then reassign the identity of the clusters to these cell types:
# Rename all identities
seurat_phase_anno <- RenameIdents(object = seurat_phase, 
                                  "0" = "Muller Glia",
                                  "1" = "Muller Glia",
                                  "2" = "RGC",
                                  "3" = "Muller Glia",
                                  "4" = "Muller Glia",
                                  "5" = "Muller Glia",
                                  "6" = "Muller Glia",
                                  "7" = "Muller Glia",
                                  "8" = "Muller Glia",
                                  "9" = "Muller Glia",
                                  "10" = "Muller Glia",
                                  "11" = "Muller Glia",
                                  "12" = "Muller Glia",
                                  "13" = "Microglia",
                                  "14" = "Bipolar",
                                  "15" = "Muller Glia",
                                  "16" = "Bipolar",
                                  "17" = "Muller Glia",
                                  "18" = "Muller Glia",
                                  "19" = "Muller Glia",
                                  "20" = "Muller Glia",
                                  "21" = "Muller Glia",
                                  "22" = "Bipolar",
                                  "23" = "Bipolar",
                                  "24" = "Muller Glia",
                                  "25" = "Muller Glia",
                                  "26" = "RGC",
                                  "27" = "Muller Glia",
                                  "28" = "RGC",
                                  "29" = "Muller Glia",
                                  "30" = "Muller Glia",
                                  "31" = "RGC")

# Plot the UMAP
DimPlot(object = seurat_phase_anno, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)

# If we wanted to remove the potentially unwanted cells, we could use the subset() function:
# Remove the unwanted clusters
seurat_phase_RemovedUnwantedClusters <- subset(seurat_phase_anno,
                                               idents = c("Muller Glia","Microglia","Bipolar"), invert = TRUE)
dim(seurat_phase_RemovedUnwantedClusters)

# Re-visualize the clusters
DimPlot(object = seurat_phase_RemovedUnwantedClusters, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3)

# ------------------------- Exclude non-RGC clusters and re-cluster with only RGC clusters
# --------------------------------------------- re-clustering (start from scratch)
#---------------------------------------------- Count normalization
seurat_phase_RemovedUnwantedClusters <- NormalizeData(seurat_phase_RemovedUnwantedClusters)
dim(seurat_phase_RemovedUnwantedClusters)

#---------------------------------------------- identification of the most variable genes
seurat_phase_RemovedUnwantedClusters <- FindVariableFeatures(seurat_phase_RemovedUnwantedClusters, 
                                                             selection.method = "vst",
                                                             nfeatures = 2000, 
                                                             verbose = FALSE)
dim(seurat_phase_RemovedUnwantedClusters)

#---------------------------------------------- Scale the counts
seurat_phase_RemovedUnwantedClusters <- ScaleData(seurat_phase_RemovedUnwantedClusters)
dim(seurat_phase_RemovedUnwantedClusters)

#---------------------------------------------- Principal Component Analysis (PCA)
# Perform PCA
seurat_phase_RemovedUnwantedClusters <- RunPCA(seurat_phase_RemovedUnwantedClusters)

# Plot the PCA colored by plate
DimPlot(seurat_phase_RemovedUnwantedClusters,
        reduction = "pca")

# Plot the PCA split by plate
PCAPlot(seurat_phase_RemovedUnwantedClusters,
        split.by = "plate") 

#---------------------------------------------- UMAP visualization
# Run UMAP
seurat_phase_RemovedUnwantedClusters <- RunUMAP(seurat_phase_RemovedUnwantedClusters, 
                                                dims = 1:40,
                                                reduction = "pca")
# Plot UMAP colored by plate                            
DimPlot(seurat_phase_RemovedUnwantedClusters) 

# Plot UMAP split by plates
DimPlot(seurat_phase_RemovedUnwantedClusters,
        split.by = "plate")  

#---------------------------------------------- Identify significant PCs 
# Explore heatmap of PCs
DimHeatmap(seurat_phase_RemovedUnwantedClusters, 
           dims = 1:18, 
           cells = 500, 
           balanced = TRUE)

# print out the top10 (or more) positive and negative genes by PCA scores driving the PCs.
# Printing out the most variable genes driving PCs
print(x = seurat_phase_RemovedUnwantedClusters[["pca"]], 
      dims = 1:18, 
      nfeatures = 5)

# Plot the elbow plot
# The elbow plot is another helpful way to determine how many PCs to use for clustering so that we are capturing majority of the variation in the data.
ElbowPlot(object = seurat_phase_RemovedUnwantedClusters, 
          ndims = 40)

#--------------------------
# Create metadata dataframe
cleanFinal_metadata <- seurat_phase_RemovedUnwantedClusters@meta.data

# Add cell IDs to metadata
cleanFinal_metadata$cells <- rownames(cleanFinal_metadata)

# Rename columns
cleanFinal_metadata <- cleanFinal_metadata %>%
  dplyr::rename(seq_folder1 = orig.ident,
                nCount = nCount_RNA,
                nFeature = nFeature_RNA)

# Add metadata back to Seurat object
seurat_phase_RemovedUnwantedClusters@meta.data <- cleanFinal_metadata

# Visualize the number of cell counts per plate
cleanFinal_metadata %>% 
  ggplot(aes(x=plate, fill=plate)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "white") +
  ggtitle("NCells")

#---------------------------------------------- Second round of clustering the cells
# To perform clustering, we determine the genes that are most different in their expression between cells. 
# Then, we use these genes to determine which correlated genes sets are responsible for the largest differences in expression between cells.
# Seurat uses a graph-based clustering approach, which embeds cells in a graph structure, using a K-nearest neighbor (KNN) graph (by default), with edges drawn between cells with similar gene expression patterns. Then, it attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.

# Determine the K-nearest neighbor graph
seurat_phase_RemovedUnwantedClusters <- FindNeighbors(object = seurat_phase_RemovedUnwantedClusters, dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_phase_RemovedUnwantedClusters <- FindClusters(object = seurat_phase_RemovedUnwantedClusters, resolution = c(0.4, 0.6, 0.8, 1, 1.4, 1.6))

# Explore resolutions
seurat_phase_RemovedUnwantedClusters@meta.data %>% 
  View()

# Assign identity of clusters
Idents(object = seurat_phase_RemovedUnwantedClusters) <- "RNA_snn_res.1"

# Plot the UMAP
DimPlot(seurat_phase_RemovedUnwantedClusters,
        reduction = "umap",
        label = TRUE,
        label.size = 3)

#---------------------------------------------- Segregation of clusters by plate
# Extract identity and plate information from Seurat object to determine the number of cells per cluster per plate
n_cells <- FetchData(seurat_phase_RemovedUnwantedClusters, 
                     vars = c("ident", "plate")) %>%
  dplyr::count(ident, plate) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

# UMAP of cells in each cluster by plates
DimPlot(seurat_phase_RemovedUnwantedClusters, 
        label = TRUE, 
        split.by = "plate", ncol=3, label.size = 3)

# Determine metrics to plot present in seurat_phase_RemovedUnwantedClusters@meta.data
metrics <-  c("nUMI", "nGene", "mitoRatio", "riboRatio")
FeaturePlot(seurat_phase_RemovedUnwantedClusters, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.3, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

#---------------------------------------------- Exploration of the PCs driving the different clusters
# Defining the information in the Seurat object of interest
columns <- c(paste0("PC_", 1:20),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the Seurat object
pc_data <- FetchData(seurat_phase_RemovedUnwantedClusters, 
                     vars = columns)

# Extract the UMAP coordinates for the first 10 cells
seurat_phase_RemovedUnwantedClusters@reductions$umap@cell.embeddings[1:10, 1:2]

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_phase_RemovedUnwantedClusters, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
# we can see that the genes driving in each PC exhibit higher expression in which cluster
map(paste0("PC_", 1:20), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# For instance, the genes driving PC_2 exhibit higher expression in clusters 6, 11, and 17 (maybe a bit higher in 15, too). 
# We could look back at our genes driving this PC to get an idea of what the cell types might be.
# Examine PCA results 
print(seurat_phase_RemovedUnwantedClusters[["pca"]], dims = 1:20, nfeatures = 5)

#---------------------------------------------- Exploring known cell type markers (Pan-RGCs)
DimPlot(object = seurat_phase_RemovedUnwantedClusters, 
        reduction = "umap", 
        label = TRUE)

# To access the expression levels of all genes, rather than just the 3000 most highly variable genes, 
# we can use the normalized count data stored in the RNA assay slot.

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_phase_RemovedUnwantedClusters) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_phase_RemovedUnwantedClusters <- NormalizeData(seurat_phase_RemovedUnwantedClusters, verbose = FALSE)

FeaturePlot(seurat_phase_RemovedUnwantedClusters, 
            reduction = "umap", 
            features = c("Rbpms", "Thy1", "Slc17a6", "Pou4f1", "Pou4f2", "Pou4f3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# Vln plot 
VlnPlot(object = seurat_phase_RemovedUnwantedClusters,
        features = c("Rbpms", "Thy1", "Slc17a6", "Pou4f1", "Pou4f2", "Pou4f3"))

# Dot Plot
DotPlot(object = seurat_phase_RemovedUnwantedClusters,
        features = c("Rbpms", "Thy1", "Slc17a6", "Pou4f1", "Pou4f2", "Pou4f3"))

# --------------------------------- Single-cell RNA-seq marker identification
# --------------------------------- Adding Gene Annotations
annotations <- read.csv("MetaData/00_Master_RGC_Y.csv")
# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat_phase_RemovedUnwantedClusters.markers <- FindAllMarkers(seurat_phase_RemovedUnwantedClusters)
seurat_phase_RemovedUnwantedClusters.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(seurat_phase_RemovedUnwantedClusters.markers, "All markers for each clusters in Round2 clustering.csv", row.names=FALSE)

# --------------------------------- Evaluating marker genes
# Extract top 30 markers per cluster
seurat_phase_RemovedUnwantedClusters.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top30

# Visualize top 30 markers per cluster
View(top30)
write.csv(top30, "top30 markers for each clusters in round2 clustering.csv", row.names=FALSE)

# --------------------------------- Visualizing marker genes
# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 30 markers for each cluster.
DoHeatmap(object = seurat_phase_RemovedUnwantedClusters, features= top30$gene)
DotPlot(object   = seurat_phase_RemovedUnwantedClusters,
        features = as.character(unique(top30$gene)),
        #group.by = "cluster",
        assay    = "RNA") +
  coord_flip()

#-------------------------
# To get a better idea of cell type identity for clusters we can explore the expression of different identified markers by cluster using the FeaturePlot() function.
# Plot interesting marker gene expression for clusters
FeaturePlot(object = seurat_phase_RemovedUnwantedClusters, 
            features = c("Sgcd","Spink10","Spink13"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)

# Vln plot - cluster 0
VlnPlot(object = seurat_phase_RemovedUnwantedClusters, 
        features = c("Sgcd","Spink10","Spink13"), split.by = 'plate')

#----------------------------- Based on Tran et. al. paper
# Dot Plot
DotPlot(object = seurat_phase_RemovedUnwantedClusters,
        features = c("Tbr1", "Foxp2", "Tusc5", "Opn4", "Spp1", "Satb2", "Cartpt"))

DotPlot(object   = seurat_phase_RemovedUnwantedClusters,
        features = c("Serpine2","Amigo2","Lypd1","Foxp2","Lrx4","Pde1a","Tbr1","Pcdh20","Zic1","Tbx20","Tagln2","Prkcq","Tac1","Slc7a11","Plpp4","Gpr88","Serpinb1b","Gm17750","Mmp17","Ntrk1","Cartpt","Vit","Apela","Col25a1","4833423E24Rik","Penk","Prdm8","Slc24a2","Gal","Calca","Cdhr1","Prokr1","Fam19a4","Slc17a7","lgfbp5","Prkcg","Cdk15","Stxbp6","Prlr","Postn","Spp1","Rhox5","Adcyap1","Opn4","Tpbg","Lgfbp4","Chrm2","Coch","Ceacam10","Anxa3","Neurod2","S100b","Nmb","Kit","Fes","Ll1rapl2","Bhlhe22","Fxyd6"))

# --------------------------------- Assigning cell type identity to clusters based on Tran et. al.
# We can then reassign the identity of the clusters to these cell types:
# Rename all identities
seurat_phase_RemovedUnwantedClusters_anno <- RenameIdents(object = seurat_phase_RemovedUnwantedClusters, 
                                                          "0" = "RGC1",
                                                          "1" = "RGC2",
                                                          "2" = "RGC3",
                                                          "3" = "RGC4",
                                                          "4" = "RGC5",
                                                          "5" = "RGC6",
                                                          "6" = "RGC7",
                                                          "7" = "RGC8",
                                                          "8" = "RGC9",
                                                          "9" = "RGC10",
                                                          "10" = "RGC11",
                                                          "11" = "RGC12",
                                                          "12" = "RGC13",
                                                          "13" = "RGC14")

# Plot the UMAP
DimPlot(object = seurat_phase_RemovedUnwantedClusters_anno, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)

