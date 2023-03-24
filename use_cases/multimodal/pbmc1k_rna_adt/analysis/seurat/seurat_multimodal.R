library(dplyr)
library(Seurat)
library(patchwork)


# Load the PBMC dataset
data.filtered <- Read10X(data.dir = "../data/filtered_gene_bc_matrices/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc.full <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc.full[["percent.mt"]] <- PercentageFeatureSet(pbmc.full, pattern = "^MT-")

pbmc <- subset(pbmc.full, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Save foms in simple flat files
Matrix::writeMM(GetAssayData(object = pbmc.full, slot = "counts", assay = "RNA"), file = "../output/fom.rna.full.full.counts.mtx")
R.utils::gzip("../output/fom.rna.full.full.counts.mtx")
Matrix::writeMM(GetAssayData(object = pbmc, slot = "counts", assay = "RNA"), file = "../output/fom.rna.clean.full.counts.mtx")
R.utils::gzip("../output/fom.rna.clean.full.counts.mtx")
Matrix::writeMM(GetAssayData(object = pbmc, slot = "data", assay = "RNA"), file = "../output/fom.rna.clean.full.data.mtx")
R.utils::gzip("../output/fom.rna.clean.full.data.mtx")
write.table(GetAssayData(object = pbmc, slot = "scale.data", assay = "RNA"), gzfile("../output/fom.rna.clean.var.scale.txt.gz"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

# Save reductions (tranposing to keep cells in columns)
write.table(t(Embeddings(Reductions(pbmc, "pca"))), gzfile("../output/fom.rna.clean.full.pca.txt.gz"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(t(Embeddings(Reductions(pbmc, "umap"))), gzfile("../output/fom.rna.clean.full.umap.txt.gz"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# Save observation annotations
write.table(pbmc.full[[]], gzfile("../output/oam.rna.full.txt.gz"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(pbmc[[]], gzfile("../output/oam.rna.clean.txt.gz"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Save IDs
write.table(colnames(pbmc.full), gzfile("../output/oid.rna.full.txt.gz"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(colnames(pbmc), gzfile("../output/oid.rna.clean.txt.gz"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(rownames(pbmc), gzfile("../output/fid.rna.full.txt.gz"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(rownames(GetAssayData(object = pbmc, slot = "scale.data", assay = "RNA")), gzfile("../output/fid.rna.variable.txt.gz"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(colnames(Embeddings(Reductions(pbmc, "pca"))), gzfile("../output/fid.rna.full.pca.txt.gz"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(colnames(Embeddings(Reductions(pbmc, "umap"))), gzfile("../output/fid.rna.full.umap.txt.gz"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# Save feature annotations
write.table(Assays(pbmc, slot = "RNA")@meta.features, gzfile("../output/fam.rna.full.txt.gz"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(data.frame("Stdev" = Stdev(Reductions(pbmc, "pca"))), gzfile("../output/fam.rna.clean.pca.txt.gz"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Save Seurat objects
saveRDS(pbmc.full, "../output/pbmc.full.rds")
saveRDS(pbmc, "../output/pbmc.rds")
sessionInfo()
