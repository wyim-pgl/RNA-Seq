library(cluster)
library(Biobase)
library(qvalue)
library(fastcluster)
options(stringsAsFactors = FALSE)
NO_REUSE = F

# try to reuse earlier-loaded data if possible
if (file.exists("RSEM_expected_count.tsv.RData") && ! NO_REUSE) {
    print('RESTORING DATA FROM EARLIER ANALYSIS')
    load("RSEM_expected_count.tsv.RData")
} else {
    print('Reading matrix file.')
    primary_data = read.table("output/counts/RSEM/RSEM_expected_count.tsv", header=T, com='', row.names=1, check.names=F, sep='\t')
    primary_data = as.matrix(primary_data)
}
source("/data/gpfs/assoc/pgl/bin/miniconda3/envs/PGRP_Snakemake/opt/trinity-2.13.2/Analysis/DifferentialExpression/R/heatmap.3.R")
source("/data/gpfs/assoc/pgl/bin/miniconda3/envs/PGRP_Snakemake/opt/trinity-2.13.2/Analysis/DifferentialExpression/R/misc_rnaseq_funcs.R")
source("/data/gpfs/assoc/pgl/bin/miniconda3/envs/PGRP_Snakemake/opt/trinity-2.13.2/Analysis/DifferentialExpression/R/pairs3.R")
source("/data/gpfs/assoc/pgl/bin/miniconda3/envs/PGRP_Snakemake/opt/trinity-2.13.2/Analysis/DifferentialExpression/R/vioplot2.R")
data = primary_data
myheatcol = colorpanel(75, 'purple','black','yellow')
samples_data = read.table("replication_relationship.txt", header=F, check.names=F, fill=T)
samples_data = samples_data[samples_data[,2] != '',]
colnames(samples_data) = c('sample_name', 'replicate_name')
sample_types = as.character(unique(samples_data[,1]))
rep_names = as.character(samples_data[,2])
data = data[, colnames(data) %in% rep_names, drop=F ]
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
names(sample_colors) = sample_types
sample_type_list = list()
for (i in 1:nsamples) {
    samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
    sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
data = data[rowSums(data)>=10,]
initial_matrix = data # store before doing various data transformations
cs = colSums(data)
data = t( t(data)/cs) * 1e6;
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix

# Centering rows
data = t(scale(t(data), scale=F))

write.table(data, file="RSEM_expected_count.tsv.minRow10.CPM.log2.centered.dat", quote=F, sep='	');
if (nrow(data) < 2) { stop("

**** Sorry, at least two rows are required for this matrix.

");}
if (ncol(data) < 2) { stop("

**** Sorry, at least two columns are required for this matrix.

");}
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
write.table(sample_cor, file="RSEM_expected_count.tsv.minRow10.CPM.log2.centered.sample_cor.dat", quote=F, sep='	')
sample_dist = dist(t(data), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')
pdf("RSEM_expected_count.tsv.minRow10.CPM.log2.centered.sample_cor_matrix.pdf")
sample_cor_for_plot = sample_cor
if (is.null(hc_samples)) { RowV=NULL; ColV=NULL} else { RowV=as.dendrogram(hc_samples); ColV=RowV }
heatmap.3(sample_cor_for_plot, dendrogram='both', Rowv=RowV, Colv=ColV, col = myheatcol, scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=1, cexRow=1, cex.main=0.75, main=paste("sample correlation matrix
", "RSEM_expected_count.tsv.minRow10.CPM.log2.centered") , ColSideColors=sampleAnnotations, RowSideColors=t(sampleAnnotations))
dev.off()
pdf("RSEM_expected_count.tsv.minRow10.CPM.log2.centered.prcomp.principal_components.pdf")
prin_comp_data = data
pca = prcomp(prin_comp_data, center = FALSE, scale. = FALSE)
pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1));
write.table(pca$rotation, file="RSEM_expected_count.tsv.minRow10.CPM.log2.centered.PCA.prcomp.scores", quote=F, sep="	")
write.table(pca$x, file="RSEM_expected_count.tsv.minRow10.CPM.log2.centered.PCA.prcomp.loadings", quote=F, sep="	")
PCA.loadings=pca$x
PCA.scores = pca$rotation
for (i in 1:(max(3,2)-1)) {
    xrange = range(PCA.scores[,i])
    yrange = range(PCA.scores[,i+1])
    samples_want = rownames(PCA.scores) %in% sample_type_list[[sample_types[1]]]
    pc_i_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i]*100)
    pc_i_1_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i+1]*100)
    plot(PCA.scores[samples_want,i], PCA.scores[samples_want,i+1], xlab=paste('PC',i, pc_i_pct_var), ylab=paste('PC',i+1, pc_i_1_pct_var), xlim=xrange, ylim=yrange, col=sample_colors[1])
    for (j in 2:nsamples) {
        samples_want = rownames(PCA.scores) %in% sample_type_list[[sample_types[j]]]
        points(PCA.scores[samples_want,i], PCA.scores[samples_want,i+1], col=sample_colors[j], pch=j)
    }
    plot.new()
    legend('topleft', as.vector(sample_types), col=sample_colors, pch=1:nsamples, ncol=2)
}

par(def.par)
pcloadings_mat_vals = PCA.loadings[,1:3]
print(dim(pcloadings_mat_vals))
pcloadings_mat = matrix_to_color_assignments(pcloadings_mat_vals, col=colorpanel(256,'purple','black','yellow'), by='col')
print(dim(pcloadings_mat))
colnames(pcloadings_mat) = paste('PC', 1:ncol(pcloadings_mat))
dev.off()
gene_cor = NULL
