library(Glimma)
library(limma)
library(DESeq2)

condition = snakemake@params[['condition']]
cat(sprintf(c('Condition: ',condition,'\n')))

title = snakemake@params[["contrast"]]

contrast = c(condition, snakemake@params[["contrast"]])
rds = snakemake@input[['rds']]
cat(sprintf(c('RDS object: ',rds,'\n')))

rds = readRDS(rds)
groups.df = as.data.frame(colData(rds))
glMDSPlot(rds, groups=groups.df)

rnaseq = as.data.frame(counts(rds, normalized=T))
genes = as.data.frame(row.names(res))
colnames(genes) = 'GeneID'

#### by contrasts
contrasts_to_plot = resultsNames(rds)
res <- results(rds, contrast=contrast)
res$padj[is.na(res$padj)] = 1
status_frame = res[,c('log2FoldChange','padj')]
status_frame['status'] = 0
status_frame$padj[is.na(status_frame$padj)] = 1
status_frame[status_frame$padj<0.05 & status_frame$log2FoldChange < 0 ,'status'] = -1
status_frame[status_frame$padj<0.05 & status_frame$log2FoldChange > 0 ,'status'] = 1

glMDPlot(res, anno=genes, status=status_frame$status, samples=colnames(rnaseq), 
         counts=log2(rnaseq + 0.0001), groups=groups.df$Condition, main=strsplit(res@elementMetadata$description[2],': ')[[1]][2], 
         transform=F, side.ylab='Log2-expression',side.main='GeneID', html = paste(title,'ma_plot',sep='.'))

## Volcano plot
glXYPlot(x=res$log2FoldChange, y=-log10(res$pvalue), xlab="logFC", ylab="logodds",
         status=status_frame$status, counts=log2(rnaseq + 0.0001), groups=groups.df$Condition, anno=genes,main=strsplit(res@elementMetadata$description[2],': ')[[1]][2],html = paste(title,'volcano_plot',sep='.'))

