library(Glimma)
library(limma)
library(DESeq2)

project_id = snakemake@params[['project_id']]

rds = snakemake@input[['rds']]
cat(sprintf(c('RDS object: ',rds,'\n')))

rds = readRDS(rds)
groups.df = as.data.frame(colData(rds))
glMDSPlot(rds, groups=groups.df, html = paste(project_id,'mds_plot',sep='.'))

