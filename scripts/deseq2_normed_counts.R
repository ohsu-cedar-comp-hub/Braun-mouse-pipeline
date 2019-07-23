library("dplyr")
library("DESeq2")

counts = snakemake@input[['counts']]

metadata <- snakemake@params[['samples']]
sampleID <- snakemake@params[['sample_id']]
Type <- snakemake@params[['linear_model']]
counts_out = snakemake@output[['counts_out']]


parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Read in metadata table and order according to sampleID
md <- read.delim(file=metadata, sep = "\t", stringsAsFactors = FALSE)
md <- md[order(md[sampleID]),]
rownames(md) <- md[[sampleID]]
md[[sampleID]] <- NULL


# Read in counts table
subdata <- read.table(counts, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
subdata <- subdata[,order(colnames(subdata))]


md <- md[colnames(subdata),,drop=FALSE]
# Check

print(md)
print(subdata)

stopifnot(rownames(md)==colnames(subdata))

# Obtain the number of genes that meet padj<0.01 for reference line in histogram
dds <- DESeqDataSetFromMatrix(countData=subdata,
                              colData=md,
                              design= as.formula(paste('~',Type)))

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file=counts_out, sep="\t", quote=F)
