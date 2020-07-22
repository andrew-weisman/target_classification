DESeqDataSetFromPython <- function(cts_file, anno_file) {
    cts <- as.matrix(read.csv(cts_file, sep="\t", row.names="gene_id"))
    coldata <- read.csv(anno_file, row.names=1)
    coldata$condition <- gsub("-", "_", coldata$condition)
    coldata$condition <- gsub(",", ".", coldata$condition)
    coldata$condition <- gsub(" ", "", coldata$condition)
    coldata$condition <- factor(coldata$condition)
    rownames(coldata) <- gsub("-", "_", rownames(coldata))
    colnames(cts) <- gsub("\\.", "_", colnames(cts))
    all(rownames(coldata) == colnames(cts))
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
    return(dds)
}


rlogUnsupervised <- function(dds) {
    rld <- rlog(dds)
    dists <- dist(t(assay(rld)))
    plot(hclust(dists))
}


vstUnsupervised <- function(dds) {
    vsd <- varianceStabilizingTransformation(dds)
    dists <- dist(t(assay(vsd)))
    plot(hclust(dists))
}


rlogSupervised <- function(ddsOld, ddsNew) {

    # run the rlog transformation on one dataset
    design(ddsOld) <- ~ 1
    ddsOld <- estimateSizeFactors(ddsOld)
    ddsOld <- estimateDispersions(ddsOld)
    rldOld <- rlog(ddsOld, blind=FALSE)

    # apply the parameters to a new sample
    mcols(ddsNew)$dispFit <- mcols(ddsOld)$dispFit
    betaPriorVar <- attr(rldOld,"betaPriorVar")
    intercept <- mcols(rldOld)$rlogIntercept
    rldNew <- rlog(ddsNew, blind=FALSE, intercept=intercept, betaPriorVar=betaPriorVar)

    return(rldNew)

}


vstSupervised <- function(ddsOld, ddsNew) {

    # learn the dispersion function of a dataset
    design(ddsOld) <- ~ 1
    ddsOld <- estimateSizeFactors(ddsOld)
    ddsOld <- estimateDispersions(ddsOld)

    # use the previous dispersion function for a new sample
    ddsNew <- estimateSizeFactors(ddsNew)
    dispersionFunction(ddsNew) <- dispersionFunction(ddsOld)
    vsdNew <- varianceStabilizingTransformation(ddsNew, blind=FALSE)

    return(vsdNew)

}


# Import relevant library
library("DESeq2")

# Original workflow
dds <-- DESeqDataSetFromPython("/data/BIDS-HPC/private/projects/dmi2/data/gene_counts.tsv", "/data/BIDS-HPC/private/projects/dmi2/data/annotation.csv")
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds) # note no blind option; this does log2(count(dds,normalized=TRUE) + 1)
plotPCA(vsd, intgroup="condition")

# rlog - unsupervised
dds <- makeExampleDESeqDataSet(m=6,betaSD=1)
rlogUnsupervised(dds)

# rlog - supervised
ddsOld <- makeExampleDESeqDataSet(m=6,betaSD=1)
ddsNew <- makeExampleDESeqDataSet(m=6,betaSD=1)
rlogSupervised(ddsOld, ddsNew)

# vst - unsupervised
dds <- makeExampleDESeqDataSet(m=6)
vstUnsupervised(dds)

# vst - supervised
ddsOld <- makeExampleDESeqDataSet(m=6)
ddsNew <- makeExampleDESeqDataSet(m=1)
vstSupervised(ddsOld, ddsNew)
