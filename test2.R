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

    #rld <- rlog(dds) # original example
    rld <- rlog(dds, blind=FALSE) # my version makes it not blind

    #plot(hclust(dist(t(assay(rld))))) # original example
    plotPCA(rld, intgroup="condition") # my version runs PCA instead of plotting a distance matrix
}


vstUnsupervised <- function(dds) {

    #vsd <- varianceStabilizingTransformation(dds) # original example
    vsd <- vst(dds, blind=FALSE) # my version uses vst() instead of varianceStabilizingTransformation() and also makes it not blind

    #plot(hclust(dist(t(assay(vsd))))) # original example
    plotPCA(vsd, intgroup="condition") # my version runs PCA instead of plotting a distance matrix

}


rlogSupervised <- function(ddsOld, ddsNew) {

    # Run the rlog transformation on one dataset
    #design(ddsOld) <- ~ 1 # from ?vst (at least the webpage) I believe this makes rlog() strictly blind to the experimental design, which we don't want, so I'm commenting it out
    ddsOld <- estimateSizeFactors(ddsOld) # determine the size factors of ddsOld
    ddsOld <- estimateDispersions(ddsOld) # determine the dispersions of ddsOld
    rldOld <- rlog(ddsOld, blind=FALSE) # determine the rlog of ddsOld

    # Apply the parameters to a new sample
    mcols(ddsNew)$dispFit <- mcols(ddsOld)$dispFit # make the dispersion function of ddsNew the same as that of ddsOld
    betaPriorVar <- attr(rldOld,"betaPriorVar") # determine the beta prior variance of the rlog of ddsOld
    intercept <- mcols(rldOld)$rlogIntercept # determine the rlog intercept of the rlog of ddsOld
    rldNew <- rlog(ddsNew, blind=FALSE, intercept=intercept, betaPriorVar=betaPriorVar) # determine the rlog of ddsNew using the beta prior variance and rlog intercept of the rlog of ddsOld

    return(rldNew)

}


vstSupervised <- function(ddsOld, ddsNew) {

    # Learn the dispersion function of a dataset
    #design(ddsOld) <- ~ 1 # from ?vst (at least the webpage) I believe this makes vst() (even though it's not explicitly run?) strictly blind to the experimental design, which we don't want, so I'm commenting it out
    ddsOld <- estimateSizeFactors(ddsOld) # determine the size factors of ddsOld
    ddsOld <- estimateDispersions(ddsOld) # determine the dispersions of ddsOld

    # Use the previous dispersion function for a new sample
    ddsNew <- estimateSizeFactors(ddsNew) # determine the size factors of ddsNew
    dispersionFunction(ddsNew) <- dispersionFunction(ddsOld) # set the dispersions of ddsNew to those of ddsOld
    vsdNew <- varianceStabilizingTransformation(ddsNew, blind=FALSE) # determine the VST of ddsNew using the dispersions of ddsOld

    return(vsdNew)

}


# Import relevant library
library("DESeq2")

# Original workflow
dds <-- DESeqDataSetFromPython("/data/BIDS-HPC/private/projects/dmi2/data/gene_counts.tsv", "/data/BIDS-HPC/private/projects/dmi2/data/annotation.csv")
dds <- DESeq(dds) # can we omit this line in this block of three lines?
ntd <- normTransform(dds) # note no blind option; this does log2(count(dds,normalized=TRUE) + 1)

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
