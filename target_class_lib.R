# Load data exported from Python (see main.ipynb) into a DESeqDataSet object
DESeqDataSetFromPython <- function(cts_file, anno_file) {
    library("DESeq2")
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


# Plot and save PCA and tSNE plots on the top-ntop-variance genes
plot_pca_and_tsne <- function(transformed_data, transformation_name, data_dir, intgroup="condition", ntop=500) {
    # Note this function is based off the function plotPCA() in the DESeq2 package

    # Don't worry about supervised normal transformation as for vst and rlog as it's probably not recommended and not easy to just transform to new data using the same transformations
    # all(assay(dds) == counts(dds))
    # all(log2(counts(dds, normalized=TRUE)+1) == assay(ntd))

    # Sample call: plot_pca_and_tsne(ntd, "normal transformation", "/data/BIDS-HPC/private/projects/dmi2/data")

    # Import relevant libraries
    library("DESeq2")
    library("ggplot2")
    library("Rtsne")

    # Set the random seed for the tSNE analysis below
    # set.seed(44)    

    # Determine the top-ntop-variance genes
    rv <- rowVars(assay(transformed_data)) # take the variances of each row (gene) probably
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))] # take the gene indexes having the top ntop variances over the samples, i.e., the genes that matter most for the sample set

    # Run PCA
    pca <- prcomp(t(assay(transformed_data)[select,])) # do a principle components analysis on these most important genes (when I was doing it from scratch I was taking $rotation instead of $x, running it on assay(transformed_data) instead of t(assay(transformed_data)), and not taking just the top 500 genes)

    # Get the grouping data
    intgroup.df <- as.data.frame(colData(transformed_data)[, intgroup, drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
        colData(transformed_data)[[intgroup]]
    }

    # Plot and save the PCA
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, intgroup.df, name = colnames(transformed_data))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    png(paste0(data_dir, "/pca_", tolower(gsub(" ", "_", transformation_name)), ".png"), height = 6, width = 6, res = 150, units = 'in')
    print(ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + coord_fixed() + ggtitle(paste0("PCA - ", transformation_name))) # plot(pca$x[,1:2] ,col=colData(dds)$condition) produces something similar I believe
    dev.off()

    # Set the perplexity hyperparameter for tSNE to a reasonable value
    # “Typical values for the perplexity range between 5 and 50”; from help page: should not be bigger than 3 * perplexity < nrow(X) - 1
    nsamples <- dim(t(assay(transformed_data)))[1]
    perpl <- round((min(floor((nsamples-1)/3), 50) + 5) / 2)

    # Run tSNE
    tsne <- Rtsne(t(assay(transformed_data)[select,]), perplexity=perpl, pca=FALSE, theta=0.0) 

    # Plot and save the tSNE
    d <- data.frame(AX1 = tsne$Y[, 1], AX2 = tsne$Y[, 2], group = group, intgroup.df, name = colnames(transformed_data))
    png(paste0(data_dir, "/tsne_", tolower(gsub(" ", "_", transformation_name)), ".png"), height = 6, width = 6, res = 150, units = 'in')
    print(ggplot(data = d, aes_string(x = "AX1", y = "AX2", color = "group")) +
        geom_point(size = 3) +
        xlab("Axis 1") +
        ylab("Axis 2") +
        coord_fixed() +
        ggtitle(paste0("tSNE (perplexity: ", perpl, ") - ", transformation_name))) # this results in the same plot as plot(tsne$Y,col=colData(dds)$condition, asp=1)
    dev.off()

}


rlogSupervised <- function(ddsOld, ddsNew) {
    library("DESeq2")

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
    library("DESeq2")

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


# # rlog - supervised
# ddsOld <- makeExampleDESeqDataSet(m=6,betaSD=1)
# ddsNew <- makeExampleDESeqDataSet(m=6,betaSD=1)
# rlogSupervised(ddsOld, ddsNew)

# # vst - supervised
# ddsOld <- makeExampleDESeqDataSet(m=6)
# ddsNew <- makeExampleDESeqDataSet(m=1)
# vstSupervised(ddsOld, ddsNew)
