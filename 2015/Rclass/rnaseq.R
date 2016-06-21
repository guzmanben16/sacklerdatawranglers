# load your libraries
library("RNAseqData.HNRNPC.bam.chr14")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("limma")
library("edgeR")
library("statmod")

# set your working directory
setwd("~/Desktop/Rclass")

# the BAM files were installed on your computer, let's find out where
# if this was for real, just make a string of the path to your BAM files
# dir <- "~/Desktop/Rclass"
dir <- system.file("extdata", package="RNAseqData.HNRNPC.bam.chr14", mustWork=TRUE)

# print what files you have
list.files(dir)

# load your sample info data
# you make this file yourself
targets <- read.delim("targets.txt", sep='\t', header=TRUE)

# attach the paths to the BAM files and make them into Rsamtools objects
filenames <- file.path(dir, targets$fileName)
bamfiles <- BamFileList(filenames, yieldSize=2000000)

# get the GTF file of known genes from your favourite database
# 
txdb = makeTxDbFromGFF("UCSC_genes_chr14.gtf", 
                               format="gtf", 
                               circ_seqs=character())

# group the exons into genes
ebg <- exonsBy(txdb, by="gene")

# collect a matrix of the number of reads per gene
se <- summarizeOverlaps(features=ebg, 
                        reads=bamfiles, 
                        mode="Union", 
                        singleEnd=FALSE, 
                        ignore.strand=TRUE, 
                        fragments=TRUE)

# store the counts matrix as a variable
# we will be using this from now on
counts <- assay(se)

# plot checks
# check distributions for batch effects
boxplot(log2(counts+1),
        las=2,
        range=0,
        ylab="log2 intensity")

# check PCA plot for batch effects and target separability
plotMDS(log2(counts+1), 
        labels=targets$batch, 
        col=ifelse(targets$batch==1,"blue","red"),
        xlab="PC 1",
        ylab="PC 2")

plotMDS(log2(counts+1), 
        labels=targets$condition, 
        col=ifelse(targets$condition=="Cntl","blue","red"),
        xlab="PC 1",
        ylab="PC 2")

# remember, we made a variable called targets with file data
# load in accession annotations
kgxref <-  read.delim("kgxref.csv", header=TRUE, row.names=1)
acc <- data.frame(accessions=rownames(counts))
anno <- merge(acc, kgxref,
              by.x="accessions",
              by.y="row.names", 
              all.x=TRUE)

# make DGEList object
y <- DGEList(counts=counts, genes=anno$geneSymbol)

# filter and reset libraries
A <- rowSums(y$counts)
isexpr <- A > 5
y <- y[isexpr, ]
y$samples$lib.size <- colSums(y$counts)

# design matrix
f <- factor(targets$condition)
design <- model.matrix(~0+f)

# change the contrast to pick the comparisons you want
# for differential expression, do one at a time 
contrast <- makeContrasts(fsiRNA1-fCntl, 
                          levels=design)

# normalization and voom
y <- calcNormFactors(y)

png(filename="voom.png")
v <- voom(y, design, plot=TRUE)
dev.off()

# estimate correlation between technical reps
cor <- duplicateCorrelation(v, 
                            design=design, 
                            block=targets$replicate)

# make the fit
fit <- lmFit(v, 
             design, 
             block=targets$replicate, 
             correlation=cor$consensus, 
             weights=v$weights)

fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
tt <- topTable(fit2, number=1000, adjust="BH", p.value=0.05)

# write data to file
write.table(tt, 
            file="topTable.csv", 
            row.names=TRUE, 
            col.names=NA, 
            sep=",")

