# Plots genes with significant change in CERES
# Input:
# @deltas: tissue agnostic significant differences between mutants and wild-types
# @tissue.specific.sig.genes: tissue specific significant differences between mutants and wild-types

arguments <- commandArgs(TRUE)
deltas <- read.table(file = arguments[1], sep = '\t', header=T, stringsAsFactors=F)
tissue.specific.sig.genes <- readLines(con=arguments[2])

# Sort the tissue agnostic deltas for waterfall plots
ordered.deltas <- deltas[order(deltas$q, decreasing=TRUE, na.last=NA),]
ordered.deltas$q <- -log10(ordered.deltas$q)

png(filename = "change_in_CERES_SCORE_plot.png")
par(mar = c(5, 5, 2, 2)) 

# Plot the change in CERES score FDR-corrected nomal p-value (-log10) for significant deltas
barplot(ordered.deltas$q, col="black", border="black", space=0.5, names.arg = ordered.deltas$Gene, cex.axis = 0.6,
        ylab="Genes", xlab="change in CERES score FDR-corrected normal p-value (-log10)", cex.lab=0.8, las=2, 
        cex.names = 0.5, horiz = TRUE, margin(10,10))

dev.off()

# Overlap genes between the tissue specific and agnostic
overlap.genes <- which(tissue.specific.sig.genes %in% ordered.deltas$Gene)
ordered.deltas.sub <- ordered.deltas[tissue.specific.sig.genes[overlap.genes],]
ordered.deltas.sub <- ordered.deltas.sub[order(ordered.deltas.sub$q, decreasing=FALSE, na.last=NA),]
ordered.deltas.sub$q <- -log10(ordered.deltas.sub$q)

# Plot the change in CERES score FDR-corrected nomal p-value for overlap genes
barplot(ordered.deltas.sub$q, col="black", border="black", space=0.5, 
        ylab="change in CERES score FDR-corrected normal p-value (-log10)", names.arg = ordered.deltas.sub$Gene, 
        cex.axis = 0.8, xlab="Genes", cex.lab=0.8, las=2, cex.names = 0.7, ylim = c(0,100))


