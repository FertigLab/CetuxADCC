# analysis of exported Illumina data using beadarray

require(lumi)
require(limma)
require(reshape)
require(lattice)
require(gplots)

##### graphical parameters
lattice.options(
	default.theme = modifyList(standard.theme(color=FALSE), list(strip.background = list(col = "transparent")))
)

##### helper functions

write.tab <- function(df, sep="\t", quote=FALSE, row.names=FALSE, ...) {
	write.table(df, sep=sep, quote=quote, row.names=row.names, ...)
}

find.gene <- function(df, pattern) {
	hits <- grep(pattern, df$SYMBOL)
	hits <- intersect(hits, grep(pattern, df$DEFINITION))
	return(df[ hits, ])
}

plot.gene <- function(df, symbol, as.summary=FALSE) {
	df.for.symbol <- df[ fData(df)$SYMBOL == symbol, ]
	if (nrow(df.for.symbol) < 1) {
		stop(paste("Unable to find data for", symbol))
	}
	df.for.symbol.melt <- melt(exprs(df.for.symbol))
	names(df.for.symbol.melt) <- c("Probe.ID", "Sample", "Expression")
	df.for.symbol.melt$Treatment <- factor(substr(df.for.symbol.melt$Sample, 1, 2))
	df.for.symbol.melt$ADCC.Status <- factor(
		ifelse(substr(df.for.symbol.melt$Sample, 4, 4) == "S", "A431", "A431/ADCCR")
	)
	if (as.summary) {
		p <- bwplot(
			Expression ~ ADCC.Status,
			df.for.symbol.melt,
			ylab=paste(symbol, " transcript level,\n", length(unique(df.for.symbol.melt$Probe.ID)), " unique probes", sep=""),
		)
	} else {
		p <- bwplot(
			Expression ~ ADCC.Status | Probe.ID,
			df.for.symbol.melt,
			ylab=paste(symbol, " transcript level"),
			layout=c(length(levels(df.for.symbol.melt$Probe.ID)), 1)
		)		
	}
	p <- update(p, 
		ylim=c(0, max(df.for.symbol.melt$Expression) + 1),
		scales=list(x=list(rot=45)),
	)
	return(p)
}

##### start analysis

do.analysis <- TRUE
if (do.analysis) {

##### read in data

sample.probe.file <- "../Data/sample_probe_profile.txt"
qc.probe.file <- "../Data/control_probe_profile.no_hyb.txt"

lumi.data <- lumiR.batch(sample.probe.file, annotationColumn=c("ENTREZ_GENE_ID", "SYMBOL", "CHROMOSOME", "DEFINITION"))
lumi.data <- addControlData2lumi(qc.probe.file, lumi.data)
sampleNames(lumi.data) <- c(paste(30:35, "S", sep="-"), paste(30:35, "R", sep="-"))

##### begin analysis

### process data using variance stabilition (log2) and normalization (quantile)
lumi.data.analyzed <- lumiExpresso(
	lumi.data,
	varianceStabilize.param=list(method="log2"), 
	normalize.param=list(method="quantile")
)
write.exprs(lumi.data.analyzed, file="output/lumi.data.vs_log2.n_quantile.analyzed.txt")

### output files for gene set enrichment analysis
lumi.data.analyzed.gsea <- lumi.data.analyzed[ fData(lumi.data.analyzed)$SYMBOL != "", ]
lumi.data.analyzed.gsea.txt <- data.frame(
	name=fData(lumi.data.analyzed.gsea)$SYMBOL, 
	description=fData(lumi.data.analyzed.gsea)$DEFINITION
)
lumi.data.analyzed.gsea.txt <- cbind(lumi.data.analyzed.gsea.txt, exprs(lumi.data.analyzed.gsea)[ , 1:12])
write.tab(lumi.data.analyzed.gsea.txt, file="output/lumi.data.vs_log2.n_quantile.analyzed.gsea.txt")

##### assess for differentially expressed genes

# used only detected ("expressed") probes
lumi.data.analyzed.selected <- lumi.data.analyzed[ detectionCall(lumi.data) > 0, ]

# assess sample differences across sensitive vs. resistant phenotypes
lumi.data.samples <- c(rep("Sensitive", 6), rep("Resistant", 6))
lumi.data.design <- model.matrix(~ factor(lumi.data.samples))
colnames(lumi.data.design) <- c("Sensitive", "Sensitive-Resistant")
lumi.data.analyzed.selected.fit <- lmFit(lumi.data.analyzed.selected, lumi.data.design)
lumi.data.analyzed.selected.fit <- eBayes(lumi.data.analyzed.selected.fit)

# find top 100 differentially-regulated genes
lumi.data.analyzed.selected.fit.top100 <- topTable(lumi.data.analyzed.selected.fit, coef="Sensitive-Resistant", number=100, adjust.method="fdr")

# find all differentially-regulated genes with FDR adjusted p-values < 0.01
lumi.data.analyzed.selected.fit.top <- topTable(lumi.data.analyzed.selected.fit, coef="Sensitive-Resistant", number=Inf, adjust.method="fdr", p.value=0.01)
#top.data <- lumi.data.analyzed.selected[ fData(lumi.data.analyzed.selected)$PROBE_ID %in% lumi.data.analyzed.selected.fit.top$PROBE_ID, ]

# ... with FDR-adjust p-values < 0.01 and fold change > 2
lumi.data.analyzed.selected.fit.top.p_lt_0.01_lfc_gte_1 <- topTable(lumi.data.analyzed.selected.fit, coef="Sensitive-Resistant", number=Inf, adjust.method="fdr", p.value=0.01, lfc=1)
top.data <- lumi.data.analyzed.selected[ fData(lumi.data.analyzed.selected)$PROBE_ID %in% lumi.data.analyzed.selected.fit.top.p_lt_0.01_lfc_gte_1$PROBE_ID, ]
write.tab(
	lumi.data.analyzed.selected.fit.top.p_lt_0.01_lfc_gte_1, 
	file="output/lumi.data.analyzed.selected.fit.top.p_lt_0.01_lfc_gte_1.txt"
)

# ... with FDR-adjusted p-value < 0.05 and fold change > 1.5
lumi.data.analyzed.selected.fit.top.p_lt_0.05_fc_gte_1.5 <- topTable(lumi.data.analyzed.selected.fit, coef="Sensitive-Resistant", number=Inf, adjust.method="fdr", p.value=0.05, lfc=log(1.5))
top.data <- lumi.data.analyzed.selected[ fData(lumi.data.analyzed.selected)$PROBE_ID %in% lumi.data.analyzed.selected.fit.top.p_lt_0.05_fc_gte_1.5$PROBE_ID, ]
write.tab(
	lumi.data.analyzed.selected.fit.top.p_lt_0.05_fc_gte_1.5, 
	file="output/lumi.data.analyzed.selected.fit.top.p_lt_0.05_fc_gte_1.5.txt"
)

}

##### plots

# visualize top hits (adjusted p<0.01 % fc > 2)
top.data.heatmap <- heatmap.2(exprs(top.data), scale="row", trace="none", labCol=NA, labRow=NA, density.info="none", rowsep=NA, col=bluered)

# microRNA plot
mir.rows <- substr(fData(lumi.data.analyzed.selected)$SYMBOL, 1, 3) == "MIR" #grep("^MIR*$", fData(lumi.data.analyzed.selected)$SYMBOL)
mir.data <- lumi.data.analyzed.selected[ mir.rows, ]
mir.heatmap <- heatmap.2(exprs(mir.data), trace="none", rowsep=NA, col=bluered, key=TRUE, labRow=fData(mir.data)$SYMBOL)

##### relationships to ADCC and network analysis

# check for differential expression among ADCC gene hits
ADCC.screen.genes <- read.table("ADCC_screen_target_genes.txt", header=TRUE)
ADCC.screen.genes.de <- intersect(lumi.data.analyzed.selected.fit.top$SYMBOL, ADCC.screen.genes$SYMBOL)

# check for differential expression among 638 EGFR siRNA library
EGFR.library.genes <- read.table("EGFR_siRNA_library_target_genes.txt", header=TRUE)
EGFR.library.genes.de <- intersect(lumi.data.analyzed.selected.fit.top$SYMBOL, EGFR.library.genes$SYMBOL)

# check for differential expression among 61 (or 62 if you include AURK?) or 93 EGFR siRNA library "hits" found on EGFR siRNA "Newplate"
EGFR.library.gene.hits <- read.table("EGFR_siRNA_library_target_gene_hits.txt", header=TRUE)
EGFR.library.gene.hits.de <- intersect(lumi.data.analyzed.selected.fit.top$SYMBOL, EGFR.library.gene.hits$SYMBOL)

# more heatmaps for EGFR network stuff
top.data.network <- top.data[ fData(top.data)$SYMBOL %in% EGFR.library.genes$SYMBOL, ]
#top.data.network.heatmap <- heatmap.2(exprs(top.data.network), scale="row", trace="none", density.info="none", col=bluered, key=FALSE, labRow=fData(top.data.network)$SYMBOL, cexRow=1.2)

top.data.network.hits <- top.data[ fData(top.data)$SYMBOL %in% ADCC.screen.genes$SYMBOL, ]
#top.data.network.hits.heatmap <- heatmap.2(exprs(top.data.network.hits), scale="row", trace="none", density.info="none", col=bluered, key=TRUE, labRow=fData(top.data.network.hits)$SYMBOL, cexRow=1.1)


##### generation of heatmap data files for circos

do.circos <- FALSE
if (do.circos) {

circos.dir <- "circos/data/"
gene.chr <- read.delim("ucsc/genes.txt")

# exclude 'other' chromosomal data
gene.chr <- gene.chr[ 
	gene.chr$hg19.knownGene.chrom %in% c(
		paste0("chr", 1:22), 
		paste0("chr", c("X", "Y"))
	),
]
# simplify some names and adjust chromosomes to indicate homo sapiens ('hs')
names(gene.chr)[ names(gene.chr) == "hg19.knownGene.txStart" ] <- "start"
names(gene.chr)[ names(gene.chr) == "hg19.knownGene.txEnd" ] <- "end"
names(gene.chr)[ names(gene.chr) == "hg19.kgXref.geneSymbol" ] <- "gene.symbol"
gene.chr$chr <- paste0("hs", substr(gene.chr$hg19.knownGene.chrom, 4, 6))

# output each sample to a separate heatmap data file
data.for.circos <- top.data
for (sample in sampleNames(phenoData(data.for.circos))) {
	heatmap.filename <- paste0(circos.dir, sample, ".txt")
	sample.data <- data.frame(
		gene.symbol=fData(data.for.circos)[ , c("SYMBOL")],
		value=exprs(data.for.circos)[ , which(sampleNames(phenoData(data.for.circos)) == sample)]
	)
	# TODO: account for multiple probe sets at given gene locus
	sample.data <- cast(
		sample.data,
		gene.symbol ~ .,
		fun.aggregate=c(median, max)
	)
	sample.data <- join(
		sample.data,
		gene.chr,
		by="gene.symbol"
	)
	sample.data <- sample.data[ !is.na(sample.data$chr), ]
	write.tab(sample.data[ , c("chr", "start", "end", "max")], file=heatmap.filename, col.names=FALSE)
}

# write gene loci
genes.for.circos <- data.frame(
	gene.symbol=fData(data.for.circos)[ , c("SYMBOL")]
)
genes.for.circos <- join(
	genes.for.circos,
	gene.chr,
	by="gene.symbol"
)
genes.for.circos.min <- cast(
	genes.for.circos,
	chr + gene.symbol ~ .,
	fun.aggregate=min,
	value="start"
)
genes.for.circos.min$start <- genes.for.circos.min$"(all)"
genes.for.circos.max <- cast(
	genes.for.circos,
	chr + gene.symbol ~ .,
	fun.aggregate=max,
	value="end"
)
genes.for.circos.max$end <- genes.for.circos.max$"(all)"
genes.for.circos <- join(
	genes.for.circos.min,
	genes.for.circos.max,
	by="gene.symbol"
)
write.tab(genes.for.circos[ , c("chr", "start", "end", "gene.symbol")], file=paste0(circos.dir, "genes.txt"), col.names=FALSE)

}
