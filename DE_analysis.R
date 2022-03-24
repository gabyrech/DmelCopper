library("DESeq2")

unSortedCounts <- read.table(file="geneCounts.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
countsFile <- unSortedCounts[ , order(names(unSortedCounts))]

# Sample Sheet
unSortedSheet <- read.table(file="sample_sheet.txt", sep="\t", stringsAsFactors=F, header=T)
targetsFile <- unSortedSheet[order(unSortedSheet$Sample),]

rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample

library("DESeq2")

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Strain + Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

rld <- rlogTransformation(dds)
vsd <- vst(dds)

#DE results:
res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)


###### DE Experiments:
#Resistant
# Gene counts
countsFile <- read.table(file="counts_Resistant.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_Resistant.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Strain + Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)


#Sensitive
# Gene counts
countsFile <- read.table(file="counts_Sensitive.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_Sensitive.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Strain + Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)



#Per strain

#MUN33
# Gene counts
countsFile <- read.table(file="counts_M3_MUN33.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_MUN33.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)


#FIN19
# Gene counts
countsFile <- read.table(file="counts_F9_FIN19.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_FIN19.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)


#BZ37
# Gene counts
countsFile <- read.table(file="counts_B7_BZ37.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_BZ37.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)


#DEN62
# Gene counts
countsFile <- read.table(file="counts_D6_DEN62.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_DEN62.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)



#DEN62
# Gene counts
countsFile <- read.table(file="counts_G4_GIM44.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_GIM44.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)



#MUN24
# Gene counts
countsFile <- read.table(file="counts_M4_MUN24.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_MUN24.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)



# Compare ONLY controls "Resistant" vs "Sensitive":

# Gene counts
countsFile <- read.table(file="counts_Controls.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_Controls.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Resistance)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Resistance <- relevel(dds$Resistance, "Sensitive")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)
