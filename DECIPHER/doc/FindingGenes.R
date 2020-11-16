### R code from vignette source 'FindingGenes.Rnw'

###################################################
### code chunk number 1: FindingGenes.Rnw:47-52
###################################################
options(continue=" ")
options(width=80)
options(SweaveHooks=list(fig=function()
	par(mar=c(4.1, 4.1, 0.3, 0.1))))
set.seed(123)


###################################################
### code chunk number 2: startup
###################################################
library(DECIPHER)


###################################################
### code chunk number 3: expr1
###################################################
# specify the path to your genome:
genome_path <- "<<path to genome FASTA file>>"
# OR use the example genome:
genome_path <- system.file("extdata",
	"Chlamydia_trachomatis_NC_000117.fas.gz",
	package="DECIPHER")

# read the sequences into memory
genome <- readDNAStringSet(genome_path)
genome


###################################################
### code chunk number 4: expr2
###################################################
getOption("SweaveHooks")[["fig"]]()
orfs <- FindGenes(genome, showPlot=TRUE, allScores=TRUE)


###################################################
### code chunk number 5: expr3
###################################################
orfs


###################################################
### code chunk number 6: expr4
###################################################
genes <- orfs[orfs[, "Gene"]==1,]


###################################################
### code chunk number 7: expr5
###################################################
colnames(genes)


###################################################
### code chunk number 8: expr6
###################################################
dna <- ExtractGenes(genes, genome)
dna


###################################################
### code chunk number 9: expr7
###################################################
table(subseq(dna[-1], 1, 3))


###################################################
### code chunk number 10: expr8
###################################################
w <- which(!subseq(dna, 1, 3) %in% c("ATG", "GTG", "TTG"))
w
w <- w[-1] # drop the first sequence because it is a fragment
w
dna[w]
genes[w,]


###################################################
### code chunk number 11: expr9
###################################################
aa <- ExtractGenes(genes, genome, type="AAStringSet")
aa


###################################################
### code chunk number 12: expr10
###################################################
subseq(aa, 2, -2)


###################################################
### code chunk number 13: expr11
###################################################
getOption("SweaveHooks")[["fig"]]()
pairs(genes[, 5:16], pch=46, col="#00000033", panel=panel.smooth)


###################################################
### code chunk number 14: expr12
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(orfs, interact=FALSE)


###################################################
### code chunk number 15: sessinfo
###################################################
toLatex(sessionInfo(), locale=FALSE)


