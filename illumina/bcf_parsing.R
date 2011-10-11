library(Rsamtools)

chr_names<- c("I", "II", "III", "IV", "V", "X", "MtDNA")
chr_lengths <- c(15072421, 15279323, 13783685, 17493784, 20924143, 17718854, 13794)
names(chr_lengths) <- chr_names
elegans_gr <-GRanges(seqnames = Rle(chr_names, rep(1, 7)),
            ranges = IRanges(rep(1,7), width = chr_lengths), 
            seqlengths = chr_lengths)


ecaRAD <- BcfFile('/kruglyak/shared/illumina/data/genomic/eca_22lib_variants.bcf')
bcf_header <- scanBcfHeader(ecaRAD)
param <- ScanBcfParam(which=elegans_gr[1])
bcf <- scanBcf(ecaRAD, param = param)

names(bcf)
