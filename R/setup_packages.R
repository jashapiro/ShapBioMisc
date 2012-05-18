#packages to install when setting up a new R installation

#Hadley's essentials
hadley <- c("reshape","stringr", "memoise", "plyr", "ggplot2")
install.packages(hadley)

#parallel computation
parallel_packs <- c("foreach", "doParallel")
install.packages(parallel_packs)

#reproducible research
install.packages("knitr")

#Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("ShortRead", "GO.db", "GOstats"))

#Saccharomyces cerevisiae annotation
biocLite("org.Sc.sgd.db")

#Caenorhabditis elegans annotation
biocLite("org.Ce.eg.db")

#Streptomyces coelicolor annotation
biocLite("org.Sco.eg.db")