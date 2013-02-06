library(SRAdata)

args <- commandArgs(trailingOnly=TRUE)
datadir <- args[1]
print(datadir)

tab <- SRATable(datadir)
show(tab)

tab@info <- SRAGetTableInfo(tab)
show(tab)

dat <- SRATableReadFastq(tab, lane=7, tile=69, max.nreads=10)
show(dat)
show(sread(dat))
show(id(dat))
show(quality(dat))

int <- SRATableReadIntensities(tab, lane=7, tile=69, max.nreads=10)
show(int)
str(intensity(int))
intensity(int)[[1,,]]


