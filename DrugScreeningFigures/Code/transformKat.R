
library(data.table)
library(dplyr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("This should contain on <output.plot.from.rule.parse.txt> <mapping.txt> <file.pdf>", call.=FALSE)
}

#args = c("Processed_data/Kat/E2plus.longformat.txt", "Processed_data/Kat/E2plus.reformatted.txt","Processed_data/Kat/E2plus.gr50.reformatted.txt")


dat <- fread(args[1],sep="\t",header=T)

rdat <- dcast(dat, Plate + Drug + Vol + Normalized + DrugSet ~ Variable, value.var = 'Value')

write.table(rdat,args[2], sep="\t", quote=F, row.names=F)

dat2 <- dat[which(dat$Normalized %in% c("Raw","RawDay0","Vehicle")),]
rdat2 <- dcast(dat2, Plate + Drug + Vol + DrugSet + Variable ~ Normalized, value.var = 'Value')

write.table(rdat2,args[3], sep="\t", quote=F, row.names=F)

