library(data.table)
library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(viridis)
library(ggthemes)

'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  stop("This should contain on <formatted.drugscreen> <output.name>", call.=FALSE)
}


#args <- c("Processed_data/Kat/HCI.gr50.scores.txt","Data/HCI_clinical.txt","Data/drugClasses.txt","Figures/OverallHeatmap.v1.pdf")
#args <- c("Processed_data/HCI.gr50.scores.txt","Data/HCI_clinical_Mar2021.txt","Data/drugClasses.txt","Figures/fig5a_OverallHeatmap.pdf")

dat <- fread(args[1])
dat$GR50c <- log10(ifelse(dat$GR50 == 0, 999,dat$GR50))
dd <- dcast(dat, HCIid ~ drug, value.var="GR50c")
dd$larotrectinib <- NULL
dd$ly3039478 <- NULL
dd$fulvestrant <- NULL
dd$dmso <- NULL
dd$epz011989 <- NULL
rownames(dd) <- dd$HCIid
dd$HCIid <- NULL
is.na(dd) <- sapply(dd,is.infinite)
toRemove <- c("Control","HCI-011.E2","HCI-017.E2")
dd <- dd[which(rownames(dd) %!in% toRemove),]
rnames <- rownames(dd)


dd2 <- data.frame(lapply(dd,function(x) ifelse(x >2.99,NA,x)))
dds <- data.frame(lapply(dd2, function(x) scale(x)))
dds2 <- data.frame(lapply(dds, function(x) ifelse(is.na(x),max(x,na.rm=T)+.5,x)))
row.names(dds2) <- rnames
m <- as.matrix(dds2)
row.names(m) <- rnames

#Clinical
clin <- fread(args[2])
#subset to the tested models 
sclin <- clin[which(clin$SampleID %in% rownames(m)),]
favclin <- sclin[,c("SampleID","CollectionAge","Anatomy","ERstatus","PRstatus","HER2status")]

#Drug input
drug <- fread(args[3])
colnames(drug) <- c("drug","group")
drug$group2 <- as.character(drug$group)


#Colors 
age_col = colorRamp2(c(40, 55, 93), c("yellow", "orange", "brown"))
anatomy_col = structure(brewer.pal(length(unique(sclin$Anatomy)), "Dark2"), names = unique(sclin$Anatomy))
er_col = structure(brewer.pal(12, "Set3")[c(1,2)], names = unique(sclin$ERstatus))
pr_col = structure(brewer.pal(12, "Set3")[c(3,4)], names = unique(sclin$PRstatus))
her2_col = structure(brewer.pal(12, "Set3")[c(5,6)], names = unique(sclin$HER2status))


#Row annotation 
#More annotations (Age and Anatomy)
#row_ha = rowAnnotation(Age=favclin$CollectionAge, Anatomy=favclin$Anatomy, ER=favclin$ERstatus, PR=favclin$PRstatus, HER2=favclin$HER2status, col=list(Age=age_col,Anatomy=anatomy_col,ER=er_col,PR=pr_col,HER2=her2_col))
row_ha = rowAnnotation(ER=favclin$ERstatus, PR=favclin$PRstatus, HER2=favclin$HER2status, col=list(ER=er_col,PR=pr_col,HER2=her2_col))

#GRAOC
ddaoc <- dcast(dat, HCIid ~ drug, value.var="GR_AOC")
ddaoc$dmso <- NULL
rownames(ddaoc) <- ddaoc$HCIid
rnames <- ddaoc$HCIid
ddaoc$HCIid <- NULL
is.na(ddaoc) <- sapply(ddaoc,is.infinite)
toRemove <- c("Control","HCI-011.E2","HCI-017.E2")
ddaoc <- ddaoc[which(rownames(ddaoc) %!in% toRemove),]
ddaoc[ddaoc>1.5] <- NA
ddaoc[ddaoc < -1.5] <- NA

m <- as.matrix(ddaoc)
hmi <- Heatmap(m,col=viridis(option="inferno",n=100,direction=-1),right_annotation = row_ha)
hmv <- Heatmap(m,col=viridis(option="viridis",n=100,direction=-1),right_annotation = row_ha)
hmm <- Heatmap(m,col=viridis(option="magma",n=100,direction=-1),right_annotation = row_ha)

pdf("Figures/ComplexHeatmap.v2.aoc.inferno.pdf",useDingbats=F,height=7,width=14)
print(hmi)
dev.off()

pdf("Figures/ComplexHeatmap.v2.aoc.viridis.pdf",useDingbats=F,height=7,width=14)
print(hmv)
dev.off()

pdf(args[4],useDingbats=F,height=7,width=14)
print(hmm)
dev.off()

#This is to make a GR50 stacked plot with the Possible scores 
ddgr50 <- dcast(dat, HCIid ~ drug, value.var="GR50")
ddgr50$dmso <- NULL
rownames(ddgr50) <- ddgr50$HCIid
rownames(ddgr50) <- ddgr50$HCIid
ddgr50$HCIid <- NULL

is.na(ddgr50) <- sapply(ddgr50,is.infinite)
toRemove <- c("Control","HCI-011.E2","HCI-017.E2")
ddgr50 <- ddgr50[which(rownames(ddgr50) %!in% toRemove),]
matgr50 <- as.matrix(ddgr50)
ddgr50[ddgr50>100] <- NA
ddgr50[ddgr50 < 6.25e-7] <- NA
ddgr50l <- log10(ddgr50)

rnames <- row.names(ddgr50l)

#This is going to make stacked Barplots for AOC based on dd 
#dd from above
#DRG = "birinapant"
drugs <- colnames(ddaoc)
for(DRG in drugs){
  ddp <- ddaoc
  drg = DRG
  d <- as.data.frame(ddp[,which(colnames(ddp) == drg)])
  d$HCIids <- rownames(dd)
  colnames(d) <- c('V1','HCIids')
  d$ordHCIids <- factor(d$HCIids,levels=(d[order(d$V1,decreasing=F,na.last=FALSE),]$HCIids))

  p <- ggplot(d,aes(x=drg,ordHCIids))
  p <- p+geom_tile(aes(fill=V1))
  p <- p+scale_fill_viridis(option = "inferno",direction =-1,limits=c(min(ddp),max(ddp)))
  p <- p+xlab("")+ylab("")
  p

  oname = paste("Figures/StackedAOC/",drg,".ordered.pdf",sep="")
  print(oname)
  pdf(oname,width=2,height=4)
  print(p)
  dev.off()
}



#This will produce the stacked plots with slightly different rankings if needed 
#
##DRG = 'birinapant'
#drugs <- colnames(ddgr50l)
#for(DRG in drugs){
#  d <- as.data.frame(ddgr50l[,which(colnames(ddgr50l) == DRG)])
#  drg = DRG
#  d$HCIids <- rnames
#  colnames(d) <- c('V1','HCIids')
#  d$ordHCIids <- factor(d$HCIids,levels=(d[order(d$V1,decreasing=T,na.last=FALSE),]$HCIids))
#  d$NotM <- ifelse(is.na(d$V1),.25,NA)
#
#  p <- ggplot(d,aes(x=drg,ordHCIids))
#  p <- p+geom_tile(aes(fill=V1))
#  p <- p+scale_fill_viridis(option = "viridis",direction =1,limits=c(min(ddgr50l,na.rm=T),max(ddgr50l,na.rm=T)))
#  p <- p+xlab("")+ylab("")
#  p <- p+geom_point(size=d$NotM)
#  p
#
#  oname = paste("Figures/StackedGR50/",drg,".ordered.pdf",sep="")
#  print(oname)
#  pdf(oname,width=2,height=4)
#  print(p)
#  dev.off()
#}
#
#
#
#
##For testing This is a SCALED GR50 score to capture outliers in group responses.
##DRG = "epirubucin"
#drugs <- colnames(dds2)
#for(DRG in drugs){
#  dds2p <- dds2
#  drg = DRG
#  d <- as.data.frame(dds2[,which(colnames(dds2p) == drg)])
#  d$HCIids <- rnames
#  colnames(d) <- c('V1','HCIids')
#  d$ordHCIids <- factor(d$HCIids,levels=(d[order(d$V1,decreasing=T,na.last=FALSE),]$HCIids))
#  d$NotM <- ifelse(d$V1 == max(d$V1),.25,NA)
#
#  p <- ggplot(d,aes(x=drg,ordHCIids))
#  p <- p+geom_tile(aes(fill=V1))
#  p <- p+scale_fill_viridis(option = "viridis",direction =1,limits=c(min(dds2p),max(dds2p)))
#  p <- p+xlab("")+ylab("")
#  p <- p+geom_point(size=d$NotM)
#  p
#
#  oname = paste("Figures/StackedScaled/",drg,".ordered.pdf",sep="")
#  print(oname)
#  pdf(oname,width=2,height=4)
#  print(p)
#  dev.off()
#}
#
