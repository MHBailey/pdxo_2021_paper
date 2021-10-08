library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(plotrix)


'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("This should contain 9 files See Snakefile for more details")
}

#args <- c('Processed_data/HCI.gr50.scores.txt','Data/HCI_clinical_Mar2021.txt','Figures/fig5b_Hormone.mTOR.TNBC.v1.pdf')

#Process the GR50 data
gr50 = fread(args[1])
dd <- dcast(gr50, HCIid ~ drug, value.var="GR_AOC")
dd50 <- dcast(gr50, HCIid ~ drug, value.var="GR50")
ddmax <- dcast(gr50, HCIid ~ drug, value.var="GRmax")
dd$dmso <- NULL
rownames(dd) <- dd$HCIid
rnames <- dd$HCIid
dd$HCIid <- NULL
is.na(dd) <- sapply(dd,is.infinite)
toRemove <- c("Control","HCI-011.E2","HCI-017.E2")
dd <- dd[which(rownames(dd) %!in% toRemove),]
rownames(dd) <- str_replace(rownames(dd),"-","")
screened <- rownames(dd)
#This is the file to move forward with
grScores = melt(as.matrix(dd,rownames=rownames(dd)))

#Process the Alana-Clinical file for ER integation
alana <- fread(args[2],sep="\t")
alana$HCIid <- gsub(pattern="-",replacement="",x=alana$"Model ID")


#Bring the data together 
grpamclin <- merge(grScores,alana,by.x="Var1",by.y="HCIid",all.x=T)

chemos = c("sn-38", "docetaxel", "paclitaxel", "gemcitabine", "doxorubicin", "epirubucin", "carboplatin", "romidepsin")

mtor = c("ink 128", "everolimus", "rapamycin", "tak-228", "sapanisertib","apitilosib", "taselisib", "azd5363", "bkm120")


grpamclinr <- data.frame(grpamclin %>% group_by(Var2) %>% mutate(rank = rank(value, ties.method = "first")))

#Now after the call today (Sep 2, 2020) we are moving to a slightly different model 
#Rank by sample for drug class 
tnbc <- c("HCI001", "HCI002", "HCI010", "HCI015", "HCI016", "HCI019", "HCI023", "HCI024", "HCI025", "HCI027")
grpamclinrs <- grpamclinr %>% group_by(Var1) %>% mutate(rankXsample = rank(value, ties.method = "first"))
grpamclinrs$DrugClass <- ifelse(grpamclinrs$Var2 %in% mtor, "mTOR","Other")
grpamclinrs$DrugClass <- ifelse(grpamclinrs$Var2 %in% chemos, "Chemo",grpamclinrs$DrugClass)
grpamclinrs$TNBC <- ifelse(grpamclinrs$Var1 %in% tnbc,"TNBC","HR-positive\n& HER2") 


p <- ggplot(grpamclinrs,aes(x=Var1,y=rankXsample,fill=DrugClass))
p <- p + geom_tile()
#p <- p + geom_text(aes(label=Var2),size=2)
p <- p + theme_minimal()
p <- p + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5),legend.position="bottom")
p <- p + scale_fill_brewer(palette="Set2")
p <- p + ggtitle("Across all drugs")
p <- p + xlab("") + ylab("Ordered GRaoc")
p <- p + facet_grid(.~ TNBC,scales="free_x",space="free_x")
p

pdf("Figures/fig5b_Hormone.mTOR.TNBC.v1.pdf",height=4.5,width=2,useDingbats=F)
print(p)
dev.off()

#MTOR
#Now for the stats of this Hormone v non-hormone 
hrp_mTOR <- grpamclinrs[which(grpamclinrs$TNBC != "TNBC" & grpamclinrs$DrugClass == "mTOR"),]$rankXsample
tnbc_mTOR <-  grpamclinrs[which(grpamclinrs$TNBC == "TNBC" & grpamclinrs$DrugClass == "mTOR"),]$rankXsample
#BY sample Rank
wilcox.test(hrp_mTOR,tnbc_mTOR,conf.int=T)

#CHEMO
hrp_aoc_chemo <- grpamclinrs[which(grpamclinrs$TNBC != "TNBC" & grpamclinrs$DrugClass == "Chemo"),]$value
tnbc_aoc_chemo <-  grpamclinrs[which(grpamclinrs$TNBC == "TNBC" & grpamclinrs$DrugClass == "Chemo"),]$value
wilcox.test(tnbc_aoc_chemo,hrp_aoc_chemo,conf.int=T)


