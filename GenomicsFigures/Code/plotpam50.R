#Updated Jun21, for the manuscript rebuttal

library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(viridis)
library(scales)
library(ComplexHeatmap)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
  stop("Check out the rule PAM50_exprression_profile, not enough files", call.=FALSE)
}

#TUMOR ONLY 
#args <- c('Processed_data/normalized.DESEQ.RSEM.batch','Data/GeneLists/PAM.genes.txt','Data/Meta/to.rna.txt','Data/GeneLists/Collated_Pam50_PDFs.20190731.txt','Data/Clinical/HCI_clinical.txt','Figures/pam50.Counts.pdf')

#LDEV
#args <- c("Processed_data/normalized.DESEQ.RSEM.batch","Data/GeneLists/PAM.genes.txt","Data/Meta/ldev.rna.txt","Data/GeneLists/Collated_Pam50_PDFs.20190731.txt","Data/Clinical/HCI_clinical_Dec2020.txt","Figures/pam50.Counts.LDEV.pdf")

#LAR
#args <- c("Processed_data/normalized.DESEQ.RSEM.batch","Data/GeneLists/LAR.genes.txt","Data/Meta/lar.rna.txt","Data/GeneLists/Collated_Pam50_PDFs.20190731.txt","Data/Clinical/HCI_clinical.txt","Figures/pam50.Counts.AR.pdf")

#REBUTALL 
#args <- c("Processed_data/Rebuttal/normalized.DESEQ.RSEM.txt","Data/GeneLists/PAM.genes.txt","Data/Meta/indepth.rna.20210120.txt","Data/GeneLists/Collated_Pam50_PDFs.20190731.txt","Data/Clinical/HCI_clinical_Mar2021.txt","Figures/pam50.Counts.Indepth.20210618.pdf")

dat <- fread(args[1])
pam50gl <- fread(args[2],header=F)
colnames(pam50gl) <- c("Gene","Classifier")
pam50glord <- pam50gl[order(pam50gl$Classifier),]
rnato <- fread(args[3],header=F)
pam50class <- fread(args[4])

#Process the pam50classifications
pam50class$HCIid <- toupper(str_split_fixed(pam50class$SampleName,",",2)[,1])
pam50class$Model1 <- str_split_fixed(pam50class$SampleName,",",2)[,2] 
pdxp50class <- pam50class[which(grepl("PDX",pam50class$Model1)),]
rnato$HCIid <- str_split_fixed(rnato$V2,"_",5)[,1]
#https://stackoverflow.com/questions/13863599/insert-a-character-at-a-specific-location-in-a-string
rnato$HCIid2 <- gsub('^([A-Z]{3})(.+)', '\\1-\\2',x=rnato$HCIid)
#Subset the pam50c Data to the RNAdata 
p50c <- pdxp50class[which(pdxp50class$HCIid %in% rnato$HCIid2 ),c("HCIid","Subtype")]
p50co <- unique(p50c[order(p50c$HCIid),])

#Collapse duplicates with a comma! 
p50coc <- p50co[, toString(Subtype), by = list(HCIid)]
pam50_col = structure(brewer.pal(9, "Set1")[c(1:(length(unique(p50coc$V1))+1))], names = c(unique(p50coc$V1),"NA"))
rnato2 <- merge(rnato,p50coc,by.x="HCIid2",by.y="HCIid",all.x=T)
ha = HeatmapAnnotation(PAM50=rnato2$V1.y,col=list(PAM50=pam50_col))


#Subset expression to 50 genes
datpam <- dat[which(dat$V1 %in% pam50gl$Gene),]
rownames(datpam) <- datpam$V1
datpam$V1 <- NULL
datpam2 <- datpam %>% select(rnato$V2) 
#turn it into a matrix
dpm <- as.matrix(datpam2,rownames=rownames(datpam2))
#This is where I need to log and then scale dpm 
dpml <- log10(dpm+1)
dpms <- t(apply(dpml,1,scale))
colnames(dpms) <- colnames(dpml)
ra = rowAnnotation(foo = anno_text(rownames(dpms), gp = gpar(fontsize = 2)))


#order the rows how I want them 
#https://stackoverflow.com/questions/36410485/how-to-order-a-data-frame-based-on-row-names-in-another-data-frame
orddpm <-  dpms[match(pam50glord$Gene, rownames(dpms)), ]
colnames(orddpm) <- colnames(dpms)
split <- pam50glord$Classifier
colnames(orddpm) <- str_split_fixed(colnames(orddpm),"_",5)[,1]
colnames(dpms) <- str_split_fixed(colnames(dpms),"_",5)[,1]


#Add the clinical data 
clin <- fread(args[5])
clin$HCIid <- str_replace(clin$SampleID,"-","")
tumsub <- clin[,c("HCIid","Anatomy","TripleNeg","Type","BMI","smokinghistory","ERstatus","PRstatus","HER2status")]

HCIsamples <- sort(unique(colnames(orddpm)))
HCIsamples[HCIsamples == 'TOW18S'] <- "HCI041"
HCIsamples[HCIsamples == 'TOW19'] <- "HCI043"
subclin <- clin[which(clin$HCIid %in% HCIsamples),]

HCIsamples2 <- data.frame("HCIid"=colnames(orddpm))
subclin2 <- plyr::join(rnato,clin,by="HCIid") 


#Make a standard set of colors for heatmap 
er_col = structure(c("white","black","grey"),names=c(unique(clin$ERstatus),"NA"))
pr_col = structure(c("white","black","grey"),names=c(unique(clin$PRstatus),"NA"))
her2_col = structure(c("white","black","grey"),names=c(unique(clin$HER2status),"NA"))
type_col = structure(magma((length(unique(subclin$Type))+1))[2:(length(unique(subclin$Type))+1)],names=c(unique(subclin$Type)))
pam50_col = structure(brewer.pal(9, "Set1")[c(1:(length(unique(p50coc$V1))+1))], names = c(unique(p50coc$V1),"NA")
)


clinanno = HeatmapAnnotation(
        ER=subclin$ERstatus,
        PR=subclin$PRstatus,
        HER2=subclin$HER2status,
        Pathology=subclin$Type,
#        PAM50=c(unique(rnato2$V1.y),NA),
#        col = list( ER=er_col, PR=pr_col,HER2=her2_col,Pathology=type_col,PAM50=pam50_col)
        col = list( ER=er_col, PR=pr_col,HER2=her2_col,Pathology=type_col)
        )


clinanno2 <- HeatmapAnnotation(
    Model=rnato$V1,
    ER=subclin2$ERstatus,
    PR=subclin2$PRstatus,
    HER2=subclin2$HER2status,
#    Pathology=subclin2$Type,
#     PAM50=c(unique(rnato2$V1.y),NA),
#     col = list( ER=er_col, PR=pr_col,HER2=her2_col,Pathology=type_col,PAM50=pam50_col)
    col = list( ER=er_col, PR=pr_col,HER2=her2_col,Pathology=type_col)
) 

#I need to make sure that I re-order the genes based off of the orders 


#Make the heatmap
hm <-  Heatmap(orddpm,
    bottom_annotation = clinanno2, 
    row_order=rownames(orddpm),
    row_names_gp = gpar(fontsize = 4),
    row_split=split,
    clustering_distance_columns = "pearson"
)


#Heatmap(orddpm, bottom_annotation = ha, row_order=rownames(orddpm) ,row_names_gp = gpar(fontsize = 4),row_split=split)



pdf(args[6],height=6,width=9)
print(hm)
dev.off()

#This is done for the LDEV samples. 
dpml <- log10(dpm+0.00001)
orddpm <-  dpml[match(pam50glord$Gene, rownames(dpml)), ]
hm <-  Heatmap(orddpm,
    #bottom_annotation = clinanno, 
    row_order=rownames(orddpm),
    row_names_gp = gpar(fontsize = 4),
    row_split=split,
    clustering_distance_columns = "pearson"
)


##Lets do some analysis of the LDEV
#normal = c("KRT17","KRT5","SFRP1","BCL2","KRT14","MLPH","MDM2","FGFR4","MYC")
#pam50glord$Other <- ifelse(pam50glord$Gene %in% normal, "Normal", pam50glord$Classifier)
#
#look10 <- abs(orddpm[,1]-orddpm[,2])
#sort(look10)
#
#look13 <- abs(orddpm[,3]-orddpm[,4])
#sort(look13)
#
#look13EI <- abs(orddpm[,5]-orddpm[,6])
#sort(look13EI)
#
#
##Custom enrichment factor HCI-010
#N = dim(pam50glord)[1] #Total number of genes in list 
#k = dim(pam50glord[which(pam50glord$Other == "Normal"),])[1] #The genes in the pathway of interst
#kgenes = pam50glord[which(pam50glord$Other == "Normal"),]$Gene #The genes in the pathway of interst
##NOTE: Options here: 
##M = length(look10[which(look10 > 1)]) #Greater than 1
##mgenes = names(look10[which(look10 > 1)]) #Greater than 1
#M = 10
#mgenes = names(tail(sort(look10),10)) #TOP 10 are from Normal genes 
#n = length(intersect(kgenes,mgenes))
##Calc E #enrichment 
#ef = (n*N)/(k*M)
##Hypergeometric Distributions
#(choose(k,n) * choose(N-k,M-n)) / choose(N,M) #probability
#
##Custom enrichment factor HCI-013
#N = dim(pam50glord)[1] #Total number of genes in list 
#k = dim(pam50glord[which(pam50glord$Other == "Normal"),])[1] #The genes in the pathway of interst
#kgenes = pam50glord[which(pam50glord$Other == "Normal"),]$Gene #The genes in the pathway of interst
##NOTE: Options here: 
##M = length(look13[which(look13 > 1)]) #or look other
##mgenes = names(look13[which(look13 > 1)]) #Greater than 1
#M = 10
#mgenes = names(tail(sort(look13),10))
#n = length(intersect(kgenes,mgenes))
##Calc E #enrichment 
#ef = (n*N)/(k*M)
##Hypergeometric Distributions
#(choose(k,n) * choose(N-k,M-n)) / choose(N,M) #probability
#
##Custom enrichment factor HCI-013EI
#N = dim(pam50glord)[1] #Total number of genes in list 
#k = dim(pam50glord[which(pam50glord$Other == "Normal"),])[1] #The genes in the pathway of interst
#kgenes = pam50glord[which(pam50glord$Other == "Normal"),]$Gene #The genes in the pathway of interst
##NOTE: Options here: 
##M = length(look13EI[which(look13EI > 1)]) #or look other
##mgenes = names(look13EI[which(look13EI > 1)]) #Greater than 1
#M = 10
#mgenes = names(tail(sort(look13EI),10))
#n = length(intersect(kgenes,mgenes))
##Calc E #enrichment 
#ef = (n*N)/(k*M)
##Hypergeometric Distributions
#(choose(k,n) * choose(N-k,M-n)) / choose(N,M) #probability
#
#
