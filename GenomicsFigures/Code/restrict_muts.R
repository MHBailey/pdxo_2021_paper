library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(viridis)

#this is source activate stats 
uniqString <- function(x){
    yo <- unique(unlist(strsplit(x, ";")))
    return(paste(yo,collapse=";"))
}
'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=7) {
  stop("Check out the rule restict mutations", call.=FALSE)
}

 
#args <- c("Data/MAFs/Indepth_2020_04_22.all.HG38.PASS.maf","Data/GeneLists/299.BRCA.glist.txt","Data/GeneLists/ChasmPlus.collin.pancan.cancer.txt","Data/GeneLists/Mutation.CTAT.3D.Scores.processed.txt","Data/Meta/indepth.dna.txt","Processed_data/Indepth_2020_04_22.all.HG38.PASS.indels.cancer_miss.maf","No")

#args <- c("Data/MAFs/Indepth_2020_04_22.all.HG38.PASS.maf","Data/GeneLists/299.cancer.glist.txt","Data/GeneLists/ChasmPlus.collin.pancan.cancer.txt","Data/GeneLists/Mutation.CTAT.3D.Scores.processed.txt","Data/Meta/indepth.dna.txt","Processed_data/Indepth_2020_04_22.all.HG38.PASS.indels.cancer_miss.maf","No")

#args <- c("Data/MAFs/TumorOnly_2020_07_31.all.HG38.maf","Data/GeneLists/299.cancer.glist.txt","Data/GeneLists/ChasmPlus.collin.pancan.cancer.txt","Data/GeneLists/Mutation.CTAT.3D.Scores.processed.txt","Data/Meta/to.dna.txt","Processed_data/Indepth_2020_04_31.all.HG38.PASS.indels.cancer_miss.maf","No")

#THIS IS A SPECIAL ONE-OFF CASE FOR CW-01
#args <- c("Data/MAFs/CW01_pdo_P4_263_WXS.maf","Data/GeneLists/299.cancer.glist.txt","Data/GeneLists/ChasmPlus.collin.pancan.cancer.txt","Data/GeneLists/Mutation.CTAT.3D.Scores.processed.txt","Data/Meta/to.dna.txt","Processed_data/Indepth_2020_04_31.all.HG38.PASS.indels.cancer_miss.maf","No")

#This is to do another one off for ALL CW-01 ascites++
#args <- c("Data/MAFs/CW01all.maf","Data/GeneLists/299.cancer.glist.txt","Data/GeneLists/ChasmPlus.collin.pancan.cancer.txt","Data/GeneLists/Mutation.CTAT.3D.Scores.processed.txt","Data/Meta/to.dna.txt","Processed_data/Indepth_2020_04_31.all.HG38.PASS.indels.cancer_miss.maf","Yes")


dat <- fread(args[1])
g299 <- fread(args[2],header=F)
chasmp <- fread(args[3])
pancanM <- fread(args[4])
meta <- fread(args[5],header=F)
istumoronly = args[7]

if(istumoronly == "Yes"){
    dat1 <- dat[which(dat$FILTER == "PASS" | dat$FILTER == "germline_risk"),]
    dat2 <- dat1[which(dat1$Variant_Classification != "Silent"),]   
    dat <- dat2
}

if(istumoronly != "Yes"){
    dat1 <- dat[which(dat$FILTER == "PASS"),]
    dat2 <- dat1[which(dat1$Variant_Classification != "Silent"),]
    dat <- dat2
}

#Add column GENE:HGVSp_Short 
dat$Mut = paste(dat$Hugo_Symbol,dat$HGVSp_Short,sep=":")
chasmp$Mut =  paste(chasmp$"HUGO Symbol",chasmp$mutation,sep=":")
pancanM$Mut = paste(pancanM$gene,pancanM$protein_change,sep=":")

#Remove any samples that are going to move forward in the analysis 
maf <- dat[which(dat$Tumor_Sample_Barcode %in% meta$V2),]
#For the one off scripts 
#maf <- dat

#Remove non-coding mutations 
noncoding <- c("3'UTR","5'UTR","3'Flank","5'Flank","Intron","Splice_Region")
maf2 <- maf[which(maf$Variant_Classification %!in% noncoding),]

#Split maf based on the Missense mutations and otheri
kuans = c("APC","AR","ATM","ATR","BARD1","BLM","BRCA1","BRCA2","BRIP1","BUB1B","CDKN2A","CHEK2","COL7A1","DOCK8","EPCAM","FANCA","FANCC","FANCM","GJB2","MAX","MUTYH","NF1","PALB2","PMS2","POLE","POT1","PRDM9","PTCH1","PTPN11","RAD51B","RAD51C","RECQL","RET","SETBP1","TP53","UROD","ESR1","FANCD2","FANCE","FANCF","FANCG","HMBS","POLD1","FANCI","FANCL","RAD51B")

missmaf <- maf2[which(maf2$Variant_Classification == "Missense_Mutation"),]
missmaf2 <- missmaf[which( missmaf$Mut %in% chasmp$Mut | missmaf$Mut %in% pancanM$Mut),]
missmaf3 <- missmaf[which( missmaf$Hugo_Symbol %in% kuans),] #Some of Kuan's genes 
missmaf4 <- missmaf3[which( grepl("deleterious",missmaf3$SIFT) | grepl("damaging",missmaf3$PolyPhen) | grepl("pathogenic",missmaf3$CLIN_SIG) | grepl("HIGH",missmaf3$IMPACT) ), ]
missmaf5 <- missmaf4[which(is.na(missmaf4$gnomAD_AF)),]
missmafAll <- rbind(missmaf2,missmaf5)


#Append the indels/splice/nonsense mutations in cancer genes 
indelmaf <- maf2[which(maf2$Variant_Classification != "Missense_Mutation"),]
indelmaf2 <- indelmaf[which(indelmaf$Hugo_Symbol %in% g299$V1 | indelmaf$Hugo_Symbol %in% kuans),]
#find the peak density of the GNOMAD 
na.omit(indelmaf$gnomAD_AF)
wind <- which.max(density(na.omit(indelmaf$gnomAD_AF))$y) #which indel
threshold <- density(na.omit(indelmaf$gnomAD_AF))$x[wind] #calc GNOMAD Threshold
#indelmaf3 <- indelmaf2[which(indelmaf2$gnomAD_AF < threshold | is.na(indelmaf2$gnomAD_AF)),]
indelmaf3 <- indelmaf2[which(indelmaf2$gnomAD_AF < 0.001 | is.na(indelmaf2$gnomAD_AF)),] #https://www.nature.com/articles/s41431-018-0169-4/metrics

#VAF remove mutions less than 1%
catmaf <- rbind(missmafAll,indelmaf3)
catmaf$VAF = catmaf$t_alt_count/catmaf$t_depth
catmaf10 <- catmaf[which(catmaf$VAF >= 0.1),]

#This is a warning to say that a filtering strategy completely removed a sample
allsamps <- sort(unique(maf$Tumor_Sample_Barcode))
filtsamps <- sort(unique(catmaf10$Tumor_Sample_Barcode))

if (length(setdiff(allsamps,filtsamps))>0) {
  print("Check out the rule restict mutations, you lost a sample")
}
#If you need this is some code to add NA's if necessary 
#Add empty column https://stackoverflow.com/questions/18214395/add-empty-columns-to-a-dataframe-with-specified-names-from-a-vector
#addbackdna <- setdiff(allsamps,filtsamps)
#glmaf3[,addbackdna] <- NA

dat <- fread(args[1])


#Write this file as is, and then go from there
write.table(catmaf10,args[6],sep="\t",quote=F,row.names=F)



