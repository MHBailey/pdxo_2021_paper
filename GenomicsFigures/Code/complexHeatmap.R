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

stde <- function(x){
    sd(x)/sqrt(length(x))
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=8) {
  stop("Check out the rule complexHeatmap", call.=FALSE)
}


#Testing DNA
#args <- c("Data/MAFs/Indepth_2020_04_22.all.HG38.PASS.maf","Data/Meta/cat_DNA_Indepth.txt","Data/GeneLists/299.BRCA.glist.txt","Data/CNV/cnv_combined.txt","Data/Meta/cat_CNV_old.txt","Processed_data/normalized.batch","Data/Meta/cat_RNA_Indepth.txt","Data/RNAseq/metadata_20200313.txt","Figures/pdxo.genomicsFig.pdf")

#BYU Meta
#args <- c("Data/MAFs/Indepth_2020_04_22.all.HG38.PASS.maf","Data/Meta/BYU.dna.txt","Data/GeneLists/299.BRCA.glist.txt","Data/CNV/cnv_combined.txt","Data/Meta/BYU.cnv.txt","Processed_data/normalized.batch","Data/Meta/BYU.rna.txt","Data/RNAseq/metadata_20200313.txt","Figures/pdxo.genomicsFig.pdf")

#Indepth
#args <- c("Processed_data/Indepth_2020_04_22.all.HG38.PASS.indels.cancer_miss.maf","Data/Meta/indepth.dna.txt","Data/CNV/combined.xingyi.20200724.txt","Data/Meta/indepth.cnv.txt","Processed_data/normalized.TPM.batch","Data/Meta/indepth.rna.txt","Figures/pdxo.genomicsFig.pdf")

#Indepth TO
#args <- c("Processed_data/Indepth_2020_04_22.all.HG38.PASS.indels.cancer_miss.maf","Data/Meta/indepth.dna.txt","Data/CNV/combined.xingyi.20200724.txt","Data/Meta/indepth.cnv.txt","Processed_data/normalized.TPM.batch","Data/Meta/indepth.rna.txt","Figures/pdxo.genomicsFig.pdf")

#Indepth TO
#args <- c("Processed_data/Indepth_2020_09_16.all.HG38.PASS.indels.cancer_miss.maf","Data/Meta/indepth.dna.txt","Data/CNV/combined.xingyi.20200724.txt","Data/Meta/indepth.cnv.txt","Processed_data/normalized.TPM.batch","Data/Meta/indepth.rna.txt","Figures/pdxo.genomicsFig.pdf")

#Indepth TO - Different Normalization strategies and TOP DE genes
#args <- c("Processed_data/Indepth_2020_09_16.all.HG38.PASS.indels.cancer_miss.maf","Data/Meta/indepth.dna.txt","Data/CNV/combined.xingyi.20200724.txt","Data/Meta/indepth.cnv.txt","Processed_data/normalized.DESEQ.RSEM.batch","Data/Meta/indepth.rna.txt","Processed_data/top1000.DEgenes.batched.txt","Figures/pdxo.genomicsFig.pdf")

#With more ER+ models
#args <- c("Processed_data/Indepth_2021_01_20.all.HG38.PASS.indels.cancer_miss_mr.maf","Data/Meta/indepth.dna.20210120.txt","Data/CNV/combined.xingyi.20210118.txt","Data/Meta/indepth.cnv.20210120.txt","Processed_data/normalized.DESEQ.RSEM.batch","Data/Meta/indepth.rna.20210120.txt","Processed_data/top1000.DEgenes.batched.txt","Figures/pdxo.genomicsFig.20210120.pdf")

#Getting numbers for the BioRxiv Manuscript 
#args <- c("Processed_data/Indepth_2021_01_20.all.HG38.PASS.indels.cancer_miss_mr.maf","Data/Meta/indepth.dna.20210120.txt","Data/CNV/combined.xingyi.20210118.txt","Data/Meta/indepth.cnv.20210120.txt","Processed_data/normalized.DESEQ.RSEM.batch","Data/Meta/indepth.rna.20210120.txt","Processed_data/top1000.DEgenes.batched.txt","Figures/pdxo.genomicsFig.20210222.pdf")

#Cleaning up the MAF a bit, removing noncoding mutations 
maf <- fread(args[1])
maf$HCIid <- str_split_fixed(maf$Tumor_Sample_Barcode,"_",5)[,1]
mutspsample <- data.frame(maf %>% group_by(Tumor_Sample_Barcode) %>% tally())

#READ META
mafmeta <- fread(args[2],header = F)
colnames(mafmeta) <- c("Category","Tumor_Sample_Barcode")
mafmeta$hawkid <- str_split_fixed(mafmeta$Tumor_Sample_Barcode,"_",5)[,1]
mafmeta$meta_elid <- paste(mafmeta$hawkid,mafmeta$Category,sep="_")

#MERGE META INTO MAF
catmaf <- merge(maf,mafmeta,by="Tumor_Sample_Barcode")
catmaf$elid <- paste(catmaf$HCIid,catmaf$Category,sep="_") 
allsamps <- unique(catmaf$elid)


glmaf2 <- as.data.table(catmaf)[, toString(Variant_Classification), by = list(elid, Hugo_Symbol)]

glmaf2$V1 <- str_replace_all(glmaf2$V1,", ",";")
glmaf3 <- dcast(glmaf2, Hugo_Symbol~elid,value.var="V1")
rownames(glmaf3) <- glmaf3$Hugo_Symbol
glmaf3$Hugo_Symbol = NULL

#ORDER COLUMS and ROWS
glmaf4 <- glmaf3[mafmeta$meta_elid]

alldna <- colnames(glmaf4)
dmat <- as.matrix(glmaf4,rownames=rownames(glmaf4))

#NOW ADD CNV
cnv <- fread(args[3])
#cnv <- fread("Processed_data/normalized.CNV.batch")
cnvmeta <- fread(args[4],header=F)
colnames(cnvmeta) <- c("Category","Tumor_Sample_Barcode")
cnv$Gloc <- paste(cnv$"Approved symbol",cnv$Chromosome,sep=":")
rownames(cnv) <- cnv$Gloc
cnv2 <- subset(cnv,select=cnvmeta$Tumor_Sample_Barcode)
rownames(cnv2) <- rownames(cnv)


cnvmat <- as.matrix(cnv2,rownames=rownames(cnv2))
cnvmatna <- na.omit(cnvmat) #look into which genes were lost 

corcnv <- cor(cnvmat, use = "pairwise.complete.obs", method = "spearman")
corcnvna <- cor(na.omit(cnvmat),method="spearman") #There in on average .1% difference between these two corr plots 
Heatmap(corcnv,row_order=rownames(corcnv),column_order=colnames(corcnv),show_column_names=F,row_names_gp = gpar(fontsize = 6),name = "Spearman CNV",col=plasma(100))



#NOW ADD RNASEQ 
#testing out a couple or RNAseqs 
rna <- fread(args[5])
#rna <- fread("Processed_data/normalized.DESEQ.RSEM.batch")
#rna <- fread("Processed_data/normalized.DESEQ.RSEM.txt")
#rna <- fread("Processed_data/combined.Ecnt.welmRNA.wide.cleaned.txt")
#rna <- fread("Processed_data/normalized.TPM.batch")
rnameta <- fread(args[6],header=F)
toprna <- fread(args[7])
topgenes <- head(toprna,100)$V1
colnames(rnameta) <- c("eid","hawkid")
rnameta$HCIid <- str_split_fixed(rnameta$hawkid,"_",2)[,1]
rnameta$meta_elid <- paste(rnameta$HCIid,rnameta$eid,sep="_")

#Special consideration for 
r <- data.frame(rna %>% select(rnameta$hawkid))
rownames(r) <- rna$V1
rg <- r[which(rownames(r) %in% topgenes),] 
rg <- r 
lr <- log10(rg+1)
lrs <- t(apply(lr,1,scale))
colnames(lrs) <- colnames(lr)
lrf  <- lrs[which(is.finite(rowSums(lrs))),]
#Turn into a matrix
mr <- as.matrix(rg)
mrl <- as.matrix(lr)
mrlf <- as.matrix(lrf)
#Capture the most variable regions 
mrsd <- apply(mr,1,sd)
mrlsd <- apply(mrl,1,sd)
#Justification 

hist(mrlsd,200)
abline(v=.05)
abline(v=.1,col="blue")

sdgenes <- names(mrlsd[which(mrlsd > 0.1)])

sdmrl <- lr[which(rownames(lr) %in% sdgenes),]
sdmrlm <- as.matrix(sdmrl)

corrna <- cor(sdmrlm,method="spearman",use = "pairwise.complete.obs")

colnames(corrna) <- rnameta$meta_elid
f1 = colorRamp2(c(0.8,.87,.875,.88,.885,.89,.895,.9,.905,.91,.915,.92,.925,.93,.935,.94,1), colorRampPalette(c("white", "red"))(17))
Heatmap(corrna,cluster_rows = FALSE,cluster_columns = FALSE,col=f1)




#COLORS FOR ONCOPRINT


rcb <- brewer.pal(name="Set1",n=9)
ccc <- brewer.pal(name="Set2",n=3)

rcb[9] <- "#DCDCDC"
col = c(Missense_Mutation=rcb[1], Translation_Start_Site=rcb[2], Nonsense_Mutation=rcb[2], Nonstop_Mutation=rcb[3], Silent=rcb[3], Splice_Site=rcb[4], Splice_Region=rcb[4], In_Frame_Ins=rcb[5],In_Frame_Del=rcb[5], Frame_Shift_Ins=rcb[7], Frame_Shift_Del=rcb[7],"3'UTR"=rcb[8],"5'UTR"=rcb[8],"3'Flank"=rcb[8],"5'Flank"=rcb[8],"Intron"=rcb[8],"background"=rcb[9],Amplification=ccc[1],Deletion=ccc[2],Neutral=ccc[3])



alter_fun = list(
    In_Frame_Del = function(x, y, w, h) {
        grid.rect(x, y, w*0.90, h*0.33, gp = gpar(fill = col['In_Frame_Del'], col = NA))
    },
    In_Frame_Ins = function(x, y, w, h) {
        grid.rect(x, y, w*0.33, h*0.95, gp = gpar(fill = col['In_Frame_Ins'], col = NA))
    },
    Frame_Shift_Ins = function(x, y, w, h) {
        grid.rect(x, y, w*0.33, h*0.95, gp = gpar(fill = col['Frame_Shift_Ins'], col = NA))
    },
    Frame_Shift_Del = function(x, y, w, h) {
        grid.rect(x, y, w*0.90, h*0.33, gp = gpar(fill = col['Frame_Shift_Del'], col = NA))
    },
    Missense_Mutation = function(x, y, w, h) {
        #grid.rect(x, y, w*0.95, h*0.33, gp = gpar(fill = col['Missense_Mutation'], col = NA))
        grid.points(x, y, pch = 1, gp = gpar(fill = col['Missense_Mutation'], col = col['Missense_Mutation']))
    },
    Translation_Start_Site = function(x, y, w, h) {
        #grid.rect(x, y, w*0.95, h*0.33, gp = gpar(fill = col['Translation_Start_Site'], col = NA))
        grid.points(x, y, pch = 1,  gp = gpar(fill = col['Translation_Start_Site'], col = col['Translation_Start_Site']))
    },
    Nonsense_Mutation = function(x, y, w, h) {
        #grid.rect(x, y, w*0.95, h*0.33, gp = gpar(fill = col['Nonsense_Mutation'], col = NA))
        grid.points(x, y, pch = 1,  gp = gpar(fill = col['Nonsense_Mutation'], col = col['Nonsense_Mutation']))
    },
    Nonstop_Mutation = function(x, y, w, h) {
        #grid.rect(x, y, w*0.95, h*0.33, gp = gpar(fill = col['Nonstop_Mutation'], col = NA))
        grid.points(x, y, pch = 1,  gp = gpar(fill = col['Nonstop_Mutation'], col = col['Nonstop_Mutation']))
    },
    Silent = function(x, y, w, h) {
        #grid.rect(x, y, w*0.95, h*0.33, gp = gpar(fill = col['Silent'], col = NA))
        grid.points(x, y, pch = 1,  gp = gpar(fill = col['Silent'], col = col['Silent']))
    },
    Splice_Site = function(x, y, w, h) {
        #grid.rect(x, y, w*0.95, h*0.33, gp = gpar(fill = col['Splice_Site'], col = NA))
        grid.points(x, y, pch = 1,  gp = gpar(fill = col['Splice_Site'], col = col['Splice_Site']))
    },
    Splice_Region = function(x, y, w, h) {
        #grid.rect(x, y, w*0.95, h*0.33, gp = gpar(fill = col['Splice_Region'], col = NA))
        grid.points(x, y, pch = 1,  gp = gpar(fill = col['Splice_Region'], col = col['Splice_Region']))
    },
    "3'UTR" = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = rcb[8], col = NA))
    },
    "5'UTR" = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = rcb[8], col = NA))
    },
    "3'Flank" = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = rcb[8], col = NA))
    },
    "5'Flank" = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = rcb[8], col = NA))
    },
    Intron = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = rcb[8], col = NA))
    },
    Amplification = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = ccc[1], col = NA))
    },
    Deletion = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = ccc[2], col = NA))
    },
    Neutral = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = ccc[3], col = NA))
    },
    background = function(x, y, w, h) {
        #grid.rect(x, y, w*0.95, h*0.1, gp = gpar(fill = col['background'], col = NA))
        grid.points(x, y, pch = 3,  gp = gpar(fill = col['background'], col = col['background']))
    }
)

##OLD COPY N"UMBE#R DATA 



op <- oncoPrint(dmat,
    get_type = function(x) strsplit(x, ";")[[1]],
    alter_fun = alter_fun,
    col = col,
    show_column_names = TRUE,
    column_order = colnames(dmat),
    split = gsplit,
    top_annotation = HeatmapAnnotation(#cbar = anno_oncoprint_barplot(),
        HCIid = str_split_fixed(colnames(dmat),"_",3)[,1],
        Model = str_split_fixed(colnames(dmat),"_",3)[,2],
        EL = str_split_fixed(colnames(dmat),"_",3)[,3]),
    row_names_gp = gpar(fontsize = 6)

) %v% Heatmap(corcnv,row_order=rownames(corcnv),column_order=colnames(corcnv),show_column_names=F,row_names_gp = gpar(fontsize = 6),name = "Spearman CNV") %v% Heatmap(corrna,row_order=rownames(corrna),column_order=colnames(corrna),show_column_names=F,row_names_gp = gpar(fontsize = 6),name = "Spearman RNA",col=f1)

#Make sure to add this once XingYi is done 
#%v% Heatmap(corcnv,row_order=rownames(corcnv),column_order=colnames(corcnv),show_column_names=F,row_names_gp = gpar(fontsize = 6),name = "Pearson CNV")



pdf(args[8],height=10,width=7)
print(op)
dev.off()

checkMAFfile <- catmaf[,c(2:16,1,17:139)]
write.table(checkMAFfile,"Processed_data/checkMAF.txt",sep="\t",quote=F,row.names=F)
 
#This is where I'm going to generate some correlation statistics bases and inter and intra correlation variablility 
diag(corrna)=NA
rna_diag <- corrna * lower.tri(corrna, diag = FALSE)
rna_diag[rna_diag==0] <- NA 


HCIids <- unique(str_split_fixed(colnames(corrna),"_",2)[,1])
#Testing s="HCI001"
RNAtable = NULL
for(s in HCIids){
    print(s)
    a = rna_diag[which(grepl(s,rownames(rna_diag))),which(grepl(s,colnames(rna_diag)))]
    b = corrna[which(!grepl(s,rownames(corrna))),which(grepl(s,colnames(corrna)))]
    svalues = na.omit(melt(a))$value
    ovalues = na.omit(melt(b))$value
    smean = mean(svalues)
    omean = mean(ovalues)
    ssd = sd(svalues)
    osd = sd(ovalues)
    sse = stde(svalues)
    ose = stde(ovalues)
    onesideT = t.test(svalues,ovalues,alternative="greater")
    twosideT = t.test(svalues,ovalues)
    normalcyS = shapiro.test(svalues)
    normalcyO = shapiro.test(ovalues)
    nonpara1side = wilcox.test(svalues,ovalues,alternative="greater",conf.int=T)
    nonpara2side = wilcox.test(svalues,ovalues,conf.int=T)

    myrow = data.frame("Sample" = s, "IntraMean" = smean, "IntraSE" = sse, "IntraSD" = ssd, "InterMean"= omean, "InterSE" = ose, "InterSD"= osd, "NormalcySample" = normalcyS$p.value , "NormalcyOther" = normalcyO$p.value, "OneSideT" =onesideT$p.value, "TwoSideT" = twosideT$p.value , "ManWhitneyU1side" = nonpara1side$p.value, "ManWhitneyU2side" = nonpara2side$p.value, "ManWhitneyU2sideEstimateSpread" = as.numeric(nonpara2side$estimate), "ManWhitneyU2sideCIlow"=nonpara2side$conf.int[1],"ManWhitneyU2sideCIhigh"=nonpara2side$conf.int[2] )
    RNAtable = rbind(RNAtable,na.omit(myrow))
}



#So I have the individual scores, what do the collective scores look like? 
diag(corrna)=NA
rna_diag <- corrna * lower.tri(corrna, diag = FALSE)
rna_diag[rna_diag==0] <- NA

mrna <- na.omit(melt(rna_diag))
mrna$a <- str_split_fixed(mrna$Var1,"_",5)[,1]
mrna$b <- str_split_fixed(mrna$Var2,"_",5)[,1]
mrna$Group <- ifelse(mrna$a == mrna$b,"Intra","Inter")

rna_intra <- mrna[which(mrna$Group == "Intra"),]$value
rna_inter <- mrna[which(mrna$Group == "Inter"),]$value

rnagroupsII <- t.test(rna_inter,rna_intra)
rnagroupsIII <- wilcox.test(rna_inter,rna_intra,conf.int=T)

#this is to check the statement that HR+ tumors are less similar (intra-model) than TNBC 
HRpos <- c("HCI003","HCI005","HCI011","HCI017")
mrna$HRposIntra <- ifelse(mrna$Group == "Intra" & mrna$a %in% HRpos, "HRpos",NA)
mrna$HRposIntra <- ifelse(mrna$Group == "Intra" & mrna$a %!in% HRpos, "TNBC",mrna$HRposIntra)

HRpos_intra <- mrna[which(mrna$HRposIntra == "HRpos"),] 
TNBC_intra <- mrna[which(mrna$HRpos=="TNBC"),]

HRpos_intra_tum <- mrna[which(mrna$HRposIntra == "HRpos" & grepl("TUMOR",mrna$Var2)),]
TNBC_intra_tum <- mrna[which(mrna$HRposIntra == "TNBC" & grepl("TUMOR",mrna$Var2)),]

mean(HRpos_intra_tum$value)
mean(TNBC_intra_tum$value)
 
wilcox.test(HRpos_intra_tum$value,TNBC_intra_tum$value,conf.int=T)



#This is where I'm going to generate some correlation statistics bases and inter and intra correlation variablility 
diag(corcnv)=NA
cnv_diag <- corcnv * lower.tri(corcnv, diag = FALSE)
cnv_diag[cnv_diag==0] <- NA


HCIids <- unique(str_split_fixed(colnames(corcnv),"_",2)[,1])
#Testing s="HCI001"
CNVtable = NULL
for(s in HCIids){
    print(s)
    a = cnv_diag[which(grepl(s,rownames(cnv_diag))),which(grepl(s,colnames(cnv_diag)))]
    b = corcnv[which(!grepl(s,rownames(corcnv))),which(grepl(s,colnames(corcnv)))]
    svalues = na.omit(melt(a))$value
    ovalues = na.omit(melt(b))$value
    smean = mean(svalues)
    omean = mean(ovalues)
    ssd = sd(svalues)
    osd = sd(ovalues)
    sse = stde(svalues)
    ose = stde(ovalues)
    onesideT = t.test(svalues,ovalues,altecnvtive="greater")
    twosideT = t.test(svalues,ovalues)
    normalcyS = shapiro.test(svalues)
    normalcyO = shapiro.test(ovalues)
    nonpara1side = wilcox.test(svalues,ovalues,altecnvtive="greater",conf.int=T)
    nonpara2side = wilcox.test(svalues,ovalues,conf.int=T)

    myrow = data.frame("Sample" = s, "IntraMean" = smean, "IntraSE" = sse, "IntraSD" = ssd, "InterMean"= omean, "InterSE" = ose, "InterSD"= osd, "NormalcySample" = normalcyS$p.value , "NormalcyOther" = normalcyO$p.value, "OneSideT" =onesideT$p.value, "TwoSideT" = twosideT$p.value , "ManWhitneyU1side" = nonpara1side$p.value, "ManWhitneyU2side" = nonpara2side$p.value, "ManWhitneyU2sideEstimateSpread" = as.numeric(nonpara2side$estimate), "ManWhitneyU2sideCIlow"=nonpara2side$conf.int[1],"ManWhitneyU2sideCIhigh"=nonpara2side$conf.int[2] )
    CNVtable = rbind(CNVtable,na.omit(myrow))
}

#So I have the individual scores, what do the collective scores look like? 
mcnv <- na.omit(melt(cnv_diag))
mcnv$a <- str_split_fixed(mcnv$Var1,"_",5)[,1]
mcnv$b <- str_split_fixed(mcnv$Var2,"_",5)[,1]
mcnv$Group <- ifelse(mcnv$a == mcnv$b,"Intra","Inter")

cnv_intra <- mcnv[which(mcnv$Group == "Intra"),]$value
cnv_inter <- mcnv[which(mcnv$Group == "Inter"),]$value

cnvgroupsII <- t.test(cnv_inter,cnv_intra)
cnvgroupsIII <- wilcox.test(cnv_intra,cnv_inter,conf.int=T)

#Now it is time to caluclate gene level differences AMP - from the PT. 
#Step one - convert LLR to CALL 
#mcnv2 <- melt(cnvmat)
#mcnv2$Plat <- str_split_fixed(mcnv2$Var2,"_",5)[,5]
#mcnva <- mcnv2[which(mcnv2$Plat == "CNVA"),]
#mcnvi <- mcnv2[which(mcnv2$Plat == "CNVI"),]
#mcnva$Call <- ifelse(mcnva$value > .9, "Amplification", ifelse(mcnva$value < -.9, "Deletion", NA))
#mcnvi$Call <- ifelse(mcnvi$value > .6, "Amplification", ifelse(mcnvi$value < -.6, "Deletion", NA))
#mcnv3 <- rbind(mcnva,mcnvi)
#mcnv3$HCI <- str_split_fixed(mcnv3$Var2,"_",5)[,1]
#cnvsamps <- unique(mcnv3$HCI)
#
#CNT_cnvs <- NULL 
#CNVcalls <- data.frame(mcnv3 %>% group_by(HCI,Var1,Call) %>% count())
#dcnv <- dcast(CNVcalls,HCI+Var1~Call,value.var="n")
#
##https://stackoverflow.com/questions/37801338/count-nas-per-row-in-dataframe
#dcnv$na_count <- apply(dcnv, 1, function(x) sum(is.na(x)))
#table(dcnv$na_count)
#dcnv[which(!is.na(dcnv$Amplification) & !is.na(dcnv$Deletion)),] #There are none that switch from AMP to DEL.between PT to PDxO.
#dcnv$gloc <- str_split_fixed(dcnv$Var1,":",2)[,2]
#cnv_HCI027 <- dcnv[which(dcnv$HCI == "HCI027" & dcnv$na_count == 1),] 
#cnv_HCI027$gloc <- str_split_fixed(cnv_HCI027$Var1,":",2)[,2]
#
#
#
#
#
##Take a deeper look 
#mcnv2[which(mcnv2$Var1 == "SNCG:10q23.2"),]
##Check the regions around the genes 
#mcnv2[which(grepl("10q23.2",mcnv2$Var1)),]
#
#mcnv2$Tech <- str_split_fixed(mcnv2$Var2,"_",5)[,5]
#plot(density(na.omit(mcnv2[which(mcnv2$Tech == "CNVA"),]$value)))
#lines(density(na.omit(mcnv2[which(mcnv2$Tech == "CNVI"),]$value)),col="red")
#
#
##Take a deeper look into HCI-001 and HCI-010  
#unique(mcnv2[,c("HCI","Tech")])
#ai <- c("HCI001","HCI002","HCI010")
#mcnv3 <- mcnv2[which(mcnv2$HCI %in% ai),]
#plot(density(na.omit(mcnv3[which(mcnv3$Tech == "CNVA"),]$value)))
#lines(density(na.omit(mcnv3[which(mcnv3$Tech == "CNVI"),]$value)),col="red")
#
#mcnv4 <- mcnv2[which(mcnv2$HCI == 'HCI001'),]
#plot(density(na.omit(mcnv4[which(mcnv4$Tech == "CNVI"),]$value)),col="red",xlim=c(-3,3))
#lines(density(na.omit(mcnv4[which(mcnv4$Tech == "CNVA"),]$value)))
#abline(v=.9)
#abline(v=-.9)
#
#p <- ggplot(mcnv4,aes(x=value,fill=Var2))
#p <- p + geom_density(alpha=.5)
#p <- p + theme_classic()
#p <- p + geom_vline(xintercept=.9, linetype="dashed", color = "blue", size=1)
#p <- p + geom_vline(xintercept=1.2, linetype="dashed", color = "green", size=1)
#p <- p + geom_vline(xintercept=-.9, linetype="dashed", color = "green", size=1)
#p <- p + geom_vline(xintercept=-.7, linetype="dashed", color = "blue", size=1)
#p <- p + xlim(c(-3,3))
#p <- p + theme(legend.position="bottom")
#p
#
#
#
##I need to figure out if there are any genes that are unique to one technology 
#mcnv2[which(is.na(mcnv2$value)),c(mcnv2$Tech,mcnv2$Var1)]
#
#
##Double check genes in common CNVs
#
#tcgaBRCAcnv <- c("RB1","PTEN","CDKN2A","KMT2C","MAP2K4","TP53","ERBB2","PIK3CA","EGFR","FOXA1","MDM2","CCND1","CSMD1","PTPRD","STK11")
#cbioBRCAcnv <- c("FGF3","CDK12","NOTCH2","H3P6","ERBB2","SIPA1L3","ADCY9","FAM72C","SDK2","ZNF217","MYC","FGF3","CDK18","STX4","CSMD1","PTEN","TNFRSF10C","NKX3-1")
#pancan2013 <- c("CCND1","EGFR","MYC","ERBB2","CCNE1","MCL1","MDM2","NSD3","FGFR1","TERC","TERT","RMRP","ATM","NOTCH1")
#pancan2010 <- c("MYC","CCND1","ERBB2","CDK4","NKX2-1","MDM2","EGFR","FGFR1","KRAS","MCL1","BCL2L1")
#gynbreast <- c("CSMD1","STK11","MECOM","ZNF217","MCL1")
#foundation <- c("LYN","JAK1","JAK2","CD274","PDCD1LG2","ERBB4")
#ambry = c("ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CHEK2", "MRE11A", "MUTYH", "NBN", "NP1", "PALB2", "PTEN", "RAD50", "RAD51C", "RA051D", "STK11", "TP53")
#
#cnvglist <- unique(c(tcgaBRCAcnv,cbioBRCAcnv,pancan2013,pancan2010,gynbreast,foundation,ambry))
#cnvglgrep <- paste(cnvglist,sep=":|",collapse=":|") 
#dcnv$Gene <- str_split_fixed(dcnv$Var1,":",2)[,1]
#gldcnv <- dcnv[which(dcnv$Gene %in% cnvglist),]
#gldcnv[which(gldcnv$na_count < 2),]
#mcnv2$Gene <- str_split_fixed(mcnv2$Var1,":",2)[,1]
#glmcnv2 <- mcnv2[which(mcnv2$Gene %in% cnvglist),]
#
#glmcnv2$gloc <- glmcnv2$Var1
#glmcnv2[which(grepl("BRCA2:",glmcnv2$Var1)),]
##Take a look  glmcnv2[order(glmcnv2$Var1),]
##Questionable genes (TERC, FAM72C, LYN, CSMD1, ERBB4)
##HCI-023 ERBB4? 
##HCI-001 PTEN
##HCI-023,002 NKX3-1
##HCI-010 -002 PTPRD
##HCI-010 CDKN2A 
##HCI-002 NSD3
##HCI-010 EGFR
##HCI-010 MCL1
##HCI-001 MCL1
##001,002 MYC
##Hci-010 kras
##001,010 stk11
##001,002 CDK18
##HCI010   SIPA1L3
##HCI010,002       NBN
##HCI010     KMT2C
##HCI010      CDH1
#
##To try to set appropriate cutoffs for LRR affy and Illumina
#
# 
##https://stackoverflow.com/questions/37801338/count-nas-per-row-in-dataframe
#cnvna$na_count <- apply(cnvna, 1, function(x) sum(x))
#
##So I'm going to simulate some cutoffs for CNVA and CNVI in cancer genes to set approriate cutoffs. 
##And I want to figure out with one will minimize the difference between affy and illumina runs
#alrrs <- seq(0.7,1.25,0.05)
#dlrrs <- seq(-1.2,-0.7,0.05)
#
#mcnv2 <- melt(cnvmat)
#mcnv2$HCI <- str_split_fixed(mcnv2$Var2,"_",5)[,1]
#mcnv2$Tech <- str_split_fixed(mcnv2$Var2,"_",5)[,5]
#mcnv2$Call <- ifelse(mcnv2$value > .9, "Amplification", ifelse(mcnv2$value < -.9, "Deletion", NA))
#SIM_CUTOFFS = NULL
#cnt = 0
#cnvi <- colnames(cnvmat)[grepl("CNVI",colnames(cnvmat))]
#cnva <- colnames(cnvmat)[grepl("CNVA",colnames(cnvmat))]
#cnvmati <- cnvmat[,as.character(cnvi)]
#cnvmata <- cnvmat[,as.character(cnva)]
#
#
##First attempt to get the logic right
#calc_cutoff <- function(mcnv,iamp,idel,aamp,adel,cnt){
#                mcnv2$Call <- ifelse(mcnv2$Tech == "CNVA" & mcnv2$value > aamp, "Amplification", ifelse(mcnv2$value < adel, "Deletion", NA))
#                mcnv2$Call <- ifelse(mcnv2$Tech == "CNVI" & mcnv2$value > iamp, "Amplification", ifelse(mcnv2$value < idel, "Deletion", mcnv2$Call))
#                CNVcalls <- data.frame(mcnv2 %>% group_by(HCI,Var1,Call) %>% count())
#                dcnv <- dcast(CNVcalls,HCI+Var1~Call,value.var="n")
#                dcnv$na_count <- apply(dcnv, 1, function(x) sum(is.na(x)))
#                dcnvtbl <- table(dcnv$na_count) 
#                o = data.frame("IlluAMP"=iamp,"IlluDEL"=idel,"AffyAMP"=aamp,"AffyDEL"=adel,"State_3"=dcnvtbl[1],"State_2"=dcnvtbl[2],"State_1"=dcnvtbl[3])
#                return(o)
#
#}
#
#for(iamp in alrrs){
#    for(idel in dlrrs){
#        for(aamp in alrrs){
#            for(adel in dlrrs){
#                cnt = cnt+1 
#                o1 <- calc_cutoff(mcnv,iamp,idel,aamp,adel,cnt)
#                #print(system.time(o2 <- calc_matrix(cnvati,cnvmatm,iamp,idel,aamp,adel,cnt)))
#                SIM_CUTOFFS = rbind(SIM_CUTOFFS,o1)                    
#                print(cnt)
#            }
#        }
#    }
#}
#
##Now trying this in parallel
#res <- foreach(iamp=alrrs, .combine='rbind') %:%
#    foreach(idel=dlrrs,.inorder=FALSE) %dopar% {
#        for(aamp in alrrs	){
#            for(adel in dlrrs){
#                cnt = cnt+1
#                print(cnt)
#                calc_cutoff(mcnv,iamp,idel,aamp,adel,cnt)
#            }
#        }
#    }
#res
#
#
#
##I may need to move this to a python loop for speed.... 
#calc_matrix <- function(cnvati,cnvmatm,iamp,idel,aamp,adel,cnt){
#    cnvmati_call <- ifelse(cnvmati > 1, "AMP",ifelse(cnvmati < -1, "DEL", "NA"))
#    cnvmata_call <- ifelse(cnvmata > 1, "AMP",ifelse(cnvmata < -1, "DEL", "NA"))
#    mci = melt(cnvmati_call)
#    mca = melt(cnvmata_call)
#    mc <- rbind(mci,mca)
#    mc$HCI <- str_split_fixed(mc$Var2,"_",5)[,1]
#    CNVcalls <- data.frame(mc %>% group_by(HCI,Var1,value) %>% count())
#    dcnv <- dcast(CNVcalls,HCI+Var1~value,value.var="n")
#    dcnv$na_count <- apply(dcnv, 1, function(x) sum(is.na(x)))
#    dcnvtbl <- table(dcnv$na_count)
#    o = data.frame("IlluAMP"=iamp,"IlluDEL"=idel,"AffyAMP"=aamp,"AffyDEL"=adel,"State_3"=dcnvtbl[1],"State_2"=dcnvtbl[2],"State_1"=dcnvtbl[3])
#    return(o)
#}
#
#
