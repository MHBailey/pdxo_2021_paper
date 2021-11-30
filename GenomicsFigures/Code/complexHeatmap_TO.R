library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(viridis)
library(data.table)

#this is source activate stats 
uniqString <- function(x){
    yo <- unique(unlist(strsplit(x, ";")))
    return(paste(yo,collapse=";"))
}
'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=8) {
  stop("Check out the rule complexHeatmap", call.=FALSE)
}


#Testing DNA
#args <- c("Data/MAFs/Indepth_2020_04_22.all.HG38.PASS.maf","Data/Meta/cat_DNA_Indepth.txt","Data/GeneLists/299.BRCA.glist.txt","Data/CNV/cnv_combined.txt","Data/Meta/cat_CNV_old.txt","Processed_data/normalized.batch","Data/Meta/cat_RNA_Indepth.txt","Data/RNAseq/metadata_20200313.txt")

#BYU Meta
#args <- c("Data/MAFs/Indepth_2020_04_22.all.HG38.PASS.maf","Data/Meta/BYU.dna.txt","Data/GeneLists/299.BRCA.glist.txt","Data/CNV/cnv_combined.txt","Data/Meta/BYU.cnv.txt","Processed_data/normalized.batch","Data/Meta/BYU.rna.txt","Data/RNAseq/metadata_20200313.txt")

#TO_DATA everything
#args <- c("Processed_data/TumorOnly_2020_09_15.all.HG38.PASS.indels.cancer_miss_mr.maf","Data/Meta/to.dna.txt","Data/CNV/combined.xingyi.20200724.txt","Data/Meta/to.cnv.txt","Processed_data/normalized.DESEQ.RSEM.batch","Data/Meta/to.rna.txt","Data/Clinical/HCI_clinical.txt","Figures/pdx.to.genomicsFig.20200923.pdf")

#TO_DATA Dec 10, 2020 
#args <- c("Processed_data/TumorOnly_2020_12_10.all.HG38.PASS.indels.cancer_miss_mr.maf","Data/Meta/to.dna.20201210.txt","Data/CNV/combined.xingyi.20200724.txt","Data/Meta/to.cnv.txt","Processed_data/normalized.DESEQ.RSEM.batch","Data/Meta/to.rna.txt","Data/Clinical/HCI_clinical.txt","Figures/pdx.to.genomicsFig.20201210.pdf")

#TO_DATA Dec 23, 2020 
#args <- c("Processed_data/TumorOnly_2020_12_23.all.HG38.PASS.indels.cancer_miss_mr.maf","Data/Meta/to.dna.20201223.txt","Data/CNV/combined.xingyi.20200724.txt","Data/Meta/to.cnv.txt","Processed_data/normalized.DESEQ.RSEM.batch","Data/Meta/to.rna.txt","Data/Clinical/HCI_clinical_Dec2020.txt","Figures/pdx.to.genomicsFig.20201223.pdf")

#LV DATA Jan 19, 2021
#args <- c("Processed_data/TumorOnly_2020_12_23.all.HG38.PASS.indels.cancer_miss_mr.maf","Data/Meta/lv.dna.20210119.txt","Data/CNV/combined.xingyi.20210118.txt","Data/Meta/lv.cnv.20210119.txt","Processed_data/normalized.DESEQ.RSEM.batch","Data/Meta/lv.rna.20210119.txt","Data/Clinical/HCI_clinical_Dec2020.txt","Figures/pdx.to.genomicsFig.20210119.LV.pdf")

#Cleaning up the MAF a bit, removing noncoding mutations 
maf <- fread(args[1])
maf <- maf[which(maf$"Variant_Classification" != "Silent"),]
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

allMUTS <- data.frame("HUGO"=catmaf$Hugo_Symbol,"ELID"=catmaf$elid,"Var_Class"=catmaf$Variant_Classification)

glmaf2 <- as.data.table(catmaf)[, toString(Variant_Classification), by = list(elid, Hugo_Symbol)]

glmaf2$V1 <- str_replace_all(glmaf2$V1,", ",";")
glmaf3 <- dcast(glmaf2, Hugo_Symbol~elid,value.var="V1")
rownames(glmaf3) <- glmaf3$Hugo_Symbol
glmaf3$Hugo_Symbol = NULL

#ORDER COLUMS and ROWS
glmaf4 <- as.data.frame(glmaf3)[mafmeta$meta_elid]
rownames(glmaf4) <- rownames(glmaf3)
alldna <- colnames(glmaf4)
dmat <- as.matrix(glmaf4,rownames=rownames(glmaf4))


#NOW ADD CNV

cnv <- fread(args[3])
cnvmeta <- fread(args[4],header=F)
colnames(cnvmeta) <- c("Category","Tumor_Sample_Barcode")
cnvmeta$hawkid <- str_split_fixed(cnvmeta$Tumor_Sample_Barcode,"_",5)[,1]
cnvmeta$meta_elid <- paste(cnvmeta$hawkid,cnvmeta$Category,sep="_")
cnvmeta$tech <- str_split_fixed(cnvmeta$Tumor_Sample_Barcode,"_",5)[,5]
cnv$Gloc <- paste(cnv$"Approved symbol",cnv$Chromosome,sep=":")
rownames(cnv) <- cnv$Gloc


#List of CNV genes 
tcgaBRCAcnv <- c("RB1","PTEN","CDKN2A","KMT2C","MAP2K4","TP53","ERBB2","PIK3CA","EGFR","FOXA1","MDM2","CCND1","CSMD1","PTPRD","STK11")
cbioBRCAcnv <- c("FGF3","CDK12","NOTCH2","H3P6","ERBB2","SIPA1L3","ADCY9","FAM72C","SDK2","ZNF217","MYC","FGF3","CDK18","STX4","CSMD1","PTEN","TNFRSF10C","NKX3-1")
pancan2013 <- c("CCND1","EGFR","MYC","ERBB2","CCNE1","MCL1","MDM2","NSD3","FGFR1","TERC","TERT","RMRP","ATM","NOTCH1")
pancan2010 <- c("MYC","CCND1","ERBB2","CDK4","NKX2-1","MDM2","EGFR","FGFR1","KRAS","MCL1","BCL2L1")
gynbreast <- c("CSMD1","STK11","MECOM","ZNF217","MCL1")
foundation <- c("LYN","JAK1","JAK2","CD274","PDCD1LG2","ERBB4")
ambry = c("ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CHEK2", "MRE11A", "MUTYH", "NBN", "NP1", "PALB2", "PTEN", "RAD50", "RAD51C", "RAD51D", "STK11", "TP53")

cnvglist <- unique(c(tcgaBRCAcnv,cbioBRCAcnv,pancan2013,pancan2010,gynbreast,foundation,ambry))
cnvgl <- cnv[which(cnv$"Approved symbol" %in% cnvglist),]
rownames(cnvgl) <- cnvgl$Gloc
#rownames(cnvgl) <- cnvgl$"Approved symbol" #This is add to catmaf if I need to... 

cnv2 <- subset(cnvgl,select=cnvmeta$Tumor_Sample_Barcode)
rownames(cnv2) <- rownames(cnvgl)
colnames(cnv2) <- cnvmeta$meta_elid

cnv3 <- as.matrix(cnv2,rownames=rownames(cnv2))
cnv4tmp <- melt(cnv3)
cnv4 <- merge(cnv4tmp,cnvmeta,by.x="Var2",by.y="meta_elid",all.x=T) #left join 
cnv4a <- cnv4[which(cnv4$tech == "CNVA"),]
cnv4i <- cnv4[which(cnv4$tech == "CNVI"),]
cnv4a$Call <- ifelse(cnv4a$value > .9, "Amplification", ifelse(cnv4a$value < -.9, "Deletion", NA))
cnv4i$Call <- ifelse(cnv4i$value > .6, "Amplification", ifelse(cnv4i$value < -.6, "Deletion", NA))
cnv4recall = rbind(cnv4a,cnv4i)
cnv5 <- na.omit(cnv4recall)
cnv6 <- data.frame("HUGO"=cnv5$Var1,"ELID"=cnv5$Var2,"Var_Class"=cnv5$Call)
cnv7 <- dcast(cnv4recall, Var1~Var2,value.var="Call")
cnv7[is.na(cnv7)] <- ""
rownames(cnv7) <- cnv7$Var1
cnv7$Chrom <- str_split_fixed(cnv7$Var1,":",2)[,2]
dupped <- cnv7$Chrom[which(duplicated(cnv7$Chrom))]
cnv7[which(cnv7$Chrom %in% dupped),]
cnv7$Var1 <- NULL 
cnv7$Chrom <- NULL

cnvmat <- as.matrix(cnv7,rownames=rownames(cnv7)) #This is for AMP/DEL calls 

#This is to calculate number of events per sample. 
dfcnv <- data.frame(cnvmat)
dfcnv2 <- dfcnv %>% mutate_all(as.character)
dfcnv2[dfcnv2 == ""] <- 0 
dfcnv2[dfcnv2 == "Amplification"] <- 1
dfcnv2[dfcnv2 == "Deletion"] <- 1
dfcnv3 <- dfcnv2 %>% mutate_all(as.numeric)
dfcnv3$sums <- rowSums(dfcnv3)



#THIS IS TO PUT IT ALL IN ONE PLOT
allMUTS$ELID <- as.character(allMUTS$ELID)
allMUTS$ELID[allMUTS$ELID == "TOW18S_PDX_ONLY"] <- "HCI041_PDX_ONLY"
allMUTS$ELID[allMUTS$ELID == "TOW19_PDX_ONLY"] <- "HCI043_PDX_ONLY"

allMUTS <- rbind(allMUTS,cnv6)
allMUTS2 <- as.data.table(allMUTS)[, toString(Var_Class), by = list(ELID, HUGO)]
allMUTS2$V1 <- str_replace_all(allMUTS2$V1,", ",";")
allMUTS3 <- dcast(allMUTS2, HUGO~ELID,value.var="V1")
rownames(allMUTS3) <- allMUTS3$HUGO
allMUTS3$HUGO = NULL

allMUTS4 <- as.matrix(allMUTS3,rownames=rownames(allMUTS3))

#NOW ADD RNASEQ 
rna <- fread(args[5])
rnameta <- fread(args[6],header=F)
colnames(rnameta) <- c("eid","hawkid")

#Special consideration for 
r <- data.frame(rna %>% select(rnameta$hawkid))

rownames(r) <- rna$V1

#THIS IS THE BIG MOVE HERE: the 23 Primary tumor is actuall going to be the 23
#THESE ARE FROM THE SAME PATIENT BUT IT IS CRAZY TO THINK ABOUT.
names(r)[names(r) == "HCI024_tumor_patient_XXX_RNA"] <- "HCI023_tumor_patient_XXX_RNA" 

lr <- log2(r+1)
lrf  <- lr[which(is.finite(rowSums(lr))),]

mr <- as.matrix(r)
mrl <- as.matrix(lrf)

corrna <- cor(mr)
f1 = colorRamp2(seq(0.5,1, length = 3), c("blue","#EEEEEE","red"))
#Heatmap(corrna,cluster_rows = FALSE,cluster_columns = FALSE,col=f1)


#COLORS FOR ONCOPRINT


rcb <- brewer.pal(name="Set1",n=9)
ccc <- brewer.pal(name="Set2",n=3)

rcb[9] <- "#DCDCDC"
col = c(Missense_Mutation=rcb[1], Translation_Start_Site=rcb[2], Nonsense_Mutation=rcb[2], Nonstop_Mutation=rcb[3], Silent=rcb[3], Splice_Site=rcb[4], Splice_Region=rcb[4], In_Frame_Ins=rcb[5],In_Frame_Del=rcb[5], Frame_Shift_Ins=rcb[7], Frame_Shift_Del=rcb[7],"3'UTR"=rcb[8],"5'UTR"=rcb[8],"3'Flank"=rcb[8],"5'Flank"=rcb[8],"Intron"=rcb[8],"background"=rcb[9],Amplification=ccc[2],Deletion=ccc[1],Neutral=ccc[3])



alter_fun = list(
    background = function(x, y, w, h) {
        #grid.rect(x, y, w*0.95, h*0.1, gp = gpar(fill = col['background'], col = NA))
        grid.points(x, y, pch = 3,  gp = gpar(fill = col['background'], col = col['background']))
    },
    In_Frame_Del = function(x, y, w, h) {
        grid.rect(x, y, w*0.95, h*0.33, gp = gpar(fill = col['In_Frame_Del'], col = NA))
    },
    In_Frame_Ins = function(x, y, w, h) {
        grid.rect(x, y, w*0.33, h*0.95, gp = gpar(fill = col['In_Frame_Ins'], col = NA))
    },
    Frame_Shift_Ins = function(x, y, w, h) {
        grid.rect(x, y, w*0.33, h*0.95, gp = gpar(fill = col['Frame_Shift_Ins'], col = NA))
    },
    Frame_Shift_Del = function(x, y, w, h) {
        grid.rect(x, y, w*0.95, h*0.33, gp = gpar(fill = col['Frame_Shift_Del'], col = NA))
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
        #grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = ccc[1], col = NA))
        grid.points(x, y, pch = 2,  gp = gpar(fill = col['Amplification'], col = col['Amplification']))
    },
    Deletion = function(x, y, w, h) {
        #grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = ccc[2], col = NA))
        grid.points(x, y, pch = 6,  gp = gpar(fill = col['Deletion'], col = col['Deletion']))
    },
    Neutral = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = ccc[3], col = NA))
    }
)

#Now I'm going to add the clinical data to the second plot and we can go from there 
clin <- fread(args[7])
clin$HCIid <- str_replace(clin$SampleID,"-","")
tumsub <- clin[,c("HCIid","Anatomy","TripleNeg","Type","BMI","smokinghistory","ERstatus","PRstatus","HER2status")]


HCIsamples <- str_replace(str_replace(colnames(dmat),"-",""),"_PDX_ONLY","")
HCIsamples[HCIsamples == 'TOW18S'] <- "HCI041"
HCIsamples[HCIsamples == 'TOW19'] <- "HCI043"
subclin <- clin[which(clin$HCIid %in% HCIsamples),]


#Quick change to DMAT- one change TOW18 and TOW19
colnames(dmat)[colnames(dmat) == 'TOW18S_PDX_ONLY'] <- 'HCI041_PDX_ONLY'  
colnames(dmat)[colnames(dmat) == 'TOW19_PDX_ONLY'] <- 'HCI043_PDX_ONLY'  

#This is just DMAT as I'm updating
#op <- oncoPrint(dmat, get_type = function(x) strsplit(x, ";")[[1]], alter_fun = alter_fun, col = col, show_column_names = TRUE,column_order = colnames(dmat))
#pdf(args[8],height=9.5, width=10)
#print(op)
#dev.off()

#Remove silent singletons
okay <- dmat
okay[is.na(okay)] <- ""
#https://stackoverflow.com/questions/37260361/selecting-only-unique-values-from-a-comma-separated-string
yo <- data.frame(apply(okay, 1, paste, collapse=","))
colnames(yo) <- "Muts"
yo$Muts2<- sapply(strsplit(as.character(yo$Muts), ",", fixed = TRUE), function(x) {paste(unique(x), collapse = ",")})

silent <- c(",Silent","Silent,")
good <- which(yo$Muts2 %!in% silent)
dmat2 <- dmat[good,]

if(length(colnames(dmat2)) != length(colnames(dmat2))){
    print("You lost a sample that only had singleton silent mutations")
}

dnames <- str_split_fixed(colnames(dmat2),"_",3)[,1]
dnames2 <- gsub('^(.{3})(.*)$', '\\1-\\2', dnames)
colnames(dmat2) <- dnames2
colnames(cnvmat) <- dnames2

#A quick check for no CNV alterations
cnvmat <- cnvmat[rowSums(cnvmat=="") != ncol(cnvmat),]



op <- oncoPrint(dmat2,
    get_type = function(x) strsplit(x, ";")[[1]],
    alter_fun = alter_fun,
    col = col,
    show_column_names = TRUE,
    column_order = colnames(dmat2),
    split = gsplit,
    bottom_annotation = HeatmapAnnotation(#cbar = anno_oncoprint_barplot(),
        ER=subclin$ERstatus,
        PR=subclin$PRstatus,
        HER2=subclin$HER2status,
        Type=subclin$Type,
        Anatomy=subclin$Anatomy,
        SmokingHistory=subclin$smokinghistory,
        Race=subclin$Race),
    row_names_gp = gpar(fontsize = 6)
) %v% oncoPrint(cnvmat,
    get_type = function(x) strsplit(x, ";")[[1]],
    alter_fun = alter_fun,
    col = col,
    show_column_names = TRUE,
    column_order = colnames(cnvmat),
    split = gsplit,
    row_names_gp = gpar(fontsize = 6)
)

pdf(args[8],height=19, width=10)
print(op)
dev.off()

######## This is just a bit of digging into the CNV study to get that going ####### 
#same <- c("HCI005_PDX_ONLY","HCI006_PDX_ONLY","HCI007_PDX_ONLY")
#which(colnames(cnv3) %in% same)
#cnv3s <- cnv3[,5:7]
#m3 <- melt(cnv3s)
#cnv3s[cnv3s > .9] <- .9 
#cnv3s[cnv3s < -.9] <- -.9
#m3s <- melt(cnv3s)
#
#p <- ggplot(m3s,aes(value,Var1,color=Var2))
#p <- p + geom_point(alpha=.5,size=3,pch=16)
#p <- p + theme_classic()
#p <- p + geom_vline(xintercept=0)
#p <- p + geom_vline(xintercept=0.9,color="red")
#p <- p + geom_vline(xintercept=-0.9,color="blue")
#p <- p + xlab("median centered logRratios") + ylab("Gene:Location")
#p 
#
#
#p <- ggplot(m3, aes(value,color=Var2))
#p <- p + geom_density(alpha=0.4)
#p <- p + theme_classic()
#p <- p + theme(legend.position="bottom")
#p <- p + xlab("median centered logRratios")
#p
#
#
#
#
#plot(density(na.omit(cnv[,which(grepl("CNVI",colnames(cnv)))])))
#cnvi <- as.matrix(cnv %>% dplyr::select(colnames(cnv)[which(grepl("CNVI",colnames(cnv)))]))
#rownames(cnvi) = rownames(cnv)
#
#cnva <- as.matrix(cnv %>% dplyr::select(colnames(cnv)[which(grepl("CNVA",colnames(cnv)))]))
#
#plot(density(na.omit(cnva)))
#lines(density(na.omit(cnvi)),col="blue")
#
#mcnva <- melt(cnva)
#mcnvi <- melt(cnvi)
#kmeans(na.omit(mcnva$value),centers=5,nstart=50)
#kmeans(na.omit(mcnvi$value),centers=5,nstart=50)
#df <- data.frame(cnv3)
#
#
#op2 <- oncoPrint(dmat2,
#    get_type = function(x) strsplit(x, ";")[[1]],
#    alter_fun = alter_fun,
#    col = col,
#    show_column_names = TRUE,
#    column_order = colnames(dmat2),
#    split = gsplit,
#    bottom_annotation = HeatmapAnnotation(#cbar = anno_oncoprint_barplot(),
#        ER=subclin$ERstatus,
#        PR=subclin$PRstatus,
#        HER2=subclin$HER2status,
#        Type=subclin$Type,
#        Anatomy=subclin$Anatomy,
#        SmokingHistory=subclin$smokinghistory,
#        Race=subclin$Race),
#    row_names_gp = gpar(fontsize = 6)
#) %v% Heatmap(cnv8,
#    show_column_names = TRUE,
#    column_order = colnames(cnvmat),
#    row_names_gp = gpar(fontsize = 4),
#    row_order = row_order(yo)
#)
#
##
