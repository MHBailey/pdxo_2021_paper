library(data.table)
library(limma)
library(stringr)
library(dplyr)
library(matrixStats)
library('DESeq2')

'%!in%' <- function(x,y)!('%in%'(x,y))

#https://stackoverflow.com/questions/20159275/removing-matrix-rows-if-values-of-a-cloumn-are-outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}



args = commandArgs(trailingOnly=TRUE)

if (length(args)!=13) {
  stop("Check out the rule normalize RNA", call.=FALSE)
}

#args <- c("Data/RNAseq/combined.TPM.welmRNA.wide.cleaned.txt","Data/RNAseq/metadata_20200714.txt","Data/Meta/RNAtechnologies.txt","Data/GeneLists/HUGO_ENSG.20200625.txt","Processed_data/normalized.TPM.batch","Figures/RNA_Batches_HCI001.pdf","Data/RNAseq/combined.Ecnt.welmRNA.wide.cleaned.txt","Processed_data/normalized.DESEQ.RSEM.txt","Processed_data/normalized.DESEQ.RSEM.batch","Data/nolist.20200921.txt","Processed_data/combined.Ecnt.welmRNA.wide.cleaned.txt","Processed_data/top1000.DEgenes.nobatch.txt","Processed_data/top1000.DEgenes.batched.txt")



df <- fread(args[1])
nolist <- fread(args[10],header=F)
dfcol <- colnames(df)
goodsamps <- dfcol[which(dfcol %!in% nolist$V1)]
dat <- df %>% dplyr::select(goodsamps)
bat <- fread(args[2], header=F)
colnames(bat) <- c("hawkid","gnomex")
bat$batch <- as.numeric(str_split_fixed(bat$gnomex,"X",2)[,1])
tech <- fread(args[3],header=F)
colnames(tech) <- c("gbatch","technology")
hugo <- fread(args[4])
pseudo <- hugo[which(grepl("pseudogene",hugo$"Approved name")),]$"Approved symbol"
RNAgenes <- hugo[which(grepl("RNA",hugo$"Approved name")),]$"Approved symbol"

#SuppFigure Building POLYA v RIBOZERO: 
h1 <- names(dat)
h2 <- h1[which(grepl("HCI001",h1))]
hci001 <- data.frame(dat %>% select(h2))
rownames(hci001) <- dat$'Approved symbol'

bt <- merge(bat,tech,by.x="batch",by.y="gbatch",all.x=T)
bthci001 <- bt[which(bt$hawkid %in% h2),]
rat001 <- data.frame(hci001 %>% select(HCI001_pdx_P2_XXX_RNA,HCI001_pdx_P3_XXX_RNA))
#toshow
data.frame(bt[order(bt$hawkid),])

#Patient POLYA, PDX RIBOZ
pdf("Figures/Supp.Poly_vs_Poly.TPM.pdf",height=5,width=5)
plot(x=log2(hci001$HCI001_tumor_patient_XXX_RNA+1),y=log2(hci001$HCI001_pdx_P2_XXX_RNA+1),pch=16)
abline(0,1,col='red')
dev.off()

pdf("Figures/Supp.Poly_vs_RiboZ.TPM.pdf",height=5,width=5)
plot(x=log2(hci001$HCI001_tumor_patient_XXX_RNA+1),y=log2(hci001$HCI001_pdx_P3_XXX_RNA+1),pch=16)
abline(0,1,col='red')
dev.off()



#GENEnames to remove
yuck <- c("MT-")
yuckstr <- paste(yuck,sep="|",collapse="|")

#BUILDING THE MATRIX
raw <- dat
rawg <- raw[which(!grepl(yuckstr,raw$"Approved symbol") & raw$"Approved symbol" %!in% RNAgenes & raw$"Approved symbol" %!in% pseudo),]
rownames(rawg) <- rawg$"Approved symbol"
rawg$"Approved symbol" <- NULL
rawg$gene_id <- NULL
lrawg <- log2(rawg+1)
ldf <- data.frame(lrawg)
rownames(ldf) <- rownames(rawg)
lrf  <- ldf[which(is.finite(rowSums(ldf))),,drop=FALSE]
mraw <- as.matrix(lrf,rownames=rownames(lrf))

#BUILDING coldata
colnames(mraw) <- colnames(lrawg)
tokeep <- colnames(mraw)
coldat <- bt[which(bt$hawkid %in% tokeep),]

oc3 <- coldat[order(coldat$hawkid),]

#Set to integers 

traw <- removeBatchEffect(mraw,batch=oc3$technology)
dt <- data.frame(traw)

plot(x=dt$HCI001_pdx_P2_XXX_RNA,y=dt$HCI001_pdx_P3_XXX_RNA,pch=16)
abline(0,1,col='red')

braw <- removeBatchEffect(traw,batch=factor(oc3$batch))
db <- data.frame(braw)
dbs <- db[c("HCI001_tumor_patient_XXX_RNA","HCI001_pdx_P2_XXX_RNA")]
plot(x=db$HCI016_pdx_P2_XXX_RNA,y=db$HCI016_tumor_patient_XXX_RNA,pch=16)
abline(0,1,col='red')


#write.table(traw,args[5],sep="\t",quote=F) #this was good too 
write.table(braw,args[5],sep="\t",quote=F) #Going with this one as raw values for correlation


#############So I found a publication that I need to implement here the aggregate  data from PolyA and RiboZero rRNA reduction strategies. #### This should work well. Bush et al. 2017 BMC Bioinfo. We'll start with RAWG 
#
#poly = bt[which(bt$technology == "polyAselction"),]$hawkid
#ribo = bt[which(bt$technology == "RiboZero" & bt$hawkid %in% tokeep),]$hawkid
#
#prawg = rawg %>% select(poly)
#rrawg = rawg %>% select(ribo)
#
#origcols <- c(colnames(prawg),colnames(rrawg))
#
#apr <- rowMedians(as.matrix(prawg,colnames=colnames(prawg)))
#arr <- rowMedians(as.matrix(rrawg,colnames=colnames(rrawg)))
#
#por <- apr/arr
#rop <- arr/apr
#
#arrawg <- rrawg*por
#aprawg <- prawg*rop
#
#pr <- data.frame(cbind(prawg,arrawg))
#pp <- data.frame(cbind(aprawg,rrawg))
#rownames(pr) <- rownames(rawg)
#prna <- na.omit(pr,drop=F)
#colnames(prna) = origcols
#ppna <- na.omit(pp,drop=F)
#colnames(ppna) = origcols 
#
#plot(x=log2(prna$HCI001_tumor_patient_XXX_RNA+1),y=log2(prna$HCI001_pdx_P2_XXX_RNA+1),pch=16)
#abline(0,1,col='red')
#
#plot(x=log2(prna$HCI001_tumor_patient_XXX_RNA+1),y=log2(prna$HCI001_pdx_P3_XXX_RNA+1),pch=16)
#abline(0,1,col='red')
#
#which(apr >10000)
#rownames(prawg)[10290]
#
##write.table(prna,args[5],sep="\t",quote=F) This one was okay if you log2 normalize again
#
##REORDER SAMPLES
#r <- data.frame(prna %>% select(colnames(rawg)))
#colnames(r) <- oc3$hawkid
#br <- removeBatchEffect(r,batch=factor(oc3$batch))
##write.table(br,args[5],sep="\t",quote=F) This was BAD
#
#plot(x=log2(br[,"HCI023_pdxo_P1_68_RNA"]+1),y=log2(br[,"HCI023_pdxo_P1_336_RNA"]+1),pch=16)
#abline(0,1,col='red')
#


rnacnt1 <- fread(args[7])
rnacnt <- rnacnt1 %>% dplyr::select(goodsamps)
rownames(rnacnt) <- rnacnt$"Approved symbol"
rnacnt$"Approved symbol" <- NULL
rnacnt$gene_id <- NULL 

mr <- as.matrix(rnacnt,rownames=rownames(rnacnt))
rr <- round(mr)
mro <- apply(rr,2,remove_outliers)

bt$HCI <- str_split_fixed(bt$hawkid,"_",5)[,1]
btm = as.matrix(bt[which(bt$hawkid %in% colnames(rr)),])

dds <- DESeqDataSetFromMatrix(countData = rr,
                              colData = btm,
                              design = ~ HCI) 
#This will run the differential expression test....
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)
res <- results(dds)
resdf <- data.frame(res)
#Get ordered list of genes
top <- data.frame(resdf[order(resdf$pvalue),])

write.table(normalized_counts, args[8], sep="\t", quote=F)
write.table(top,args[12], sep="\t", quote=F)


ddsB <- DESeqDataSetFromMatrix(countData = rr,
                              colData = btm,
                              design = ~ HCI+technology)   
ddsB <- DESeq(ddsB)
normalized_counts <- counts(ddsB, normalized=TRUE)
resB <- results(ddsB)
resdfB <- data.frame(resB)

#Get ordered list of genes 
topB <- data.frame(resdfB[order(resdfB$pvalue),])
write.table(normalized_counts, args[9], sep="\t", quote=F)
write.table(topB,args[13], sep="\t", quote=F)

write.table(dat,args[11],sep="\t", quote=F,row.names=F)

pdf(args[6],height=3,width=3)
plot(1,1)
dev.off()

