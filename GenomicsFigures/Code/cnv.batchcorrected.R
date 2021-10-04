library(data.table)
library(limma)
library(stringr)
library(dplyr)
library(ggplot2)



#this is source activate rnaseq2
uniqString <- function(x){
    yo <- unique(unlist(strsplit(x, ";")))
    return(paste(yo,collapse=";"))
}
'%!in%' <- function(x,y)!('%in%'(x,y))

stde <- function(x){
    sd(x)/sqrt(length(x))
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Check out the rule normalize CNV", call.=FALSE)
}
#This is just to test in real time. 

#args <- c("Data/CNV/combined.xingyi.20200724.txt","Processed_data/normalized.CNV.batch","Figures/CNV.threshold.adj.pdf") 
#args <- c("Data/CNV/combined.xingyi.20201231.txt","Processed_data/normalized.CNV.batch","Figures/CNV.threshold.adj.pdf") 
cnv <- fread(args[1])
cnames <- colnames(cnv)
hcimeta <- data.frame("HCIid"=cnames[which(grepl("HCI",cnames))])
hcimeta$Tech <- str_split_fixed(hcimeta$HCIid,"_",5)[,5]

annoCols <- as.character(cnames[which(!grepl("CNV",cnames))])
anno <- data.frame(cnv %>% select(annoCols))

colsIwant = c("Approved symbol",as.character(hcimeta$HCIid))
cnvm <- data.frame(cnv %>% select(colsIwant))
rownames(cnvm) <- cnvm$Approved.symbol
cnvm$Approved.symbol <- NULL 
cnvmat <- as.matrix(cnvm,rownames=rownames(cnvm))


craw <- removeBatchEffect(cnvmat,batch=hcimeta$Tech)
#Bring all of the annotations back 
crawdf <- data.frame(craw)
crawdf$"Approved.symbol" = rownames(crawdf)
crawdfanno <- merge(crawdf,anno,by="Approved.symbol")



write.table(crawdfanno,args[2],sep="\t",row.names=F,quote=F)


#this is just me exporing the variants of CNVA and CNVI 
cnvmelt <- melt(cnvmat)
cnvmelt$Tech <- str_split_fixed(cnvmelt$Var2,"_",5)[,5]


p <- ggplot(cnvmelt,aes(x=value,fill=Tech))
p <- p+geom_density(alpha=.5)
p <- p+theme_classic()
p <- p+geom_vline(xintercept=.9,col="red") + geom_vline(xintercept=-.9,col="red")
p <- p+geom_vline(xintercept=.6,col="blue") + geom_vline(xintercept=-.6,col="blue")
p <- p+xlim(-3,3)
p <- p+xlab("Log R ratios for Affyimetrix and Illumina SNPchips")
p <- p+ylab("Kernel density")
#p

pdf(args[3],height=3,width=5)
print(p)
dev.off()


#dev.new()
hci027 <- cnvmelt[which(grepl("HCI027",cnvmelt$Var2)),]
p <- ggplot(hci027,aes(x=value,fill=Var2))
p <- p+geom_density(alpha=.5)
p <- p+theme_classic()
p <- p+geom_vline(xintercept=.9,col="red") + geom_vline(xintercept=-.9,col="red")
p <- p+geom_vline(xintercept=.6,col="blue") + geom_vline(xintercept=-.6,col="blue")
p <- p+xlim(-3,3)
p <- p+xlab("Log R ratios for Affyimetrix and Illumina SNPchips")
p <- p+ylab("Kernel density")
#p


cnvadj <- melt(craw)
cnvadj$Tech <- str_split_fixed(cnvadj$Var2,"_",5)[,5]


p <- ggplot(cnvadj,aes(x=value,fill=Tech))
p <- p+geom_density(alpha=.5)
p <- p+theme_classic()
p <- p+geom_vline(xintercept=.9,col="red") + geom_vline(xintercept=-.9,col="red")
p <- p+geom_vline(xintercept=.6,col="blue") + geom_vline(xintercept=-.6,col="blue")
p <- p+xlim(-3,3)
p <- p+xlab("Log R ratios for Affyimetrix and Illumina SNPchips")
p <- p+ylab("Kernel density")
#p


