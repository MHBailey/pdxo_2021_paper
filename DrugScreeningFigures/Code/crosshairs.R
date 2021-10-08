library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(reshape2)
library(stringr)

'%!in%' <- function(x,y)!('%in%'(x,y))

std <- function(y){
    x <- na.omit(y)
    return(sd(x)/sqrt(length(x)))
}


args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  stop("This should contain 5 files, check out snakefile", call.=FALSE)
}


gibio <- fread(args[1])
grbio <- fread(args[2])


#gigr <- merge(grall,giall)
ggbio <- merge(grbio,gibio,by.x=c("cell_line","drug"),by.y=c("Plate","drug"),all.x=T)
#Convert INF to NA 
ggbio2 <- do.call(data.frame,lapply(ggbio, function(x) replace(x, is.infinite(x),NA)))
badsamps <- c("Control","HCI-017.E2","HCI-011.E2")
ggbio3 <- ggbio2[which(ggbio2$HCIid %!in% badsamps),]
ggbio3$ECstatic <- ifelse(ggbio3$ECstatic==-9,NA,ggbio3$ECstatic)

#To remove outlier samples replicates 
toremove = c("HCI-010_Day_102_(2)_plate_2.0","HCI-010_Day_102_(2)_plate_1.0","HCI-024_Day_202_(1)_plate_2.0","HCI-024_Day_202_(1)_plate_1.0","HCI-016_Day_327_(2)_plate_2.0")
#toremove = ""
ggbio4 <- ggbio3[which(ggbio3$cell_line %!in% toremove),]

ggstat <- ggbio4 %>% group_by(HCIid,drug) %>% summarize("Mean_GR50"=mean(log10(GR50),na.rm=T),"Mean_IC50"=mean(log10(IC50),na.rm=T),"Mean_GI50"=mean(log10(ECstatic),na.rm=T),"SE_GR50"=std(log10(GR50)),"SE_IC50"=std(log10(IC50)),"SE_GI50"=std(log10(ECstatic)))


#This is to calculate all of the missing values  
ggstat2 <- ggbio3 %>% group_by(HCIid,drug) %>% summarize("Mean_GR50"=mean(log10(GR50)),"Mean_IC50"=mean(log10(IC50)),"Mean_GI50"=mean(log10(ECstatic)),"SE_GR50"=std(log10(GR50)),"SE_IC50"=std(log10(IC50)),"SE_GI50"=std(log10(ECstatic)))



#REMOVE DMSO
ggdmso <- ggstat2[which(ggstat2$drug != "dmso"),]

#Get rid of INF
#ggclean <- do.call(data.frame,lapply(ggstat, function(x) replace(x, is.infinite(x),NA)))
ggclean <- do.call(data.frame,lapply(ggdmso, function(x) replace(x, is.infinite(x),NA)))


#GR GI 
plot(ggclean$Mean_GR50,ggclean$Mean_GI50,pch=16)

nagg <- na.omit(ggclean)
gigr_cortest <- cor.test(ggclean$Mean_GI50,ggclean$Mean_GR50)

sum(is.na(ggclean$Mean_GR50))
#[1] 350
sum(is.na(ggclean$Mean_IC50))
#[1] 216
sum(is.na(ggclean$Mean_GI50))
#[1] 271


p <- ggplot(ggclean,aes(x=Mean_GR50,y=Mean_GI50))
p <- p + geom_point()
p <- p + geom_errorbarh(aes(xmin=Mean_GR50-SE_GR50,xmax=Mean_GR50+SE_GR50),col="blue")
p <- p + geom_errorbar(aes(ymin=Mean_GI50-SE_GI50,ymax=Mean_GI50+SE_GI50),col="red")
p <- p + theme_bw()
#p <- p + ggtitle("Correlation of GR50 and GI50 by biological replicates")
p <- p + geom_abline(intercept=0,slope=1,col="orange")
p <- p + annotate("text", x = -2, y = 2, label = paste("Pearson corr coef: ",as.character(round(gigr_cortest$estimate,3)),sep=""))
p 

pdf(args[3],height=5,width=5)
print(p)
dev.off()


#This is looking at the biggest differences between samples 
ggclean$Direction <- ifelse(ggclean$Mean_GR50 < ggclean$Mean_GI50, "Inflated","Deflated")


yo <- lm(ggclean$Mean_GI50~ggclean$Mean_GR50)
yo2 <- data.frame(Yo=yo$residuals)
gg <- merge(ggclean,yo2,by="row.names")

pdf("Figures/grgi50.density.residutals.pdf",height=4,width=4)
plot(density(gg$Yo))
dev.off()

faster <- c("HCI-002","HCI-012","HCI-017","HCI-010","HCI-019","HCI-003","HCI-025","HCI-016")
slower <- c("HCI-023","HCI-001","HCI-024","HCI-008","HCI-011","HCI-015","HCI-027","HCI-005")

gg2 <- gg[which(gg$Yo > 0.5 | gg$Yo < -0.5),]
tgg <- data.frame(table(gg2$HCIid))
tggs <- tgg[which(tgg$Var1 %in% slower),]
tggf <- tgg[which(tgg$Var1 %in% faster),]

gg2$Group <- ifelse(gg2$HCIid %in% slower,"Slower",NA)
gg2$Group <- ifelse(gg2$HCIid %in% faster,"Faster",gg2$Group)

gg2 <- gg2[which(!is.na(gg2$Group)),]
gg2$Freq <- 1

p <- ggplot(gg2,aes(x=Group,y=Freq,fill=Direction))
p <- p + geom_bar(stat="identity")
p <- p + theme_classic()
p

pdf("Figures/CountsOfBiggestDifferences.pdf",height=4,width=4,useDingbats=F)
print(p)
dev.off()


#This shows that the biggest difference between GR50 and GI50 are not attributable to inflation of slow growing models. 


