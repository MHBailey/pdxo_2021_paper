library(data.table)
library(dplyr)
library(stringr)
library(viridis)
library(ggplot2)

'%!in%' <- function(x,y)!('%in%'(x,y))

std <- function(y){
    x <- na.omit(y)
    return(sd(x)/sqrt(length(x)))
}

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("This should contain two files", call.=FALSE)
}

#args <- c("Processed_data/supp2.gr50.reformatted.txt","Figures/fig58_variable_technical.pdf")
dat <- fread(args[1])

dat$Drug = ifelse(dat$Drug == "tak-228", "sapanisertib", dat$Drug)
dat$replicate = ifelse(dat$replicate == "V1" & dat$Drug == "sapanisertib", "V5", dat$replicate)
dat$replicate = ifelse(dat$replicate == "V2" & dat$Drug == "sapanisertib", "V6", dat$replicate)
dat$replicate = ifelse(dat$replicate == "V3" & dat$Drug == "sapanisertib", "V7", dat$replicate)
dat$replicate = ifelse(dat$replicate == "V4" & dat$Drug == "sapanisertib", "V8", dat$replicate)
dat$Drug = ifelse(dat$Drug == "ink 128", "sapanisertib", dat$Drug)
dat$Day0norm <- dat$Raw/dat$RawDay0

maybe <- c("HCI-010_Day_102_(2)_plate_2.0","HCI-010_Day_102_(2)_plate_1.0")
maybe <- c("")


dat$HCIid <- str_split_fixed(dat$Plate,"_",5)[,1]
datp <- dat[which(dat$Plate %!in% maybe),]
dat2 <- data.frame(datp %>% group_by(HCIid,Vol,Drug) %>% summarize("Rmean" = mean(Raw),Rse=std(Raw)))
toremove <- c("Control","HCI-017.E2","HCI-011.E2")
dat3 <- dat2[which(dat2$HCIid %!in% toremove),]

dat4 <- data.frame(datp %>% group_by(HCIid,Vol,Drug) %>% summarize("D0mean" = mean(Day0norm),"D0se"=std(Day0norm)))
dat5 <- dat4[which(dat4$HCIid %!in% toremove),]

p <- ggplot(dat3,aes(factor(Vol),Drug,size=Rmean,color=Rse))
p <- p + geom_point(alpha=0.7,pch=16)
p <- p + theme_classic()
p <- p + facet_wrap(~HCIid)
p <- p + scale_color_viridis_c()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position="bottom")
p



p <- ggplot(dat5,aes(factor(Vol),Drug,size=D0mean,color=D0se))
p <- p + geom_point(alpha=0.7,pch=16)
p <- p + theme_classic()
p <- p + facet_wrap(~HCIid)
p <- p + scale_color_viridis_c()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position="bottom")
p


pdf(args[2],height=22,width=11,useDingbats=F)
print(p)
dev.off()


