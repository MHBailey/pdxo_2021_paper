library(data.table)
library(dplyr)
library(reshape2)
library(GRmetrics)
library(stringr)
library(ggplot2)
library(viridis)
library(drc)
library(grid)
library(gridExtra)


'%!in%' <- function(x,y)!('%in%'(x,y))

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  stop("This should contain on <input.rule.transform_kat> <output.rule.compute_gr50.scores>", call.=FALSE)
}

 
#args <- c("Processed_data/supp2.gr50.reformatted.txt","Processed_data/HCIday.gr50.scores.txt","Figures/Lines/growthRates.pdf","Figures/Lines/example.screeningScore.pdf")



dat <- fread(args[1],header=T,sep="\t")
colnames(dat) <- c("cell_line","drug","concentration","drugset","replicate","cell_count","cell_count__time0","cell_count__ctrl")
dat$time <- 96
look2 <- (str_split(dat$cell_line,"_"))
tmpday <- str_split(dat$cell_line,"_plate_")
dat$plate <- unlist(lapply(look2,tail,n=1))
dat$HCIid <- unlist(lapply(look2,head,n=1))
dat$HCIday <- unlist(lapply(tmpday,head,n=1))

#So this is the change that I need to make 
dat$drug = ifelse(dat$drug == "tak-228", "sapanisertib", dat$drug)
dat$replicate = ifelse(dat$replicate == "V1" & dat$drug == "sapanisertib", "V5", dat$replicate)
dat$replicate = ifelse(dat$replicate == "V2" & dat$drug == "sapanisertib", "V6", dat$replicate)
dat$replicate = ifelse(dat$replicate == "V3" & dat$drug == "sapanisertib", "V7", dat$replicate)
dat$replicate = ifelse(dat$replicate == "V4" & dat$drug == "sapanisertib", "V8", dat$replicate)
dat$drug = ifelse(dat$drug == "ink 128", "sapanisertib", dat$drug)


#DATA PROVENANCE
dim(dat[which(dat$HCIid != "Control" & dat$drug != "dmso"),])
length(unique(dat$HCIid))
length(unique(dat$drug))


#DRCurve generation 
drc_output = GRfit(dat, groupingVariables = c("HCIday","drug"))

#Calculate scores
do <-  GRgetMetrics(drc_output)

#Clean up to sample 
tmp <- str_split(do$HCIday,"_")
do$HCIid <- unlist(lapply(tmp,head,n=1))


#Write out the results tables I don't need to repeat this step
write.table(do,args[2], sep="\t", quote=F, row.names=F) #This is Line

toremove = c("Control","HCI-017.E2","HCI-011.E2")
doc <- do[which(do$HCIid %!in% toremove),]

#2x Cell growth plot 
ccd <- data.frame(doc %>% group_by(HCIday) %>% summarize("cell2x"=mean(ctrl_cell_doublings)))
tmp <- str_split(ccd$HCIday,"_")
ccd$HCIid <- unlist(lapply(tmp,head,n=1))

ccdagg <- data.frame(ccd %>% group_by(HCIid) %>% summarize("mean"=mean(cell2x),"stddev"=sqrt(var(cell2x)),"coefvar"=sd(cell2x)/mean(cell2x),"stderr" = sd(cell2x)/sqrt(sum(cell2x))))

ordbycv <- ccdagg[order(ccdagg$coefvar),]$HCIid 
ordbygr <- ccdagg[order(ccdagg$mean),]$HCIid
ccd$ordHCIidCV <- factor(ccd$HCIid, levels=ordbycv)
ccd$ordHCIidGR <- factor(ccd$HCIid, levels=ordbygr)

#DOC
doc$cell_line = rownames(doc)
ccd3 <- data.frame(doc %>% group_by(cell_line) %>% summarize("cell2x"=mean(ctrl_cell_doublings)))

p <- ggplot(ccd,aes(x=HCIday,y=cell2x),fill=cell2x)
#p <- p + stat_summary(fun.data = mean_se, geom = "errorbar")
p <- p + stat_summary(fun.y = mean, geom = "bar") 
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
p <- p + ggtitle("Variable growth")
p <- p + ylab("Estimated doubling rate") + xlab("Model ID")

pdf(args[3],useDingbats=F,height=4,width=4)
print(p)
dev.off()

#DIVE into HCI027
samp <- "HCI-027|HCI-023|HCI-015|HCI-016|HCI-003|HCI-005|HCI-002|HCI-010"
drug <- "ro4929097 "
onesamp <- doc[which(grepl(samp,doc$HCIid)),]
onedrug <- onesamp[which(grepl(drug,onesamp$drug)),c("HCIday","GR50","GR_AOC")]


#This makes sure that I have access to the values for different heatplots
dov <- GRgetValues(drc_output)
#write.table(dov3,args[4],sep="\t",quote=F,row.names=F) r
dovc <- dov[which(dov$HCIid %!in% toremove),]


#This little section is for alternative orders and adding code for adj GR50 for different drug "efficacies"
d <- do

is.na(d)<-sapply(d, is.infinite)
d[is.na(d)]<-0

d$GR50tmp <- ifelse(d$GR50 == 0, 999,d$GR50)
d$GR50c <- log10(d$GR50tmp)

dd <- dcast(d, HCIday ~ drug, value.var="GR50c")
#dd$larotrectinib <- NULL
#dd$ly3039478 <- NULL
#dd$fulvestrant <- NULL
#dd$dmso <- NULL
#dd$epz011989 <- NULL
rownames(dd) <- dd$HCIday
rnames <- dd$HCIday
dd$HCIday <- NULL

dd2 <- data.frame(lapply(dd,function(x) ifelse(x >2.99,NA,x)))
dds <- data.frame(lapply(dd2, function(x) scale(x)))
#dds2 <- data.frame(lapply(dds, function(x) ifelse(is.na(x),max(x,na.rm=T)+.5,x)))
dds2 <- dds #This is because the weird max value will bite me for the GR50adj see line above for why that would be the case
is.na(dds2)<-sapply(dds2, is.infinite)

rownames(dds2) <- rnames
remsamps <- c("Control","HCI-011.E2","HCI-017.E2")
dds3 <- dds2[which(rownames(dds2) %!in% remsamps),]
samps <- rownames(dds3)
colnames(dds3) <- colnames(dd)
ords <- lapply(dds3,order)

#nord <- reorder(samps,ords$vistusertib)


##This is going to plot GR stats for each day (biological replicate)
drugs <- unique(dovc$drug)
samples <- unique(dovc$HCIday)
d = "ly3039478"
s = "HCI-027_Day_201_(3)"
remBAD <- c("")



#####################ORDER BY AOCadj###################################
d <- do

for(s in samples){
    print(s)
    ns <- dovc[which(dovc$HCIday == s),]
    #this is a normalization step for above 1 and below -1
    ns$GRvalue2 <- ifelse(ns$GRvalue < -1, -1, ns$GRvalue)
    ns$GRvalue3 <- ifelse(ns$GRvalue2 > 3, 3, ns$GRvalue2)
    ns$FCvalue <- ns$cell_count/ns$cell_count__time0
    #Orderin the drugs
    daoc <- dcast(d, HCIday ~ drug, value.var="GR_AOC")
    daoc2 <- daoc[which(daoc$HCIday %!in% remsamps),]
    rownames(daoc2) <- daoc2$HCIday
    daoc2$HCIday <- NULL
    no <- data.frame(t(daoc2[which(rownames(daoc2) == s),]))
    no$drug <- rownames(no)
    colnames(no) <- c("AOC","drug")
    nona <- na.omit(no)
    nord <- nona[order(-nona$AOC),]$drug
    nsdrugs <- unique(ns$drug)
    sd <- setdiff(nsdrugs,nord)
    sda <- paste(sd,"*",sep="")
    ndo <- c(nord,sd)
    nda <- c(nord,sda)


    ns$ordDrug <- factor(ns$drug,levels=ndo)
    nsx <- data.frame(ns %>% group_by(ordDrug,concentration) %>% summarize("FCmean"=mean(FCvalue)))
    
    p <- ggplot(na.omit(nsx),aes(x=ordDrug,y=factor(concentration)))
    p <- p + geom_tile(aes(fill=FCmean))
    p <- p + scale_fill_viridis(option="magma")
    #p <- p + scale_fill_gradient2(low = "#FDE725", mid = "white", high = "#440154", midpoint = 1,limits=c(0,4))
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
    p <- p + theme(legend.position="top")
    p <- p + ylab("Drug concentration") + xlab("")
    p <- p + ggtitle(toupper(s))
    #p <- p + scale_x_discrete(labels = nda)

    pdf(paste("Figures/Lines/OutputAOC/",s,".grAOC.pdf",sep=""),useDingbats=F,height=5,width=6)
    print(p)
    dev.off()

}
#HERE BECAUSE OF THE NATURE OF THIS STUDY I NEED TO MAKE AN ALL REPLICATES or one with the proper scaling  by samples 
dovc$HCIid <- str_split_fixed(dovc$HCIday,"_",2)[,1]
d <- do

hcis <- unique(dovc$HCIid) 
for(s in hcis){
    print(s)
    ns <- dovc[which(dovc$HCIid == s),]
    #this is a normalization step for above 1 and below -1
    ns$GRvalue2 <- ifelse(ns$GRvalue < -1, -1, ns$GRvalue)
    ns$GRvalue3 <- ifelse(ns$GRvalue2 > 3, 3, ns$GRvalue2)
    ns$FCvalue <- ns$cell_count/ns$cell_count__time0
    #Orderin the drugs
    daoc <- dcast(d, HCIday + HCIid ~ drug, value.var="GR_AOC")
    daoc2 <- daoc[which(daoc$HCIid %!in% remsamps),]
    daoc3 <- daoc2[which(daoc2$HCIid == s),]
    rownames(daoc3) <- daoc3$HCIday
    #no <- melt(daoc3) This was to order the drugs and I don't think that we have to do that anymore 
    

    p <- ggplot(ns,aes(x=HCIday,y=factor(concentration)))
    p <- p + geom_tile(aes(fill=FCvalue))
    p <- p + scale_fill_viridis(option="magma")
    #p <- p + scale_fill_gradient2(low = "#FDE725", mid = "white", high = "#440154", midpoint = 1,limits=c(0,4))
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
    p <- p + theme(legend.position="top")
    p <- p + ylab("Drug concentration") + xlab("")
    p <- p + ggtitle(toupper(s))
    p <- p + facet_wrap(~drug,scale=)
    p

    pdf(paste("Figures/Lines/OutputBYsample/",s,".sample.pdf",sep=""),useDingbats=F,height=20,width=16)
    print(p)
    dev.off()
    
}
    
#d = 'birinapant'
for(d in drugs){
    print(d)
    nd <- dovc[which(dovc$drug == d),]
    #this is a normalization step for above 1 and below -1
    nd$GRvalue2 <- ifelse(nd$GRvalue < -1, -1, nd$GRvalue)
    nd$GRvalue3 <- ifelse(nd$GRvalue2 > 3, 3, nd$GRvalue2)
    nd$FCvalue <- nd$cell_count/nd$cell_count__time0

    ndsamp <- unique(nd$HCIday)
    nds <- ndsamp[! ndsamp %in% remBAD ]
    no <- doc[which(doc$drug == d & is.finite(doc$GR_AOC)),]
    nord <- no[order(-no$GR_AOC),]$HCIday
    sd <- setdiff(nds,nord)
    sda <- paste(sd,"*",sep="")
    nso <- c(nord,sd)
    nsa <- c(nord,sda)

    nd$ordHCIday <- factor(nd$HCIday,levels=nso)
    ndx <- data.frame(nd %>% group_by(HCIday,ordHCIday,concentration) %>% summarize("FCmean"=mean(FCvalue)))
    p <- ggplot(na.omit(ndx),aes(x=HCIday,y=factor(concentration)))
    p <- p + geom_tile(aes(fill=FCmean))
    p <- p + scale_fill_viridis(option = "magma",direction = 1)
    #p <- p + scale_fill_gradient2(low = "#FDE725", mid = "white", high = "#440154",  midpoint = 1,limits=c(0,4))
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
    p <- p + ylab("Drug concentration") + xlab("")
    p <- p + ggtitle(toupper(d))
    #p <- p + scale_x_discrete(labels = nsa)
    pdf(paste("Figures/Lines/OutputAOC/",d,".aoc.pdf",sep=""),useDingbats=F,height=4,width=6)
    print(p)
    dev.off()
}


#####################NOW FOR A SIMPLE DRC OF GRVALUES############### 

drcSample = "HCI-010_Day_65_(2)" 
drcDrg = "navitoclax"
myexp <- paste(drcSample,drcDrg,collapse=" ")
thisGR50 <- do[which(do$HCIid == drcSample & do$drug == drcDrg),]$GR50
thisEC50 <- do[which(do$HCIid == drcSample & do$drug == drcDrg),]$EC50

p <- GRdrawDRC(drc_output, experiments = myexp, plotly = F, min = 10^(-6), max = 10^2)
p <- p + theme_bw()
p <- p + theme(legend.position="bottom")
p <- p + geom_vline(xintercept = log10(thisGR50), linetype="dotted", color = "blue", size=1.5) 
p <- p + geom_vline(xintercept = log10(thisEC50), linetype="dashed", color = "red", size=1.5)
p 

pdf(args[4],useDingbats=F,height=4,width=4)
print(p)
dev.off()

#For all points
#explist <- unique(do$experiment[grepl(drcDrg,do$experiment)])
#explist2 <- explist[which(!grepl("Control",explist))]


#GRdrawDRC(drc_output, experiments = explist2, plotly = TRUE, min = 10^(-4), max = 10^2)



#################THESE BELOW ARE JUST NOTES###################################
#plop <- unique(dat[,c("drug","drugset")])
#plop[which(plop$drugset == "Pilot"),]
#
#
#dat %>% group_by(HCIid) %>% summarise(Unique_Compounds = n_distinct(drug))
#plates <- data.frame(dat %>% group_by(HCIid,drug) %>% summarise(Unique_Plates = n_distinct(replicate)))
#
#plates <- data.frame(dat %>% group_by(HCIid,drug) %>% tally(replicates)
#
#Unique_Plates = n_distinct(replicate)))
#
#	
##### LOOKING into DO3#### 
#do <- GRgetMetrics(drc_output)
#do3 <- GRgetMetrics(drc_out3)
#
#data.frame(do3a %>% group_by(drug) %>% summarise(Unique_GR50 = n_distinct(GR50))
#goodGR50 <- data.frame(do3a %>% group_by(drug) %>% summarise(Unique_GR50 = n_distinct(HCIi)))
#
#
###### LOOK into DO 
#look <- str_split(do$cell_line,"_")
#do$HCIid <- unlist(lapply(look,head,n=1))
#
#doa <- do[which(is.finite(do$GR50)),]
#goodGR50 <- data.frame(doa %>% group_by(drug) %>% summarise(Unique_GR50 = n_distinct(cell_line)))
#
#p <- ggplot(goodGR50,aes(x=reorder(drug,-Unique_GR50),y=Unique_GR50))
#p <- p + geom_bar(stat="identity")
#p <- p + theme_bw()
#p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
#p <- p + ylab("n Samples with GR50") + xlab("Compounds")
#p
#
#
##Now to compare to the data that I have here. 
#mydrc <- fread("Data/E2plus.combined_noAVG.txt")
#
#md <- mydrc[which(mydrc$Flex != "NA"),]
#goodDRC <- data.frame(md %>% group_by(Drug) %>% summarise(Unique_GR50 = n_distinct(Plate)))
#
#
#goodGR50$type = "GR50"
#goodDRC$type = "myDRC"
#colnames(goodGR50) <- c("drug","n","t")
#colnames(goodDRC) <- c("drug","n","t")
#
#
#p <- ggplot(goodDRC[which(goodDRC$drug != "dmso" ),],aes(x=reorder(drug,-n),y=n))
#p <- p + geom_bar(stat="identity")
#p <- p + theme_bw()
#p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
#p <- p + ylab("n Samples with DRC") + xlab("Compounds")
#p
#
#
#yo <- rbind(goodGR50,goodDRC)
#yo2 <- yo[which(yo$drug != "dmso"),]
#
#p <- ggplot(yo2,aes(x=drug,y=n,fill=t))
#p <- p + geom_bar(stat="identity",position = "dodge")
#p <- p + theme_bw()
#p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
#p <- p + ylab("n Samples with GR50") + xlab("Compounds")
#p
#
##playing around with everything:
#mydrug = "birinapant"
#explist <- unique(do$experiment[grepl(mydrug,do$experiment)])
#explist2 <- explist[which(!grepl("Control",explist))]
#
#GRdrawDRC(drc_output, experiments = explist2, plotly = FALSE, min = 10^(-4), max = 10^2)
#
#
#ccd <- data.frame(do %>% group_by(cell_line) %>% summarize("cell2x"=mean(ctrl_cell_doublings)))
#
#ccda <- data.frame(doa %>% group_by(HCIid) %>% summarize("cell2x"=mean(ctrl_cell_doublings)))
#
#GRdrawDRC(drc_output, experiments = explist2, plotly = FALSE, min = 10^(-4), max = 10^2)
#
#
##DO3 
#mydrug = "birinapant"
#explist <- unique(do3$experiment[grepl(mydrug,do3$experiment)])
#explist2 <- explist[which(!grepl("Control",explist))]
#
#GRdrawDRC(drc_out3, experiments = explist2, plotly = TRUE, min = 10^(-4), max = 10^2)
#
#
#ccd <- data.frame(do3 %>% group_by(cell_line) %>% summarize("cell2x"=mean(ctrl_cell_do3ublings)))
#
#ccda <- data.frame(do3a %>% group_by(HCIid) %>% summarize("cell2x"=mean(ctrl_cell_do3ublings)))
#
#GRdrawDRC(drc_output, experiments = explist2, plotly = FALSE, min = 10^(-4), max = 10^2)
#
###### SNOOPING FOR ISSUES ##### 
#
#drcSample = "HCI-010"
#drcDrg = "aslan002"
#myexp <- paste(drcSample,drcDrg,collapse=" ")
#thisdata <- do3[which(do3$HCIid == drcSample & do3$drug == drcDrg),]
#thisGR50 <- do3[which(do3$HCIid == drcSample & do3$drug == drcDrg),]$GR50
#thisEC50 <- do3[which(do3$HCIid == drcSample & do3$drug == drcDrg),]$EC50
#
#grval <- GRmetrics::GRgetValues(fitData=drc_out3)
#thisval <- grval[which(grval$experiment == myexp),]
#
#
#p <- GRdrawDRC(drc_out3, experiments = myexp, plotly = F, min = 10^(-6), max = 10^2)
#p <- p + theme_bw()
#p <- p + theme(legend.position="bottom")
#p <- p + geom_vline(xintercept = log10(thisGR50), linetype="dotted", color = "blue", size=1.5)
#p <- p + geom_vline(xintercept = log10(thisEC50), linetype="dashed", color = "red", size=1.5)
##p
#
#look <- data.frame(dat[which(dat$HCIid == drcSample & dat$drug == drcDrg ),])
#
#
#
