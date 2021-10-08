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

if (length(args) != 6) {
  stop("This should contain on <input.rule.transform_kat> <output.rule.compute_gr50.scores>", call.=FALSE)
}

 
#args <- c("Processed_data/supp2.gr50.reformatted.txt","Processed_data/HCI.gr50.scores.txt","Processed_data/Lines.gr50.scores.txti","Processed_data/Lines.gr50.values.txt","Figures/fig67b_growthRates.pdf","Figures/example.screeningScore.pdf")



dat <- fread(args[1],header=T,sep="\t")
colnames(dat) <- c("cell_line","drug","concentration","drugset","replicate","cell_count","cell_count__time0","cell_count__ctrl")
dat$time <- 96
look2 <- (str_split(dat$cell_line,"_"))
dat$plate <- unlist(lapply(look2,tail,n=1))
dat$HCIid <- unlist(lapply(look2,head,n=1))

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
drc_output = GRfit(dat, groupingVariables = c("cell_line","drug"))
drc_out3 <- GRfit(dat, groupingVariables = c("HCIid","drug"),cap=T) #Forces some data...We'll see what it does. 


#Calculate scores
do <-  GRgetMetrics(drc_output)
do3 <- GRgetMetrics(drc_out3)

#do <- fread(args[3])
#do3 <- fread(args[2])

#Clean up to sample 
tmp <- str_split(do$cell_line,"_")
do$HCIid <- unlist(lapply(tmp,head,n=1))

#Write out the results tables
write.table(do3,args[2],sep="\t", quote=F, row.names=F) #This is HCI
write.table(do,args[3], sep="\t", quote=F, row.names=F) #This is Line

toremove = c("Control","HCI-017.E2","HCI-011.E2")
drugs_toremove = c("")

doc <- do[which(do$HCIid %!in% toremove & do$drug %!in% drugs_toremove),]
doc3 <- do3[which(do3$HCIid %!in% toremove & do3$drug %!in% drugs_toremove),]

#2x Cell growth plot 
ccd <- data.frame(doc %>% group_by(cell_line) %>% summarize("cell2x"=mean(ctrl_cell_doublings)))
tmp <- str_split(ccd$cell_line,"_")
ccd$HCIid <- unlist(lapply(tmp,head,n=1))

ccdagg <- data.frame(ccd %>% group_by(HCIid) %>% summarize("mean"=mean(cell2x),"stddev"=sqrt(var(cell2x)),"coefvar"=sd(cell2x)/mean(cell2x),"stderr" = sd(cell2x)/sqrt(sum(cell2x))))

ordbycv <- ccdagg[order(ccdagg$coefvar),]$HCIid 
ordbygr <- ccdagg[order(ccdagg$mean),]$HCIid
ccd$ordHCIidCV <- factor(ccd$HCIid, levels=ordbycv)
ccd$ordHCIidGR <- factor(ccd$HCIid, levels=ordbygr)

#DOC3
ccd3 <- data.frame(doc3 %>% group_by(HCIid) %>% summarize("cell2x"=mean(ctrl_cell_doublings)))

p <- ggplot(ccd,aes(x=ordHCIidGR,y=cell2x),fill=cell2x)
p <- p + stat_summary(fun.data = mean_se, geom = "errorbar",color="grey")
p <- p + stat_summary(fun.y = mean, geom = "bar") 
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
p <- p + ggtitle("Variable growth")
p <- p + ylab("Estimated doubling over screen") + xlab("Model ID")
p <- p + scale_x_discrete(labels = ccd %>%
                     group_by(ordHCIidGR) %>%
                     summarize(n=length(cell2x)) %>%
                     mutate(lab = paste0(ordHCIidGR, " (n=", n, ")")) %>%
                     pull(lab))
p <- p + geom_jitter(alpha=.7,height=0)
p


pdf(args[5],useDingbats=F,height=4,width=4)
print(p)
dev.off()

#DIVE into HCI027
samp <- "HCI-027|HCI-023|HCI-015|HCI-016|HCI-003|HCI-005|HCI-002|HCI-010"
drug <- "ro4929097 "
onesamp <- doc[which(grepl(samp,doc$HCIid)),]
onedrug <- onesamp[which(grepl(drug,onesamp$drug)),c("cell_line","GR50","GR_AOC")]


onesamp3 <- doc3[which(grepl(samp,doc3$HCIid)),]
onedrug <- onesamp3[which(grepl(drug,onesamp3$drug)),c("HCIid","GR50","GR_AOC","EC50","IC50")]


#This makes sure that I have access to the values for different heatplots
dov3 <- GRgetValues(drc_out3)
write.table(dov3,args[4],sep="\t",quote=F,row.names=F)
dov3c <- dov3[which(dov3$HCIid %!in% toremove & dov3$drug %!in% drugs_toremove),]


#This little section if for alternative orders and adding code for adj GR50 for different drug "efficacies"
d <- do3

is.na(d)<-sapply(d, is.infinite)
d[is.na(d)]<-0

d$GR50tmp <- ifelse(d$GR50 == 0, 999,d$GR50)
d$GR50c <- log10(d$GR50tmp)
d$GEC50c <- log10(d$GEC50)

dd <- dcast(d, HCIid ~ drug, value.var="GR50c")
#dd$larotrectinib <- NULL
#dd$ly3039478 <- NULL
#dd$fulvestrant <- NULL
#dd$dmso <- NULL
#dd$epz011989 <- NULL
rownames(dd) <- dd$HCIid
rnames <- dd$HCIid
dd$HCIid <- NULL

#Here are some interesting caveats to deal with missing data. 
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
drugs <- unique(dov3c$drug)
samples <- unique(dov3c$HCIid)
d = "birinapant"
s = "HCI-027"
remBAD <- c("")

#####################ORDER BY AOCadj###################################
d <- doc3

for(s in samples){
    print(s)
    ns <- dov3c[which(dov3c$HCIid == s),]
    #this is a normalization step for above 1 and below -1
    ns$GRvalue2 <- ifelse(ns$GRvalue < -1, -1, ns$GRvalue)
    ns$GRvalue3 <- ifelse(ns$GRvalue2 > 3, 3, ns$GRvalue2)
    ns$FCvalue <- ns$cell_count/ns$cell_count__time0
    #Orderin the drugs
    daoc <- dcast(d, HCIid ~ drug, value.var="GR_AOC")
    daoc2 <- daoc[which(daoc$HCIid %!in% remsamps),]
    rownames(daoc2) <- daoc2$HCIid
    daoc2$HCIid <- NULL
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
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
    p <- p + theme(legend.position="top")
    p <- p + ylab("Drug concentration") + xlab("")
    p <- p + ggtitle(toupper(s))
    p <- p + scale_x_discrete(labels = nda)

    graocs <- as.data.frame(t(daoc2[which(rownames(daoc2) == s),]))
    graocs$ordDrug <- factor(rownames(graocs),levels=ndo)
    colnames(graocs) <- c("grAOC","ordDrug")

    #Add second part of the plot
    p2 <- ggplot(graocs,aes(x=ordDrug))
    p2 <- p2 + geom_point(aes(y=grAOC),color="limegreen",size=2)
    p2 <- p2 + theme_minimal()
    p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position="bottom")
    p2 <- p2 + xlab("")+ylab("GR_AOC")
    p2

    gA <- ggplotGrob(p)
    gB <- ggplotGrob(p2)

    maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
    gA$widths[2:5] <- as.list(maxWidth)
    gB$widths[2:5] <- as.list(maxWidth)


    pdf(paste("Figures/GRAOC/AOC.",s,".grAOC.pdf",sep=""),useDingbats=F,height=5,width=6)
    grid.arrange(gA, gB, ncol = 1, heights = c(5, 2))
    dev.off()

}

#d = 'birinapant'
for(d in drugs){
    print(d)
    nd <- dov3c[which(dov3c$drug == d),]
    #this is a normalization step for above 1 and below -1
    nd$GRvalue2 <- ifelse(nd$GRvalue < -1, -1, nd$GRvalue)
    nd$GRvalue3 <- ifelse(nd$GRvalue2 > 3, 3, nd$GRvalue2)
    nd$FCvalue <- nd$cell_count/nd$cell_count__time0

    ndsamp <- unique(nd$HCIid)
    nds <- ndsamp[! ndsamp %in% remBAD ]
    no <- doc3[which(doc3$drug == d & is.finite(doc3$GR_AOC)),]
    nord <- no[order(-no$GR_AOC),]$HCIid
    sd <- setdiff(nds,nord)
    sda <- paste(sd,"*",sep="")
    nso <- c(nord,sd)
    nsa <- c(nord,sda)

    nd$ordHCIid <- factor(nd$HCIid,levels=nso)
    ndx <- data.frame(nd %>% group_by(ordHCIid,concentration) %>% summarize("FCmean"=mean(FCvalue)))
    p <- ggplot(na.omit(ndx),aes(x=ordHCIid,y=factor(concentration)))
    p <- p + geom_tile(aes(fill=FCmean))
    p <- p + scale_fill_viridis(option = "magma",direction = 1,limit=c(-0.25,4.5))
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
    p <- p + ylab("Drug concentration") + xlab("")
    p <- p + ggtitle(toupper(d))
    p <- p + scale_x_discrete(labels = nsa)
    
    pdf(paste("Figures/GRAOC/",d,".aoc.pdf",sep=""),useDingbats=F,height=4,width=6)
    print(p)
    dev.off()


    p <- ggplot(na.omit(ndx),aes(x=ordHCIid,y=factor(concentration)))
    p <- p + geom_tile(aes(fill=FCmean))
    p <- p + scale_fill_viridis(option = "magma",direction = 1)
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
    p <- p + ylab("Drug concentration") + xlab("")
    p <- p + ggtitle(toupper(d))
    p <- p + scale_x_discrete(labels = nsa) 

    pdf(paste("Figures/GRAOC_indcolor/",d,".aoc.indcolor.pdf",sep=""),useDingbats=F,height=4,width=6)
    print(p)
    dev.off()
}





#####################NOW FOR A SIMPLE DRC OF GRVALUES############### 

drcSample = "HCI-015" 
drcDrg = "birinapant"
myexp <- paste(drcSample,drcDrg,collapse=" ")
thisGR50 <- do3[which(do3$HCIid == drcSample & do3$drug == drcDrg),]$GR50
thisEC50 <- do3[which(do3$HCIid == drcSample & do3$drug == drcDrg),]$EC50

p <- GRdrawDRC(drc_out3, experiments = myexp, plotly = F, min = 10^(-6), max = 10^2)
p <- p + theme_bw()
p <- p + theme(legend.position="bottom")
p <- p + geom_vline(xintercept = log10(thisGR50), linetype="dotted", color = "blue", size=1.5) 
p <- p + geom_vline(xintercept = log10(thisEC50), linetype="dashed", color = "red", size=1.5)
p 

pdf(args[6],useDingbats=F,height=4,width=4)
print(p)
dev.off()






#EXTRA NOTES 
#For all points - just exta stuff
#explist <- unique(do3$experiment[grepl(drcDrg,do3$experiment)])
#explist2 <- explist[which(!grepl("Control",explist))]
#GRdrawDRC(drc_out3, experiments = explist2, plotly = TRUE, min = 10^(-4), max = 10^2)


##########################THIS NEXT SECTION OF CODE WILL DO GR50 plots if you'd like.
#
#
#for(d in drugs){
#    print(d)
#    nd <- dov3c[which(dov3c$drug == d),]
#    #this is a normalization step for above 1 and below -1
#    nd$GRvalue2 <- ifelse(nd$GRvalue < -1, -1, nd$GRvalue)
#    nd$GRvalue3 <- ifelse(nd$GRvalue2 > 3, 3, nd$GRvalue2)
#    nd$FCvalue <- nd$cell_count/nd$cell_count__time0
#
#    ndsamp <- unique(nd$HCIid)
#    nds <- ndsamp[! ndsamp %in% remBAD ]
#    no <- doc3[which(doc3$drug == d & is.finite(doc3$GR50)),]
#    nord <- no[order(no$GR50),]$HCIid
#    sd <- setdiff(nds,nord)
#    sda <- paste(sd,"*",sep="")
#    nso <- c(nord,sd)
#    nsa <- c(nord,sda)
#
#    nd$ordHCIid <- factor(nd$HCIid,levels=nso)
#    ndx <- data.frame(nd %>% group_by(ordHCIid,concentration) %>% summarize("FCmean"=mean(FCvalue)))
#    p <- ggplot(na.omit(ndx),aes(x=ordHCIid,y=factor(concentration)))
#    p <- p + geom_tile(aes(fill=FCmean))
#    p <- p + scale_fill_viridis(option = "magma",direction = 1)
#    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
#    p <- p + ylab("Drug concentration") + xlab("")
#    p <- p + ggtitle(toupper(d))
#    p <- p + scale_x_discrete(labels = nsa)
#    pdf(paste("Figures/GR50/",d,".gr50.pdf",sep=""),useDingbats=F,height=4,width=6)
#    print(p)
#    dev.off()
#}
#
##THIS plot is now also going to show the GR50 scores that are associated
#for(s in samples){
#    print(s)
#    ns <- dov3c[which(dov3c$HCIid == s),]
#    #this is a normalization step for above 1 and below -1
#    ns$GRvalue2 <- ifelse(ns$GRvalue < -1, -1, ns$GRvalue)
#    ns$GRvalue3 <- ifelse(ns$GRvalue2 > 3, 3, ns$GRvalue2)
#    ns$FCvalue <- ns$cell_count/ns$cell_count__time0
#    #Orderin the drugs
#    no <- doc3[which(doc3$HCIid == s & is.finite(doc3$GR50)),]
#    nord <- no[order(no$GR50),]$drug
#    nsdrugs <- unique(ns$drug)
#    sd <- setdiff(nsdrugs,nord)
#    sda <- paste(sd,"*",sep="")
#    ndo <- c(nord,sd)
#    nda <- c(nord,sda)
#
#    ns$ordDrug <- factor(ns$drug,levels=ndo)
#    nsx <- data.frame(ns %>% group_by(ordDrug,concentration) %>% summarize("FCmean"=mean(FCvalue)))
#    p <- ggplot(na.omit(nsx),aes(x=ordDrug,y=factor(concentration)))
#    p <- p + geom_tile(aes(fill=FCmean))
#    p <- p + scale_fill_viridis(option="magma")
#    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
#    p <- p + theme(legend.position="top")
#    p <- p + ylab("Drug concentration") + xlab("")
#    p <- p + ggtitle(toupper(s))
#    p <- p + scale_x_discrete(labels = nda)
#    
#    gr50s <- as.data.frame(t(dd[which(rownames(dd) == s),]))
#    gr50s$ordDrug <- factor(rownames(gr50s),levels=ndo)
#    colnames(gr50s) <- c("gr50log","ordDrug")
#	
#    #Add second part of the plot
#    p2 <- ggplot(gr50s,aes(x=ordDrug))
#    p2 <- p2 + geom_point(aes(y=gr50log),color="blue",size=2)
#    p2 <- p2 + theme_minimal()
#    p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position="bottom")
#    p2 <- p2 + xlab("")+ylab("log10(GR50)")
#    p2
#
#
#    gA <- ggplotGrob(p)
#    gB <- ggplotGrob(p2)
#
#    maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
#    gA$widths[2:5] <- as.list(maxWidth)
#    gB$widths[2:5] <- as.list(maxWidth)
#
#
#    pdf(paste("Figures/GR50/",s,".gr50.pdf",sep=""),useDingbats=F,height=5,width=6)
#    grid.arrange(gA, gB, ncol = 1, heights = c(5, 2))
#
#    dev.off()
#
#}
#
##This is going to be for the adjusted with GR50adj scores to accopany these plots
#
#for(s in samples){
#    print(s)
#    ns <- dov3c[which(dov3c$HCIid == s),]
#    #this is a normalization step for above 1 and below -1
#    ns$GRvalue2 <- ifelse(ns$GRvalue < -1, -1, ns$GRvalue)
#    ns$GRvalue3 <- ifelse(ns$GRvalue2 > 3, 3, ns$GRvalue2)
#    ns$FCvalue <- ns$cell_count/ns$cell_count__time0
#    #Orderin the drugs
#    no <- data.frame(t(dds3[which(rownames(dds3) == s),]))
#    no$drug <- rownames(no)
#    colnames(no) <- c("GR50adj","drug") 
#    nona <- na.omit(no)
#    nord <- nona[order(nona$GR50adj),]$drug
#    nsdrugs <- unique(ns$drug)
#    sd <- setdiff(nsdrugs,nord)
#    sda <- paste(sd,"*",sep="")
#    ndo <- c(nord,sd)
#    nda <- c(nord,sda)
#
#    ns$ordDrug <- factor(ns$drug,levels=ndo)
#    nsx <- data.frame(ns %>% group_by(ordDrug,concentration) %>% summarize("FCmean"=mean(FCvalue)))
#    p <- ggplot(na.omit(nsx),aes(x=ordDrug,y=factor(concentration)))
#    p <- p + geom_tile(aes(fill=FCmean))
#    p <- p + scale_fill_viridis(option="magma")
#    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
#    p <- p + theme(legend.position="top")
#    p <- p + ylab("Drug concentration") + xlab("")
#    p <- p + ggtitle(toupper(s))
#    p <- p + scale_x_discrete(labels = nda)
#    
#    gr50s <- as.data.frame(t(dds3[which(rownames(dds3) == s),]))
#    gr50s$ordDrug <- factor(rownames(gr50s),levels=ndo)
#    colnames(gr50s) <- c("gr50adj","ordDrug")
#
#    #Add second part of the plot
#    p2 <- ggplot(gr50s,aes(x=ordDrug))
#    p2 <- p2 + geom_point(aes(y=gr50adj),color="limegreen",size=2)
#    p2 <- p2 + theme_minimal()
#    p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position="bottom")
#    p2 <- p2 + xlab("")+ylab("scaled(log10(GR50))")
#    p2
#
#
#    gA <- ggplotGrob(p)
#    gB <- ggplotGrob(p2)
#
#    maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
#    gA$widths[2:5] <- as.list(maxWidth)
#    gB$widths[2:5] <- as.list(maxWidth)
#
#
#    pdf(paste("Figures/GR50/",s,".gr50adj.pdf",sep=""),useDingbats=F,height=5,width=6)
#    grid.arrange(gA, gB, ncol = 1, heights = c(5, 2))
#    dev.off()
#
#}
#
#
#
####### ADJ is just an adjust way to prioritize the data 
##DIVE into HCI027
#samp <- "HCI-027|HCI-023|HCI-015|HCI-016|HCI-003|HCI-005|HCI-002|HCI-010"
#drug <- "ro4929097 "
#onesamp <- doc[which(grepl(samp,doc$HCIid)),]
#onedrug <- onesamp[which(grepl(drug,onesamp$drug)),c("cell_line","GR50","GR_AOC")]
#
#
#onesamp3 <- doc3[which(grepl(samp,doc3$HCIid)),]
#onedrug <- onesamp3[which(grepl(drug,onesamp3$drug)),c("HCIid","GR50","GR_AOC","EC50","IC50")]
#
#
##This makes sure that I have access to the values for different heatplots
#dov3 <- GRgetValues(drc_out3)
#write.table(dov3,args[4],sep="\t",quote=F,row.names=F)
#dov3c <- dov3[which(dov3$HCIid %!in% toremove & dov3$drug %!in% drugs_toremove),]
#
#
##This little section if for alternative orders and adding code for adj GR50 for different drug "efficacies"
#d <- do3
#
#is.na(d)<-sapply(d, is.infinite)
#d[is.na(d)]<-0
#
#d$GR50tmp <- ifelse(d$GR50 == 0, 999,d$GR50)
#d$GR50c <- log10(d$GR50tmp)
#d$GEC50c <- log10(d$GEC50)
#
#dd <- dcast(d, HCIid ~ drug, value.var="GR50c")
##dd$larotrectinib <- NULL
##dd$ly3039478 <- NULL
##dd$fulvestrant <- NULL
##dd$dmso <- NULL
##dd$epz011989 <- NULL
#rownames(dd) <- dd$HCIid
#rnames <- dd$HCIid
#dd$HCIid <- NULL
#
##Here are some interesting caveats to deal with missing data. 
#dd2 <- data.frame(lapply(dd,function(x) ifelse(x >2.99,NA,x)))
#dds <- data.frame(lapply(dd2, function(x) scale(x)))
##dds2 <- data.frame(lapply(dds, function(x) ifelse(is.na(x),max(x,na.rm=T)+.5,x)))
#dds2 <- dds #This is because the weird max value will bite me for the GR50adj see line above for why that would be the case
#is.na(dds2)<-sapply(dds2, is.infinite)
#
#rownames(dds2) <- rnames
#remsamps <- c("Control","HCI-011.E2","HCI-017.E2")
#dds3 <- dds2[which(rownames(dds2) %!in% remsamps),]
#samps <- rownames(dds3)
#colnames(dds3) <- colnames(dd)
#ords <- lapply(dds3,order)
#
##nord <- reorder(samps,ords$vistusertib)
#drugs <- unique(dov3c$drug)
#samples <- unique(dov3c$HCIid)
#d = "birinapant"
#s = "HCI-027"
#remBAD <- c("")
