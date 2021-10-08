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

if (length(args) != 5) {
  stop("This should contain 5 files, check the snakemake file", call.=FALSE)
}


#args <- c("Processed_data/supp2.gr50.reformatted.txt",'Processed_data/HCIday.GI50.scores.txt')




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

#I want to tackle this with with variability for all the things 
dat$HCIid_tech <- paste(dat$cell_line,dat$replicate,sep="_")
drc_output = GRfit(dat, groupingVariables = c("HCIid_tech","drug"))
drc_bio = GRfit(dat,groupingVariables=c("cell_line","drug"))

#Calculate scores
do <-  GRgetMetrics(drc_output)
do_bio <- GRgetMetrics(drc_bio)
do$HCIid <- str_split_fixed(do$HCIid_tech,"_",4)[,1]
do_bio$HCIid <- str_split_fixed(do_bio$cell_line,"_",4)[,1]

#So I'm going to branch this off and go from there 
#dat$GIadj <- 100 × (T - T0)/(C - T0) = 50
dat$GIadj <- 100 * (dat$cell_count - dat$cell_count__time0)/(dat$cell_count__ctrl-dat$cell_count__time0)


#NOTE: mind the s's and t's here. 

f <- function(p) as.numeric(PR(mod2,p))
ecstat <- function(dat,mod,plate_samp) {

    beta = coef(mod)[1] #Beta
    bottom = coef(mod)[2] #Bottom
    top = coef(mod)[3] #Top
    flex = coef(mod)[4] #Flex

    ecstatic = NULL
    modelled = NULL
    #Calculate and chech ECSTATIC

    if(bottom < 50 && PR(mod,100000) < 50 && PR(mod,6.25e-07) > 50){ #this is to account for max and min of these 
        z <- function(x, y, mod) y - predict(mod, data.frame(conc = x))[1]
        ecstatic <- as.numeric(unlist(uniroot(z, c(0, 100000), y = 50, mod)[1])) #this is the score to pull this off 50 is the GR50 mark
        if(ecstatic < 100){
            modelled = "Modelled"
        }
        else{
            modelled = "Predicted"
        }
    }
    else{
        modelled = "NonConverged"
        ecstatic = -9
    }

    #Model measures
    GoF = cor(dat$GIadj,predict(mod)) #Goodness of fit
    e2 = dat$GIadj-predict(mod) #resids
    nlm_error = sqrt(mean(e2^2))

    #Plot the figures. 
    pdf(paste("Figures/GI50s/",e2flag,"_",s,"_",d,".pdf",sep=""),height=4,width=4,useDingbats=F)
    plot(mod,type="all")
    #lines(dat$Vol,predict(mod),col="red",lty=2,lwd=3)
    abline(h=1,v=ecstatic)
    dev.off()

    o = NULL
    #Keep track of the data
    if(modelled == "NonConverged"){ #Don't integrate with NA ecstatic
        o = data.frame("Plate" = plate_samp, "drug" = d, "drugset" = rp, "ECstatic" = ecstatic, "respInt"=NA, "nonrInt" = NA, "nonrUnd" = NA, "GoodFit" = GoF, "Error"=nlm_error, "Beta"=beta, "Flex"=flex, "Bottom"=bottom, "Top"=top, "Modelled" = modelled)
    }
    else if(modelled == "Predicted"){
        o = data.frame("Plate" = plate_samp, "drug" = d, "drugset" = rp, "ECstatic" = ecstatic, "respInt"=NA, "nonrInt" = NA, "nonrUnd" = NA, "GoodFit" = GoF, "Error"=nlm_error, "Beta"=beta, "Flex"=flex, "Bottom"=bottom, "Top"=top, "Modelled" = modelled)
    }
    else{
        z <- function(x, y, mod2) y - predict(mod2, data.frame(conc = x))[1]
        mod2 <- checkMod(dat$GIadj, dat$Con)
        if(is.null(mod2)){
            o = data.frame("Plate" = plate_samp, "drug" = d, "drugset" = rp, "ECstatic" = ecstatic, "respInt"= NA, "nonrInt" = NA, "nonrUnd" = NA, "GoodFit" = GoF, "Error"=nlm_error, "Beta"=beta, "Flex"=flex, "Bottom"=bottom, "Top"=top, "Modelled" = modelled)
        }
        else if(PR(mod2,19) < 1 && PR(mod2,0) > 1){
            ecstatic2 <- as.numeric(unlist(uniroot(z, c(0, 10000), y = 1, mod2)[1]))
            f <- function(p) as.numeric(PR(mod2,p))
            n = integrate(f, min(dat$Con), ecstatic2)$value - ecstatic2
            r = max(dat$Con)-ecstatic2 - integrate(f, ecstatic2, max(dat$Con))$value
            u = integrate(f, ecstatic2, max(dat$Con))$value
            o = data.frame("Plate" = plate_samp, "drug" = d, "drugset" = rp, "ECstatic" = ecstatic, "respInt"=r, "nonrInt" = n, "nonrUnd" = u, "GoodFit" = GoF, "Error"=nlm_error, "Beta"=beta, "Flex"=flex, "Bottom"=bottom, "Top"=top, "Modelled" = modelled)
        }
        else{
            o = data.frame("Plate" = plate_samp, "drug" = d, "drugset" = rp, "ECstatic" = ecstatic, "respInt"= NA, "nonrInt" = NA, "nonrUnd" = NA, "GoodFit" = GoF, "Error"=nlm_error, "Beta"=beta, "Flex"=flex, "Bottom"=bottom, "Top"=top, "Modelled" = modelled)
        }
    }
    return(o)
}



checkMod <- function(Value, Vol) {
    mod <- tryCatch(
        {
           out <- NULL
           out <- drm(Value ~ Vol, fct=LL.4())
        }, error = function(e) {
           print(paste("MY_ERROR: DID NOT CONVERGE with drc", e))
           return(NULL)
        }, finally = {
           return(out)
        }
    )
    return(mod)
}


#Capture the run
e2flag <- "GI"

#So those are just a couple of old functions, and I'll see what else I can pull out to find the 50 line, not the ECstatic line. 
d0 <- dat[which(!grepl("Control",dat$cell_line)),] 
d0$HCIid_tech <- paste(d0$cell_line,d0$replicate,sep="_")
samples = unique(d0$cell_line)
techs = unique(d0$HCIid_tech)
#samples() 
ECSTATIC = NULL


for(s in samples){ #sample and plate
    sp = d0[which(d0$cell_line == s),]
    drugs = unique(sp$drug)
    #drugs = "dmso"
    #d = drugs
    for(d in drugs){#drug
        thisdrug = NULL
        dd <- sp[which(sp$drug == d),]
        drugset = unique(dd$drugset)
        for(rp in drugset){
            ddd <- dd[which(dd$drugset == rp),]
            print(c(s,d))
            #mymin <- min(ddd[which(ddd$Value > 0),]$Value)
            #mymin <- ifelse(mymin < 0.05, mymin, 0.05) #set min
            #ddd$Value <- ifelse(ddd$Value < 0, mymin, ddd$Value)#This will remove negative values... 
            print(paste(ddd$GIadj,collapse=","))
            print(paste(ddd$concentration,collapse=","))
            ddd$Con = log(ddd$concentration) + (-1*min(log(ddd$concentration)))
            mod <- checkMod(ddd$GIadj, ddd$concentration)
            if(is.null(mod)){
                print("HERE")
                thisdrug = data.frame("Plate" = s, "drug" = d, "drugset" = rp, "ECstatic"=NA, "respInt"=NA, "nonrInt" = NA, "nonrUnd" = NA, "GoodFit" = NA, "Error"=NA, "Beta"=NA, "Flex"=NA, "Bottom"=NA, "Top"=NA, "Modelled" = NA)
                print("HERE2")
                    pdf(paste("Figures/GI50s/",e2flag,"_",s,"_",d,".pdf",sep=""),height=4,width=4,useDingbats=F)
                    plot(ddd$Con,ddd$GIadj)
                    abline(h=50)
                    dev.off()
            }else{
                print("HERE3)")
                thisdrug <- ecstat(ddd,mod,s)
                print("HERE4)")
            }
            ECSTATIC <- rbind(ECSTATIC,thisdrug)
        }
    }
}


ECSTATIC2 = NULL
for(t in techs){ #sample and plate
    sp = d0[which(d0$HCIid_tech == t),]
    drugs = unique(sp$drug)
    #drugs = "dmso"
    #d = drugs
    for(d in drugs){#drug
        thisdrug = NULL
        dd <- sp[which(sp$drug == d),]
        drugset = unique(dd$drugset)
        #rp = "Retro"
        for(rp in drugset){
            ddd <- dd[which(dd$drugset == rp),]
            print(c(s,d))
            #mymin <- min(ddd[which(ddd$Value > 0),]$Value)
            #mymin <- ifelse(mymin < 0.05, mymin, 0.05) #set min
            #ddd$Value <- ifelse(ddd$Value < 0, mymin, ddd$Value)#This will remove negative values... 
            print(paste(ddd$GIadj,collapse=","))
            print(paste(ddd$concentration,collapse=","))
            ddd$Con = log(ddd$concentration) + (-1*min(log(ddd$concentration)))
            mod <- checkMod(ddd$GIadj, ddd$concentration)
            if(is.null(mod)){
                thisdrug = data.frame("Plate" = t, "drug" = d, "drugset" = rp, "ECstatic"=NA, "respInt"=NA, "nonrInt" = NA, "nonrUnd" = NA, "GoodFit" = NA, "Error"=NA, "Beta"=NA, "Flex"=NA, "Bottom"=NA, "Top"=NA, "Modelled" = NA)
                    pdf(paste("Figures/GI50s/",e2flag,"_",s,"_",d,".pdf",sep=""),height=4,width=4,useDingbats=F)
                    plot(ddd$Con,ddd$GIadj)
                    abline(h=50)
                    dev.off()
            }else{
                thisdrug <- ecstat(ddd,mod,t)
            }
            ECSTATIC2 <- rbind(ECSTATIC2,thisdrug)
        }
    }
}

#So the ECSTATICS and ECSTATIC2 
#I'm just going to write this files and work on this project in a different window. 

#GI50 scores that I calculated 
write.table(ECSTATIC,args[2],sep="\t",row.names=F,quote=F)
write.table(ECSTATIC2,args[3],sep="\t",row.names=F,quote=F)

#GR50 scores for all replicates and bio remplicates 
write.table(do,args[4],sep="\t",row.names=F,quote=F)
write.table(do_bio,args[5],sep="\t",row.names=F,quote=F)


