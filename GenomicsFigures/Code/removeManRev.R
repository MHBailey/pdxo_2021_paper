library(data.table)

#this is source activate stats 
uniqString <- function(x){
    yo <- unique(unlist(strsplit(x, ";")))
    return(paste(yo,collapse=";"))
}
'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Check out the rule restict mutations", call.=FALSE)
}

#Indepth man review removal 
#args = c('Processed_data/Indepth_2020_04_22.all.HG38.PASS.indels.cancer_miss.maf','Processed_data/ManReview_Indepth_False.1-16.txt','Processed_data/Indepth_2020_04_22.all.HG38.PASS.indels.cancer_miss_mr.maf')


maf <- fread(args[1])
manr <- fread(args[2])


maf$KEY <- do.call(paste0, maf[,1:16])
manr$KEY <- do.call(paste0, manr[,1:16])


maf2 <- maf[which(maf$KEY %!in% manr$KEY),]
maf2$KEY <- NULL 

write.table(maf2,args[3],sep="\t",row.names=F,quote=F)

