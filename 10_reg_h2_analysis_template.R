library(asreml)
library(dplyr)
source("makeGRM.R")
source("ASReml.EstEffects.R")

#Load antler data

antlers <- read.table("antler_data_comp.txt", header = T, stringsAsFactors = F)
#antlers <- read.table("antler_resid_PCs.txt", header = T, stringsAsFactors = F)
colnames(antlers)[1]<-"ID"
antlers$ID        <- as.factor(antlers$ID)
antlers$MeCaYear  <- as.factor(antlers$MeCaYear)
antlers$BirthYear <- as.factor(antlers$BirthYear)

antlers$ID2 <- antlers$ID
 
#make the matrix for the region
com_reg<-paste0("./gcta64 --bfile Deer31.v2 --autosome --autosome-num 29 --extract " , w , ".txt --keep idlist.txt --make-grm-gz --out ", w , "_region_GRM")
system(com_reg)
com_reg_adj<-paste0("./gcta64 --grm-gz ", w , "_region_GRM --grm-adj 0 --make-grm-gz --out ", w , "_region_GRM_adj")
system(com_reg_adj)
  
#make the matrix for rest of genome minus region
com_genome<-paste0("./gcta64 --bfile Deer31.v2 --autosome --autosome-num 29 --exclude ", w , ".txt --keep idlist.txt --make-grm-gz --out ", w , "_genome_GRM")
system(com_genome)
com_genome_adj<-paste0("./gcta64 --grm-gz ", w , "_genome_GRM --grm-adj 0 --make-grm-gz --out ", w , "_genome_GRM_adj")
system(com_genome_adj)
  
#read in the matrices and make the GRM for ASREML
  
grm.region <- read.table(paste0(w , "_region_GRM_adj.grm.gz"))
ids.region <- read.table(paste0(w , "_region_GRM_adj.grm.id"))
  
grm.genome <- read.table(paste0(w , "_genome_GRM_adj.grm.gz"))
ids.genome <- read.table(paste0(w , "_genome_GRM_adj.grm.id"))
  
region.inv <- makeGRM(grm.region, ids.region, id.vector = antlers$ID)
genome.inv <- makeGRM(grm.genome, ids.genome, id.vector = antlers$ID)
  
#remove any IDs from antlers that aren't in the GRM
  
antlers <- droplevels(subset(antlers, ID %in% ids.genome$V2))

#loop

#create list of antler traits (response variable) to do regional heritability analysis for
antler_traits<-list()
k=1

for ( i in c(7:16)) {
  unit<- colnames(antlers)[i]
  antler_traits[[k]]<-unit
  k=k+1
}

# for ( i in c(5:15)) {
#   unit<- colnames(antlers)[i]
#   antler_traits[[k]]<-unit
#   k=k+1
# }


for (trait in antler_traits){
  
  fmla=as.formula(paste0(trait, "~ I(MeCaAge)+I(MeCaAge^2)"))
  #fmla=as.formula(paste0(trait, "~ 1"))  
  
  repeat {
    
    #run a model with the genome-wide GRM
    
    model1 <- asreml(fixed = fmla,
                     random = ~ giv(ID) + ide(ID) + MeCaYear + BirthYear,
                     data = antlers,
                     ginverse =  list(ID = genome.inv),
                     na.method.X = "omit", na.omit.Y = "na.omit",
                     workspace = 500e+6, pworkspace = 500e+6)
    
    #run a model with the genome-wide and region GRMs
    
    model2 <- asreml(fixed = fmla,
                     random = ~ giv(ID) + giv(ID2) + ide(ID) + MeCaYear + BirthYear,
                     data = antlers,
                     ginverse =  list(ID = genome.inv, ID2 = region.inv),
                     na.method.X = "omit", na.omit.Y = "na.omit",
                     workspace = 500e+6, pworkspace = 500e+6)
    
    #likelihood ratio test
    
    logli1 <- summary(model1, all = T)$loglik
    logli2 <- summary(model2, all = T)$loglik
    
    if (logli1 !=0 && logli2 !=0) break;
    
  }
  
  p.value<-pchisq(q = 2*(logli2 - logli1), df = 1, lower.tail = F) #pvalue
    
  df<-as.data.frame(ASReml.EstEffects(model1))
  df$model<-"partial"
  df$model_variable<-rownames(df)
  rownames(df)<-NULL
  df$loglik<-logli1
    
  df2<-as.data.frame(ASReml.EstEffects(model2))
  df2$model<-"full"
  df2$model_variable<-rownames(df2)
  rownames(df2)<-NULL
  df2$loglik<-logli2
    
  df_full<-rbind(df, df2)
  df_full$p_value<-p.value
  df_full$Trait<-trait
  df_full$window_name<-w
  df_full<-df_full[, c(8,9,1,2,3,4,5,6,7,10,11,12,13)]
    
  saveRDS(df_full , file=paste0(trait, "_results.", w, ".RDS"))
  
}


