source("makeGRM.R")
library(asreml)

grm.auto <- read.table("Deer33.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.auto <- read.table("Deer33.grm.id")  # CONTAINS ID LIST

#read in data set
antler_data<-read.table("antler_data_comp.txt", header = TRUE)
colnames(antler_data)[1]<-"ID"
#ASreml likes categories to be factors
antler_data$ID<-as.factor(antler_data$ID)
antler_data$MeCaYear<-as.factor(antler_data$MeCaYear)
antler_data$BirthYear<-as.factor(antler_data$BirthYear)
nrow(antler_data)
#all IDs in antler_data must be represented in GRM, so subset data to remove individuals not in GRM
antler_data<-subset(antler_data, ID %in% ids.auto$V2)
nrow(antler_data)
no.ids<-length(unique(antler_data$ID))
#make variable for quadratic Age effect
antler_data$MeCaAge.2<-(antler_data$MeCaAge)^2
#drop levels!
antler_data<-droplevels(antler_data)
#make inverse GRM
grminv <- makeGRM(grm.auto, ids.auto, antler_data$ID) # vector of IDs from the datasset that you use for the asreml model
attr(grminv, which = "INVERSE") <- TRUE

#loop for pairwise genetic correlation calculation

source("ASReml.EstEffects.R")

var_table_biv_all<-data.frame()
gen.cov_results<- data.frame()

for (m in c(7:16)){
  
  antler_traits_sub_ist<-as.list(colnames(antler_data)[-c(1:m, 17)])  
 
  for (n in antler_traits_sub_ist){
    
    fixed_fmla<-as.formula(paste0("cbind(",colnames(antler_data)[m],",", n,")" , "~ trait + trait:MeCaAge + trait:MeCaAge.2"))
   
    gen.cov_model <- asreml(fixed  = fixed_fmla, 
                            random = ~ us(trait):vm(ID, grminv) + idh(trait):ide(ID)+ 
                              idh(trait):BirthYear + idh(trait):MeCaYear,
                            residual   = ~ units:idh(trait, init = NA),
                            data = antler_data,
                            maxit = 100)
    
    #summary.asreml(gen.cov_model, coef=T)$coef.fixed
    var_table<-summary.asreml(gen.cov_model, coef = T)$varcomp
    var_table_biv_all<-rbind(var_table_biv_all, var_table)
    
    #constrain rA to be zero
    
    myinit.0 <- c(0, 1, 1)
    names(myinit.0) <- c("F", "U", "U")
    
    gen.cov_model.0 <- asreml(fixed  = fixed_fmla, 
                            random = ~ corgh(trait, init = myinit.0):vm(ID, grminv) + idh(trait):ide(ID) + 
                              idh(trait):BirthYear + idh(trait):MeCaYear,
                            residual = ~ units:idh(trait),
                            data = antler_data,
                            maxiter = 100)
    
    #LRT to test if genetic correlation is different from 0
    
    chi2.0<- 2*(gen.cov_model$loglik - gen.cov_model.0$loglik)
    p.0<-pchisq(chi2.0, df = 1, lower.tail = F)
    
    
    #constrain rA to be one (or very close to one)
    
    myinit.1        <- c(0.99999, 1, 1)
    names(myinit.1) <- c("F", "U", "U")
    
    gen.cov_model.1 <- asreml(fixed  = fixed_fmla, 
                              random = ~ corgh(trait, init = myinit.1):vm(ID, grminv) + idh(trait):ide(ID) + 
                                idh(trait):BirthYear + idh(trait):MeCaYear,
                              residual = ~ units:idh(trait),
                              data = antler_data,
                              maxiter = 100)
    
    
    #LRT to test if genetic correlation is different from 1
    
    chi2.1<- 2*(gen.cov_model$loglik - gen.cov_model.1$loglik)
    p.1<-pchisq(chi2.1, df = 1, lower.tail = F)
  
    #calculate genetic correlation coefficients
    ra<-vpredict(gen.cov_model, ra~V6/(sqrt(V5*V7)))$Estimate[1]
    ra.se<-vpredict(gen.cov_model, ra~V6/(sqrt(V5*V7)))$SE[1]
  
    df<-data.frame(traits=paste0(colnames(antler_data)[m],"-", n), cor=ra, 
                   SE=ra.se, p_value.0=p.0, p_value.1=p.1)
    
    gen.cov_results<-rbind(gen.cov_results, df)
     
  }
  
}


write.table(gen.cov_results, file="Genetic_correlations_antler_measures.txt", sep="\t", col.names = T, row.names = F)


gen.cov_results<-read.table("Genetic_correlations_antler_measures.txt", header = T)
head(gen.cov_results)

#split trait column into the two separate traits for heat map
gen.cov_results$traits<-as.character(gen.cov_results$traits)
trait_split<-strsplit(gen.cov_results$traits, split = "-")

#vector containing trait 1
trait1<-character()
len<-length(trait_split)

for (i in c(1:len)){
  unit<-trait_split[[i]][1]
  trait1<-append(trait1, unit)
}

#vector containing trait 2

trait2<-character()
len<-length(trait_split)

for (i in c(1:len)){
  unit<-trait_split[[i]][2]
  trait2<-append(trait2, unit)
}

library(plyr)
library(dplyr)
library(ggplot2)

gen.cov_results<-mutate(gen.cov_results, trait_1=trait1, trait_2=trait2)
head(gen.cov_results)

#make data frame with all combinations of traits in order to make symmetric matrix
gen.cov_results2<-data.frame(traits=gen.cov_results$traits, cor=gen.cov_results$cor, SE=gen.cov_results$SE,
                             p_value.0=gen.cov_results$p_value.0, p_value.1=gen.cov_results$p_value.1,
                             trait_1=gen.cov_results$trait_2, trait_2=gen.cov_results$trait_1)

gen.cov_results_double<-rbind(gen.cov_results, gen.cov_results2)

#make symmetric matrix and 
gen.cov_results_matrix<-dcast(gen.cov_results_double, trait_1~trait_2, value.var = "cor")
rownames(gen.cov_results_matrix)<-gen.cov_results_matrix$trait_1
gen.cov_results_matrix<-as.matrix(gen.cov_results_matrix[, -1])

#select upper triangle and melt for plotting
gen.cov_results_uppertri<-gen.cov_results_matrix
gen.cov_results_uppertri[lower.tri(gen.cov_results_uppertri)]<-NA
gen.cov_results_uppertri_melt<-melt(gen.cov_results_uppertri)
colnames(gen.cov_results_uppertri_melt)<-c("trait_1", "trait_2", "cor")

#get rid off rows where the same measure in trait_1 or trait_2 only has NA correlation (trait_1 = UpperBeam, trait-2=AntlerWt)
match_trait_1<-grep("UpperBeam", gen.cov_results_uppertri_melt$trait_1)
gen.cov_results_uppertri_melt<-gen.cov_results_uppertri_melt[-match_trait_1, ]
match_trait_2<-grep("AntlerWt", gen.cov_results_uppertri_melt$trait_2)
gen.cov_results_uppertri_melt<-gen.cov_results_uppertri_melt[-match_trait_2, ]

p<-ggplot(gen.cov_results_uppertri_melt, aes(trait_1, trait_2))+geom_raster(aes(fill=cor))+
  scale_fill_gradientn(colours = c("orange", "white", "blue"), na.value="white", name = parse(text=paste("r^2")))+
  theme_classic()+
  scale_y_discrete(labels=c("Weight", "Brow Length", "Coronet-Brow-Junction", "Coronet Circumference",  "Coronet-Tray-Junction",
                            "Form", "Length", "Lower Beam", "Tray Length","Upper Beam"), expand = c(0.08, 0.08))+
  scale_x_discrete(labels=c("Weight", "Brow Length", "Coronet-Brow-Junction", "Coronet Circumference",  "Coronet-Tray-Junction",
                            "Form", "Length", "Lower Beam", "Tray Length","Upper Beam"),expand = c(0.08, 0.08), position = "top")+
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black", angle=90, 
                                   hjust=0, vjust=1), 
        axis.text.y = element_text(size=20, colour = "black", vjust=0.4),
        axis.line = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))+ylab(NULL)+xlab(NULL)
  
  
gen.cov_results_double_non.sig<-subset(gen.cov_results_double, p_value.0>=0.05)
gen.cov_results_uppertri_melt_label<-join(gen.cov_results_double_non.sig[, c("trait_1", "trait_2")], 
                                          gen.cov_results_uppertri_melt[, c("trait_1", "trait_2", "cor")])

gen.cov_results_uppertri_melt_label<-na.omit(gen.cov_results_uppertri_melt_label)

p2<-p+geom_text(data=gen.cov_results_uppertri_melt_label, label="x", size=6)

ggsave("genetic_correlations_antler_measures_heatmap.png", p2, scale = 1, width = 25, height=20, units="cm")

#save correlation matrix
write.table(gen.cov_results_matrix, file = "Genetic_correlations_antler_measures_matrix.txt", sep = "\t",
            col.names = T, row.names = T)
