#-------------------------------------------------------------------------------------#
# Animal models of antler principle components using pedigree and genomic relatedness #
#-------------------------------------------------------------------------------------#
#

library(lattice)
library(asreml)
library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(data.table)
library(reshape2)
library(forcats)

#pedigree models

#read in data
antler_PCs_residuals<-read.table("antler_resid_PCs.txt", header = T)
colnames(antler_PCs_residuals)[1]<-"ID"
antler_PCs_residuals$ID<-as.factor(antler_PCs_residuals$ID)
antler_PCs_residuals$MeCaYear<-as.factor(antler_PCs_residuals$MeCaYear)
antler_PCs_residuals$BirthYear<-as.factor(antler_PCs_residuals$BirthYear)

#pedigree
#subset pedigree so it only contains ID, Dam and sire
pedigree<-read.table("Pedigree_Deer_2017-12.txt", header=T)
pedigree<- pedigree[, c(1:3)]
colnames(pedigree)[1:3]<-c("ID", "Dam","Sire")

#ASReml likes categories to be as factors

pedigree$ID <- as.factor(pedigree$ID)
pedigree$Dam <- as.factor(pedigree$Dam)
pedigree$Sire <- as.factor(pedigree$Sire)

#exclude individuals not in pedigree
antler_PCs_residuals<-subset(antler_PCs_residuals, ID %in% pedigree$ID)

#drop levels!
antler_PCs_residuals<-droplevels(antler_PCs_residuals)

#create inverse relatedness matrix

ainv <- ainverse(pedigree)

#animal model loop


PC_pedigree_h2_results <- data.frame()

for (p in c(5:15)){
  fixed_fmla<-as.formula(paste0(colnames(antler_PCs_residuals)[p],"~ 1"))
  
  model1<- asreml(fixed = fixed_fmla,
                  random = ~ vm(ID, ainv) + ide(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                  data = antler_PCs_residuals,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  var_table<-summary.asreml(model1, coef = T)$varcomp[c(1:5), ]    # variance components table
  
  #get proportion and associated std error for all variance components
  #add all of them apart from units!R 
  
  var_table_prop<-data.frame(variable=rownames(var_table), Effect=c(
    vpredict(model1, RY ~ V1/(V1+V2+V3+V4+V5))$Estimate[1],
    vpredict(model1, BY ~ V2/(V1+V2+V3+V4+V5))$Estimate[1],
    vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$Estimate[1],
    vpredict(model1, ID ~ V4/(V1+V2+V3+V4+V5))$Estimate[1],
    vpredict(model1, RS ~ V5/(V1+V2+V3+V4+V5))$Estimate[1]),
    SE=c(
      vpredict(model1, RY ~ V1/(V1+V2+V3+V4+V5))$SE[1],
      vpredict(model1, BY ~ V2/(V1+V2+V3+V4+V5))$SE[1],
      vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$SE[1],
      vpredict(model1, ID ~ V4/(V1+V2+V3+V4+V5))$SE[1],
      vpredict(model1, RS ~ V5/(V1+V2+V3+V4+V5))$SE[1]))
  
  var_table<-cbind(var_table, var_table_prop[-1])
  
  model2<- asreml(fixed = fixed_fmla,
                  random = ~ idv(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                  data = antler_PCs_residuals,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  effect<-vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$Estimate[1]
  se<-vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$SE[1]
  chi2<- 2*(model1$loglik - model2$loglik)
  p.h2<-pchisq(chi2, df = 1, lower.tail = F)
  df_h2_temp<-data.frame(Trait= colnames(antler_PCs_residuals)[p], Heritability=effect,
                         SE=se, p.value=p.h2)
  
  PC_pedigree_h2_results<-rbind(PC_pedigree_h2_results, df_h2_temp)
  write.table(var_table, file=paste0("variance_table_", colnames(antler_PCs_residuals)[p],
                                                     "_ped.txt"), sep="\t", row.names = T)
  
}


write.table(PC_pedigree_h2_results, file = "antler_h2_PCs.txt", sep = "\t", row.names = F)


antler_h2_PCs<-read.table("antler_h2_PCs.txt", header = T)

#label data frames to signify significance

#label df for PCs significant at <0.001 ***
label.df1<-data.frame(Trait=c("PC1","PC3", "PC4","PC5", "PC6", "PC7"), 
                      Heritability=c(antler_h2_PCs$Heritability[c(1, 3:4, 5:7)]+antler_h2_PCs$SE[c(1, 3:4, 5:7)]+0.02))

#label df for PCs significant at <0.01 **
label.df2<-data.frame(Trait=c("PC2","PC8", "PC11"), 
                      Heritability=c(antler_h2_PCs$Heritability[c(2, 8, 11)]+antler_h2_PCs$SE[c(2, 8, 11)]+0.02))

#label df for PCs significant at <0.05 *
label.df3<-data.frame(Trait=c("PC9","PC10"), 
                      Heritability=c(antler_h2_PCs$Heritability[c(9:10)]+antler_h2_PCs$SE[c(9:10)]+0.02))

pg<-ggplot(antler_h2_PCs, aes(fct_inorder(Trait), Heritability)) + geom_bar(aes(fill=fct_inorder(Trait)), stat = "identity") + 
  scale_fill_viridis(discrete=T, option = "C") + 
  geom_errorbar(aes(ymin=Heritability-SE, ymax=Heritability+SE)) + 
  geom_text(aes(label= sprintf("%0.4f", round(Heritability, digits = 5))), vjust=-0.25, size=4) + 
  theme_classic() + coord_cartesian(ylim = c(0, 0.65))

pg2<-pg + geom_text(data=label.df1, label="***")

pg3<-pg2 + geom_text(data=label.df2, label="**")

pg4<-pg3 + geom_text(data=label.df3, label="*")+
  theme(axis.text.x = element_text(size = 12, vjust = -0.01, colour = "black"), 
        axis.text.y = element_text(size=12, colour = "black"), axis.title = element_text(size=16), 
        legend.text = element_text(size=12))+
  xlab("Trait")+theme(legend.position="none")

ggsave("antler_h2_PCs_plot.png", pg4, scale = 1, width = 25, height=25, units="cm")


#get variance components

antler_PCs_residuals<-read.table("antler_resid_PCs.txt", header = T)
colnames(antler_PCs_residuals)[1]<-"ID"
antler_traits<-as.list(colnames(antler_PCs_residuals)[5:15])

var_table_comp.ped<-data.frame()

for (m in antler_traits) {
  var_table<-read.table(paste0("variance_table_", m, "_ped.txt"), header=T)
  var_table$measure<-paste0(m)
  var_table_comp.ped<-rbind(var_table_comp.ped, var_table)
  
}

var_table_comp.ped$variable.names<-rep(c("RutYear", "BirthYear", "genetic", "permanent", "residual"), 11)
rownames(var_table_comp.ped)<-NULL
var_table_comp.ped<-var_table_comp.ped[, c(9, 1:8)]


#get repatability for each antler measure - permanent plus genetic plus birthyear component
head(var_table_comp.ped)

var_table_rep.ped<-data.frame()

for (m in antler_traits){
  var_table_temp<-var_table_comp.ped %>%
    filter(variable.names == "genetic" | variable.names == "permanent" | variable.names=="BirthYear")%>%
    filter(measure==m )%>%
    mutate(repeatability=sum(Effect))
  
  var_table_sub<-var_table_comp.ped%>%
    filter(measure==m)%>%
    mutate(repeatability=var_table_temp[["repeatability"]][1])
  
  var_table_rep.ped<-rbind(var_table_rep.ped, var_table_sub)
}

write.table(var_table_rep.ped, file="Variance_components_h2_animal_model_PCs_ped.txt" , sep="\t",
            row.names = F, col.names = T)

antler_variables<-as.list(unique(as.character(var_table_comp.ped[["variable.names"]])))
v<-antler_variables[[1]]

#which measures have the highest/lowest contributions from which varience components?

var_table_min_max<-data.frame()

for (v in antler_variables) {
  var_table_temp_max<-var_table_comp.ped %>%
    filter(variable.names==v ) %>%
    filter(Effect==max(Effect))%>%
    mutate(maximum=Effect, measure_max=measure)%>%
    select(c(variable.names, maximum, measure_max))
  
  var_table_temp_min<-var_table_comp.ped %>%
    filter(variable.names==v ) %>%
    filter(Effect==min(Effect))%>%
    mutate(minimum=Effect, measure_min=measure)%>%
    select(c(variable.names, minimum, measure_min))
  
  var_table_temp<-join(var_table_temp_max, var_table_temp_min)
  var_table_min_max<-rbind(var_table_min_max, var_table_temp)
}



#GRM models ========================================================================================================================================================================================================================================================================================================================================

#read in data
antler_PCs_residuals<-read.table("antler_resid_PCs.txt", header = T)
colnames(antler_PCs_residuals)[1]<-"ID"
antler_PCs_residuals$ID<-as.factor(antler_PCs_residuals$ID)
antler_PCs_residuals$MeCaYear<-as.factor(antler_PCs_residuals$MeCaYear)
antler_PCs_residuals$BirthYear<-as.factor(antler_PCs_residuals$BirthYear)

grm.auto <- read.table("Deer33.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.auto <- read.table("Deer33.grm.id")  # CONTAINS ID LIST

#all IDs in antler_data must be represented in GRM, so subset data to remove individuals not in GRM
antler_PCs_residuals<-subset(antler_PCs_residuals, ID %in% ids.auto$V2)

#drop levels!
antler_PCs_residuals<-droplevels(antler_PCs_residuals)

#make inverse GRM
source("makeGRM.R")
grminv <- makeGRM(grm.auto, ids.auto, antler_PCs_residuals$ID) # vector of IDs from the datasset that you use for the asreml model
attr(grminv, which = "INVERSE") <- TRUE

#loop for grm PC animal models

PC_grm_h2_results= data.frame()

for (p in c(5:15)){
  fixed_fmla<-as.formula(paste0(colnames(antler_PCs_residuals)[p],"~ 1"))
  
  model1<- asreml(fixed = fixed_fmla,
                  random = ~ vm(ID, grminv) + ide(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                  data = antler_PCs_residuals,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  var_table<-summary.asreml(model1, coef = T)$varcomp[c(1:5), ]    # variance components table
  
  #get proportion and associated std error for all variance components
  #add all of them apart from units!R 
  
  var_table_prop<-data.frame(variable=rownames(var_table), Effect=c(
    vpredict(model1, RY ~ V1/(V1+V2+V3+V4+V5))$Estimate[1],
    vpredict(model1, BY ~ V2/(V1+V2+V3+V4+V5))$Estimate[1],
    vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$Estimate[1],
    vpredict(model1, ID ~ V4/(V1+V2+V3+V4+V5))$Estimate[1],
    vpredict(model1, RS ~ V5/(V1+V2+V3+V4+V5))$Estimate[1]),
    SE=c(
      vpredict(model1, RY ~ V1/(V1+V2+V3+V4+V5))$SE[1],
      vpredict(model1, BY ~ V2/(V1+V2+V3+V4+V5))$SE[1],
      vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$SE[1],
      vpredict(model1, ID ~ V4/(V1+V2+V3+V4+V5))$SE[1],
      vpredict(model1, RS ~ V5/(V1+V2+V3+V4+V5))$SE[1]))
  
  var_table<-cbind(var_table, var_table_prop[-1])
  
  model2<- asreml(fixed = fixed_fmla,
                  random = ~ idv(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                  data = antler_PCs_residuals,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  effect<-vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$Estimate[1]
  se<-vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$SE[1]
  chi2<- 2*(model1$loglik - model2$loglik)
  p.h2<-pchisq(chi2, df = 1, lower.tail = F)
  df_h2_temp<-data.frame(Trait= colnames(antler_PCs_residuals)[p], Heritability=effect,
                         SE=se, p.value=p.h2)
  
  PC_grm_h2_results<-rbind(PC_grm_h2_results, df_h2_temp)
  write.table(var_table, file=paste0("variance_table_", colnames(antler_PCs_residuals)[p],
                                                     "_grm.txt"), sep="\t", row.names = T)
  
}


write.table(PC_grm_h2_results, file = "antler_h2_grm_PCs.txt", sep = "\t", row.names = F)
antler_h2_grm_PCs<-read.table("antler_h2_grm_PCs.txt", header=T)

#plot heritabilities

#label data frames

#label df for PC significant at <0.001 ***
label.df1<-data.frame(Trait=c("PC1", "PC2", "PC3", "PC5", "PC6", "PC7", "PC10", "PC11"), 
                      Heritability=c(antler_h2_grm_PCs$Heritability[c(1:3, 5:7, 10:11)]+antler_h2_grm_PCs$SE[c(1:3, 5:7, 10:11)]+0.02))

#label df for PCs significant at <0.01 **
label.df2<-data.frame(Trait=c("PC4", "PC8", "PC9"), Heritability=c(antler_h2_grm_PCs$Heritability[c(4,8:9)]+antler_h2_grm_PCs$SE[c(4,8:9)]+0.02))


pg<-ggplot(antler_h2_grm_PCs, aes(fct_inorder(Trait), Heritability)) + geom_bar(aes(fill=fct_inorder(Trait)), stat = "identity") + 
  scale_fill_viridis(discrete=T, option = "C") + geom_errorbar(aes(ymin=Heritability-SE, ymax=Heritability+SE)) + 
  geom_text(aes(label= sprintf("%0.4f", round(Heritability, digits = 5))), vjust=-0.25, size=4) + 
  theme_classic() + coord_cartesian(ylim = c(0, 0.65))

pg2<-pg + geom_text(data=label.df1, label="***")
  
pg3<-pg2 + geom_text(data=label.df2, label="**")

pg4<-pg3 +theme(axis.text.x = element_text(size = 12, vjust = -0.01, colour = "black"), 
        axis.text.y = element_text(size=12, colour = "black"), axis.title = element_text(size=16), 
        legend.text = element_text(size=12))+
  xlab("Trait")+ theme(legend.position="none")


ggsave("antler_h2_grm_PCs_plot.png", pg4, scale = 1, width = 25, height=25, units="cm")


#plot variance components 

antler_data<-read.table("antler_resid_PCs.txt", header = TRUE)
colnames(antler_data)[1]<-"ID"
antler_traits<-as.list(colnames(antler_data)[5:15])

var_table_comp.grm<-data.frame()

for (m in antler_traits) {
  var_table<-read.table(paste0("variance_table_", m, "_grm.txt"), header=T)
  var_table$measure<-paste0(m)
  var_table_comp.grm<-rbind(var_table_comp.grm, var_table)
  
}

var_table_comp.grm$variable.names<-rep(c("RutYear", "BirthYear", "genetic", "permanent", "residual"), 11)
rownames(var_table_comp.grm)<-NULL
var_table_comp.grm<-var_table_comp.grm[, c(9, 1:8)]
head(var_table_comp.grm)


#get repatability for each antler measure - permanent plus genetic plus birth year component
var_table_rep.grm<-data.frame()

for (m in antler_traits){
  #get value for total phenotypic variance
  var_table.grm<-var_table_comp.grm%>%
    filter(measure==m)%>%
    mutate(total_var=sum(component))
  
  var_table_temp<-var_table_comp.grm %>%
    filter(variable.names == "genetic" | variable.names == "permanent" | variable.names=="BirthYear")%>%
    filter(measure==m )%>%
    mutate(repeatability_comp=sum(component))
  
  var_table_sub<-var_table.grm%>%
    mutate(repeatability=var_table_temp[["repeatability_comp"]][1]/total_var)
  
  var_table_rep.grm<-rbind(var_table_rep.grm, var_table_sub)
}

write.table(var_table_rep.grm, file="Variance_components_h2_animal_model_PCs_grm.txt" , sep="\t",
            row.names = F, col.names = T)
var_table_rep.grm<-read.table("Variance_components_h2_animal_model_PCs_grm.txt", header=T)

antler_variables<-as.list(unique(as.character(var_table_rep.grm[["variable.names"]])))
v<-antler_variables[[1]]


var_table_flipped_effect.grm<-dcast(var_table_rep.grm[, c(1, 7:11)],
                                    measure+repeatability+total_var~variable.names, 
                                    value.var = c("Effect"))

var_table_flipped_SE.grm<-dcast(var_table_rep.grm[, c(1, 7:9)], measure~variable.names, 
                                value.var = c("SE"))

colnames(var_table_flipped_SE.grm)[-1]<-paste0(colnames(var_table_flipped_SE.grm)[-1], "_SE")

var_table_flipped.grm<-join(var_table_flipped_effect.grm, var_table_flipped_SE.grm, by="measure")
var_table_flipped.grm<-var_table_flipped.grm[c(1, 4:11, 2,3), c("measure", "total_var", "genetic", "genetic_SE",
                                                                "permanent", "permanent_SE", "BirthYear", "BirthYear_SE",
                                                                "RutYear", "RutYear_SE","residual", "residual_SE","repeatability")]

write.table(var_table_flipped.grm, file="Variance_composition_h2_animal_model_PCs_grm.txt" , sep="\t",
            row.names = F, col.names = T)


#which measures have the highest/lowest contributions from which varient components?

var_table_min_max.grm<-data.frame()

for (v in antler_variables) {
  var_table_temp_max<-var_table_comp.grm %>%
    filter(variable.names==v ) %>%
    filter(Effect==max(Effect))%>%
    mutate(maximum=Effect, measure_max=measure)%>%
    select(c(variable.names, maximum, measure_max))
  
  var_table_temp_min<-var_table_comp.grm %>%
    filter(variable.names==v ) %>%
    filter(Effect==min(Effect))%>%
    mutate(minimum=Effect, measure_min=measure)%>%
    select(c(variable.names, minimum, measure_min))
  
  var_table_temp<-join(var_table_temp_max, var_table_temp_min)
  var_table_min_max.grm<-rbind(var_table_min_max.grm, var_table_temp)
}



#stacked plot

#calculate total variance for each PC and proportion that each component makes up of that variance

var_table_total.grm<-ddply(var_table_rep.grm, .(measure), summarise,
                           total.variance=sum(component))

var_table_comp.grm<-join(var_table_rep.grm,var_table_total.grm)

var_table_comp.grm<-mutate(var_table_comp.grm, prop.variance=component/total.variance)

head(var_table_comp.grm)

p.var.pc<-ggplot(var_table_comp.grm, aes(fct_inorder(measure), prop.variance, group=fct_inorder(variable.names), fill=fct_inorder(variable.names), width = 0.8))+
  geom_bar(stat = "identity") + 
  scale_fill_brewer(type="qual", palette = "Accent", labels = c("rut year", "birth year", "additive genetic", "permanent", "residual")) + theme_classic()+ 
  theme(axis.text.x = element_text(size = 14, colour = "black", angle = 90, vjust=0.3, hjust=1), 
        axis.text.y = element_text(size=14, colour = "black"),
        axis.title = element_text(size=18),legend.text = element_text(size=13),
        legend.title = element_blank())+xlab(NULL)+ 
  ylab("proportion of variance")+theme(strip.text.x = element_blank())

ggsave("antler_variance_components_grm_plot_PCs.png", p.var.pc, scale = 1, width = 25, height=25, units="cm")



#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#PCs comparative plot (pedigree vs, grm heritabilities)

#read in h2 data and rename heritability to identify origin of data e.g. grm or pedigree
antler_h2_ped<-read.table("antler_h2_PCs.txt", header=T)
colnames(antler_h2_ped)[2]<-c("pedigree")
antler_h2_grm<-read.table("antler_h2_grm_PCs.txt", header = T)
colnames(antler_h2_grm)[2]<-c("grm")

# join and melt heritability part of data frames 
antler_h2<-join(antler_h2_ped[,c(1:2)], antler_h2_grm[,c(1:2)], by="Trait")
antler_h2<-melt(antler_h2)
colnames(antler_h2)[3]<-"Heritability"
#join and melt standard error data
antler_SE_ped<-antler_h2_ped[, c(1,3)]
colnames(antler_SE_ped)[2]<-"pedigree"
antler_SE_grm<-antler_h2_grm[, c(1,3)]
colnames(antler_SE_grm)[2]<-"grm"
antler_SE<-join(antler_SE_ped, antler_SE_grm, by="Trait")
antler_SE<-melt(antler_SE)
colnames(antler_SE)[3]<-"SE"

#join melted data frames with heritabilities and SE
antler_h2_full<-join(antler_h2, antler_SE)

variable_levels<-c("pedigree", "grm")
antler_h2_full$variable<-factor(antler_h2_full$variable, levels = variable_levels)


p<-ggplot(antler_h2_full, aes(x = fct_inorder(Trait), y= Heritability, fill = variable)) +
  geom_bar(stat="identity", width=.8, position = "dodge")+
  geom_errorbar(aes(ymin=Heritability-SE, ymax=Heritability+SE), 
                position = position_dodge(width=0.8))+
  scale_fill_brewer(type="qual", palette = "Accent")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 20, colour = "black", angle = 90, vjust=0.6, hjust=1), 
        axis.text.y = element_text(size=16, colour = "black"), 
        axis.title = element_text(size=20),legend.text = element_text(size=16),
        legend.title = element_blank()) +xlab("Trait")
  
p

antler_h2_ped<-read.table("antler_h2_PCs.txt", header=T)

antler_h2_grm<-read.table("antler_h2_grm_PCs.txt", header = T)

label.df<-data.frame(Trait=antler_h2_ped$Trait,
                     variable="pedigree",
                     Heritability=c(antler_h2_ped$Heritability+antler_h2_ped$SE+0.02))


label.df1<-label.df[c(1, 3:7), ] #***
label.df1.na<-label.df[c(2, 8:11), ]
#label.df1.na$variable<-NA
label.df1.na$Heritability<-NA
label.df1.full<-rbind(label.df1, label.df1.na)

label.df1.1<-label.df[c(2, 8, 11), ] #**
label.df1.1.na<-label.df[c(1, 3:7, 9:10), ]
#label.df1.1.na$variable<-NA
label.df1.1.na$Heritability<-NA
label.df1.1.full<-rbind(label.df1.1, label.df1.1.na)

label.df1.2<-label.df[c(9:10), ] #*
label.df1.2.na<-label.df[c(1:9, 11), ]
#label.df1.2.na$variable<-NA
label.df1.2.na$Heritability<-NA
label.df1.2.full<-rbind(label.df1.2, label.df1.2.na)

label.df2<-data.frame(Trait=antler_h2_grm$Trait,
                     variable="grm",
                     Heritability=c(antler_h2_grm$Heritability+antler_h2_grm$SE+0.02))


label.df2.1<-label.df2[c(1:3, 5:7, 10:11), ] #***
label.df2.1.na<-label.df2[c(4, 8:9), ]
#label.df2.1.na$variable<-NA
label.df2.1.na$Heritability<-NA
label.df2.1.full<-rbind(label.df2.1, label.df2.1.na)

label.df2.2<-label.df2[c(4, 8:9), ] #**
label.df2.2.na<-label.df2[c(1:3, 5:7, 10:11), ]
#label.df2.2.na$variable<-NA
label.df2.2.na$Heritability<-NA
label.df2.2.full<-rbind(label.df2.2, label.df2.2.na)


p1<-p+geom_text(data=label.df2.1.full, label="***", hjust=-0.3) #

p2<-p1+geom_text(data=label.df2.2.full, label="**", hjust=-1)

p3<-p2 + geom_text(data=label.df1.full, label="***", hjust=1.4, stat = "identity")

p4<-p3 + geom_text(data=label.df1.1.full, label="**", hjust=1.8, stat = "identity")

p5<-p4 + geom_text(data=label.df1.2.full, label="*", hjust=3.4, stat = "identity")



ggsave("antler_h2_biplot_PCs.png", p5, scale = 1, width = 27, height=25, units="cm")


#make table with all h2 estimates from pedigree and GRM for both antler measures and antler PCs

#antler measure heritabilities
antler_h2_ped<-read.table("antler_h2_ped.txt", header=T)
antler_h2_grm<-read.table("antler_h2_grm.txt", header=T)

names(antler_h2_ped)[names(antler_h2_ped)=="Heritability"]<-"Pedigree"
names(antler_h2_grm)[names(antler_h2_grm)=="Heritability"]<-"GRM"

antler_h2<-cbind(antler_h2_ped, antler_h2_grm[, -1])

#antler PCs heritabilities
antler_h2_PC_ped<-read.table("antler_h2_PCs.txt", header=T)
antler_h2_PC_grm<-read.table("antler_h2_grm_PCs.txt", header=T)

names(antler_h2_PC_ped)[names(antler_h2_PC_ped)=="Heritability"]<-"Pedigree"
names(antler_h2_PC_grm)[names(antler_h2_PC_grm)=="Heritability"]<-"GRM"

antler_h2_PCs<-cbind(antler_h2_PC_ped, antler_h2_PC_grm[, -1])

antler_h2_all<-rbind(antler_h2, antler_h2_PCs)

write.table(antler_h2_all, file="antler_h2_measures_and_PCs_ped_grm.txt", sep="\t", 
            row.names=F, col.names=T, quote=F)

