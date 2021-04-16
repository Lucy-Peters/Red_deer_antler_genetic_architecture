#-------------------------------------------------------------------------#
# Animal models of antler measures using pedigree and genomic relatedness #
#-------------------------------------------------------------------------#

library(lattice)
library(asreml)
library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(data.table)
library(reshape2)
library(forcats)

#read in data

antler_data1<-read.table("antler_data_comp.txt", header = TRUE)
colnames(antler_data1)[1]<-"ID"
#pedigree analysis
#subset pedigree so it only contains ID, Dam and sire
pedigree<-read.table("Pedigree_Deer_2017-12.txt", header=T, stringsAsFactors=FALSE)
pedigree<- pedigree[, c(1:3)]
colnames(pedigree)[1:3]<-c("ID", "Dam","Sire")

#ASReml likes categories to be as factors
pedigree$ID <- as.factor(pedigree$ID)
pedigree$Dam <- as.factor(pedigree$Dam)
pedigree$Sire <- as.factor(pedigree$Sire)

#join pedigree information and antler data
antler_data1<- join(antler_data1, pedigree)
antler_data1$ID<-as.factor(antler_data1$ID)
antler_data1$MeCaYear<-as.factor(antler_data1$MeCaYear)
antler_data1$BirthYear<-as.factor(antler_data1$BirthYear)
antler_data1$MeCaAge<-as.numeric(antler_data1$MeCaAge)

#all IDs in antler_data must be in the pedigree, so subset data to remove individuals not in pedigree
antler_data1<-subset(antler_data1, ID %in% pedigree$ID)
nrow(antler_data1)
#drop levels!
antler_data1<-droplevels(antler_data1)

#create inverse relatedness matrix
ainv <- ainverse(pedigree)

#read in antler summaries table to mean standardise effect sizes of fixed effects
antler_summaries<-read.table("antler_summaries.txt", header=T)


#animal model loop

pedigree_h2_results<- data.frame()
pedigree_fixed_effects_results<-data.frame()

#make variable for quadratic Age effect
antler_data1$MeCaAge.2<-(antler_data1$MeCaAge)^2

for (m in c(7:16)){
  fixed_fmla<-as.formula(paste0(colnames(antler_data1)[m],"~ MeCaAge+MeCaAge.2"))
  model1<- asreml(fixed = fixed_fmla,
                  random = ~ vm(ID, ainv) + ide(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                  data = antler_data1,
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
                  data = antler_data1,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  #fixed effects full model
  fixed_effects_temp<-as.data.frame(summary.asreml(model1, coef = T)$coef.fixed)
  
  wald.test<-as.data.frame(wald.asreml(model1))
  wald_test_p_values<-na.omit(data.frame(variable=rownames(wald.test), 
                                         p_value=wald.test$`Pr(Chisq)`))
  antler_summaries_sub<-subset(antler_summaries, Trait==colnames(antler_data1)[m])
  
  fixed_effects<-data.frame(variable = rownames(fixed_effects_temp),
                            effect_size = fixed_effects_temp$solution,
                            std.error=fixed_effects_temp$`std error`, 
                            Z_ratio=fixed_effects_temp$z.ratio,
                            Trait = colnames(antler_data1)[m])
  
  fixed_effects<-join(fixed_effects, wald_test_p_values, by="variable")
  fixed_effects<-fixed_effects%>%
    mutate(effect_mean_std=fixed_effects$effect_size/antler_summaries_sub$Mean)
  
  fixed_effects<-fixed_effects[, c(1:4,6,7,5)]
  pedigree_fixed_effects_results<-rbind(pedigree_fixed_effects_results, fixed_effects)
  
  effect<-vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$Estimate[1]
  se<-vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$SE[1]
  chi2<- 2*(model1$loglik - model2$loglik)
  p.h2<-pchisq(chi2, df = 1, lower.tail = F)
  df_h2_temp<-data.frame(Trait= colnames(antler_data1)[m], Heritability=effect,
                         SE=se, p.value=p.h2)
  
  pedigree_h2_results<-rbind(pedigree_h2_results, df_h2_temp)
  write.table(var_table, file=paste0("variance_table_", colnames(antler_data1)[m],
                                                     "_ped.txt"), sep="\t", row.names = T)
  
}


write.table(pedigree_h2_results, file = "antler_h2_ped.txt", sep = "\t", row.names = F)
write.table(pedigree_fixed_effects_results, file="antler_h2_model_fixed_ped.txt", sep = "\t", row.names = F)


#compile variance component results

antler_data<-read.table("antler_data_comp.txt", header = TRUE)
colnames(antler_data)[1]<-"ID"
antler_traits<-as.list(colnames(antler_data)[7:16])

var_table_comp.ped<-data.frame()

for (m in antler_traits) {
  var_table<-read.table(paste0("variance_table_", m, "_ped.txt"), header=T)
  var_table$measure<-paste0(m)
  var_table_comp.ped<-rbind(var_table_comp.ped, var_table)
  
}

var_table_comp.ped$variable.names<-rep(c("RutYear", "BirthYear", "genetic", "permanent", "residual"), 10)
rownames(var_table_comp.ped)<-NULL
var_table_comp.ped<-var_table_comp.ped[, c(9, 1:8)]

#get repatability for each antler measure - permanent plus genetic plus birth year component

var_table_rep.ped<-data.frame()

for (m in antler_traits){
  var_table.ped<-var_table_comp.ped%>%
    filter(measure==m)%>%
    mutate(total_var=sum(component))
  
  var_table_temp<-var_table_comp.ped %>%
    filter(variable.names == "genetic" | variable.names == "permanent" | variable.names=="BirthYear")%>%
    filter(measure==m )%>%
    mutate(repeatability_comp=sum(component))
  
  var_table_sub<-var_table.ped%>%
    mutate(repeatability=var_table_temp[["repeatability_comp"]][1]/total_var)
  
  var_table_rep.ped<-rbind(var_table_rep.ped, var_table_sub)
}

write.table(var_table_rep.ped, file="Variance_components_h2_animal_model_ped.txt" , sep="\t",
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


#plot heritabilities

antler_h2_ped<-read.table("antler_h2_ped.txt", header=T)

#label data frames to signify significance

#label df for antler measures significant at <0.001 ***
label.df<-data.frame(Trait=antler_h2_ped$Trait[c(1:7, 10)],
                     Heritability=c(antler_h2_ped$Heritability[c(1:7, 10)]+antler_h2_ped$SE[c(1:7, 10)]+0.02))

#label df for antler measures significant at <0.01 **
label.df2<-data.frame(Trait=antler_h2_ped$Trait[c(8:9)],
                      Heritability=c(antler_h2_ped$Heritability[c(8:9)]+antler_h2_ped$SE[c(8:9)]+0.02))
#label for h2 value
label.df3<-data.frame(Trait=antler_h2_ped$Trait,
                      Heritability=c(antler_h2_ped$Heritability+antler_h2_ped$SE+0.04))


p<-ggplot(antler_h2_ped, aes(Trait, Heritability)) + 
  geom_bar(aes(fill=Trait), stat = "identity") + 
  scale_fill_viridis(discrete=T, option = "C") + 
  geom_errorbar(aes(ymin=Heritability-SE, ymax=Heritability+SE)) + 
  theme_classic()

p2<-p + theme(axis.text.x = element_text(angle = 90, vjust=1, hjust=1))

p3<-p2 + geom_text(data=label.df3, 
                   label=sprintf("%0.3f", round(antler_h2_ped$Heritability, digits = 4))) +
        theme(axis.text.x = element_text(size = 14, vjust = -0.01, colour = "black"), 
        axis.text.y = element_text(size=14, colour = "black"),
        axis.title = element_text(size=16), legend.position="none")


p4<-p3 + geom_text(data=label.df, label="***") 
p5<-p4+ geom_text(data=label.df2, label="**")

ggsave ("antler_h2_ped_plot.png", p5, scale = 1, width = 25, height=25, units="cm")


#GRM models ========================================================================================================================================================================================================================================================================================================================================

source("makeGRM.R")

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
no.ids<-unique(antler_data$ID)
#drop levels!
antler_data<-droplevels(antler_data)

#make inverse GRM
grminv <- makeGRM(grm.auto, ids.auto, antler_data$ID) # vector of IDs from the datasset that you use for the asreml model
attr(grminv, which = "INVERSE") <- TRUE

#animal model loop

#read in antler summaries table to mean standardise effect sizes of fixed effects
antler_summaries<-read.table("antler_summaries.txt", header=T)

#make variable for quadratic Age effect
antler_data$MeCaAge.2<-(antler_data$MeCaAge)^2

grm_h2_results <- data.frame()

grm_fixed_effects_results<- data.frame()

for (m in c(7:16)){
  fixed_fmla<-as.formula(paste0(colnames(antler_data)[m],"~ MeCaAge+MeCaAge.2"))
  model1<- asreml(fixed = fixed_fmla,
                  random = ~ vm(ID, grminv) + ide(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                  data = antler_data,
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
                  data = antler_data,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
  
  #fixed effects full model
  fixed_effects_temp<-as.data.frame(summary.asreml(model1, coef = T)$coef.fixed)
  
  wald.test<-as.data.frame(wald.asreml(model1))
  wald_test_p_values<-na.omit(data.frame(variable=rownames(wald.test), 
                                         p_value=wald.test$`Pr(Chisq)`))
  antler_summaries_sub<-subset(antler_summaries, Trait==colnames(antler_data)[m])
  
  fixed_effects<-data.frame(variable = rownames(fixed_effects_temp),
                            effect_size = fixed_effects_temp$solution,
                            std.error=fixed_effects_temp$`std error`, 
                            Z_ratio=fixed_effects_temp$z.ratio,
                            Trait = colnames(antler_data)[m])
  
  fixed_effects<-join(fixed_effects, wald_test_p_values, by="variable")
  fixed_effects<-fixed_effects%>%
    mutate(effect_mean_std=fixed_effects$effect_size/antler_summaries_sub$Mean)
  
  fixed_effects<-fixed_effects[, c(1:4,6,7,5)]
  grm_fixed_effects_results<-rbind(grm_fixed_effects_results, fixed_effects)
  
  effect<-vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$Estimate[1]
  se<-vpredict(model1, h2 ~ V3/(V1+V2+V3+V4+V5))$SE[1]
  chi2<- 2*(model1$loglik - model2$loglik)
  p.h2<-pchisq(chi2, df = 1, lower.tail = F)
  df_h2_temp<-data.frame(Trait= colnames(antler_data)[m], Heritability=effect,
                         SE=se, p.value=p.h2)
  
  grm_h2_results<-rbind(grm_h2_results, df_h2_temp)
  write.table(var_table, file=paste0("variance_table_", colnames(antler_data)[m],
                                     "_grm.txt"), sep="\t", row.names = T)
  
}


write.table(grm_h2_results, file = "antler_h2_grm.txt", sep = "\t", row.names = F)
write.table(grm_fixed_effects_results, file="antler_h2_model_fixed_grm.txt", sep = "\t", row.names = F)


#plot grm heritabilities

antler_h2_grm<-read.table("antler_h2_grm.txt", header = T)

label.df<-data.frame(Trait=antler_h2_grm$Trait[c(1:8, 10)],
                      Heritability=c(antler_h2_grm$Heritability[c(1:8, 10)]+antler_h2_grm$SE[c(1:8, 10)]+0.02))

label.df2<-data.frame(Trait=antler_h2_grm$Trait[9],
                     Heritability=c(antler_h2_grm$Heritability[9]+antler_h2_grm$SE[9]+0.02))

label.df3<-data.frame(Trait=antler_h2_grm$Trait,
                      Heritability=c(antler_h2_grm$Heritability+antler_h2_grm$SE+0.04))



pg<-ggplot(antler_h2_grm, aes(Trait, Heritability)) + 
  geom_bar(aes(fill=Trait), stat = "identity") + scale_fill_viridis(discrete=T, option = "C") + 
  geom_errorbar(aes(ymin=Heritability-SE, ymax=Heritability+SE)) + 
  theme_classic() + coord_cartesian(ylim = c(0, 0.55))

pg2<-pg + theme(axis.text.x = element_text(angle = 90, vjust=1, hjust=1))

pg3<-pg2 + geom_text(data=label.df, label="***")+ 
  theme(axis.text.x = element_text(size = 14, vjust = -0.01, colour = "black"), 
        axis.text.y = element_text(size=14, colour = "black"), 
        axis.title = element_text(size=16), legend.position="none")

pg4<-pg3 + geom_text(data=label.df2, label="**")

pg5<-pg4 + geom_text(data=label.df3, label=sprintf("%0.3f", round(antler_h2_grm$Heritability, digits = 4)))

ggsave("antler_h2_grm_plot.png", pg5, scale = 1, width = 25, height=25, units="cm")



#plot variance components

antler_data<-read.table("antler_data_comp.txt", header = TRUE)
colnames(antler_data)[1]<-"ID"
antler_traits<-as.list(colnames(antler_data)[7:16])

var_table_comp.grm<-data.frame()

for (m in antler_traits) {
  var_table<-read.table(paste0("variance_table_", m, "_grm.txt"), header=T)
  var_table$measure<-paste0(m)
  var_table_comp.grm<-rbind(var_table_comp.grm, var_table)
  
}

var_table_comp.grm$variable.names<-rep(c("RutYear", "BirthYear", "genetic", "permanent", "residual"), 10)
rownames(var_table_comp.grm)<-NULL
var_table_comp.grm<-var_table_comp.grm[, c(9, 1:8)]

#get repatability for each antler measure - permanent plus genetic component
var_table_comp.grm<-read.table("Variance_components_h2_animal_model_grm.txt", header=T)
head(var_table_comp.grm)

antler_traits<-as.list(unique(as.character(var_table_comp.grm[["measure"]])))

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

write.table(var_table_rep.grm, file="Variance_components_h2_animal_model_grm.txt" , sep="\t",
            row.names = F, col.names = T)


var_table_flipped_effect.grm<-dcast(var_table_rep.grm[, c(1, 7:11)],
                                    measure+repeatability+total_var~variable.names, 
                                    value.var = c("Effect"))

var_table_flipped_SE.grm<-dcast(var_table_rep.grm[, c(1, 7:9)], measure~variable.names, 
                                value.var = c("SE"))

colnames(var_table_flipped_SE.grm)[-1]<-paste0(colnames(var_table_flipped_SE.grm)[-1], "_SE")

var_table_flipped.grm<-join(var_table_flipped_effect.grm, var_table_flipped_SE.grm, by="measure")
var_table_flipped.grm<-var_table_flipped.grm[c(7,4,8,10,3,5,2,9,1,6), c("measure", "total_var", "genetic", "genetic_SE",
                                                                        "permanent", "permanent_SE", "BirthYear", "BirthYear_SE",
                                                                        "RutYear", "RutYear_SE","residual", "residual_SE","repeatability")]

write.table(var_table_flipped.grm, file="Variance_composition_h2_animal_model_grm.txt" , sep="\t",
            row.names = F, col.names = T)


antler_variables<-as.list(unique(as.character(var_table_comp.grm[["variable.names"]])))
v<-antler_variables[[1]]

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


#make variable with measure group for facet wrap

var_table_comp.grm$measure.group<-ifelse(var_table_comp.grm$measure=="BrowLength", "group5",
                                         ifelse(var_table_comp.grm$measure=="CoronetTrayJunc", "group5",
                                                ifelse(var_table_comp.grm$measure=="Length", "group2",
                                                       ifelse(var_table_comp.grm$measure=="CoronetCirc", "group4",
                                                              ifelse(var_table_comp.grm$measure=="AntlerWt", "group1",
                                                                     ifelse(var_table_comp.grm$measure=="TrayLength","group5","group3"))))))

#stacked plots

#calculate total variance for each measure and proportion that each component makes up of that variance

var_table_total.grm<-ddply(var_table_comp.grm, .(measure), summarise,
                           total.variance=sum(component))
var_table_comp.grm<-join(var_table_comp.grm,var_table_total.grm)

var_table_comp.grm<-mutate(var_table_comp.grm, prop.variance=component/total.variance)


var_table_comp.grm$measures.named<-rep(c("Length", "Coronet Circumference", "Lower Beam", "Upper Beam",
                                         "Coronet-Brow-Junction", "Coronet-Tray-Junction", "Brow Length",
                                         "Tray Length", "Weight", "Form"), each=5)

p.var<-ggplot(var_table_comp.grm, aes(fct_inorder(measures.named), prop.variance, group=fct_inorder(variable.names), fill=fct_inorder(variable.names), width = 0.8))+
  geom_bar(stat = "identity") + 
  scale_fill_brewer(type="qual", palette = "Accent", labels = c("rut year", "birth year", "additive genetic", "permanent", "residual")) + theme_classic()+ 
  theme(axis.text.x = element_text(size = 14, colour = "black", angle = 90, vjust=0.3, hjust=1), 
        axis.text.y= element_text(size=14, colour="black"),
        axis.title = element_text(size=18), legend.text = element_text(size=13),
        legend.title = element_blank())+xlab(NULL)+
  ylab("proportion of variance")+theme(strip.text.x = element_blank())

ggsave("antler_variance_components_grm_plot.png", p.var, scale = 1, width = 25, height=25, units="cm")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#make comparative plots contrasting pedigree and genomic relatedness estimates of heritability and add Loeske Kruuk's (2014) previous estimates

#read in h2 data and rename heritability to identify origin of data e.g. grm or pedigree
antler_h2_ped<-read.table("antler_h2_ped.txt", header=T)
colnames(antler_h2_ped)[2]<-c("pedigree")
antler_h2_grm<-read.table("antler_h2_grm.txt", header = T)
colnames(antler_h2_grm)[2]<-c("grm")
#create data frame for Loeske's previous rsults on h2 of antler weight & form (2014)
se.w<-0.078
h2.w<-0.378
se.b<-0.071
h2.b<-0.373 
se.l<-0.047
h2.l<-0.384
se.c<-0.082
h2.c<-0.364
se.f<- 0.053
h2.f<-0.240
antler_h2_previous<-data.frame(Trait=c("AntlerWt","BrowLength", "Length", "CoronetCirc","Form"), 
                               previous=c(h2.w, h2.b, h2.l, h2.c, h2.f))
# join and melt heritability part of data frames (including Loeske's result)
antler_h2<-join(antler_h2_ped[,c(1:2)], antler_h2_grm[,c(1:2)], by="Trait")
antler_h2<-join(antler_h2, antler_h2_previous, by="Trait")
antler_h2<-melt(antler_h2)
colnames(antler_h2)[3]<-"Heritability"
#join and melt standard error data- including se from Loeske's previous paper
antler_SE_ped<-antler_h2_ped[, c(1,3)]
colnames(antler_SE_ped)[2]<-"pedigree"
antler_SE_grm<-antler_h2_grm[, c(1,3)]
colnames(antler_SE_grm)[2]<-"grm"
antler_SE_previous<-data.frame(Trait=c("AntlerWt","BrowLength", "Length", "CoronetCirc","Form"), 
                               previous=c(se.w, se.b, se.l, se.c, se.f))
antler_SE<-join(antler_SE_ped, antler_SE_grm, by="Trait")
antler_SE<-join(antler_SE, antler_SE_previous, by="Trait")
antler_SE<-melt(antler_SE)
colnames(antler_SE)[3]<-"SE"

#join melted data frames with heritabilities and SE
antler_h2_full<-join(antler_h2, antler_SE)
variable_levels<-c("pedigree","grm","previous")
antler_h2_full$variable<-factor(antler_h2_full$variable, levels = variable_levels)


#make plot with pedigree, grm and previous (Loeske's) heritabilities

p<-ggplot(antler_h2_full, aes(x = fct_inorder(Trait), y= Heritability, fill = variable)) +
  geom_bar(stat="identity", width=.8, position = "dodge")+
  geom_errorbar(aes(ymin=Heritability-SE, ymax=Heritability+SE), 
                position = position_dodge(width=0.8))+
  scale_fill_brewer(type="qual", palette = "Accent")+
  theme_classic()+ ylim(0,0.6)+
  scale_x_discrete(labels=c("Antler Length", "Coronet Circumference", "Lower Beam Circ.", 
                            "Upper Beam Circ.", "Coronet-Brow Junc.", "Coronet-Tray Junc.", "Brow Length", "Tray Length", 
                            "Antler Weight", "Form"))+
  theme(axis.text.x = element_text(size = 20, colour = "black", angle = 90, vjust=0.4, hjust=1), 
        axis.text.y = element_text(size=15, colour = "black"), 
        axis.title = element_text(size=20),legend.text = element_text(size=16),
        legend.title = element_blank())+xlab(element_blank())

antler_h2_ped<-read.table("antler_h2_ped.txt", header=T)

antler_h2_grm<-read.table("antler_h2_grm.txt", header = T)


#pedigree labels
label.df.p<-data.frame(Trait=antler_h2_ped$Trait,
                        variable="pedigree",
                        Heritability=c(antler_h2_ped$Heritability+antler_h2_ped$SE+0.02))
#ped sig p<0.001 ***
label.df.p.1<-label.df.p[c(1:7, 10), ] #***
label.df.p.1.na<-label.df.p[c(8:9), ]
#label.df.p.1.na$variable<-NA
label.df.p.1.na$Heritability<-NA
label.df.p.1.full<-rbind(label.df.p.1, label.df.p.1.na)

#ped sig p<0.01 **
label.df.p.2<-label.df.p[c(8:9), ] #**
label.df.p.2.na<-label.df.p[c(1:7, 10), ]
#label.df.p.2.na$variable<-NA
label.df.p.2.na$Heritability<-NA
label.df.p.2.full<-rbind(label.df.p.2, label.df.p.2.na)


#grm labels
label.df.g<-data.frame(Trait=antler_h2_grm$Trait,
                      variable="grm",
                      Heritability=c(antler_h2_grm$Heritability+antler_h2_grm$SE+0.02))

#grm sig p<0.001 ***
label.df.g.1<-label.df.g[c(1:8, 10), ] #***
label.df.g.1.na<-label.df.g[9, ]
#label.df.g.1.na$variable<-NA
label.df.g.1.na$Heritability<-NA
label.df.g.1.full<-rbind(label.df.g.1, label.df.g.1.na)

#grm sig p<0.01 **
label.df.g.2<-label.df.g[9, ] #**
label.df.g.2.na<-label.df.g[c(1:8, 10), ]
#label.df.g.2.na$variable<-NA
label.df.g.2.na$Heritability<-NA
label.df.g.2.full<-rbind(label.df.g.2, label.df.g.2.na)

label.df.pv<-data.frame(Trait=antler_h2_previous$Trait, variable="previous", 
                      Heritability=antler_h2_previous$previous+antler_SE_previous$previous+0.02)

p1<-p + geom_text(data=label.df.p.1.full, label="***", hjust=1.6) 
p2<-p1 + geom_text(data=label.df.g.1.full, label="***", hjust=0.5) 
p3<-p2+ geom_text(data=label.df.g.2.full, label="**", hjust=0.5)
p4<-p3+geom_text(data=label.df.p.2.full, label="**", hjust=1.6)
p5<-p4+geom_text(data=label.df.pv, label="*", hjust=-3) 


ggsave("antler_h2_biplot_with_previous.png", p5, scale = 1, width = 25, height=25, units="cm")
