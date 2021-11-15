#SNP effect size estimation in animal model. Test for effects on antler phenotypes 
#of SNPs with highest effect sizes from GWAS/Bayesian FDR (ashR)

source("makeGRM.R")
library(asreml)
library(msm)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)
library(forcats)

#read in genotypes and map (containing SNP names)
SNP_geno<-read.table("Top_SNPs_PCs_genotypes.ped", header = F)
SNP_map<-read.table("Top_SNPs_PCs_genotypes.map", header = F)

#match genotypes to SNP IDs - take marker IDs from map table and make them column names of genotype table

SNP_IDs<-as.character(SNP_map$V2)
colnames(SNP_geno)[c(7:ncol(SNP_geno))]<-SNP_IDs
colnames(SNP_geno)[2]<-"ID"

#read in GRM
grm.auto <- read.table("Deer33.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.auto <- read.table("Deer33.grm.id")  # CONTAINS ID LIST

#read in antler data set
antler_data<-read.table("antler_resid_PCs.txt", header = TRUE)
colnames(antler_data)[1]<-"ID"
head(antler_data)

#ASreml likes categories to be factors
antler_data$ID<-as.factor(antler_data$ID)
antler_data$MeCaYear<-as.factor(antler_data$MeCaYear)
antler_data$BirthYear<-as.factor(antler_data$BirthYear)

#all IDs in antler_data must be represented in GRM, so subset data to remove individuals not in GRM
antler_data<-subset(antler_data, ID %in% ids.auto$V2)
#remove individuals with no genotype data
antler_data<-subset(antler_data, ID %in% SNP_geno$ID)
#drop levels!
antler_data<-droplevels(antler_data)

#make inverse GRM
grminv <- makeGRM(grm.auto, ids.auto, antler_data$ID) # vector of IDs from the datasset that you use for the asreml model
attr(grminv, which = "INVERSE") <- TRUE

#read in tables with significant SNP results from Bayesian FDR to associate SNPs with antler trait
antler_PCs_sig.SNPs<-read.table("PC.SNP_effect_size_estimation_results_sig.txt", header=T)
head(antler_PCs_sig.SNPs)
#lists of antler traits and SNPs to test in loop
antler_traits<-as.list(unique(antler_PCs_sig.SNPs$measure))
SNPs<-as.list(SNP_map$V2)


#test variables ----------------------------------------------------------------------------------------------------------------------------
# levels(SNP_geno$cela1_red_4_33583475)
m<-"PC2"
s<-"cela1_red_x_33757283"

#-------------------------------------------------------------------------------------------------------------------------------------------

#results data frames
model_summary_all.cat<-data.frame()
model_wald.test_all.cat<-data.frame()
model_summary_all.add<-data.frame()
model_wald.test_all.add<-data.frame()

var_table.comp<-data.frame()

co_var.table.fixed_all<-data.frame()

#loop to test the effect of all the top SNPs (highest effect sizes in GWAS/Bayesian FDR) 
#on their respective antler measures in models treating genotypes either as categories or additive

for (m in antler_traits){
  m<-as.character(m)
  #subset significant SNP data to only contain one antler measure
  antler_PCs_sig.SNPs_sub<-subset(antler_PCs_sig.SNPs, measure==m)
  
  for (s in SNPs){
    s<-as.character(s)
    match<-grep(s, antler_PCs_sig.SNPs_sub$SNP.name, value = T)
    
    if (length(match)>0){
      
      #make sure there are no "00" genotypes for SNPs
      zero.rows<-grep("00", SNP_geno[, paste0(s)])
      if (length(zero.rows)>0){
        SNP_geno<-SNP_geno[c(-zero.rows), ]
        SNP_geno<-droplevels(SNP_geno) 
      }  
      
      
      #exclude entries for which focal antler measure is missing to avoid 
      #null solution for genotype with no associated antler measure
      antler_data_sub<-antler_data[, c("ID", "MeCaYear", "BirthYear", "MeCaAge", paste0(m))]
      antler_data_sub<-na.omit(antler_data_sub)
      
      #join SNP genotype to antler data
      antler_data_sub<-join(antler_data_sub, SNP_geno[, c("ID", paste0(s))])
      
      
      #add copy number variable - G or C as focal allele
      antler_data_sub$SNP_copy_no<-ifelse(antler_data_sub[, paste0(s)] =="AA", 0,
                                          ifelse(antler_data_sub[, paste0(s)] == "GA", 1,
                                                 ifelse(antler_data_sub[, paste0(s)] == "AG", 1,
                                                        ifelse(antler_data_sub[, paste0(s)] == "GG", 2,
                                                               ifelse(antler_data_sub[, paste0(s)] == "CA", 1,
                                                                      ifelse(antler_data_sub[, paste0(s)] == "AC", 1,
                                                                             ifelse(antler_data_sub[, paste0(s)] == "CC", 2, NA)))))))
      #add second copy number variable - A as focal allele
      antler_data_sub$SNP_copy_no.2<- ifelse(antler_data_sub[, paste0(s)] =="AA", 2,
                                             ifelse(antler_data_sub[, paste0(s)] == "GA", 1,
                                                    ifelse(antler_data_sub[, paste0(s)] == "AG", 1,
                                                           ifelse(antler_data_sub[, paste0(s)] == "GG", 0,
                                                                  ifelse(antler_data_sub[, paste0(s)] == "CA", 1,
                                                                         ifelse(antler_data_sub[, paste0(s)] == "AC", 1,
                                                                                ifelse(antler_data_sub[, paste0(s)] == "CC", 0, NA)))))))
      
      
      
      
      antler_data_sub[, paste0(s)]<-as.factor(antler_data_sub[, paste0(s)])
      #drop levels!
      antler_data_sub<-droplevels(antler_data_sub)
      
      #SNP as category - AA, GA, GG
      fixed_formula1<-as.formula(paste0(m, "~", s))
      
      model1<- asreml(fixed = fixed_formula1,
                      random = ~ vm(ID, grminv) + ide(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                      data = antler_data_sub,
                      na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
      
      
      fixed_formula1.1<-as.formula(paste0(m, "~1"))
      
      model1.1<- asreml(fixed = fixed_formula1.1,
                        random = ~ vm(ID, grminv) + ide(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                        data = antler_data_sub,
                        na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
      
      
      #is the SNP significant? LRT to test SNP as a whole, not differences between levels
      lrt.df<-model1$aov[2,2]
      chi2<- 2*(model1$loglik - model1.1$loglik)
      p.value<-pchisq(chi2, df = lrt.df, lower.tail = F)
      
      #summarise fixed effects model results - sum = summary with effect estimates, wald= significance testing
      n.variables<-nrow(summary.asreml(model1, coef = T)$coef.fixed)
      
      model.category.sum<-data.frame(variable=row.names(summary.asreml(model1, coef = T)$coef.fixed),
                                     solution=summary.asreml(model1, coef = T)$coef.fixed[c(1:n.variables)],
                                     SE=summary.asreml(model1, coef = T)$coef.fixed[c((1+n.variables):(2*n.variables))],
                                     z_ratio=summary.asreml(model1, coef = T)$coef.fixed[c((2*n.variables+1):(3*n.variables))],
                                     model_type="categorical", measure=m)
      
      model.category.sum<-mutate(model.category.sum, p.z_ratio=2*pnorm(-abs(z_ratio)))
      
      model.category.sum<-model.category.sum[, c(1,2,3,4,7,5,6)]
      
      model.category.wald<-as.data.frame(wald.asreml(model1))
      model.category.wald$variable<-rownames(model.category.wald)
      rownames(model.category.wald)<-NULL
      model.category.wald<-model.category.wald[, c(5,1,2,3,4)]
      model.category.wald<-mutate(model.category.wald, p.LRT.SNP=p.value, model_type="categorical", measure=m)
      
      #variance table for categorical model
      var_table.cat<-summary.asreml(model1, coef = T)$varcomp[c(1:5), ]    # variance components table
      
      #get proportion and associated std error for all variance components
      #add all of them apart from units!R 
      
      var_table.cat_prop<-data.frame(variable=rownames(var_table.cat), Effect=c(
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
      
      var_table.cat_comp<-cbind(var_table.cat, var_table.cat_prop[-1])
      var_table.cat_comp<-mutate(var_table.cat_comp, measure=m, SNP=s, model="categorical")
      
      
      #run model to get variance-covariance matrix for fixed effects (Cfixed=T); needed to calculate error for variacne explained by SNP
      #needs to be run in categorical model, otherwise no covariance or specific variance for different genotypes;
      #change fixed formula to omit intercept- more intuitive to get variances or covariances of fixed effects if none of the levels is
      #null estimate (in case for example AA has the lowest effect, then A would be focal allele (because convention is to look at reducing allele)
      #but no co-/variance estimate for AA genotype in Cfixed matrix (would be in intercept column))
      
      
      fixed_formula1.2<-as.formula(paste0(m, "~ -1+", s))
      asreml.options(Cfixed=TRUE)
      
      model1.2<- asreml(fixed = fixed_formula1.2,
                        random = ~ vm(ID, grminv) + ide(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                        data = antler_data_sub,
                        na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
      
      
      
      co_var.table.fixed<-data.frame(variable=rep(rownames(model1.2$Cfixed), nrow(model1.2$Cfixed)), 
                                     value=model1.2$Cfixed[c(1:(nrow(model1.2$Cfixed)*ncol(model1.2$Cfixed)))], 
                                     variable2=rep(colnames(model1.2$Cfixed), each=nrow(model1.2$Cfixed)),
                                     measure=m, SNP.name=s)
      
      co_var.table.fixed_all<-rbind(co_var.table.fixed_all, co_var.table.fixed)
      
      
      #SNP as continuous variable - 0, 1, 2 copy numbers of allele 1(G or C)
      
      fixed_formula2<-as.formula(paste0(m, "~ SNP_copy_no"))
      
      model2<- asreml(fixed = fixed_formula2,
                      random = ~ vm(ID, grminv) + ide(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                      data = antler_data_sub,
                      na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
      
      n.variables<-nrow(summary.asreml(model2, coef=T)$coef.fixed)
      
      
      #is there a C allele instead of G?
      C.allele<-grep("C", SNP_geno[, paste0(s)], value = T)
      
      if (length(C.allele)>0){
        focal.allele<-"C"
      } else {
        focal.allele<-"G"
      }
      
      
      #summarise fixed effects model results - sum = summary with effect estimates, wald= significance testing
      model.additive1.sum<-data.frame(variable=row.names(summary.asreml(model2, coef=T)$coef.fixed),
                                      solution=summary.asreml(model2, coef=T)$coef.fixed[c(1:n.variables)],
                                      SE=summary.asreml(model2, coef=T)$coef.fixed[c((1+n.variables):(2*n.variables))],
                                      z_ratio=summary.asreml(model2, coef=T)$coef.fixed[c((2*n.variables+1):(3*n.variables))],
                                      model_type="additive", main_allele=focal.allele, SNP=paste0(s), measure=m)
      
      
      model.additive1.wald<-as.data.frame(wald.asreml(model2))
      model.additive1.wald$variable<-rownames(model.additive1.wald)
      rownames(model.additive1.wald)<-NULL
      model.additive1.wald<-model.additive1.wald[, c(5,1,2,3,4)]
      model.additive1.wald<-mutate(model.additive1.wald, model_type="additive", main_allele=focal.allele, SNP=paste0(s), measure=m)
      
      #variance table for additive model
      var_table.add<-summary.asreml(model2, coef = T)$varcomp[c(1:5), ]    # variance components table
      
      #get proportion and associated std error for all variance components
      #add all of them apart from units!R 
      
      var_table.add_prop<-data.frame(variable=rownames(var_table.add), Effect=c(
        vpredict(model2, RY ~ V1/(V1+V2+V3+V4+V5))$Estimate[1],
        vpredict(model2, BY ~ V2/(V1+V2+V3+V4+V5))$Estimate[1],
        vpredict(model2, h2 ~ V3/(V1+V2+V3+V4+V5))$Estimate[1],
        vpredict(model2, ID ~ V4/(V1+V2+V3+V4+V5))$Estimate[1],
        vpredict(model2, RS ~ V5/(V1+V2+V3+V4+V5))$Estimate[1]),
        SE=c(
          vpredict(model2, RY ~ V1/(V1+V2+V3+V4+V5))$SE[1],
          vpredict(model2, BY ~ V2/(V1+V2+V3+V4+V5))$SE[1],
          vpredict(model2, h2 ~ V3/(V1+V2+V3+V4+V5))$SE[1],
          vpredict(model2, ID ~ V4/(V1+V2+V3+V4+V5))$SE[1],
          vpredict(model2, RS ~ V5/(V1+V2+V3+V4+V5))$SE[1]))
      
      var_table.add_comp<-cbind(var_table.add, var_table.add_prop[-1])
      var_table.add_comp<-mutate(var_table.add_comp, measure=m, SNP=s, model="additive")
      
      #SNP as continuous variable - 0, 1, 2 copy numbers of allele 1(A)
      
      fixed_formula3<-as.formula(paste0(m, "~ SNP_copy_no.2"))
      
      model3<- asreml(fixed = fixed_formula3,
                      random = ~ vm(ID, grminv) + ide(ID) + BirthYear + MeCaYear, #vm(ID, ainv) is the relatedness matrix, ide(ID) is the individual identity
                      data = antler_data_sub,
                      na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units))
      
      
      #summarise fixed effects model results - sum = summary with effect estimates, wald= significance testing
      n.variables<-nrow(summary.asreml(model3, coef=T)$coef.fixed)
      
      model.additive2.sum<-data.frame(variable=row.names(summary.asreml(model3, coef=T)$coef.fixed),
                                      solution=summary.asreml(model3, coef=T)$coef.fixed[c(1:n.variables)],
                                      SE=summary.asreml(model3, coef=T)$coef.fixed[c((1+n.variables):(2*n.variables))],
                                      z_ratio=summary.asreml(model3, coef=T)$coef.fixed[c((2*n.variables+1):(3*n.variables))],
                                      model_type="additive", main_allele="A", SNP=paste0(s), measure=m)
      
      
      model.additive2.wald<-as.data.frame(wald.asreml(model3))
      model.additive2.wald$variable<-rownames(model.additive2.wald)
      rownames(model.additive2.wald)<-NULL
      model.additive2.wald<-model.additive2.wald[, c(5,1,2,3,4)]
      model.additive2.wald<-mutate(model.additive2.wald, model_type="additive", main_allele="A", SNP=paste0(s),measure=m)
      
      
      #bind addditve model results and then combine all results to data frame outside loop for both additive and categorical models
      all_models.sum.add<-rbind(model.additive1.sum, model.additive2.sum)
      all_models.wald.add<-rbind(model.additive1.wald, model.additive2.wald)
      
      model_summary_all.cat<-rbind(model_summary_all.cat, model.category.sum)
      model_wald.test_all.cat<-rbind(model_wald.test_all.cat, model.category.wald)
      
      model_summary_all.add<-rbind(model_summary_all.add, all_models.sum.add)
      model_wald.test_all.add<-rbind(model_wald.test_all.add, all_models.wald.add)
      
      #bind variance tables of categorical and additive models
      var_table.all<-rbind(var_table.cat_comp, var_table.add_comp)
      variable.names.var<-rep(rownames(summary.asreml(model2, coef = T)$varcomp[c(1:5), ]), (nrow(var_table.all)/5))
      variable.names.var<-as.data.frame(variable.names.var)
      colnames(variable.names.var)[1]<-"variable.names"
      var_table.all<-cbind(variable.names.var, var_table.all)
      
      var_table.comp<-rbind(var_table.comp, var_table.all)
      
      #read in genotypes again (need new one for each SNP because of removal of "00" genotypes)
      SNP_geno<-read.table("Top_SNPs_PCs_genotypes.ped", header = F)
      
      #match genotypes to SNP IDs
      colnames(SNP_geno)[c(7:ncol(SNP_geno))]<-SNP_IDs
      colnames(SNP_geno)[2]<-"ID"
      
    }
  }
}




write.table(model_summary_all.add, file = "Summary_table_additive_models.top_SNPs_PCs.txt", sep = "\t", 
            col.names = T, row.names = F)

write.table(model_wald.test_all.add, file = "Sig.test_table_additive_models.top_SNPs_PCs.txt", sep = "\t", 
            col.names = T, row.names = F)


write.table(model_summary_all.cat, file = "Summary_table_categorical_models.top_SNPs_PCs.txt", sep = "\t", 
            col.names = T, row.names = F)


write.table(model_wald.test_all.cat, file = "Sig.test_table_categorical_models.top_SNPs_PCs.txt", sep = "\t", 
            col.names = T, row.names = F)

write.table(co_var.table.fixed_all, file="Variance_covariance_table_fixed_effects_SNP_models_PCs.txt",
            sep = "\t", col.names = T, row.names = F)

write.table(var_table.comp, file = "variance_table_SNP_animal_models_PCs.txt", sep = "\t", row.names = F)

#filter for significant SNPs -------------------------------------------------------------------------------------------------------------------------------
#categorical models, SNPs significant according to Wald test

#first filter for SNP entries
match<-grep("^cela1_red_", model_wald.test_all.cat$variable)

SNPs.sig.cat.wald<-model_wald.test_all.cat[c(match), ]
SNPs.sig.cat.wald<-subset(SNPs.sig.cat.wald, `Pr(Chisq)`<0.05)


#filter model summary (categorical models) for significant SNPs
sig.SNP_names<-as.list(SNPs.sig.cat.wald$variable)

summary.match<-character()

for (sp in sig.SNP_names){
  sp<-as.character(sp)
  sum.match<-grep(sp, model_summary_all.cat$variable)
  summary.match<-append(summary.match, sum.match)
}
summary.match<-as.data.frame(summary.match)

SNPs.sig.cat.sum<-subset(model_summary_all.cat, rownames(model_summary_all.cat)%in%summary.match$summary.match)


class(SNPs.sig.cat.wald$measure)
SNPs.sig.cat.sum$measure<-as.character(SNPs.sig.cat.sum$measure)

SNPs.sig.cat.wald.ordered<-SNPs.sig.cat.wald[order(SNPs.sig.cat.wald$measure), ]
SNPs.sig.cat.sum.ordered<-SNPs.sig.cat.sum[order(SNPs.sig.cat.sum$measure), ]

SNPs.sig.cat.wald.multi<-SNPs.sig.cat.wald.ordered[c(rep(1:23, each=3), rep(24:27, each=2), rep(28:35, each=3), 
                                                     rep(36:37, each=2), rep(38:44, each=3),rep(45:46, each=2),
                                                     rep(47:53, each=3), rep(54, each=2), rep(55, each=3)), ]
colnames(SNPs.sig.cat.wald.multi)[1]<-"SNP.name"
SNPs.sig.cat<-cbind(SNPs.sig.cat.sum.ordered, SNPs.sig.cat.wald.multi[, c(1:6)])

write.table(SNPs.sig.cat, file="Sig.SNPs.PCs.categorical.model.txt", sep = "\t", row.names = F, col.names = T)



#significant SNPs in additive model------------------------------------------------------------------------------------------------------------------
#filter for SNP entries
match<-grep("SNP_copy_", model_wald.test_all.add$variable)

SNPs.sig.add.wald<-model_wald.test_all.add[c(match), ]
SNPs.sig.add.wald<-subset(SNPs.sig.add.wald, `Pr(Chisq)`<0.05)

#exlude non-SNP variable entries from summary table
match<-grep("SNP_copy_", model_summary_all.add$variable)
SNPs.sig.add.sum<-model_summary_all.add[c(match), ]

#filter model summary (additive) for significant SNPs
SNPs.sig.add.sum<-subset(SNPs.sig.add.sum, SNP %in% SNPs.sig.add.wald$SNP)


#join significant SNPs wald and summary tables for additive models (same order)
SNPs.sig.add<-join(SNPs.sig.add.sum, SNPs.sig.add.wald)
colnames(SNPs.sig.add)[7]<-"SNP.name"

write.table(SNPs.sig.add, file="Sig.SNPs.PCs.additive.model.txt", sep = "\t", row.names = F, col.names = T)
SNPs.sig.add<-read.table("Sig.SNPs.PCs.additive.model.txt", header=T)

SNPs.sig.add$measure<-as.character(SNPs.sig.add$measure)

#get SNP trait combination significant in additive but not categorical model
sig.traits.add<-as.list(c(unique(SNPs.sig.add$measure)))
SNPs.sig.cat<-read.table("Sig.SNPs.PCs.categorical.model.txt", header=T)
SNPs.sig.add.only<-data.frame()

for (st in sig.traits.add){
  SNPs.sig.add_sub<-subset(SNPs.sig.add, measure==st)
  SNPs.sig.cat_sub<-subset(SNPs.sig.cat, measure==st)
  sig.SNPs.sub<-as.list(SNPs.sig.add_sub$SNP.name)
  no_dup.SNPs<-seq(1, length(sig.SNPs.sub), by=2)
  sig.SNPs.sub<-sig.SNPs.sub[no_dup.SNPs]
  for (S in sig.SNPs.sub){
    S<-droplevels(S)
    S<-as.character(S)
    match<-grep(S, SNPs.sig.cat_sub$SNP.name)
    if(length(match)==0){
      df.temp<-subset(SNPs.sig.add_sub, SNP.name==paste0(S))
      df.temp<-subset(df.temp, !duplicated(SNP.name))
      SNPs.sig.add.only<-rbind(SNPs.sig.add.only, df.temp)
    }
  }
}


# #get categorical value estimates for SNPs only significant in additive model (to be able to calculate SNP variance)
# SNPs.sig.add.only_summary<- data.frame()
# SNPs.sig.add.only$SNP.name<-as.character(SNPs.sig.add.only$SNP.name)
# snps.sig.add.only<-as.list(unique(SNPs.sig.add.only$SNP.name))
# 
# traits.snps.sig.add.only<-as.list(unique(SNPs.sig.add.only$measure))
# 
# for (t in traits.snps.sig.add.only) {
#   model_summary_all.cat_sub<-subset(model_summary_all.cat, measure==paste0(t))
#   for (s in snps.sig.add.only){
#     match<-grep(paste0("^",s,"_.*$"), model_summary_all.cat_sub$variable)
#     if(length(match)>0){
#       df.temp<-model_summary_all.cat_sub[match, ]
#       SNPs.sig.add.only_summary<-rbind(SNPs.sig.add.only_summary, df.temp)
#     }
#   }
# }
# 
# 
# 
# SNPs.sig.add.only_summary<-droplevels(SNPs.sig.add.only_summary)
# SNPs.sig.add.only_summary<-mutate(SNPs.sig.add.only_summary, SNP.name=SNPs.sig.add.only_summary$variable)
# SNPs.sig.add.only_summary$SNP.name<-gsub(".{3}$", "", SNPs.sig.add.only_summary$SNP.name)
# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######### Siginificant SNP variance analysis ##########

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#how much variation is explained by significant SNPs?
#based on categorical model output (has right data - i.e. estimated for different genotypes)

#read in output data of significant SNPs 
SNPs.sig.cat<-read.table("Sig.SNPs.PCs.categorical.model.txt", header=T)

# #join SNPs only significant in additive model to all other significant SNPs in categorical model
# SNPs.sig.all<-rbind(SNPs.sig.cat[, c(1:8)], SNPs.sig.add.only_summary)

#all SNps-trait combinations significant in both categorical and additve model
SNPs.sig.all<-SNPs.sig.cat[, c(1:8)]
head(SNPs.sig.all)
#read in allele frequency data
allele.freq<-read.table("Allele_frequencies-GWAS_SNPs.txt", header = T)

#join allele frequencies to significant SNP data 
#(all significant SNP trait combinations- including the ones only significant in additive models)
SNPs.sig.cat.freq<-join(SNPs.sig.all, allele.freq)
head(SNPs.sig.cat.freq)

#~~ Specify values

#We adopt the convention that the A allele is the allele which increases
#the trait value

#specify fixed effect sizes 
antler_traits<-as.list(unique(SNPs.sig.cat.freq$measure))
SNPs<-as.list(unique(SNPs.sig.cat.freq$SNP.name))
SNP.sig.cat.AB_effects<-data.frame()
dominance.effect.SNPs<-character()


for (m in antler_traits) {
  m<-as.character(m)
  SNPs.sig.cat.freq_sub<-subset(SNPs.sig.cat.freq, measure==paste0(m))
  for (s in SNPs) {
    s<-as.character(s)
    match<-grep(s, SNPs.sig.cat.freq_sub$SNP.name, value = T)
    
    if(length(match)>0){
      
      SNPs.sig.cat.freq_sub2<-subset(SNPs.sig.cat.freq_sub, SNP.name==paste0(s))
      effect.sum<-sum(SNPs.sig.cat.freq_sub2$solution)
      
      if (effect.sum>0){
        SNPs.sig.cat.freq_sub2$effect.solution<-ifelse(SNPs.sig.cat.freq_sub2$solution==0, -(max(SNPs.sig.cat.freq_sub2$solution)),
                                                       ifelse(SNPs.sig.cat.freq_sub2$solution==max(SNPs.sig.cat.freq_sub2$solution), 0, -(max(SNPs.sig.cat.freq_sub2$solution)-(SNPs.sig.cat.freq_sub2$solution))))
        
        SNPs.sig.cat.freq_sub2$effect<-ifelse(SNPs.sig.cat.freq_sub2$effect.solution==0, "effectAA",
                                              ifelse(SNPs.sig.cat.freq_sub2$effect.solution==min(SNPs.sig.cat.freq_sub2$effect.solution), "effectBB", "effectAB"))
        
      } else {
        SNPs.sig.cat.freq_sub2$effect.solution<-SNPs.sig.cat.freq_sub2$solution
        
        SNPs.sig.cat.freq_sub2$effect<-ifelse(SNPs.sig.cat.freq_sub2$effect.solution==0, "effectAA",
                                              ifelse(SNPs.sig.cat.freq_sub2$effect.solution==min(SNPs.sig.cat.freq_sub2$effect.solution), "effectBB", "effectAB"))
        
      }
      
      max.effect<-max(abs(SNPs.sig.cat.freq_sub2$solution))
      SNPs.sig.cat.freq_sub.max<-subset(SNPs.sig.cat.freq_sub2, abs(solution)==max.effect)
      
      if (SNPs.sig.cat.freq_sub.max$effect=="effectAB"){
        dominance.effect.SNPs<-append(dominance.effect.SNPs, SNPs.sig.cat.freq_sub.max$SNP.name)
      }
      
      SNP.sig.cat.AB_effects<-rbind(SNP.sig.cat.AB_effects, SNPs.sig.cat.freq_sub2)
      
    }
  }
}




SNP.sig.cat.AB_effects<-SNP.sig.cat.AB_effects[, c(1:3, 13:14, 4:12)]
head(SNP.sig.cat.AB_effects)

SNP.sig.cat.AB_effects_BB<-subset(SNP.sig.cat.AB_effects, effect=="effectBB")
#need to fix effect classification for SNP were heterozygote has lowest solution estimate (and het and non 0 homozygote have different signs)
SNP.sig.cat.AB_effects[c(109:110), 5]<-c("effectBB", "effectAB")

#add info about which allele is reducing allele ("B" allele)
SNP.sig.cat.allele.info<-data.frame()

for (m in antler_traits){
  m<-as.character(m)
  SNP.sig.cat.AB_effects_sub<-subset(SNP.sig.cat.AB_effects, measure==paste0(m))
  
  for (s in SNPs) {
    s<-as.character(s)
    match<-grep(s, SNP.sig.cat.AB_effects_sub$SNP.name, value = T)
    
    if(length(match)>0){
      
      SNP.sig.cat.AB_effects_sub2<-subset(SNP.sig.cat.AB_effects_sub, SNP.name==paste0(s))
      SNP.sig.cat.AB_effects_sub3<-subset(SNP.sig.cat.AB_effects_sub2, effect=="effectBB")
      
      allele.match1<-length(grep("^.*_GG$", SNP.sig.cat.AB_effects_sub3$variable))
      allele.match2<-length(grep("^.*_CC$", SNP.sig.cat.AB_effects_sub3$variable))
      
      SNP.sig.cat.AB_effects_sub3$neg.allele<-ifelse(allele.match1 > 0, "G",
                                                     ifelse(allele.match2 > 0, "C", "A"))
      
      SNP.sig.cat.allele.info<-rbind(SNP.sig.cat.allele.info, SNP.sig.cat.AB_effects_sub3 )
      
    }
  }
}



#join info about reducing allele to data frame with recoded effect sizes (effectAA, effectAB, effectBB)
SNPs.sig.cat.effects_info<-join(SNP.sig.cat.AB_effects, 
                                SNP.sig.cat.allele.info[, c("measure", "SNP.name", "neg.allele")])

var_table.comp<-read.table("variance_table_SNP_animal_models_PCs.txt", header=T)
var_table.cat<-subset(var_table.comp, model=="categorical")

co_var.table.fixed_all<-read.table("Variance_covariance_table_fixed_effects_SNP_models_PCs.txt", header=T)

sig.SNP_variances<-data.frame()

#test variabels
# m<-antler_traits[[1]]
# s<-"cela1_red_2_34227169"

for (m in antler_traits) {
  m<-as.character(m)
  SNPs.sig.cat.effects_info_sub<-subset(SNPs.sig.cat.effects_info, measure==paste0(m))
  var_table.cat_sub<-subset(var_table.cat, measure==paste0(m))
  co_var.table.fixed_sub<-subset(co_var.table.fixed_all, measure==paste0(m))
  
  for (s in SNPs) {
    s<-as.character(s)
    match<-grep(s, SNPs.sig.cat.effects_info_sub$SNP.name, value = T)
    
    if(length(match)>0){
      
      SNPs.sig.cat.effects_info_sub2<-subset(SNPs.sig.cat.effects_info_sub, SNP.name==paste0(s))
      var_table.cat_sub2<-subset(var_table.cat_sub, SNP==paste0(s))
      co_var.table.fixed_sub2<-subset(co_var.table.fixed_sub, SNP.name==paste0(s))
      
      #Frequency of B allele
      q<-ifelse(SNPs.sig.cat.effects_info_sub2$minor.allele[1]==SNPs.sig.cat.effects_info_sub2$neg.allele[1],              
                SNPs.sig.cat.effects_info_sub2$minor.allele.freq[1], SNPs.sig.cat.effects_info_sub2$major.allele.freq[1])
      
      #Frequency of A allele
      p<-1-q
      
      #additive genetic variance 
      var_table.cat_Va<-subset(var_table.cat_sub2, variable.names=='vm(ID, grminv)')
      Va<-var_table.cat_Va$component
      SE_Va<-var_table.cat_Va$std.error
      
      #specify values for effectsAA, BB and AB and corresponding variable names in model
      SNP.info.AA<-subset(SNPs.sig.cat.effects_info_sub2, effect=="effectAA")
      SNP.info.AA<-droplevels(SNP.info.AA)
      variabel.name.AA<-as.character(SNP.info.AA$variable)
      effectAA<-SNP.info.AA$effect.solution
      
      SNP.info.AB<-subset(SNPs.sig.cat.effects_info_sub2, effect=="effectAB")
      SNP.info.AB<-droplevels(SNP.info.AB)
      variabel.name.AB<-as.character(SNP.info.AB$variable)
      effectAB<-SNP.info.AB$effect.solution
      
      SNP.info.BB<-subset(SNPs.sig.cat.effects_info_sub2, effect=="effectBB")
      SNP.info.BB<-droplevels(SNP.info.BB)
      variabel.name.BB<-as.character(SNP.info.BB$variable)
      effectBB<-SNP.info.BB$effect.solution
      
      
      if(length(match)==2){
        
        co_var.table.fixed_AA<-subset(co_var.table.fixed_sub2, variable==variabel.name.AA & variable2==variabel.name.AA)      
        VarAA<-co_var.table.fixed_AA$value
        
        co_var.table.fixed_BB<-subset(co_var.table.fixed_sub2, variable==variabel.name.BB & variable2==variabel.name.BB)      
        VarBB<-co_var.table.fixed_BB$value
        
        co_var.table.fixed_AA.BB<-subset(co_var.table.fixed_sub2, variable==variabel.name.AA & variable2==variabel.name.BB)      
        CovAA.BB<-co_var.table.fixed_AA.BB$value
        
        a <- (effectAA - effectBB)/2         # Calculate a
        
        Vq <- 2*p*q*(a)^2                    # Calculate Vq
        VarExplained <- Vq/(Vq + Va)         # Calculate additive genetic variance explained by QTL
        
        total.var.model<-sum(var_table.cat_sub2$component)
        h2.SNP<- Vq/(total.var.model+Vq)
        
        #covariance matrix of SNP genotypes to calculate error around genetic variance explained by SNP
        vcovMatrix <- matrix(data=c(VarAA, CovAA.BB, CovAA.BB, VarBB), nrow=2)   # VarAA, CovAABB, CovAABB, VarBB
        
        
        x1 <- vcovMatrix[1,1]   #  variance of the het effect
        x2 <- vcovMatrix[2,2]   #  variance of the hom effect
        
        beta <- c(effectAA, effectBB)   # create vactor beta for deltamethod
        
        X <- 2*p*q
        Y <- q^2
        
        Vq.se <- deltamethod(~X*(-x2/2 + (-x2/2 + x1)*Y)^2, beta, vcovMatrix)  # standard error
        
      } else {
        
        #specify the variance covariance matrix components of the fixed effects
        
        co_var.table.fixed_AB<-subset(co_var.table.fixed_sub2, variable==variabel.name.AB & variable2==variabel.name.AB)      
        VarAB<-co_var.table.fixed_AB$value
        
        co_var.table.fixed_BB<-subset(co_var.table.fixed_sub2, variable==variabel.name.BB & variable2==variabel.name.BB)      
        VarBB<-co_var.table.fixed_BB$value
        
        co_var.table.fixed_AB.BB<-subset(co_var.table.fixed_sub2, variable==variabel.name.AB & variable2==variabel.name.BB)      
        CovAB.BB<-co_var.table.fixed_AB.BB$value
        
        #calculate genetic variance explained by SNP (determined by additive component 1 and dominance component d)
        
        a <- (effectAA - effectBB)/2     # Calculate a
        d <-  effectAB - a               # Calculate d
        
        Vq <- 2*p*q*(a + d*(q - p))^2        # Calculate Vq
        VarExplained <- Vq/(Vq + Va)         # Calculate additive genetic variance explained by QTL
        
        total.var.model<-sum(var_table.cat_sub2$component)
        h2.SNP<- Vq/(total.var.model+Vq)
        
        #covariance matrix of SNP genotypes to calculate error around genetic variance explained by SNP
        vcovMatrix <- matrix(data=c(VarAB, CovAB.BB, CovAB.BB, VarBB), nrow=2)   # VarAB, CovABBB, CovABBB, VarBB
        
        
        x1 <- vcovMatrix[1,1]   #  variance of the het effect
        x2 <- vcovMatrix[2,2]   #  variance of the hom effect
        
        beta <- c(effectAB, effectBB)   # create vactor beta for deltamethod
        
        X <- 2*p*q
        Y <- q^2
        
        Vq.se <- deltamethod(~X*(-x2/2 + (-x2/2 + x1)*Y)^2, beta, vcovMatrix)  # standard error
        
      }
      
      df.temp<-data.frame(measure=m, SNP.name=s, additive_effect=a, domiance_effect=d, SNP.var=Vq, SNP.var.SE=Vq.se, SE.prop.SNP.var=Vq.se/Vq, SNP.prop.gen=VarExplained, SNP.prop.phen=h2.SNP, 
                          Var.add=Va, SE.var.add=SE_Va, Var.total.model=total.var.model)
      
      sig.SNP_variances<-rbind(sig.SNP_variances, df.temp)
      
    }
  }  
}


SNPs.sig.all.comp<-join(SNPs.sig.cat.freq, sig.SNP_variances)
head(SNPs.sig.all.comp)
SNPs.sig.all.comp<-join(SNPs.sig.all.comp, SNPs.sig.cat.effects_info)
#add wald test data
model_wald.test_all.cat<-read.table("Sig.test_table_categorical_models.top_SNPs_PCs.txt", header=T)
match<-grep("^cela1_red_", model_wald.test_all.cat$variable)
SNPs.cat.wald<-model_wald.test_all.cat[c(match), ]
colnames(SNPs.cat.wald)[1]<-"SNP.name"
SNPs.sig.all.comp<-join(SNPs.sig.all.comp, SNPs.cat.wald )

write.table(SNPs.sig.all.comp, file="SNPs.sig.add.and.cat.complete.PCs.info.txt", sep = "\t", row.names = F, col.names = T)

#join information about SNP variance 
#(and all composite information, like effect sizes used, reducing allele, allele freqiencies etc.)
#to data specific to significant SNPs in categorical and additve models 

SNP.sig.cat.comp<-join(SNPs.sig.cat, SNPs.sig.all.comp)
head(SNP.sig.cat.comp)
length(unique(SNP.sig.cat.comp$SNP.name))
write.table(SNP.sig.cat.comp, file="SNPs.sig.cat.complete.info.PCs.txt", sep = "\t", row.names = F, col.names = T)

head(SNPs.sig.all.comp)
head(SNPs.sig.add)
SNP.sig.add.comp<-join(SNPs.sig.add, SNPs.sig.all.comp[, c(7:22, 25)])
head(SNP.sig.add.comp)

#remove all duplicated rows (for both allele models in additive model rows are replicated 3 times)
SNP.sig.add.comp_allele1<-subset(SNP.sig.add.comp, variable=="SNP_copy_no")
SNP.sig.add.comp_allele1<-subset(SNP.sig.add.comp_allele1, !duplicated(SNP.name))

SNP.sig.add.comp_allele2<-subset(SNP.sig.add.comp, variable=="SNP_copy_no.2")
SNP.sig.add.comp_allele2<-subset(SNP.sig.add.comp_allele2, !duplicated(SNP.name))

SNP.sig.add.comp<-rbind(SNP.sig.add.comp_allele1, SNP.sig.add.comp_allele2)
SNP.sig.add.comp<-arrange(SNP.sig.add.comp, SNP.name)

write.table(SNP.sig.add.comp, file="SNPs.sig.add.complete.info.PCs.txt", sep = "\t", row.names = F, col.names = T)


#quantify and plot maximum effect sizes
SNP.sig.cat.comp<-read.table("SNPs.sig.cat.complete.info.PCs.txt", header=T)
head(SNP.sig.cat.comp)
library(plyr)
#order SNPs according to effect size and pick top SNPs 
antler_traits<-as.list(unique(SNP.sig.cat.comp$measure))
SNPs.high.eff<- data.frame()
#SNPs.min.eff<- data.frame()

for (m in antler_traits) {
  m<-as.character(paste0(m))
  SNP.sig.cat.comp_sub<-subset(SNP.sig.cat.comp, measure==paste0(m))
  length(SNP.sig.cat.comp_sub)
  SNP.sig.cat.comp_sub.ordered<-arrange(SNP.sig.cat.comp_sub, desc(SNP.sig.cat.comp_sub$additive_effect))
  SNPs.high<-SNP.sig.cat.comp_sub.ordered[1, ]
  SNPs.high<-droplevels(SNPs.high)
  SNPs.high.eff<-rbind(SNPs.high.eff, SNPs.high)
  
  #SNP.sig.cat.comp_sub.ordered2<-arrange(SNP.sig.cat.comp_sub, SNP.sig.cat.comp_sub$additive_effect)
  #SNPs.low<-SNP.sig.cat.comp_sub.ordered2[1, ]
  #SNPs.low<-droplevels(SNPs.low)
  #SNPs.min.eff<-rbind(SNPs.min.eff, SNPs.low)
}

SNPs.high.eff_sub<-SNPs.high.eff[, c(-(2:6))]

#get intercepts for top SNPs (to use as baseline for effect size plot)
model_summary_all.cat<-read.table("Summary_table_categorical_models.top_SNPs_PCs.txt", header=T)
head(model_summary_all.cat)

top_SNPs.intercepts<-data.frame()

for (m in antler_traits) {
  m<-as.character(paste0(m))
  model_summary_all.cat_sub<-subset(model_summary_all.cat, measure==paste0(m))
  SNPs.high.eff_sub<-subset(SNPs.high.eff, measure==paste0(m))
  match_keep<-grep(SNPs.high.eff_sub$SNP.name, model_summary_all.cat_sub$variable)
  model_summary_all.cat_sub<-model_summary_all.cat_sub[(max(match_keep)+1), ]
  df<-data.frame(SNP.name=SNPs.high.eff_sub$SNP.name, measure=m, intercept=model_summary_all.cat_sub$solution)
  top_SNPs.intercepts<-rbind(top_SNPs.intercepts, df)
  
}

# min_SNPs.intercepts<-data.frame()
# 
# for (m in antler_traits) {
#   m<-as.character(paste0(m))
#   model_summary_all.cat_sub<-subset(model_summary_all.cat, measure==paste0(m))
#   SNPs.min.eff_sub<-subset(SNPs.min.eff, measure==paste0(m))
#   match_keep<-grep(SNPs.min.eff_sub$SNP.name, model_summary_all.cat_sub$variable)
#   model_summary_all.cat_sub<-model_summary_all.cat_sub[(max(match_keep)+3), ]
#   df<-data.frame(SNP.name=SNPs.min.eff_sub$SNP.name, measure=m, intercept=model_summary_all.cat_sub$solution)
#   min_SNPs.intercepts<-rbind(min_SNPs.intercepts, df)
#   
# }
# 


SNP.sig.cat.comp_top<-data.frame()

for (m in antler_traits){
  m<-as.character(paste0(m))
  SNP.sig.cat.comp_sub<-subset(SNP.sig.cat.comp, measure==paste0(m))
  top_SNPs.intercepts_sub<-subset(top_SNPs.intercepts, measure==paste0(m))
  df<- subset(SNP.sig.cat.comp_sub, SNP.name %in% top_SNPs.intercepts_sub$SNP.name)
  SNP.sig.cat.comp_top<-rbind(SNP.sig.cat.comp_top, df)
}
head(SNP.sig.cat.comp_top)

# SNP.sig.cat.comp_min<-data.frame()
# 
# for (m in antler_traits){
#   m<-as.character(paste0(m))
#   SNP.sig.cat.comp_sub<-subset(SNP.sig.cat.comp, measure==paste0(m))
#   min_SNPs.intercepts_sub<-subset(min_SNPs.intercepts, measure==paste0(m))
#   df<- subset(SNP.sig.cat.comp_sub, SNP.name %in% min_SNPs.intercepts_sub$SNP.name)
#   SNP.sig.cat.comp_min<-rbind(SNP.sig.cat.comp_min, df)
# }
# 
# head(SNP.sig.cat.comp_min)

SNP.sig.cat.comp.int<-join(SNP.sig.cat.comp_top, top_SNPs.intercepts)
SNP.sig.cat.comp.int<-mutate(SNP.sig.cat.comp.int, effect_size="maximum")
head(SNP.sig.cat.comp.int)

# SNP.sig.cat.comp.int2<-join(SNP.sig.cat.comp_min, min_SNPs.intercepts)
# SNP.sig.cat.comp.int2<-mutate(SNP.sig.cat.comp.int2, effect_size="minimum")
# head(SNP.sig.cat.comp.int2)
# 
# SNP.sig.cat.comp.int<-rbind(SNP.sig.cat.comp.int, SNP.sig.cat.comp.int2)

AA_match<-grep("_AA$", SNP.sig.cat.comp.int$variable)
SNP.sig.cat.comp.int<-SNP.sig.cat.comp.int[-AA_match, ]
head(SNP.sig.cat.comp.int)
GG_match<-grep("_GG$", SNP.sig.cat.comp.int$variable)
GGs<-as.data.frame(SNP.sig.cat.comp.int$variable[GG_match])
GGs<-droplevels(GGs)
colnames(GGs)<-"SNP.name"
CC_match<-grep("_CC$", SNP.sig.cat.comp.int$variable)
CCs<-as.data.frame(SNP.sig.cat.comp.int$variable[CC_match])
CCs<-droplevels(CCs)
colnames(CCs)<-"SNP.name"
head(SNP.sig.cat.comp.int)


SNP.sig.cat.comp.int$SNP.status<-ifelse(SNP.sig.cat.comp.int$variable %in% GGs$SNP.name, "homozygote",
                                        ifelse(SNP.sig.cat.comp.int$variable %in% CCs$SNP.name, "homozygote", "heterozygote"))

head(SNP.sig.cat.comp.int)
SNP.sig.cat.comp.hom<-subset(SNP.sig.cat.comp.int, SNP.status=="homozygote")
head(SNP.sig.cat.comp.hom)

SNP.sig.cat.comp.hom$intercept<-ifelse(SNP.sig.cat.comp.hom$major.allele=="A", SNP.sig.cat.comp.hom$intercept, 
                                       SNP.sig.cat.comp.hom$intercept+SNP.sig.cat.comp.hom$solution)

SNP.sig.cat.comp.hom<-mutate(SNP.sig.cat.comp.hom, prop.change=abs(additive_effect/intercept))
head(SNP.sig.cat.comp.hom)
hist(SNP.sig.cat.comp.hom$intercept)
SNP.sig.cat.comp.hom$measure<-as.factor(SNP.sig.cat.comp.hom$measure)
SNP.sig.cat.comp.hom<-SNP.sig.cat.comp.hom[c(1, 4:9, 2:3),  ]

write.table(SNP.sig.cat.comp.hom, file="SNPs.sig.cat.max.add.effect.PCs.txt", sep="\t", col.names = T, row.names = F)

#SNP.sig.cat.comp.hom_max<-subset(SNP.sig.cat.comp.hom, effect_size=="maximum")

p<-ggplot(SNP.sig.cat.comp.hom, aes(fct_inorder(measure), prop.change))+geom_bar(stat = "identity", aes(fill=measure), 
         position = "dodge", width = .8)  + scale_fill_viridis(discrete=T, option = "C") +
  theme(axis.text.x = element_text(size = 14, colour = "black", angle = 90, vjust=0.3, hjust=1), 
        axis.text.y = element_text(size=14, colour = "black"), 
        axis.title = element_text(size=16),legend.text = element_text(size=13),
        legend.title = element_blank())+xlab(NULL) + ylab("proportional change")+ theme(legend.position = "none") 

ggsave("SNPs.sig.animal.model.max.PCs.png", p, scale = 1, width = 25, height=25, units="cm")
