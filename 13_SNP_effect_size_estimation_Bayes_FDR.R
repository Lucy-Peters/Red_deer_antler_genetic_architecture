library(plyr)
library(dplyr)
library(ashr)
library(ggplot2)
library(viridis)

GWAS_results<-read.table("GWAS_results.comp_linkage_group.txt", header=T)
antler_data<-read.table("antler_data_comp.txt", header=T)


#loop to apply empirical Bayes FDR to GWAS results (effect size & sd) of all 10 antler measures 

antler_measures_list<-list()
k=1

for ( i in c(7:16)) {
  unit<-colnames(antler_data[i])
  antler_measures_list[[k]]<-unit
  k=k+1
}

ash_results_all<-data.frame()
ash_results.sig_all<-data.frame()

for (m in antler_measures_list){
GWAS_results_sub<-subset(GWAS_results, measure==paste0(m)) #subset GWAS results according to antler measures

lambda<-median(GWAS_results_sub$effB, na.rm = T)+1             #calculate a lambda equivalent to account for inflation - 
GWAS_results_sub$effB<-(GWAS_results_sub$effB)/lambda          #expected median 0, so need to +1 to both expected and observed values 
GWAS_results_sub$se_effB<-(GWAS_results_sub$se_effB)/lambda    #(so: (observed median +1)/1 = obs median + 1)

betahat<-GWAS_results_sub$effB
sebetahat<-GWAS_results_sub$se_effB

measure.ash<-ash(betahat, sebetahat, mixcompdist = "uniform", method= "fdr", optmethod="mixEM", 
                 control = list(maxiter = 10000), outputlevel= 3)
ash_results<-measure.ash$result
ash_results<-cbind(GWAS_results_sub[, c("SNP", "measure")], ash_results)
colnames(ash_results)[c(1,2)]<-c("SNP.name", "measure")
ash_results$lambda<-lambda

ash_results.sig<-subset(ash_results, lfsr < 0.05)
ash_results.sig_all<-rbind(ash_results.sig_all, ash_results.sig)

ash_results_all<-rbind(ash_results_all, ash_results)

}


write.table(ash_results_all, file="SNP_effect_size_estimation_results_all.txt", sep="\t", 
            col.names = T, row.names = F)


#count frequency of significant SNPs and join counts to data frame (to see which SNPs influence multiple traits)
SNP.count<-count(ash_results.sig_all, vars = "SNP.name")
ash_results.sig_all<-join(ash_results.sig_all, SNP.count)
colnames(ash_results.sig_all)[14]<-"traits_per_SNP"

#are the significant SNPs shared with SNPs significant for PCs or unique to antler measures?
ash_results.sig_all<-read.table("SNP_effect_size_estimation_results_sig.txt", header=T)
PC.ash_results.sig_all<-read.table("PC.SNP_effect_size_estimation_results_sig.txt", header = T)
ash_results.sig_unique<-subset(ash_results.sig_all, !SNP.name %in% PC.ash_results.sig_all$SNP.name)
length(unique(ash_results.sig_unique$SNP.name))
length(unique(ash_results.sig_all$SNP.name))

ash_results.sig_all$SNP.status<-ifelse(ash_results.sig_all$SNP.name %in% ash_results.sig_unique$SNP.name, "unique", "shared")
head(ash_results.sig_all)

write.table(ash_results.sig_all, file="SNP_effect_size_estimation_results_sig.txt", sep="\t", 
            col.names = T, row.names = F)

#significant SNPs----------------------------------------------------------------------------------------------------------------------------------------------


#make summary table

#are effect sizes of SNPs poitive or negative?
ash_results.sig_all$sign<-ifelse(ash_results.sig_all$PosteriorMean<0, "negative", "positive")
#how many SNPs are positive or negative (per measure)?
ash_results_sig_summary<-ddply(ash_results.sig_all, .(measure, sign), nrow)
sum(ash_results_sig_summary$V1)
colnames(ash_results_sig_summary)[3]<-"no.SNPs"
#how many SNPs are signifcantly associated with each antler measure (per sign of effect)?
ash_results_sig_summary<-join(ash_results_sig_summary, ddply(ash_results_sig_summary, 
                                                       .(measure), summarise,
                                                       total=sum(no.SNPs)))
#what are the minima and maxima of the effect sizes for each antler trait?
ash_results_sig_summary<-join(ash_results_sig_summary, ddply(ash_results.sig_all, .(measure), summarise,
                               max.effect=max(PosteriorMean),
                               min.effect=min(PosteriorMean)))

#what are the values for the quantiles (for 50% and 75%) for effect sizes for each antler measure?
ash_results.sig_quantiles<-data.frame()
antler.traits<-as.list(as.character(unique(ash_results.sig_all$measure)))

for (i in antler.traits){
  ash_results.sig_mes<-subset(ash_results.sig_all, measure==i)
  q<-quantile(ash_results.sig_mes$PosteriorMean, probs=c(0.05, 0.95))
  df.temp<-data.frame(measure=i, lower.quant=q[1], upper.quant=q[2])
  ash_results.sig_quantiles<-rbind(ash_results.sig_quantiles, df.temp)
}


ash_results_sig_summary<-join(ash_results_sig_summary, ash_results.sig_quantiles)

#how many SNPs (per trait) are unique to antler measures and how many are shared with a PC?

ash_results_SNP.status<-ddply(ash_results.sig_all, .(measure, SNP.status), nrow)
colnames(ash_results_SNP.status)[3]<-"no.SNPs.status"
#Length only has unique SNPs - need to add row with 0 count for shared to data frame
Length_shared<-data.frame(measure="Length", SNP.status="shared", no.SNPs.status=0)
ash_results_SNP.status<-rbind(ash_results_SNP.status, Length_shared)
ash_results_SNP.status<-arrange(ash_results_SNP.status, measure)

ash_results_sig_summary<-cbind(ash_results_sig_summary, ash_results_SNP.status[, -1])

#summary continued below (add no of pleiotropic SNPs)


#filter signifcant SNPs for pleiotropic markers (influence more than one trait)
ash_results.sig_sub<-subset(ash_results.sig_all, traits_per_SNP>1)
length(unique(ash_results.sig_sub$SNP.name))
no_dups<-subset(ash_results.sig_all, !duplicated(SNP.name))
#calciulate a 'proportion' for 
ash_results.sig_sub$prop<-1/(ash_results.sig_sub$traits_per_SNP)

#make an 'ordered variable of antler measure so that this order will be used in plot
measure_levels<-c("Form", "Length" , "BrowLength", "TrayLength", "CoronetBrowJunc", 
                  "CoronetTrayJunc", "CoronetCirc", "LowerBeam")
ash_results.sig_sub$measure.ordered<-factor(ash_results.sig_sub$measure, levels = measure_levels)

#plot SNPs against 'proportion' of traits they influence (only SNPs asociated with 2 or more traits as above)
p.prop<-ggplot(ash_results.sig_sub, aes(SNP.name, prop))+geom_bar(aes(fill=measure.ordered), stat = "identity")+theme_classic()+
  scale_fill_viridis(discrete=T, option = "D", name="Antler measure", labels=c("Form","Length", "Brow Length", "Tray Length",
                                                                               "Coronet-Brow-Junction", "Coronet-Tray-Junction", 
                                                                               "Coronet Circumference", "Lower Beam"))+
  theme(axis.text.x = element_text(size = 6, colour = "black", angle = 90, vjust=0.1, hjust=1), 
        axis.text.y = element_blank(), 
        axis.title = element_blank())
  

ggsave("SNP_contribution_to_trait.png", p.prop, scale=1, width=25, height=25, unit="cm")  

#make category to group and order antler traits into types of measure
ash_results.sig_sub$Trait.cat<-ifelse(ash_results.sig_sub$measure == "Browlength", "length",
                                  ifelse(ash_results.sig_sub$measure == "CoronetBrowJunc", "length",
                                    ifelse(ash_results.sig_sub$measure == "CoronetCirc", "circumference", 
                                      ifelse(ash_results.sig_sub$measure == "CoronetBrowJunc", "length",
                                        ifelse(ash_results.sig_sub$measure == "Form", "form",
                                          ifelse(ash_results.sig_sub$measure == "Length", "length",
                                            ifelse(ash_results.sig_sub$measure == "LowerBeam", "circumference",
                                              ifelse(ash_results.sig_sub$measure == "TrayLength", "length", "circumference"))))))))


measure.count<-ddply(ash_results.sig_sub, .(measure), nrow)
colnames(measure.count)[2]<-"pleiotropic.SNPs"
ash_results_sig_summary<-join(ash_results_sig_summary, measure.count)
ash_results_sig_summary$prop.plei.SNPs<-(ash_results_sig_summary$pleiotropic.SNPs)/(ash_results_sig_summary$total)
plot(ash_results_sig_summary$total, ash_results_sig_summary$prop.plei.SNPs)

write.table(ash_results_sig_summary, file="ash_sig_SNP_summary.txt", sep="\t", row.names = F,
            col.names = T)

#mean standardise effect size values (max, min, quantiles) in ash results summary data frame

antler_summaries<-read.table("antler_summaries.txt", header=T)
colnames(antler_summaries)[1]<-"measure"

ash_results_sig_summary_mean_std<-ash_results_sig_summary %>%
  select(c(measure, max.effect, min.effect, upper.quant, lower.quant))%>%
  join(., antler_summaries[, c("measure", "Mean")], by="measure")%>%
  mutate(max.effect=max.effect/Mean, min.effect=min.effect/Mean,
         upper.quant=upper.quant/Mean, lower.quant=lower.quant/Mean)%>%
  select(-Mean)%>%
  join(., ash_results_sig_summary[, c("measure", "sign", "no.SNPs", "total", "SNP.status", 
                                     "no.SNPs.status", "pleiotropic.SNPs", "prop.plei.SNPs" )], by="measure")
  
ash_results_sig_summary_mean_std<-ash_results_sig_summary_mean_std[c(1:2, 5:6, 9:10, 13:14, 17:18, 21:22,
                                                                     25:26, 29:30), ]

write.table(ash_results_sig_summary_mean_std, file="ash_sig_SNP_summary_mean.stand.txt", 
            sep="\t", row.names = F, col.names = T)

ash_results.sig_sub<-join(ash_results.sig_sub, measure.count)
ash_results.sig_sub$prop.measure<-1/(ash_results.sig_sub$pleiotropic.SNPs)


# ggplot(ash_results.sig_sub, aes(measure.ordered, SNP.name))+geom_raster(aes(fill=PosteriorMean))+
#   scale_fill_gradient(low = "red", high = "blue")+ theme_classic()+
#   theme(axis.text.x = element_text(size = 12, colour = "black", angle = 90, vjust=0.3, hjust=1), 
#         axis.title = element_blank())
# 
# 
# #make new measure variable that is ordered in way it should be displayed on y-axis (to force ggplot to stick to that order)
# measure_levels<-c("AntlerWt", "Form", "Length" , "BrowLength", "TrayLength", "CoronetBrowJunc", 
#                   "CoronetTrayJunc", "CoronetCirc", "LowerBeam", "UpperBeam")
# ash_results.sig_all$measure.ordered<-factor(ash_results.sig_all$measure, levels = measure_levels)
# 
# ggplot(ash_results.sig_all, aes(SNP.name, measure.ordered))+geom_raster(aes(fill=PosteriorMean))+
#   scale_fill_gradient2(low = "red", mid="yellow", high = "blue", midpoint = 0, guide = "legend")+ theme_classic()+
#   theme(axis.text.x = element_text(size = 1, colour = "black", angle = 90, vjust=0.3, hjust=1), 
#         axis.title = element_blank())
# 

#split data into types of measures to display effect sizes because scale is so different (and plots cannot 
#distinguish effect sizes properly)
head(ash_results.sig_all)

ash_results.sig_lengths<-subset(ash_results.sig_all, !measure== "Form" & 
                              !measure== "CoronetCirc" & !measure=="LowerBeam")

head(ash_results.sig_lengths)
#make an 'ordered variable of antler measure so that this order will be used in plot
measure_levels<-c("Form", "Length" , "BrowLength", "TrayLength", "CoronetBrowJunc", 
                  "CoronetTrayJunc", "CoronetCirc", "LowerBeam")
ash_results.sig_lengths$measure.ordered<-factor(ash_results.sig_lengths$measure, levels = measure_levels)


library(forcats)
max(ash_results.sig_lengths$PosteriorMean)

p.eff.lengths<-ggplot(ash_results.sig_lengths, aes(SNP.name, measure.ordered))+geom_tile(aes(fill=PosteriorMean))+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", name="Posterior mean")+
  theme(axis.text.x = element_text(size = 1, colour = "black", angle = 90, vjust=0.2, hjust=1), 
        axis.text.y = element_text(size=18, colour = "black"), axis.ticks.y = element_line(size = 1),axis.title = element_blank(), 
        legend.title = element_text(size=14), legend.text = element_text(size=12))+
  scale_y_discrete(labels=c("Length", "Brow Length", "Tray Length", "Coronet-Brow-Junction","Coronet-Tray-Junction"))

ggsave("Sig_SNPs_post.effect.sizes_lengths.png", p.eff.lengths, scale=1, width = 25, height = 25, units = "cm")


ash_results.sig_circ<-subset(ash_results.sig_all, !measure== "Form" & !measure== "Length" & !measure=="BrowLength"& !measure== "TrayLength" &
                               !measure=="CoronetBrowJunc" & !measure== "CoronetTrayJunc")

max(ash_results.sig_circ$PosteriorMean)

p.eff.circ<-ggplot(ash_results.sig_circ, aes(SNP.name, measure))+geom_tile(aes(fill=PosteriorMean))+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", name="Posterior mean")+
  theme(axis.text.x = element_text(size = 6, colour = "black", angle = 90, vjust=0.2, hjust=1), 
        axis.text.y = element_text(size=18, colour = "black"), axis.ticks.y = element_line(size = 1),axis.title = element_blank(), 
        legend.title = element_text(size=14), legend.text = element_text(size=12))+
  scale_y_discrete(labels=c("Coronet Circumference", "Lower Beam"))
  

ggsave("Sig_SNPs_post.effect.sizes_circ.png", p.eff.circ, scale=1, width = 25, height = 25, units = "cm")


#for form make barplots (looks better because only 1 measure)

ash_results.sig_form<-subset(ash_results.sig_all, measure=="Form")

p.eff.form<-ggplot(ash_results.sig_form, aes(SNP.name, PosteriorMean))+
  geom_bar(aes(fill=ifelse(PosteriorMean>0, "positive", "negative")), stat = "identity")+ 
  scale_fill_manual(name=NULL, values=c("red", "mediumblue"))+
  ylab("Posterior mean") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 5, colour = "black", 
                                                                   angle = 90, vjust=0.2, hjust=1))


ggsave("Sig_SNPs_post.effect.sizes_form.png", p.eff.form, scale=1, width = 25, height = 25, units = "cm")


#all SNPs (significant or not)----------------------------------------------------------------------------------------------------------------------------------------------------------

#make density plots of posterior effect size distribution - split by category because measures have very different ranges/units

library(ggplot2)
library(viridis)

#plot all effect sizes of antler length measures

ash_results_lengths<-subset(ash_results_all, !measure=="AntlerWt" & !measure== "Form" & 
                              !measure== "CoronetCirc" & !measure=="LowerBeam"& !measure== "UpperBeam")

max(ash_results_lengths$PosteriorMean)
min(ash_results_lengths$PosteriorMean)

p.lengths<-ggplot(data = ash_results_lengths, aes(PosteriorMean)) + geom_density(aes(fill=measure))+xlim(-1,1) + 
  scale_fill_viridis(discrete=T, option = "C", alpha=0.5)+xlab("Posterior effect size")+
  ylab("Density")+theme(legend.title=element_blank())


ggsave("SNP_effect_size_distribution_lengths.png", p.lengths, scale = 1,
       width = 25, height = 25, units = "cm")


#plot all effect sizes of antler circumference measures

ash_results_circ<-subset(ash_results_all, !measure=="AntlerWt" & !measure== "Form" & !measure== "TrayLength" &
                           !measure== "CoronetTrayJunc" & !measure=="CoronetBrowJunc"& !measure== "Length" & 
                           !measure== "BrowLength")

max(ash_results_circ$PosteriorMean)
min(ash_results_circ$PosteriorMean)

p.circ<-ggplot(data = ash_results_circ, aes(PosteriorMean)) + geom_density(aes(fill=measure))+xlim(-0.5,0.5) + 
  scale_fill_viridis(discrete=T, option = "C", alpha=0.5)+xlab("Posterior effect size")+
  ylab("Density")+theme(legend.title=element_blank())

ggsave("SNP_effect_size_distribution_circumferences.png", p.circ, scale = 1,
       width = 25, height = 25, units = "cm")


#plot all effect sizes of antler weight

ash_results_wt<-subset(ash_results_all, measure=="AntlerWt")

p.wt<-ggplot(data = ash_results_wt, aes(PosteriorMean)) + geom_density(aes(fill=measure), alpha=0.6) +scale_fill_manual(values="skyblue")+
  xlab("Posterior effect size")+
  ylab("Density")+theme(legend.position="none")

max(ash_results_wt$PosteriorMean)

ggsave("SNP_effect_size_distribution_weigth.png", p.wt, scale = 1,
       width = 25, height = 25, units = "cm")


#plot all effect sizes of antler form

ash_results_form<-subset(ash_results_all, measure=="Form")

p.form<-ggplot(data = ash_results_form, aes(PosteriorMean)) + geom_density(aes(fill=measure), alpha=0.6) +xlab("Posterior effect size")+
  scale_fill_manual(values = "royalblue")+ylab("Density")+xlim(-0.2,0.2)+theme(legend.position = "none")


ggsave("SNP_effect_size_distribution_form.png", p.form, scale = 1,
       width = 25, height = 25, units = "cm")

#loop to make data frames for plots that show shrinkage of effect sizes after FDR/re-estimation of 
#effect sizes using unimodal prior in ashr

beta_post.est_all<-data.frame()

for (m in antler_measures_list){
  ash_results_sub<-subset(ash_results_all, measure==paste0(m))

  betahat_vec<-as.data.frame(ash_results_sub[, c("betahat")])
  betahat_vec$variable<-"GWAS.est"
  colnames(betahat_vec)[1]<-"effect.size"
  
  posterior_vec<-as.data.frame(ash_results_sub[, c("PosteriorMean")])
  posterior_vec$variable<-"posterior.est"
  colnames(posterior_vec)[1]<-"effect.size"
  
  beta_post.est<-rbind(betahat_vec, posterior_vec)
  beta_post.est$measure<-paste0(m)
  
  assign(paste0("beta_post.est_", m), beta_post.est)
  #beta_post.est_all<-rbind(beta_post.est_all, beta_post.est)

}  


library(ggplot2)
library(viridis)

# have to make plots outside loop - won't plot or save inside loop  
p<-ggplot(data=beta_post.est_TrayLength, aes(effect.size))+geom_density(aes(fill=variable))+
    scale_fill_viridis(discrete=T, option = "C", alpha=0.5)+xlim(-2, 2)+xlab("Effect size")+
    ylab("Density")+theme(legend.title=element_blank())
p  
max(beta_post.est_CoronetBrowJunc$effect.size, na.rm = T)

ggsave("SNP_effect_size_shrinkage_TrayLength.png", p, scale = 1, width = 25, height = 25, units = "cm")

#PC FDR -------------------------------------------------------------------------------------------------------------

library(plyr)
library(dplyr)
library(ashr)

PC.GWAS_results<-read.table("PC.GWAS_results.comp.txt", header=T)
antler_data<-read.table("antler_resid_PCs.txt", header=T)


#loop to apply empirical Bayes FDR to GWAS results (effect size & sd) of all 11 antler PCs

antler_measures_list<-list()
k=1

for ( i in c(5:15)) {
  unit<-colnames(antler_data[i])
  antler_measures_list[[k]]<-unit
  k=k+1
}

PC.ash_results_all<-data.frame()
PC.ash_results.sig_all<-data.frame()

for (m in antler_measures_list){
  PC.GWAS_results_sub<-subset(PC.GWAS_results, measure==paste0(m)) #subset GWAS results according to antler measures
  
  lambda<-median(PC.GWAS_results_sub$effB, na.rm = T)+1             #calculate a lambda equivalent to account for inflation - 
  PC.GWAS_results_sub$effB<-(PC.GWAS_results_sub$effB)/lambda          #expected median 0, so need to +1 to both expected and observed values 
  PC.GWAS_results_sub$se_effB<-(PC.GWAS_results_sub$se_effB)/lambda    #(so: (observed median +1)/1 = obs median + 1)
  
  betahat<-PC.GWAS_results_sub$effB
  sebetahat<-PC.GWAS_results_sub$se_effB
  
  measure.ash<-ash(betahat, sebetahat, mixcompdist = "uniform", method= "fdr", optmethod="mixEM", 
                   control = list(maxiter = 100000), outputlevel= 3)
  ash_results<-measure.ash$result
  ash_results<-cbind(PC.GWAS_results_sub[, c("SNP", "measure")], ash_results)
  colnames(ash_results)[c(1,2)]<-c("SNP.name", "measure")
  ash_results$lambda<-lambda
  
  ash_results.sig<-subset(ash_results, lfsr < 0.05)
  PC.ash_results.sig_all<-rbind(PC.ash_results.sig_all, ash_results.sig)
  
  PC.ash_results_all<-rbind(PC.ash_results_all, ash_results)
  
}


write.table(PC.ash_results_all, file="PC.SNP_effect_size_estimation_results_all.txt", sep="\t", 
            col.names = T, row.names = F)

write.table(PC.ash_results.sig_all, file="PC.SNP_effect_size_estimation_results_sig.txt", sep="\t", 
            col.names = T, row.names = F)

#count frequency of significant SNPs and join counts to data frame (to see which SNPs influence multiple traits)
SNP.count<-PC.ash_results.sig_all %>%
  group_by(SNP.name)%>%
  summarise(freq=n())
  
PC.ash_results.sig_all<-join(PC.ash_results.sig_all, SNP.count)
colnames(PC.ash_results.sig_all)[14]<-"traits_per_SNP"


#check how many significant SNPs are shared between antler measures and PCs
ash_results.sig<-read.table("SNP_effect_size_estimation_results_sig.txt", header=T)
ash_results.sig_shared<-subset(PC.ash_results.sig_all, SNP.name %in% ash_results.sig$SNP.name)
ash_results.sig_unique<-subset(PC.ash_results.sig_all, !SNP.name %in% ash_results.sig$SNP.name)

PC.ash_results.sig_all$SNP.status<-ifelse(PC.ash_results.sig_all$SNP.name %in% ash_results.sig_unique$SNP.name, 
                                          "unique", "shared")

head(PC.ash_results.sig_all)



write.table(PC.ash_results.sig_all, file="PC.SNP_effect_size_estimation_results_sig.txt", sep="\t", 
            col.names = T, row.names = F)

#significant SNPs----------------------------------------------------------------------------------------------------------------------------------------------


#make summary table

PC.ash_results.sig_all$sign<-ifelse(PC.ash_results.sig_all$PosteriorMean<0, "negative", "positive")

PC.ash_results_sig_summary<-ddply(PC.ash_results.sig_all, .(measure, sign), nrow)
sum(PC.ash_results_sig_summary$V1)
colnames(PC.ash_results_sig_summary)[3]<-"no.SNPs"
PC.ash_results_sig_summary<-join(PC.ash_results_sig_summary, ddply(PC.ash_results_sig_summary, 
                                                             .(measure), summarise,
                                                             total=sum(no.SNPs)))

PC.ash_results_sig_summary<-join(PC.ash_results_sig_summary, ddply(PC.ash_results.sig_all, .(measure), summarise,
                                                             max.effect=max(PosteriorMean),
                                                             min.effect=min(PosteriorMean)))

#add 5% and 95% quantiles
PC.ash_results.sig_quantiles<-data.frame()
antler.traits<-as.list(as.character(unique(PC.ash_results.sig_all$measure)))

for (i in antler.traits){
  ash_results.sig_mes<-subset(PC.ash_results.sig_all, measure==i)
  q<-quantile(ash_results.sig_mes$PosteriorMean, probs = c(0.05, 0.95))
  df.temp<-data.frame(measure=i, lower.quant=q[1], upper.quant=q[2])
  PC.ash_results.sig_quantiles<-rbind(PC.ash_results.sig_quantiles, df.temp)
}


PC.ash_results_sig_summary<-join(PC.ash_results_sig_summary, PC.ash_results.sig_quantiles)

#how many SNPs (per trait) are unique to a PC and how many are shared with antler measures?

PC.ash_results_SNP.status<-ddply(PC.ash_results.sig_all, .(measure, SNP.status), nrow)
colnames(PC.ash_results_SNP.status)[3]<-"no.SNPs.status"

#PC6, PC7, PC8, PC10, PC11 only have has unique SNPs but positive and negative effect sizes , so need to add 0 count for shared SNPs 
#before data can be joined to summary table 
missing_status_df<-data.frame(measure=c("PC7","PC10", "PC11"), SNP.status="shared", no.SNPs.status=0)
PC.ash_results_SNP.status<-rbind(PC.ash_results_SNP.status, missing_status_df)
PC.ash_results_SNP.status<-arrange(PC.ash_results_SNP.status, measure)

#join SNP status data to summary data
PC.ash_results_sig_summary<-cbind(PC.ash_results_sig_summary, PC.ash_results_SNP.status[, -1])
#PC.ash_results_sig_summary<-PC.ash_results_sig_summary[, -9]

#filter signifcant SNPs for pleiotropic markers (influence more than one trait)
PC.ash_results.sig_sub<-subset(PC.ash_results.sig_all, traits_per_SNP>1) #no pleiotropic SNPs

write.table(PC.ash_results_sig_summary, file="PC.ash_sig_SNP_summary.txt", sep="\t", row.names = F, col.names = T)

#make new measure variable that is ordered in way it should be displayed on y-axis (to force ggplot to stick to that order)
measure_levels<-c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11")
PC.ash_results.sig_all$measure.ordered<-factor(PC.ash_results.sig_all$measure, levels = measure_levels)

library(ggplot2)

p.sig.pc<-ggplot(PC.ash_results.sig_all, aes(SNP.name, measure.ordered))+geom_tile(aes(fill=PosteriorMean))+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", name="Posterior mean")+
  theme(axis.text.x = element_text(size = 1, colour = "black", angle = 90, vjust=0.2, hjust=1), axis.title = element_blank())

ggsave("Sig_SNPs_post.effect.sizes_PCs.png", p.sig.pc, scale = 1,
       width = 25, height = 25, units = "cm")

max(antler_data$PC1)
median(antler_data$PC1)
min(antler_data$PC1)



max(antler_data$PC11)
median(antler_data$PC11)
min(antler_data$PC11)


#------------------------------------------------------------------------------------------------------------------------------------------------

#make file with top 10 SNPs with highest effect sizes

PC.ash_results.sig_all<-read.table("PC.SNP_effect_size_estimation_results_sig.txt", header = T)
ash_results.sig_all<-read.table("SNP_effect_size_estimation_results_sig.txt", header=T)

library(plyr)
library(dplyr)

top_SNPs<-data.frame()

antler_traits_list<-as.list(unique(ash_results.sig_all$measure))

#antler measure top SNPs

for (t in antler_traits_list){
  t<-as.character(t)
  ash_results.sig_sub<-subset(ash_results.sig_all, measure==t)
  abs.post.mean<-order(abs(ash_results.sig_sub$PosteriorMean),decreasing=T)
  rows.top.snps<-abs.post.mean[c(1:10)]
  snps<-as.data.frame(ash_results.sig_sub[rows.top.snps, 1])
  top_SNPs<-rbind(top_SNPs, snps)
}

top_SNPs<-na.omit(top_SNPs)
colnames(top_SNPs)<-"SNP.name"
top_SNPs<-subset(top_SNPs, !duplicated(SNP.name))

write.table(top_SNPs, file="SNPs_high_effects_traits.txt", sep="\t", row.names=F, col.names=F, quote = F)

#for PCs, only look at unique SNPs
PC.ash_results.sig_all<-subset(PC.ash_results.sig_all, SNP.status=="unique")
antler_traits_list<-as.list(unique(PC.ash_results.sig_all$measure))

top_SNPs<-data.frame()

for (t in antler_traits_list){
  t<-as.character(t)
  PC.ash_results.sig_sub<-subset(PC.ash_results.sig_all, measure==t)
  rows.ordered<-order(abs(PC.ash_results.sig_sub$PosteriorMean), decreasing=T)
  top.rows<-rows.ordered[c(1:10)]
  top.rows<-na.omit(top.rows)
  snps<-as.data.frame(PC.ash_results.sig_sub$SNP.name[c(top.rows)])
  top_SNPs<-rbind(top_SNPs, snps)
}

top_SNPs<-na.omit(top_SNPs)
colnames(top_SNPs)<-"SNP.name"
top_SNPs<-subset(top_SNPs, !duplicated(SNP.name))

write.table(top_SNPs, file="SNPs_high_effects_PCs.txt", sep="\t", row.names=F, col.names=F, quote = F)


#get allele frequencies of SNPs used in GWAS
GWAS_results<-read.table("GWAS_results.comp_linkage_group.txt", header=T)
head(GWAS_results)
SNP.frq<-GWAS_results[, c("SNP", "A1", "A2", "Q.2")]
SNP.frq<-subset(SNP.frq, !duplicated(SNP))
library(plyr)
library(dplyr)
SNP.frq<-mutate(SNP.frq, major.allele.freq=1-Q.2)
colnames(SNP.frq)[c(1:4)]<-c("SNP.name", "major.allele", "minor.allele", "minor.allele.freq")
head(SNP.frq)
write.table(SNP.frq, file="Allele_frequencies-GWAS_SNPs.txt", sep = "\t", row.names = F)
