
##Antler GWAS
library(GenABEL)
library(RepeatABEL)
library(plyr)
library(dplyr)

#load data
summary(antler.gen)
#use .gen file that has cow chromosomes instead of linkage group (and right order)
antler.gen.cowChr<-load.gwaa.data(phenofile = "antler_pheno_link_group.txt", 
                           genofile = "antlerAbel.link_group.cowChr.gen")

#use .gen file that has linkage group but has same SNP order (in map file) as cow chromosome .gen file
antler.gen.ordered<-load.gwaa.data(phenofile = "antler_pheno_link_group.txt", 
                           genofile = "antlerAbel.link_group.ordered.gen")

snp_summary.cow<-summary(gtdata(antler.gen.cowChr))
snp_summary.ordered<-summary(gtdata(antler.gen.ordered))

#qc-run marker check

qc<-check.marker(antler.gen.ordered, p.level=0)
qc<-check.marker(antler.gen.cowChr, p.level=0)

summary(qc)
str(qc)

nsnps(antler.gen.ordered) # Number of SNPs
nsnps(antler.gen.cowChr) # Number of SNPs
SNP_list_noqc<-snp.names(antler.gen.ordered)
SNP_list_noqc<-as.data.frame(SNP_list_noqc)
nids(antler.gen.ordered)  # Number of IDs with genotypes

ok.markers<-c(qc$snpok, qc$Xmrkfail)
ok.ids<-qc$idok
antler.gen.qced<-antler.gen.ordered[ok.ids, ok.markers]
qc1<-check.marker(antler.gen.qced, p.level = 0)
antler.gen.qced<-Xfix(antler.gen.qced)
qc2<-check.marker(antler.gen.qced, p.level = 0)
summary(qc2)
ok.markers2<-c(qc2$snpok)
ok.ids2<-qc2$idok
antler.gen.qced2<-antler.gen.ordered[ok.ids2, ok.markers2]

qc3<-check.marker(antler.gen.qced2, p.level = 0)
summary(qc3)
antler.gen.qced2<-Xfix(antler.gen.qced2)

qc4<-check.marker(antler.gen.qced2, p.level = 0)
summary(qc4)
ok.markers3<-c(qc4$snpok)
ok.ids3<-qc4$idok
antler.gen.qced3<-antler.gen.qced2[ok.ids3, ok.markers3]
qc5<-check.marker(antler.gen.qced3, p.level = 0)
summary(qc5)#no errors

fam_pheno<- read.table("antler_pheno_link_group.txt", header=T)
pheno_file_qced<-subset(fam_pheno, id %in% idnames(antler.gen.qced3))
write.table(pheno_file_qced, file="antler_qced_pheno_linkage_group.txt", sep = "\t", row.names = F)
write.table(pheno_file_qced[, c(1,2)], file="antler_qced_linkgroup_ids.txt", sep = "\t", row.names = F, col.names = F, quote = F)

save.gwaa.data(antler.gen.qced3,
               pheno="antler_qced_pheno_linkage_group.txt",
               genofile = "antler_qced_linkage_group.gen")



antler_qced.gen<-load.gwaa.data(phenofile = "antler_qced_pheno_linkage_group.txt", 
                                genofile = "antler_qced_linkage_group.gen")
nsnps(antler_qced.gen)
nids(antler_qced.gen)
SNP_list<-snp.names(antler_qced.gen)
SNP_list<-as.data.frame(SNP_list)
write.table(SNP_list, file = "SNP_list_qced.txt", sep = "\t", row.names = F, col.names = F, quote = F)

nids(antler_qced.gen)
antler.snpsummary <- summary.snp.data(gtdata(antler_qced.gen))
head(antler.snpsummary)


#make kinship matrix (exclude X-linked markers)
selectedSNPs<-sample(autosomal(antler_qced.gen))
antler.gkin<-ibs(antler_qced.gen, snps=selectedSNPs, weight = "freq")
save(antler.gkin, file="antler_gwas_grm_linkage_group_SNPs.RData")
sqrt(4571044)



#load genabel data that passed quality control (only linkage group mapped SNPs)

antler_qced.gen<-load.gwaa.data(phenofile = "antler_qced_pheno_linkage_group.txt", 
                                genofile = "antler_qced_linkage_group.gen")

antler_qced.gen@gtdata@nsnps

#kinship matrix needs to be symmetrical
#first load saved grm: antler_gwas_grm_linkage_group_SNPs.RData
antler.gkin.sym <- antler.gkin
antler.gkin.sym[upper.tri(antler.gkin.sym)] = t(antler.gkin.sym)[upper.tri(antler.gkin.sym)]
antler.gkin.sym <- antler.gkin.sym * 2

# The rGLS() function that RepeatABEL uses can't deal with missing data. For 
# each model, a dataset that is clean of NA values for that trait can be created
# (NB but MUST contain the field "id"):
antler_data<-read.table("antler_data_comp.txt", header=T)
colnames(antler_data)[1]<-"id"
IDs<-antler_data[, 1]
IDs<-unique(IDs)
P.threshold<-1.42*10^-6

keep_ids<-idnames(antler_qced.gen)
keep_ids<-as.data.frame(keep_ids)
antler_ids_kept<-subset(antler_data, id %in% keep_ids$keep_ids)
length(unique(antler_ids_kept$id))
table(is.na(antler_ids_kept$AntlerWt))

#loop for gwas 

antler_data_list<-list()
k=1
cleaned_data_list<-list()
l=1

for ( i in c(7:16)) {
  unit<-cbind(antler_data[1:4], antler_data[i])
  antler_data_list[[k]]<-unit
  k=k+1
}
  for (data in antler_data_list) {
    element<-na.omit(data)
    cleaned_data_list[[l]]<-element
    l=l+1
  }

#code to save each data frame before putting them into list
#unit<-assign(paste0("antler_data_", colnames(antler_data)[5]), cbind(antler_data[1:4], antler_data[i]))
#element<-(paste0("cleaned_data_", colnames(data)[5]), na.omit(data))


#test for prefit function of repeatable - select one data frame from cleaned_data_list
Length_data<-cleaned_data_list[[1]]

prefit<-preFitModel(fixed = Length ~ I(MeCaAge)+ I(MeCaAge^2), random = ~1|id + 1|MeCaYear + 1|BirthYear, 
                    genabel.data = antler_qced.gen, phenotype.data = Length_data, 
                    corStruc = list(id=list("GRM","Ind"), MeCaYear=list("Ind"), BirthYear=list("Ind")),
                    GRM = antler.gkin.sym)

summary(prefit$fitted.hglm)

fmla<-as.formula("Length ~ I(MeCaAge)+ I(MeCaAge^2)")

gwas<- rGLS(fmla,
            genabel.data = antler_qced.gen,
            phenotype.data = Length_data,
            V=prefit$V)

head(results(gwas))
results_gwas_length<-results(gwas)

class(antler_data_AntlerWt$AntlerWt)
cleaned_data_sublist<-list(cleaned_data_list[[1]],cleaned_data_list[[2]] )
class(cleaned_data_list[[1]][,5]) 
#-------------------------------------------------------------------------------------------------------------------
#continue gwas loop

gwas.results<-data.frame()
gwas.results.sig<-data.frame()
gwas.objects.list=list()

m=1

for (df in cleaned_data_list){
  fmla<-as.formula(paste0(colnames(df[5]) , "~ I(MeCaAge)+ I(MeCaAge^2)"))
  
  prefit<-preFitModel(fmla, random = ~1|id + 1|MeCaYear + 1|BirthYear, 
                      genabel.data = antler_qced.gen, phenotype.data = df, 
                      corStruc = list(id=list("GRM","Ind"), MeCaYear=list("Ind"), BirthYear=list("Ind")),
                      GRM = antler.gkin.sym)
  
  gwas<- rGLS(fmla,
              genabel.data = antler_qced.gen,
              phenotype.data = df,
              V=prefit$V)
  
  #get results from gwas (correct for lambda inflation)
  results.measure <- results(gwas)
  print(head(results.measure))
  lambda.measure<-lambda(gwas)$estimate
  results.measure$chi2.1df<-qchisq(results.measure$P1df, 1, lower.tail = F)
  results.measure$chi2.1df_corrected<-results.measure$chi2.1df/lambda.measure
  results.measure$Pc1df<-pchisq(results.measure$chi2.1df_corrected, 1, lower.tail = F)
  results.ordered<-results.measure[order(results.measure$Pc1df),]
  results.measure.sig <- subset(results.measure, Pc1df < P.threshold)
  
  #create temporary data frames with results for each gwas run/loop
  df.temp1<-data.frame(measure=colnames(df)[5], SNP=rownames(results.ordered))
  df.temp2<-results.ordered[, -3 ]
  df.temp2<-as.data.frame(df.temp2, row.names=NULL)
  df.temp2$SNP<-rownames(results.ordered)
  df.temp<-merge(df.temp1, df.temp2)
  df.temp.sig<-data.frame(measure=rep(colnames(df)[5], nrow(results.measure.sig)),SNP_sig=rownames(results.measure.sig), stringsAsFactors = F)
  #combine results from each gwas run in data frame 
  gwas.results<-rbind(gwas.results, df.temp)
  gwas.results.sig<-rbind(gwas.results.sig, df.temp.sig)
  
  #create new gwas object with corrected p-values and make manhatten plots
  gwas.measure_corrected<-Create_gwaa_scan(gwas, results.measure$Pc1df, results.measure$effB)
  assign(paste0("gwas_corrected", ".", colnames(df)[5]), gwas.measure_corrected)
  #plot(gwas.measure_corrected, col=c("blue", "red"), ylim=range(0:7), main=colnames(df)[5])
  #compile new gwas objects in list
  gwas.objects.list[[m]]<-gwas.measure_corrected
  m=m+1
}

head(results.measure)
-log10(P.threshold)

write.table(gwas.results, file="Antler_GWAS_results_linkage_group.txt", sep = "\t", row.names = F)
write.table(gwas.results.sig, file="Antler_GWAS_significant_SNPs.txt", sep = "\t", row.names = F)
save(gwas.objects.list, file="GWAS_objects.RData") 

GWAS_CorCirc<-gwas.objects.list[[2]]



#PC analysis
#load genabel data that passed quality control (only linkage group mapped SNPs)

antler_qced.gen<-load.gwaa.data(phenofile = "antler_qced_pheno_linkage_group.txt", 
                                genofile = "antler_qced_linkage_group.gen")



#kinship matrix needs to be symmetrical
#first load saved grm
load("antler_gwas_grm_linkage_group_SNPs.RData")
antler.gkin.sym <- antler.gkin
antler.gkin.sym[upper.tri(antler.gkin.sym)] = t(antler.gkin.sym)[upper.tri(antler.gkin.sym)]
antler.gkin.sym <- antler.gkin.sym * 2

antler_PCs<-read.table("antler_resid_PCs.txt", header = T)
colnames(antler_PCs)[1]<-"id"
clean.data.PC <- na.omit(antler_PCs)
P.threshold<-1.42*10^-6

PC.gwas.results<-data.frame()
PC.gwas.results.sig<-data.frame()
PC.gwas.objects.list<-list()
n=1

for (pc in c(5:15)){
  fmla<-as.formula(paste0(colnames(clean.data.PC)[pc] , "~ 1"))
  
  prefit<-preFitModel(fmla, random = ~1|id + 1|MeCaYear + 1|BirthYear, 
                      genabel.data = antler_qced.gen, phenotype.data = clean.data.PC, 
                      corStruc = list(id=list("GRM","Ind"), MeCaYear=list("Ind"), BirthYear=list("Ind")),
                      GRM = antler.gkin.sym)
  
  gwas<- rGLS(fmla,
              genabel.data = antler_qced.gen,
              phenotype.data = clean.data.PC,
              V=prefit$V)
  
  #get results from gwas (correct for lambda inflation)
  results.measure <- results(gwas)
  print(head(results.measure))
  lambda.measure<-lambda(gwas)$estimate
  results.measure$chi2.1df<-qchisq(results.measure$P1df, 1, lower.tail = F)
  results.measure$chi2.1df_corrected<-results.measure$chi2.1df/lambda.measure
  results.measure$Pc1df<-pchisq(results.measure$chi2.1df_corrected, 1, lower.tail = F)
  results.ordered<-results.measure[order(results.measure$Pc1df),]
  results.measure.sig <- subset(results.measure, Pc1df < P.threshold)
  
  #create temporary data frames with results for each gwas run/loop
  df.temp1<-data.frame(measure=colnames(clean.data.PC)[pc], SNP=rownames(results.ordered))
  df.temp2<-results.ordered[, -3 ]
  df.temp2<-as.data.frame(df.temp2, row.names=NULL)
  df.temp2$SNP<-rownames(results.ordered)
  df.temp<-merge(df.temp1, df.temp2)
  df.temp.sig<-data.frame(measure=rep(colnames(clean.data.PC)[pc], nrow(results.measure.sig)),SNP_sig=rownames(results.measure.sig), stringsAsFactors = F)
  #combine results from each gwas run in data frame (only top 10 SNPs)
  PC.gwas.results<-rbind(PC.gwas.results, df.temp)
  PC.gwas.results.sig<-rbind(PC.gwas.results.sig, df.temp.sig)
  
  #create new gwas object with corrected p-values and make manhatten plots
  gwas.measure_corrected<-Create_gwaa_scan(gwas, results.measure$Pc1df, results.measure$effB)
  #assign(paste0("gwas_corrected", ".", colnames(clean.data.PC)[pc]), gwas.measure_corrected)
  #plot(gwas.measure_corrected, col=c("blue", "red"), ylim=range(0:7), main=colnames(clean.data.PC)[pc])
  #create list with corrected gwas objects
  PC.gwas.objects.list[[n]]<-gwas.measure_corrected
  n=n+1
  
}

write.table(PC.gwas.results, file="Antler.PC_GWAS_results.txt", sep = "\t", row.names = F)
save(PC.gwas.objects.list, file="PC.GWAS_objects.RData")

plot(gwas_corrected.PC3, col=c("blue", "red"), ylim=range(0:6), main="PC3")
abline(h=-log10(P.threshold), lty="dashed")
abline(h=-log10(1.42e-06), lty="dashed")

#get complete data for SNPs (MAF etc) and combine with GWAS results (antler measures and PCs0)
snp_summary<-gtdata(antler_qced.gen)
snp_table<-summary(snp_summary)
snp_frame<-snp_table[, c(1:14)]
rownames(snp_frame)<-NULL
snp_frame$SNP<-rownames(snp_table)

gwas.results.comp<-merge(gwas.results, snp_frame[, -3])
PC.gwas.results.comp<-merge(PC.gwas.results, snp_frame[, -3])

#add lambda measure to GWAs results tables
gwas.results.comp$lambda<-gwas.results.comp$chi2.1df/gwas.results.comp$chi2.1df_corrected
head(gwas.results.comp)
gwas_results_length<-subset(gwas.results.comp, measure=="Length")

write.table(gwas.results.comp, file="GWAS_results.comp_linkage_group.txt", sep = "\t", row.names = F)

PC.gwas.results.comp$lambda<-PC.gwas.results.comp$chi2.1df/PC.gwas.results.comp$chi2.1df_corrected
head(PC.gwas.results.comp)

write.table(PC.gwas.results.comp, file="PC.GWAS_results.comp.txt", sep = "\t", row.names = F)

#load genable data that only has stags in data set (to calculate allele frequencies in data set)

antler_qced.sub.gen<-load.gwaa.data(phenofile = "antler_qced_pheno_linkage_group.txt", 
                                genofile = "antlerAbel.link_group.sub.gen")
snp_summary<-summary(antler_qced.sub.gen@gtdata)

#loop to make manhattan plots for GWAS results of all traits/PCs
antler_data<-read.table("antler_data_comp.txt", header=T)
GWAS_results <- read.table("GWAS_results.comp_linkage_group.txt", header=T)

#read in linkage map and join with GWAS results data frame to plot SNPs by linkage group
linkage_map<- read.table("Cervus_elaphus_linkage_map_data.ordered.txt", header = T)
colnames(linkage_map)[1]<-"SNP"
library(plyr)
Gwas_results_linkmap<-join(GWAS_results, linkage_map[,c(1,4:5)], by="SNP")
antler_measures<-colnames(antler_data[, c(7:16)])

for (trait in antler_measures){
  df=subset(Gwas_results_linkmap, measure==trait)
  #df$Chromosome<-gsub("X", "30", df$Chromosome) #make pseudoautosomal SNPs part of X chromosome again (for plotting)
  #df$Chromosome<-as.numeric(df$Chromosome)
  df$CEL.LG<-as.numeric(df$CEL.LG)
  source("manhattan_plot_code.R")
  png(filename = paste0("manhattan_plot_", trait, ".png"), width=800, height=600, type = "windows")
  par(mgp = c(2.5, 0.7, 0))
  p<- manhattan_p(df, chr="CEL.LG", bp="CEL.order", p="Pc1df", snp="SNP", 
                  col=c("blue", "red"), chrlabs=c(1:33), suggestiveline = FALSE, 
                  genomewideline = -log10(1.42e-06), ylim=range(0:7), xlab="CEL Linkage Group",
                  cex.axis= 1.5, cex.lab=2, cex.main = 2, las=2)
  
  dev.off()
}


#title(ylab=expression('-log'[10]*italic((p))), line=0, cex.lab=1.2, family="Calibri Light")

antler_data_PC<-read.table("antler_resid_PCs.txt", header=T)
PC.GWAS_results <- read.table("PC.GWAS_results.comp.txt", header=T)
linkage_map<- read.table("Cervus_elaphus_linkage_map_data.ordered.txt", header = T)
colnames(linkage_map)[1]<-"SNP"
library(plyr)
PC.Gwas_results_linkmap<-join(PC.GWAS_results, linkage_map[,c(1,4:5)], by="SNP")

PCs<-colnames(antler_data_PC[, c(5:15)])
source("manhattan_plot_code.R")

for (pc in PCs){
  df=subset(PC.Gwas_results_linkmap, measure==pc)
  #df$Chromosome<-gsub("X", "30", df$Chromosome) #make pseudoautosomal SNPs part of X chromosome again (for plotting)
  #df$Chromosome<-as.numeric(df$Chromosome)
  df$CEL.LG<-as.numeric(df$CEL.LG)
  png(filename = paste0("manhattan_plot_", pc, ".png"), width=800, height=600, type = "windows")
  par(mgp = c(2.5, 0.7, 0))
  manhattan_p(df, chr="CEL.LG", bp="CEL.order", p="Pc1df", snp="SNP", 
              col=c("blue", "red"), chrlabs=c(1:29, "X"), suggestiveline = FALSE, 
              genomewideline = -log10(1.42e-06), ylim=range(0:7), cex.axis=0.6,
              cex.axis= 1.5, cex.lab=2, cex.main = 2, xlab="CEL Linkage Group", las=2)
  dev.off()
}


##### make facet manhattan plots for all antler measures and all PCs  #######
library(plyr)
library(dplyr)
library(ggplot2)

antler_data<-read.table("antler_data_comp.txt", header=T)
GWAS_results <- read.table("GWAS_results.comp_linkage_group.txt", header=T)

#read in linkage map and join with GWAS results data frame to plot SNPs by linkage group
linkage_map<- read.table("Cervus_elaphus_linkage_map_data.ordered.txt", header = T)
colnames(linkage_map)[1]<-"SNP"

Gwas_results_linkmap<-join(GWAS_results, linkage_map[,c(1,4:5)], by="SNP")

Gwas_results_linkmap$CEL.LG<-as.numeric(Gwas_results_linkmap$CEL.LG)
Gwas_results_linkmap<-arrange(Gwas_results_linkmap, CEL.LG, CEL.order)
cum<-nrow(Gwas_results_linkmap)
Gwas_results_linkmap$pos_cuml<-c(1:cum)
Gwas_results_linkmap<-Gwas_results_linkmap[, c("measure", "SNP", "Pc1df", "CEL.LG", "pos_cuml", "CEL.order")]
colnames(Gwas_results_linkmap)[3]<-"p_corrected"

axisdf_gwas <- ddply(Gwas_results_linkmap,.(CEL.LG), summarize, center=(max(pos_cuml) + min(pos_cuml))/2)

label_seq<-axisdf_gwas$CEL.LG[seq(1,length(axisdf_gwas$CEL.LG),by=2)]
label_seq<-append(label_seq, "X")
breaks_seq<-axisdf_gwas$center[seq_along(axisdf_gwas$center)%% 2 > 0]
breaks_seq<-append(breaks_seq, max(axisdf_gwas$center))

#make custom facet labels
measure.labs<-c("a. Antler Length", "b. Coronet Circumference", "c. Lower Beam Circ.", 
                "d. Upper Beam Circ.", "e. Coronet-Brow Junc.", "f. Coronet-Tray Junc.", "g. Brow Length", "h. Tray Length", 
                "i. Antler Weight", "j. Form")
names(measure.labs)<-c("Length", "CoronetCirc", "LowerBeam", "UpperBeam", "CoronetBrowJunc",
                       "CoronetTrayJunc", "BrowLength", "TrayLength", "AntlerWt", "Form" )

#force order of measures to be the same as labels above
measure_levels<-c("Length", "CoronetCirc", "LowerBeam", "UpperBeam", "CoronetBrowJunc",
                  "CoronetTrayJunc", "BrowLength", "TrayLength", "AntlerWt", "Form" )
Gwas_results_linkmap$measure.ordered<-factor(Gwas_results_linkmap$measure, levels = measure_levels)


library(forcats)

p_GWAS_wrap<-ggplot(Gwas_results_linkmap, aes(pos_cuml, -log10(p_corrected))) +
  geom_point(aes(color=as.factor(CEL.LG)), size=3, alpha=0.5)+
  scale_colour_manual(values = rep(c("steelblue3", "red3"), 34))+
  scale_x_continuous(label= label_seq, breaks = breaks_seq)+
  scale_y_continuous(expand = c(0,0))+
  ylim(0,7)+ 
  geom_hline(yintercept = -log10(1.42e-06), linetype="dashed", color="black", size=1)+
  xlab("CEL Linkage Group")+ylab(expression(-log[10]*italic((p))))+
  theme(legend.position = "none", 
        axis.text.x = element_text(size=10, colour = "black"),
        axis.text.y = element_text(size=15, colour = "black"), 
        axis.title.y = element_text(size=35, margin = margin(t = 0, r = 30, b = 0, l = 0)), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=35, margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.ticks.y =element_line(size=1), axis.ticks.x = element_line(size=1),
        strip.text = element_text(size=20, hjust=0.01), strip.background = element_rect(colour="black", fill="lightgray")) +
  facet_wrap(vars(measure.ordered), nrow = 5, ncol = 2, labeller = labeller(measure.ordered=measure.labs))


ggsave("manhattan_plot_GWAS_wrap_ms.png", p_GWAS_wrap, width = 30, height = 40, units = "cm")

#wraped manhattan plot for PCs

PC.GWAS_results <- read.table("PC.GWAS_results.comp.txt", header=T)

#read in linkage map and join with GWAS results data frame to plot SNPs by linkage group
linkage_map<- read.table("Cervus_elaphus_linkage_map_data.ordered.txt", header = T)
colnames(linkage_map)[1]<-"SNP"

PC.Gwas_results_linkmap<-join(PC.GWAS_results, linkage_map[,c(1,4:5)], by="SNP")

PC.Gwas_results_linkmap$CEL.LG<-as.numeric(PC.Gwas_results_linkmap$CEL.LG)
PC.Gwas_results_linkmap<-arrange(PC.Gwas_results_linkmap, CEL.LG, CEL.order)
cum<-nrow(PC.Gwas_results_linkmap)
PC.Gwas_results_linkmap$pos_cuml<-c(1:cum)
PC.Gwas_results_linkmap<-PC.Gwas_results_linkmap[, c("measure", "SNP", "Pc1df", "CEL.LG", "pos_cuml", "CEL.order")]
colnames(PC.Gwas_results_linkmap)[3]<-"p_corrected"

axisdf_gwas <- ddply(PC.Gwas_results_linkmap,.(CEL.LG), summarize, center=(max(pos_cuml) + min(pos_cuml))/2)

label_seq<-axisdf_gwas$CEL.LG[seq(1,length(axisdf_gwas$CEL.LG),by=2)]
label_seq<-append(label_seq, "X")
breaks_seq<-axisdf_gwas$center[seq_along(axisdf_gwas$center)%% 2 > 0]
breaks_seq<-append(breaks_seq, max(axisdf_gwas$center))

#make custom facet labels
measure.labs<-c("a. PC1", "b. PC2", "c. PC3", "d. PC4", "e. PC5", "f. PC6", "g. PC7", "h. PC8", 
                "i. PC9", "j. PC10", "k. PC11")
names(measure.labs)<-c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", 
                       "PC9", "PC10", "PC11")

#force order of measures to be the same as labels above
measure_levels<-c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", 
                  "PC9", "PC10", "PC11")
PC.Gwas_results_linkmap$measure.ordered<-factor(PC.Gwas_results_linkmap$measure, levels = measure_levels)


p_GWAS_wrap.PCs<-ggplot(PC.Gwas_results_linkmap, aes(pos_cuml, -log10(p_corrected))) +
  geom_point(aes(color=as.factor(CEL.LG)), size=3, alpha=0.5)+
  scale_colour_manual(values = rep(c("steelblue3", "red3"), 34))+
  scale_x_continuous(label= label_seq, breaks = breaks_seq)+
  scale_y_continuous(expand = c(0,0))+
  ylim(0,7)+ 
  geom_hline(yintercept = -log10(1.42e-06), linetype="dashed", color="black", size=1)+
  xlab("CEL Linkage Group")+ylab(expression(-log[10]*italic((p))))+
  theme(legend.position = "none", 
        axis.text.x = element_text(size=10, colour = "black"),
        axis.text.y = element_text(size=15, colour = "black"), 
        axis.title.y = element_text(size=35, margin = margin(t = 0, r = 30, b = 0, l = 0)), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=35, margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.ticks.y =element_line(size=1), axis.ticks.x = element_line(size=1),
        strip.text = element_text(size=20, hjust=0.01), strip.background = element_rect(colour="black", fill="lightgray")) +
  facet_wrap(vars(measure.ordered), nrow = 6, ncol = 2, labeller = labeller(measure.ordered=measure.labs))


ggsave("manhattan_plot_GWAS_wrap_PCs_ms.png", p_GWAS_wrap.PCs, width = 30, height = 40, units = "cm")















