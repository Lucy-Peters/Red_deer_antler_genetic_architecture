#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Red deer antler GWAS power analysis   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(plyr)
library(dplyr)
library(GenABEL)
library(RepeatABEL)
library(purrr)
library(furrr)


#define output directory (scratch)
args <- commandArgs(trailingOnly = TRUE)
resultpath <- as.character(args[1])


antler_data<-read.table("data/antler_data_comp.txt", header=T)
var_comp<-read.table("data/Variance_components_h2_animal_model_grm.txt", header=T)
allele_frq<-read.table("data/Allele_frequencies-GWAS_SNPs.txt", header=T)
geno_tab<-read.delim("data/Deer31_linkage_group.ped", sep = "", header=F)

geno_tab<-geno_tab[, c(2,7:ncol(geno_tab))]
geno_tab<-data.frame( geno_tab[1], mapply( paste0, geno_tab[-1][c(T,F)], geno_tab[-1][c(F,T)] ) )
geno_map<-read.table("data/Deer31_linkage_group_names.ordered.map", header=F)
colnames(geno_tab)[-1]<-as.character(geno_map$V2)
colnames(geno_tab)[1]<-"Code"

antler_data<-subset(antler_data, Code %in% geno_tab$Code)
allele_frq<-subset(allele_frq, minor.allele.freq>0.05)


#load genabel data that passed quality control (only linkage group mapped SNPs)
antler_qced.gen<-load.gwaa.data(phenofile = "data/antler_qced_pheno_linkage_group.txt", 
                                genofile = "data/antler_qced_linkage_group.gen")

#kinship matrix needs to be symmetrical
#first load saved grm: 
load("data/antler_gwas_grm_linkage_group_SNPs.RData")
antler.gkin.sym <- antler.gkin
antler.gkin.sym[upper.tri(antler.gkin.sym)] = t(antler.gkin.sym)[upper.tri(antler.gkin.sym)]
antler.gkin.sym <- antler.gkin.sym * 2

P.threshold<-1.42*10^-6
#trait<-"AntlerWt"
#h2SNPrange<-c(0.2, 0.4)   c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
nsim<-100

source("scripts/GWASpwrFunc.R")
#run function in safe mode
getGWASpwr_safely<-purrr::safely(getGWASpwr)
#specify no.of cores
plan(multisession, workers=20)

options(future.globals.maxSize = +Inf)

GWASpwrRes<-future_map(c(0.01, 0.05, 0.1, 0.4), getGWASpwr_safely, trait = trait, 
                       antler_data = antler_data, geno_tab = geno_tab, 
                       allele_frq = allele_frq, antler_qced.gen = antler_qced.gen,
                       antler.gkin.sym = antler.gkin.sym, var_comp = var_comp, 
                       P.threshold = P.threshold, nsim = nsim,
                       .options = furrr_options(seed = TRUE))



save(GWASpwrRes, file=paste0(resultpath, "GWASpwr", trait, ".RData"))

