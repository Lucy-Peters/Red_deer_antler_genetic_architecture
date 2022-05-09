#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Red deer antler GWAS power analysis   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(plyr)
library(dplyr)
library(GenABEL)
library(RepeatABEL)
library(purrr)
library(furrr)
library(xtable)


#define output directory (scratch)
args <- commandArgs(trailingOnly = TRUE)
resultpath <- as.character(args[1])


antler_data<-read.table("antler_data_comp.txt", header=T)
var_comp<-read.table("Variance_components_h2_animal_model_grm.txt", header=T)
allele_frq<-read.table("Allele_frequencies-GWAS_SNPs.txt", header=T)
geno_tab<-read.delim("Deer31_linkage_group.ped", sep = "", header=F)

geno_tab<-geno_tab[, c(2,7:ncol(geno_tab))]
geno_tab<-data.frame( geno_tab[1], mapply( paste0, geno_tab[-1][c(T,F)], geno_tab[-1][c(F,T)] ) )
geno_map<-read.table("Deer31_linkage_group_names.ordered.map", header=F)
colnames(geno_tab)[-1]<-as.character(geno_map$V2)
colnames(geno_tab)[1]<-"Code"

antler_data<-subset(antler_data, Code %in% geno_tab$Code)
allele_frq<-subset(allele_frq, minor.allele.freq>0.05)


#load genabel data that passed quality control (only linkage group mapped SNPs)
antler_qced.gen<-load.gwaa.data(phenofile = "antler_qced_pheno_linkage_group.txt", 
                                genofile = "antler_qced_linkage_group.gen")

#kinship matrix needs to be symmetrical
#first load saved grm: 
load("antler_gwas_grm_linkage_group_SNPs.RData")
antler.gkin.sym <- antler.gkin
antler.gkin.sym[upper.tri(antler.gkin.sym)] = t(antler.gkin.sym)[upper.tri(antler.gkin.sym)]
antler.gkin.sym <- antler.gkin.sym * 2

P.threshold<-1.42*10^-6
trait<-"AntlerWt"
#h2SNPrange<-c(0.01, 0.05, 0.1, 0.4)
nsim<-100

source("GWASpwrFunc.R")
#run function in safe mode
getGWASpwr_safely<-purrr::safely(getGWASpwr)
#specify no.of cores
plan(multisession, workers=4)

options(future.globals.maxSize = +Inf)

GWASpwrRes<-future_map(c(0.01, 0.05, 0.1, 0.4), getGWASpwr_safely, trait = trait, 
                       antler_data = antler_data, geno_tab = geno_tab, 
                       allele_frq = allele_frq, antler_qced.gen = antler_qced.gen,
                       antler.gkin.sym = antler.gkin.sym, var_comp = var_comp, 
                       P.threshold = P.threshold, nsim = nsim,
                       .options = furrr_options(seed = TRUE))

save(GWASpwrRes, file=paste0(resultpath, "GWASpwrLength.RData"))

#make separate R scripts for each trait to run power analysis in parallele on ashworth server
trait_list<-colnames(antler_data[7:ncol(antler_data)])

x<-readLines("GWASpwrTemplate.R")
f=1
for (l in trait_list){
  write(paste0("trait <- ", "'", l, "'"), paste0("GWASpwr_run", f, ".R"))
  write(x, paste0("GWASpwr_run", f, ".R"), append = T)
  f=f+1
}


GWASpwrAll<-data.frame()
GWASsigAll<-data.frame()
for (trait in trait_list){
  load(paste0("GWASpwr", trait,".RData"))
  
  GWASpwrTrait<-data.frame()
  for (i in c(1:length(GWASpwrRes))) {
    GWASpwrdf<-GWASpwrRes[[i]]$result[[1]]
    GWASpwrTrait<-rbind(GWASpwrTrait, GWASpwrdf)
  }
  GWASpwrAll<-rbind(GWASpwrAll, GWASpwrTrait)
  
  GWASsigTrait<-data.frame()
  for (i in c(1:length(GWASpwrRes))) {
    GWASsigdf<-GWASpwrRes[[i]]$result[[2]]
    GWASsigTrait<-rbind(GWASsigTrait, GWASsigdf)
  }
  GWASsigAll<-rbind(GWASsigAll, GWASsigTrait)
}

write.table(GWASpwrAll, file="GWASpwr_analysis_resultsTbl_all.txt", sep = "\t",
            col.names = T, row.names = F, quote = F)

write.table(GWASsigAll, file="GWASpwr_analysis_sigSNPsAll.txt", sep = "\t", col.names = T,
            row.names = F, quote = F)

print(xtable(GWASpwrAll, type = "latex", 
             display = c( "d", "s", "f", "f", "fg")), 
      file = "GWASpwr_tbl.tex")




