# #read in SNP map and pseudoautosomal SNPs
# map<- read.table("Deer31_ok.map", header = F)
# class(map$V2)
# pseudo_snps<-read.table("Pseudoautosomal_SNPs.txt", header = F)
# #exclude pseudoautosomal SNPs from map file
# map_sub<-subset(map, !V2 %in% pseudo_snps$V1)
# #subset map file for pseudoautosomal snps & change x to 30
# pseudo_snps_full<-subset(map, V2 %in% pseudo_snps$V1)
# pseudo_snps_full$V2<-gsub("x", "30", pseudo_snps_full$V2)
# map_sub$V1<-gsub("30", "X", map_sub$V1)
# library(dplyr)
# class(pseudo_snps_full$V1)
# class(map_sub$V1)
# pseudo_snps_full$V1<-as.character(pseudo_snps_full$V1)
# map_new<-full_join(map_sub, pseudo_snps_full)
# class(map_new$V1)
# map_sorted<-arrange(map_new, map_new$V1)
# class(map_sorted$V1)
# class(map_sorted$V2)
# write.table(map_sorted, "Deer31_v2.map", sep="\t", row.names = F, col.names =F, quote=F)


#re-code SNP map file to deer linkage groups
#read in linkage group file
linkage_group_map<-read.table("Cervus_elaphus_linkage_map_data.ordered.txt", header=T)
#make file with SNP list to use in GWAs from linkage group file- need to make ped, map and fam file based on that
SNP_list<-linkage_group_map$SNP.Name
SNP_list[1]
write.table(SNP_list, file="SNP_list.txt", sep="\t", row.names = F, col.names = F, quote = F)
#read in new map file (only SNPs in linkage group file) and recode to linkage groups -for manhattan plots after GWAS
map<-read.table("Deer31_linkage_group.map", header = F)
colnames(map)[2]<-"SNP.Name"
library(plyr)
new_map<-join(map, linkage_group_map[, c(1,4,5)])
new_map<-arrange(new_map,  CEL.LG, CEL.order)
new_map_final<-new_map[,c(5,2,4)]
new_map_final$CEL.LG<-gsub("34", "X", new_map_final$CEL.LG)
write.table(new_map_final, file = "Deer31_linkage_group_names.map", sep="\t", 
            row.names = F, col.names = F, quote=F)



#genabel phenotype data prep- read in antler data
antler_data<-read.table("antler_data_comp.txt", header = T)
antler_data<-read.table("antler_data.txt", header = T)
#make sure only IDs genotyped on chip are in data -read in ID on chip (fam file)
fam<-read.table("Deer31_linkage_group.fam", header = F)
antler_data_sub<-subset(antler_data, Code %in% fam$V2)
table(is.na(antler_data_sub$AntlerWt))
write.table(antler_data_sub, file="antler_data_gwas.txt", sep = "\t", row.names = F, col.names = T)
#make pheno file for gwaas object- family id, id and sex
fam_pheno<-fam[, c(1,2,5)]
colnames(fam_pheno)<-c("fam.id", "id", "sex")
#sex has to be recoded- 0=female, 1=male
fam_pheno$sex<-gsub("2", "0", fam_pheno$sex)
write.table(fam_pheno, "antler_pheno_link_group.txt", sep = "\t", row.names = F)


#convert plink snp files (ped and map) to genable file
library(GenABEL)
#make sure SNP names in map file are in same order for map file with linkage group and map
#file with cow chromosomes
linkage_group_map<-read.table("Deer31_linkage_group_names.map", header=F)
chromosome_map_file<-read.table("Deer31_linkage_group.map", header=F)
#replace chromosome with linkage group but keep SNP order
linkage_group_map_ordered<-join(chromosome_map_file[, -1], linkage_group_map[, c(1:2)])
linkage_group_map_ordered<-linkage_group_map_ordered[, c(4,1:3)]

write.table(linkage_group_map_ordered, file = "Deer31_linkage_group_names.ordered.map", sep="\t", 
            row.names = F, col.names = F, quote=F)

convert.snp.ped("Deer31_linkage_group.ped", "Deer31_linkage_group_names.ordered.map", #don't use linkage group recoded map file, doesn' work
                "antlerAbel.link_group.ordered.gen", mapHasHeaderLine = F)

#.gen file with cow chromosomes (and SNP order)
convert.snp.ped("Deer31_linkage_group.ped", "Deer31_linkage_group.map", #don't use linkage group recoded map file, doesn' work
                "antlerAbel.link_group.cowChr.gen", mapHasHeaderLine = F)

#make genable file that only includes stags in antler data (that passed qc from already run GWAS)
convert.snp.ped("Deer31.sub.ped", "Deer31.sub.map", 
                "antlerAbel.link_group.sub.gen", mapHasHeaderLine = F)

