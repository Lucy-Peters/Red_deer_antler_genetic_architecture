#load libraries
library(OmicKriging)
library(plyr)
library(devtools)
library(pedigreeR)


#read in GRM files from GCTA output (in .bin format) and convert into matrix 
grmFile<-"C:/Users/s1767711/Documents/Red_deer/Antler_genetic_architechture/Antler_model/Deer33.grm"
grmFileBase<-substr(grmFile,1, nchar(grmFile))
grm.matrix<-read_GRMBin(grmFileBase)

#read in pedigree
pedframe<-read.table("Pedigree_Deer_2017-12.txt", header=T)
class(pedframe)
#sort pedigree (ancestors before offspring) and make sure all sire/dam IDs appear in id column
pedEdit<- editPed(sire=pedframe$sire, dam= pedframe$dam, label=pedframe$id) 
#convert pedigree into pedigree class object and calculate relatedness matrix
ped<-with(pedEdit, pedigree(label=label, sire=sire, dam=dam))
pedA<-getA(ped)
pedA<-as.matrix(pedA)


#filter and sort both relatedness matrices so that only common IDs are included and they are ordered in the same way
ped1_cn = colnames(pedA) 
grm_cn = colnames(grm.matrix) 
common_t = sort(intersect(ped1_cn, grm_cn)) 
ped1 = pedA[,common_t] 
grm = grm.matrix[,common_t] 
ped1_d = rownames(pedA) 
grm_d = rownames(grm.matrix) 
common_d = sort(intersect(ped1_d, grm_d)) 
ped1 = ped1[common_d,] 
grm = grm[common_d,] 


#flatten the two matrices
library(reshape)

ped1.flat<-melt(ped1)
names(ped1.flat)[3]<-"ped.r" #denotes origin of relatedness value


grm.flat<-melt(grm)
names(grm.flat)[3]<-"grm.r" #denotes origin of relatedness value


#join the matrices
mat.compare <- join(ped1.flat, grm.flat)

head(mat.compare)

#remove lines where X1 == X2 (so self-self relatedness)

mat.compare <- subset(mat.compare, X1 != X2)
mean(mat.compare$ped.r)
mean(mat.compare$grm.r)
var(mat.compare$ped.r)
var(mat.compare$grm.r)
sd(mat.compare$ped.r)
sd(mat.compare$grm.r)
#plot with a "lm" smooth function for linear regression

library(ggplot2)

p<-ggplot(mat.compare, aes(ped.r, grm.r)) + geom_point(shape=1) +
  xlab("Pedigree relatedness")+ylab("Genomic relatedness")+theme_classic() + 
  theme(axis.text.x = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=18, colour = "black"), 
        axis.title = element_text(size=20))+
  geom_abline(slope = 0.8688  , intercept = -2.144e-02  , color="red", linetype=5)


p2<-p+geom_abline(slope = 1, intercept = 0, color="red")

p3<-p2+geom_abline(slope = 0.8688  , intercept = -2.144e-02  , color="red", linetype=5)

p3

+ stat_smooth(method = "lm", linetype=5, color="black")

ggsave("ped_grm_plot.png", scale = 1, width = 25, height=25, units="cm")

#get output for linear regression (slope, intercept etc)
summary(lm(grm.r ~ ped.r, data = mat.compare))


