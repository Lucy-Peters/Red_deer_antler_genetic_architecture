#---------------------------------# 
# PCA analysis of antler measures #
#---------------------------------#

library(devtools)
library(plyr)
library(dplyr)
library(scales)
library(grid)
library(ggbiplot)
library(pcaMethods)
library(lme4)
library(ggplot2)
library(reshape)
library(grid)

#read in data 
antler_data<-read.table("antler_data_compPCA.txt", header=TRUE)

#subset antler data into antler measures only
antler_measures<-antler_data[,c(1:4, 5:15)]
table((is.na(antler_measures)))

#PCA imputation

#subset antler_measures to exclude entries with too many NAs
antler_measures<-antler_measures[rowSums(is.na(antler_measures)) < 5, ]
class(antler_measures)
#estimate best number of PCs using same baysian method 
no_of_PC<-kEstimate(antler_measures[,c(-1, -2, -3, -4)], method="bpca", evalPcs=1:6, segs = 3, nruncv = 5, em="q2")
checkData(antler_measures, verbose=T)     
no_of_PC#5

#fit PCA to impute missing values with best estimated number of PCs
antler_impute<-pca(antler_measures[, c(-1,-2,-3,-4)], method = "bpca", nPcs=5, scale=c("uv"), centre=TRUE, completeObs = TRUE, cv=c("q2"))
print(antler_impute)
slplot(antler_impute, scoresLoadings=c(TRUE, FALSE), sl=NULL)
cvstat(antler_impute)
summary(antler_impute)
antler_impute_data<-completeObs(antler_impute)
antler_impute_data<-as.data.frame(antler_impute_data)
antler_imputed<- cbind(antler_measures[,c (1:4)], antler_impute_data)
write.table(antler_imputed, file="antler_data_imputed.txt", sep="\t", row.names = F)

#make test data set with no missing values and a subset of that with randomly missing values
#to evaluate error of imputation

#calculate proportion of missing values in data set
table((is.na(antler_measures)))
1104*11 #total count of measure data points
1044/12144 #count of NA over toal no of data points - 0.085
antler_measures_comp<-na.omit(antler_measures[,c(-1,-2,-3,-4)])
antler_measures_miss<-as.data.frame(lapply(antler_measures_comp, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.91, 0.09), size = length(cc), replace = TRUE) ]))
sum(is.na(antler_measures_miss))
#checke no. of missing values in antler_measures_miss - should be about 0.09 of total data point count (11*459)
11*459
0.09*5049
#impute values for fake missing data set
antler_miss_impute<-pca(antler_measures_miss, method = "bpca", nPcs=5, scale=c("uv"), centre=TRUE, completeObs = TRUE, cv=c("q2"))
summary(antler_miss_impute)

imputed<-completeObs(antler_miss_impute)
impute_se<-sum((antler_measures_comp[is.na(antler_measures_miss)] - imputed[is.na(antler_measures_miss)])^2) / sum(antler_measures_comp[is.na(antler_measures_miss)]^2)
#impute_se =  0.01531616

#linear models to account for age and environmetnal effects before PCA
antler_imputed<-read.table("antler_data_imputed.txt", header = T)
antler_imputed$MeCaYear<-as.factor(antler_imputed$MeCaYear)
class(antler_imputed$MeCaYear)

#Form1
hist(antler_imputed$Form1)
antler_imputed$Form1<-as.factor(antler_imputed$Form1)
antler_imputed$Form1<-as.numeric(antler_imputed$Form1)
class(antler_imputed$Form1)
model_form1<-glm(Form1~I(MeCaAge)+I(MeCaAge^2), family=poisson, data=antler_imputed)
summary(model_form1)
plot(model_form1)
hist(model_form1$residuals)
hist(exp(model_form1$residuals))

Tines<-data.frame(Form1=exp(model_form1$residuals))


#Form2
antler_imputed$Form2<-as.numeric(antler_imputed$Form2)
#class(antler_imputed$Form2)
model_form2<-lm(Form2~MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_form2)
#anova(model_form2)
hist(model_form2$residuals)
Points<-data.frame(Form2=model_form2$residuals)

#Form (total no of points)
antler_imputed$Form<-antler_imputed$Form1+antler_imputed$Form2
head(antler_imputed)
class(antler_imputed$Form)
model_form<-lm(Form~MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_form)
hist(model_form$residuals)
Form<-data.frame(Form=model_form$residuals)

#antler weight
model_wt<-lm(AntlerWt~MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_wt)
#anova(model_wt)
hist(model_wt$residuals)
AntlerWt<-data.frame(Weight=model_wt$residuals)

#length
model_le<-lm(Length ~ MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_le)
#anova(model_le)
hist(model_le$residuals)
Length<-data.frame(Length=model_le$residuals)

#CoronetCirc
model_cc<-lm(CoronetCirc~ MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_cc)
#anova(model_cc)
hist(model_cc$residuals)
CoronetCirc<-data.frame(CoronetCirc=model_cc$residuals)

#LowerBeam
model_lb<-lm(LowerBeam~ MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_lb)
hist(model_lb$residuals)
#anova(model_lb)
LowerBeam<-data.frame(LowerBeam=model_lb$residuals)

#UpperBeam
model_ub<-lm(UpperBeam ~ MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_ub)
#anova(model_ub)
hist(model_ub$residuals)
UpperBeam<-data.frame(UpperBeam=model_ub$residuals)

#CoronetBrowJunc
model_cb<-lm(CoronetBrowJunc~ MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_cb)
hist(model_cb$residuals)
#anova(model_cb)
CoronetBrowJunc<-data.frame(CoronetBrowJunc=model_cb$residuals)

#CoronetTrayJunc
model_ct<-lm(CoronetTrayJunc~MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_ct)
hist(model_ct$residuals)
#anova(model_ct)
CoronetTrayJunc<-data.frame(CoronetTrayJunc=model_ct$residuals)

#TrayLength
model_tl<-lm(TrayLength~ MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_tl)
#anova(model_tl)
hist(model_tl$residuals)
TrayLength<-data.frame(TrayLength=model_tl$residuals)

#BrowLength
model_bl<-lm(BrowLength~MeCaAge+I(MeCaAge^2), data=antler_imputed)
summary(model_bl)
hist(model_bl$residuals)
#anova(model_bl)
BrowLength<-data.frame(BrowLength=model_bl$residuals)

#join data sets with residuals into one 
antler_data_resid<-cbind(antler_imputed[, c(1:4)], Tines, Points, AntlerWt, Length, LowerBeam, UpperBeam, CoronetCirc, CoronetBrowJunc, CoronetTrayJunc, TrayLength, BrowLength)
head(antler_data_resid)
write.table(antler_data_resid, file="antler_data_residuals.txt", sep = "\t", row.names = F)
antler_data_resid<-read.table("antler_data_residuals.txt", header=T)


#PCA with  no max PCs on residuals (exclude form but include form1 and form2)
residual_PCA<-prcomp(antler_data_resid[, c(5:15)], centre=T, scale. = T)
summary(residual_PCA) #11 PCs = max no of possible PCs
print(residual_PCA)
#estimate best number of PCs
no_of_PC_resid<-kEstimate(antler_data_resid[,c(5:15)], method="svd", evalPcs=1:10, segs = 3, nruncv = 5, em="q2", allVariables = T)
no_of_PC_resid #10 (i.e max no that can be evaluated)

residual_PCA_summary<-as.data.frame(summary(residual_PCA)$importance)
residual_PCA_summary$variable<-rownames(residual_PCA_summary)
rownames(residual_PCA_summary)<-NULL
residual_PCA_summary<-residual_PCA_summary[, c(12, 1:11)]
write.table(residual_PCA_summary, file="PCA_summary_var_explained.txt", sep = "\t", col.names = T, 
            row.names = F, quote = F)

#calculate how much each variable contributes to a PC
aload <- abs(residual_PCA$rotation)
PC_composition<-sweep(aload, 2, colSums(aload), "/")
PC_composition<-as.data.frame(PC_composition)
colSums(PC_composition)
write.table(PC_composition, "PC_resid_composition.txt", sep="\t", row.names = T, col.names = T)
#create biplot with scores (scatterplot) overlaid by eigenvectors (loadings)
biplot<-ggbiplot(residual_PCA, scale=1, obs.scale = 1, var.scale = 0,  ellipse = T, circle = T)
ggsave("antler_PCs_biplot.png", biplot, scale = 1, width = 25, height=25, units="cm")

#extract PC scores
residual_scores<-as.data.frame(residual_PCA$x)
head(residual_scores)

#save PC scores as data frame
write.table(cbind(antler_data_resid[, c(1:4)], residual_scores), file = "antler_resid_PCs.txt", sep = "\t", row.names = F)
antler_resid_PCs<-read.table("antler_resid_PCs.txt", header = T)

#check normality
hist(residual_scores$PC1)
shapiro.test(residual_scores$PC1)
hist(residual_scores$PC2)
hist(residual_scores$PC3)
hist(residual_scores$PC4)
hist(residual_scores$PC5)
hist(residual_scores$PC6)
hist(residual_scores$PC7)
hist(residual_scores$PC8)
hist(residual_scores$PC9)
hist(residual_scores$PC10)
hist(residual_scores$PC11)

#PC heatmap
PC.comp<-read.table("PC_resid_composition.txt", header = T)
PC.comp<-PC.comp[c("Form1", "Form2", "Weight", "TrayLength", "BrowLength", "CoronetTrayJunc",
                   "CoronetBrowJunc", "UpperBeam", "LowerBeam", "CoronetCirc", "Length"), ]

PC.com_flat<-melt(PC.comp)
PC.com_flat$measure<-c(rep(rownames(PC.comp), times=11))
colnames(PC.com_flat)[c(1:2)]<-c("PC", "contribution" )

class(PC.com_flat$measure)
measure_levels<-PC.com_flat$measure[1:11]
PC.com_flat$measure.ordered<-factor(PC.com_flat$measure, levels = measure_levels)

P<-ggplot(PC.com_flat, aes(PC, measure.ordered))+geom_raster(aes(fill=contribution))+
  scale_fill_gradientn(colours = c("white", "red"))+
  geom_text(aes(label = round(contribution, 2)), size=12) + theme_bw()+
  ylab("Antler measures")+xlab("Principal component")+
  theme(axis.text.x = element_text(size = 34, colour = "black", margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=34, colour = "black", margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.title.x = element_text(size=44, margin = margin(t = 150, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=44, margin = margin(t = 0, r = 30, b = 0, l = 0)),
        legend.title = element_text(size=30, margin = margin(t = -80, r = -240, b = 50, l = 20)),
        legend.text = element_text(size=28), legend.direction = "horizontal", legend.key.width = unit(2, "cm"),
        legend.position = "top") +
  scale_y_discrete(labels=c("Form1", "Form2", "Antler Weight", "Tray Length","Brow Length", "Coronet-Tray Junc.", 
                            "Coronet-Brow Junc.", "Upper Beam Circ.", "Lower Beam Circ.", "Coronet Circumference", "Antler Length"))+
  labs(fill="Contribution")+theme(plot.margin = unit(c(1.5,1.5,1,1), "cm"))



residual_PCA_summary<-read.table("PCA_summary_var_explained.txt", header=T)
variance_explained<-subset(residual_PCA_summary, variable=="Proportion of Variance")
variance_explained<-variance_explained[-1]
rounded_variance<-round(variance_explained, 2)
percent_variance<-rounded_variance*100
print(percent_variance)
percent_variance<-c("41%", "13%", "10%","8%","6%","5%","4%", "4%","4%", "2%", "1%")

png(filename = "PC_heatmap_values.png", width=2000, height=1400, type = "windows")
P2<-P+annotate("text",x=c("PC1","PC2", "PC3", "PC4", "PC5", "PC6",
                          "PC7", "PC8", "PC9", "PC10", "PC11"), 
               y="CoronetTrayJunc",
               label=percent_variance, vjust=24.5, size=12)

gt <- ggplot_gtable(ggplot_build(P2))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
dev.off()



#calculate correlations of PCs with antler measures

antler_resid_PCs<-read.table("antler_resid_PCs.txt", header = T)
antler_data<-read.table("antler_data_compPCA.txt", header=T)

corr.test.results<-data.frame()

for (m in c(5:15)){
  antler_data_me_sub<-na.omit(antler_data[, c(1:4, m)])
  keep_data<-antler_data_me_sub[, c(1:2)]
  antler_resid_PCs_sub<-join(keep_data, antler_resid_PCs)
  antler_resid_PCs_sub<-na.omit(antler_resid_PCs_sub)
  keep_data<-antler_resid_PCs_sub[, c(1:2)]
  antler_data_me_sub<-join(keep_data, antler_data_me_sub)
  paste0(colnames(antler_data_me_sub)[5])
  for (p in c(5:15)){
    c.test<-cor.test(antler_data_me_sub[, 5], antler_resid_PCs_sub[, p] , alternative = "two.sided", method = "pearson")
    df_temp<-data.frame(comparison=paste0(colnames(antler_data_me_sub)[5], "_", colnames(antler_resid_PCs_sub)[p]), t=c.test[[1]], 
                        df=c.test[[2]], p_value=c.test[[3]], c.coeff= c.test[[4]])
    corr.test.results<-rbind(corr.test.results, df_temp)
    plot(antler_resid_PCs_sub[, p], antler_data_me_sub[, 5], 
         main = paste0(colnames(antler_data_me_sub)[5], "vs", colnames(antler_resid_PCs_sub)[p]), 
         ylab = paste0(colnames(antler_data_me_sub)[5]), xlab=paste0(colnames(antler_resid_PCs_sub)[p]))
  }
}

write.table(corr.test.results, file="Correlations_antler_traits_and_PCs.txt", sep = "\t", col.names = T, row.names = F)




class(paste0("antler_data_me_sub$", m))

paste0("antler_data_me_sub$", mes)

#create loadings plot
theta <- seq(0,2*pi,length.out = 100)

circle <- data.frame(x = cos(theta), y = sin(theta))

p <- ggplot(circle,aes(x,y)) + geom_path()



loadings <- data.frame(residual_PCA3$rotation, 
                       
                       .names = row.names(residual_PCA3$rotation))

p + geom_text(data=loadings, 
              
              mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) +
  
  coord_fixed(ratio=1) +
  
  labs(x = "PC1", y = "PC2")


#check how much missing data there is per column and sample

?apply

antler_measures_sub3<-antler_measures_sub[rowSums(is.na(antler_measures_sub))<3, ]

pMiss <- function(x){sum(is.na(x))/length(x)*100}

apply(antler_measures_sub3, 1, pMiss)

apply(antler_measures_sub3, 2, pMiss)


