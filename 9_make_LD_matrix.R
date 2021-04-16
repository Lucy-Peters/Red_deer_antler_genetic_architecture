#make LD matrix for all SNPs used in GWAS and reg h2 analysis
library(GenABEL)
library(plyr)
library(reshape)

antler_qced.gen<-load.gwaa.data(phenofile = "antler_qced_pheno_linkage_group.txt", 
                                genofile = "antler_qced_linkage_group.gen")

#use autosomal SNPs only
selectedSNPs<-sample(autosomal(antler_qced.gen))
selectedSNPs<-as.data.frame(selectedSNPs)

#deer linkage map to calculate estimayed Mb distances for each pairwise SNP LD comparison 
linkage_map<-read.table("Cervus_elaphus_linkage_map_data.ordered.txt", header=T)
linkage_map_autosomal_qced<-subset(linkage_map, SNP.Name %in% selectedSNPs$selectedSNPs)


#loop to calculate LD matrix per linkage group and calculate distance from linkage map 
#(after flattening of matrix); join all LD matrices at the end (in flattened format)

LD_matrix_r2_full<-data.frame()

for (i in c(1:33)){
  LG_map<-subset(linkage_map_autosomal_qced, CEL.LG==i)
  SNPs_LG<-as.character(LG_map$SNP.Name)
  LD_matrix_r2<-r2fast(gtdata(antler_qced.gen), snpsubset = SNPs_LG)
  LD_matrix_r2_flat<-melt(LD_matrix_r2)
  colnames(LD_matrix_r2_flat)<-c("SNP1", "SNP2", "r2")
  LD_matrix_r2_flat<-na.omit(LD_matrix_r2_flat)
  LD_matrix_r2_flat<-subset(LD_matrix_r2_flat, r2<=1)
  
  LG_map_SNP1<-LG_map[, c("SNP.Name", "CEL.LG", "Estimated.Mb.Position")]
  colnames(LG_map_SNP1)<-c("SNP1", "CEL.LG","Est.Mb.Pos.SNP1")
  LD_matrix_r2_flat<-join(LD_matrix_r2_flat, LG_map_SNP1, by="SNP1")
  
  LG_map_SNP2<-LG_map[, c("SNP.Name", "Estimated.Mb.Position")]
  colnames(LG_map_SNP2)<-c("SNP2", "Est.Mb.Pos.SNP2")
  LD_matrix_r2_flat<-join(LD_matrix_r2_flat, LG_map_SNP2, by="SNP2")
  
  LD_matrix_r2_flat$Est.Dist.Mb<-abs(LD_matrix_r2_flat$Est.Mb.Pos.SNP1-LD_matrix_r2_flat$Est.Mb.Pos.SNP2)
  LD_matrix_r2_flat$Est.Dist.Mb<-(LD_matrix_r2_flat$Est.Dist.Mb)*1e-06
  LD_matrix_r2_flat<-subset(LD_matrix_r2_flat, !Est.Dist.Mb==0)
  
  LD_matrix_r2_full<-rbind(LD_matrix_r2_full, LD_matrix_r2_flat)
  
}


save(LD_matrix_r2_full, file="LD_matrix_r2_per_LG.RData")


#load LD matrix (per linkage group) - run in eddie

load("LD_matrix_r2_per_LG.RData")

head(LD_matrix_r2_full)
LD_matrix_r2_sub<-LD_matrix_r2_full[c(1:1000), ]
LD_matrix_r2_sub<-subset(LD_matrix_r2_full,Est.Dist.Mb <=1)
LD_matrix_r2_sub<-subset(LD_matrix_r2_sub,!Est.Dist.Mb ==0)


LD_matrix_r2_sub2<-subset(LD_matrix_r2_full, !Est.Dist.Mb > 20)
LD_matrix_r2_sub2<-subset(LD_matrix_r2_sub2,!Est.Dist.Mb ==0)

library(ggplot2)

ggplot(LD_matrix_r2_sub, aes(Est.Dist.Mb, r2))+geom_smooth()+
  xlab("Distance (Mb)")+ylab("LD (r2)")+theme_classic()+
  theme(axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title = element_text(size=16))

geom_point()+stat_smooth(method = "lm", linetype=5, color="red")


#plot in base R and fit  trendline

#make exponential function
x<-LD_matrix_r2_sub$Est.Dist.Mb
y<-LD_matrix_r2_sub$r2
f<-function(x,a,b) {a * exp(b*x)}
#estimate coefficinets
fit<-nls(y~ f(x,a,b), start = c(a=0.1, b=-5))
co<-coef(fit)

plot(LD_matrix_r2_sub2$Est.Dist.Mb, LD_matrix_r2_sub2$r2, xlab = "Distance in Mb", ylab = "LD (r2)", frame.plot = F)

curve(f(x, a=co[1], b=co[2]), add = T, col="red", lwd=2)

#log fit
f2<-function(x,a,b) {a *log(x)+b}                                  
fit2<-nls(y~f2(x,a,b), start=c(a=5 ,b=-5 ))
co2<-coef(fit2)
curve(f2(x, a=co2[1], b=co2[2]), add = T, col="blue", lwd=2)

#log and exponential essentialy identical


# table(LD_matrix_r2_sub$r2)


#smooth scatter/density plot of LD decay

low_density_points<-smoothScatter(LD_matrix_r2_sub$Est.Dist.Mb, LD_matrix_r2_sub$r2, ret.selection = TRUE,
                                  nrpoints=Inf, pch=NA,  xlab = "Distance in Mb", ylab = "LD (r2)", cex=2.5, bty="n")

outlier<-low_density_points[1:100]

text(LD_matrix_r2_sub$Est.Dist.Mb[c(outlier)], LD_matrix_r2_sub$r2[c(outlier)],
     labels = "o" , cex=0.5, col="red", font=2)

LD_matrix_r2_sub$SNP1[c(outlier)]

#make categories in LD_matrix that puts r2 values in groups according to distance
#find cut off values for distance categories - step of 0.05Mb
max<-round(max(LD_matrix_r2_sub$Est.Dist.Mb), 2)
dist.seq<-seq(0.05, max,0.05)
#create data frames for distance groups and the corresponding distance value of each pair-wise comparison
Dist.table<-data.frame()
Dist.table.max<-data.frame()
n=1

for(i in c(dist.seq)){
  for (r in seq(1:nrow(LD_matrix_r2_sub))){
    LD_matrix_r2_temp<-LD_matrix_r2_sub[r, ]
    d<-LD_matrix_r2_temp$Est.Dist.Mb
    if(d>(i-0.05)&&d<=i){
      dg<-n
      Dist.table_sub<-data.frame(Dist.Group=dg, Est.Dist.Mb=d, 
                      SNP1=LD_matrix_r2_temp$SNP1, SNP2=LD_matrix_r2_temp$SNP2)#put in estimated distance value (given by d in loop) as control
      Dist.table<-rbind(Dist.table, Dist.table_sub)
    } 
    if (d>max(dist.seq)){ #if distance value is higher than rounded maximum
      max.d<-length(dist.seq)
      dg.max<-data.frame(Dist.Group=max.d, Est.Dist.Mb=d)
      Dist.table.max<-rbind(Dist.table.max, dg.max)
    }
    
  }
  n= n+1
} 


#if statement to join the two distance tables - the main one and the one including distances 
#above set rounded maximum for distance groups
if(length(Dist.table.max)>0){
Dist.table.max<-subset(Dist.table.max, !duplicated(Dist.table.max))
Dist.table.max$Dist.Group<-gsub(max.d, max(Dist.table$Dist.Group), Dist.table.max$Dist.Group)
Dist.table<-rbind(Dist.table, Dist.table.max)
}

save(Dist.table, file = "LD_distance_group_table_1Mb_0.05_steps.RData")
load("LD_distance_group_table_1Mb_0.05_steps.RData")

#join LD matrix and distance table accordiing to SNP pair
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)


LD_matrix_r2_sub<-join(LD_matrix_r2_sub[, c("SNP1", "SNP2", "r2", "CEL.LG")], Dist.table, by=c("SNP1", "SNP2"))

#calculate mean r2 for each distance group
head(LD_matrix_r2_sub)

LD_matrix_r2_sub_avg<- LD_matrix_r2_sub %>%
  group_by(Dist.Group, CEL.LG)%>%
  mutate(mean_r2_dist_group=round(mean(r2), 4))

  join(LD_matrix_r2_sub, ddply(LD_matrix_r2_sub, .(Dist.Group, CEL.LG), summarise, 
                                                   mean_r2_dist_group=round(mean(r2),4)))

#make data frame to plot average r2 per distance group for each linkage group/chromosome
LD_matrix_r2_LG_plot<-LD_matrix_r2_sub_avg[, c("CEL.LG", "Dist.Group", "mean_r2_dist_group")]

#melt data frame to get entries of r2 for each distance group AND linkage group
LD_matrix_r2_LG_plot<-melt(LD_matrix_r2_LG_plot, id.vars=c("CEL.LG", "Dist.Group"))
#remove dupicated rows
LD_matrix_r2_LG_plot<-LD_matrix_r2_LG_plot[!duplicated(LD_matrix_r2_LG_plot), ]

#rename variables, get rid of generic 'variable'column created by melt
LD_matrix_r2_LG_plot<-LD_matrix_r2_LG_plot%>%
  mutate(mean_r2_dist_group=value, Distance= Dist.Group*0.05)%>%
  select(CEL.LG, Dist.Group, mean_r2_dist_group, Distance)

LD_matrix_r2_LG_plot$CEL.LG<-as.factor(LD_matrix_r2_LG_plot$CEL.LG)


#plot average r2 per distance group (over 1 Mb) with separate point for each linkage group
p_LD_decay_avg<-ggplot(LD_matrix_r2_LG_plot, aes(Distance, mean_r2_dist_group))+geom_point(size=1.5)+
stat_smooth(method="lm", formula = y~ log(x))+ylim(0, 0.3)+xlab("Distance (Mb)")+ylab(expression('mean LD (r'^2*')'))+
  theme_classic()+
  theme(axis.text = element_text(size=20, color="black"),
        axis.title = element_text(size=25, color="black"),
        axis.line = element_line(color = "black"))

ggsave("LD_decay_plot_1Mb_avg_per_LG.png", p_LD_decay_avg, scale = 1, width = 25, height=20, units="cm")

log.model<-lm(mean_r2_dist_group~log(Distance), data=LD_matrix_r2_LG_plot)
summary(log.model)



#plot averages over distance groups
max(LD_matrix_r2_sub_avg$mean_r2_dist_group)
max_group<-max(LD_matrix_r2_sub_avg$Dist.Group)
plot(LD_matrix_r2_sub_avg$Dist.Group, LD_matrix_r2_sub_avg$mean_r2_dist_group, xlab = "Distance group (0.05 Mb interval)", 
     ylab = "LD (r2)", frame.plot = F, xaxt="n", ylim = c(0, 1))
axis(1, at=c(1, round(max_group/8), round(max_group/4), 3*round(max_group/8), 
             round(max_group/2),5*round(max_group/8), 3*round(max_group/4), 7*round(max_group/8), max_group))

#make exponential function for trendline
x<-LD_matrix_r2_sub_avg_single$Dist.Group
y<-LD_matrix_r2_sub_avg_single$mean_r2_dist_group
f<-function(x,a,b) {a * exp(b*x)}
#estimate coefficients
fit<-nls(y~ f(x,a,b), start = c(a=0.1, b=-1))
co<-coef(fit)

curve(f(x, a=co[1], b=co[2]), add = T, col="blue", lwd=2)

barplot(LD_matrix_r2_sub_avg$mean_r2_dist_group, names.arg = LD_matrix_r2_sub_avg$Dist.Group )
hist(LD_matrix_r2_sub_avg$mean_r2_dist_group)
plot(LD_matrix_r2_sub_avg$Dist.Group, LD_matrix_r2_sub_avg$r2)

#calculate LD for all SNP windows used in reg h2 analysis
#get list of SNP window names

load("SNP_windows_list.RData")

#get ID vector from genable object (keep only those also in antler data)
antler_data<-read.table("antler_data_comp.txt", header=T)
colnames(antler_data)[1]<-"id"
IDs<-antler_qced.gen@phdata$id
IDs<-as.data.frame(IDs)
IDs<-subset(IDs, IDs %in% antler_data$id)
IDs<-as.character(IDs$IDs)

#loop to calculate pairwise LD for each window separately and get averages
linkage_map<-read.table("Cervus_elaphus_linkage_map_data.ordered.txt", header=T)
avg_LD_windows<-data.frame()
SNP_list<-snp.names(antler_qced.gen)
SNP_list<-as.data.frame(SNP_list)

for (w in windows_list){
  window<-read.table(paste0("C:/Users/s1767711/Documents/Red_deer/Antler_genetic_architechture/Regional_h2/SNP_windows/", w, ".txt"))
  SNPs<-window$V1
  SNPs.df<-as.data.frame(SNPs)
  SNPs.r2<-subset(SNPs.df, SNPs %in% SNP_list$SNP_list)
  SNPs.r2<-as.character(SNPs.r2$SNPs)
  
  LD_matrix_r2<-r2fast(gtdata(antler_qced.gen), idsubset = IDs, 
                       snpsubset = SNPs.r2)
  
  LD_matrix_r2.sym<-LD_matrix_r2
  LD_matrix_r2.sym[lower.tri(LD_matrix_r2.sym)] = t(LD_matrix_r2.sym)[lower.tri(LD_matrix_r2.sym)]
  
  Dist.frame<-as.data.frame(window$V1[c(1,20)])
  colnames(Dist.frame)<-"SNP.Name"
  Dist.frame<-join(Dist.frame, linkage_map[, c("SNP.Name", "CEL.LG", "Estimated.Mb.Position")])
  Dist<-abs(Dist.frame$Estimated.Mb.Position[1]-Dist.frame$Estimated.Mb.Position[2])*1e-06
  LD_values<-as.vector(LD_matrix_r2.sym)
  LD_avg<-mean(LD_values, na.rm = T)
  df<-data.frame(window=w, mean_LD=LD_avg, LG=Dist.frame$CEL.LG[1], Est.Dist.Mb=Dist)
  avg_LD_windows<-rbind(avg_LD_windows, df)
  
}

SNPs
mean(avg_LD_windows$Est.Dist.Mb)
mean(avg_LD_windows$mean_LD)
median(avg_LD_windows$mean_LD)

quantile(avg_LD_windows$mean_LD, probs = seq(0,1, 0.25))

write.table(avg_LD_windows, file = "average_LD_all_windows.txt", sep = "\t", col.names = T, row.names = F)
avg_LD_windows<-read.table("average_LD_all_windows.txt", header=T)

hist(avg_LD_windows$mean_LD, breaks = 50)
median(avg_LD_windows$mean_LD)
mean(avg_LD_windows$mean_LD)

highlight <- function(x, value1, col.value1, value2, col.value2, col=NA, ...){
  hst <- hist(x, ...)
  idx1 <- findInterval(value1, hst$breaks)
  idx2 <- findInterval(value2, hst$breaks)
  cols <- rep(col, length(hst$counts))
  cols[idx1] <- col.value1
  cols[idx2] <- col.value2
  hist(x, col=cols, ...)
}

png(filename = "average_LD_all_windows_hist.png", width=800, height=683, type = "windows")
par(mgp = c(2.5, 0.7, 0))
highlight(avg_LD_windows$mean_LD, 0.1157019, "dimgray", 0.1601955, "dimgray", breaks=50, main=NULL, 
          xlab="Pairwise LD r2", cex.lab=2, cex.axis=1.7)
dev.off()

#LD for windows in reg h2 analysis that contained GWAS hit for CoronetCirc

#first need to get windows that contain gwas hit

load("SNP_windows_list.RData")
w<-windows_list[[1]]
sig.SNP.windows<-character()

for (w in windows_list){
  assign("window", read.table(paste0("C:/Users/s1767711/Documents/Red_deer/Antler_genetic_architechture/Regional_h2/SNP_windows/", w,".txt")))
  match<-grep("^cela1_red_14_12034664$", window$V1)
  if (length(match)>0){
    sig.SNP.windows<-append(sig.SNP.windows, w)
  }
}


window1<-read.table(paste0("C:/Users/s1767711/Documents/Red_deer/Antler_genetic_architechture/Regional_h2/SNP_windows/", sig.SNP.windows[1], ".txt"))
SNPs_window1<-window1$V1
SNPs_window1<-as.character(SNPs_window1)

LD_matrix_r2<-r2fast(gtdata(antler_qced.gen), idsubset = IDs, 
                     snpsubset = SNPs_window1)



#make r2 matrix symmetric
LD_matrix_r2.sym<-LD_matrix_r2
LD_matrix_r2.sym[lower.tri(LD_matrix_r2.sym)] = t(LD_matrix_r2.sym)[lower.tri(LD_matrix_r2.sym)]

library(LDheatmap)
library(colorspace)
library(grid)
palette<-sequential_hcl(20)

png(filename = "LD_heatmap_window_link21_pos100.png", width=800, height=683, type = "windows")
window1_heatmap<-LDheatmap(LD_matrix_r2.sym, SNP.name = "cela1_red_14_12034664", color = palette, add.map = F, 
                           title = NULL)
grid.edit(gPath("ldheatmap", "SNPnames"), gp=gpar(cex=1, col="black"))
grid.edit(gPath("ldheatmap", "Key", "title"), gp=gpar(cex=1.3))
grid.edit(gPath("ldheatmap", "Key", "labels"), gp=gpar(cex=1.2))
LDheatmap.marks(window1_heatmap, i= c(rep(10, times=19)), j=c(1:9, 11:20), pch = 0, col="red", cex=2.8)
dev.off()

sig_SNP_LD<-c(LD_matrix_r2[1:9, 10], LD_matrix_r2[10, 11:20])
class(sig_SNP_LD)
LD_values_window1<-as.vector(LD_matrix_r2.sym)
mean(sig_SNP_LD)
mean(LD_values_window1, na.rm = T)


window2<-read.table(paste0("C:/Users/s1767711/Documents/Red_deer/Antler_genetic_architechture/Regional_h2/SNP_windows/",sig.SNP.windows[2],".txt"))
SNPs_window2<-window2$V1
SNPs_window2<-as.character(SNPs_window2)

LD_matrix_r2<-r2fast(gtdata(antler_qced.gen), idsubset = IDs, 
                     snpsubset = SNPs_window2)



#make r2 matrix symmetric
LD_matrix_r2.sym<-LD_matrix_r2
LD_matrix_r2.sym[lower.tri(LD_matrix_r2.sym)] = t(LD_matrix_r2.sym)[lower.tri(LD_matrix_r2.sym)]

library(LDheatmap)
library(colorspace)
library(grid)
palette<-sequential_hcl(20)

png(filename = "LD_heatmap_window_link21_pos90.png", width=800, height=683, type = "windows")
window2_heatmap<-LDheatmap(LD_matrix_r2.sym, SNP.name = "cela1_red_14_12034664", color = palette, add.map = F, 
                           title = NULL)
grid.edit(gPath("ldheatmap", "SNPnames"), gp=gpar(cex=1, col="black"))
grid.edit(gPath("ldheatmap", "Key", "title"), gp=gpar(cex=1.3))
grid.edit(gPath("ldheatmap", "Key", "labels"), gp=gpar(cex=1.2))
LDheatmap.marks(window1_heatmap, i= c(rep(20, times=19)), j=c(1:19), pch = 0, col="red", cex=2.8)
dev.off()

sig_SNP_LD<-c(LD_matrix_r2.sym[1:19, 20])
class(sig_SNP_LD)
LD_values_window2<-as.vector(LD_matrix_r2.sym)
mean(sig_SNP_LD)
mean(LD_values_window2, na.rm = T)
