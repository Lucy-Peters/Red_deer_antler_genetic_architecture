
antlers <- read.table("antler_data_comp.txt", header = T, stringsAsFactors = F)

antler_traits<-list()
k=1

for ( i in c(7:15, 21)) {
  unit<- colnames(antlers)[i]
  antler_traits[[k]]<-unit
  k=k+1
}

#antler_traits_sub<-antler_traits[c(5:6,9)]

df_list<-list()
l=1
for (trait in antler_traits){
  df<-data.frame()
  df_list[[l]]<-df
  l=l+1
}

#df<-assign(paste0(trait, "_all_windows_results"), df)

n=1

for (trait in antler_traits){
  results<-list.files(pattern=paste0("^",trait, "_results.*RDS"))
  results<-gsub(".RDS", "", results)
  results_list<-as.list(results)
  for (r in results_list){
    assign("results_df", readRDS(paste0(r,".RDS")))
    df_full<-rbind(df_list[[n]], results_df)
    df_list[[n]]<-df_full
  }
 # df_full<-assign(paste0(trait, "_all_windows_results"), df_full)
  n=n+1
} 

save(df_list, file="antler_traits_all_reg_h2_results.RData")

library(GenABEL)
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(forcats)

#load data frame list (run in eddie)
#load("antler_traits_all_reg_h2_results7.RData")
load("antler_traits_all_PCs_reg_h2_results5.RData")

#Bonferroni p-value correction, no of tests
df1<-df_list[[1]]
no_tests<-length(unique(df1$window_name))/2 #no of windows/2
p_value<-0.05/no_tests #new p-value cut-off at alpha=0.05
3608/2
-log10(p_value)


#use loop to filter out significant results
#sig_windows_all<-data.frame()
sig_windows_all_PCs<-data.frame()


for (df in df_list){
  exclude_windows<-subset(df$window_name, df$loglik==0)
  for (win in exclude_windows){
    pa<-paste0("^", win, "$")
    match<-grep(pa, df$window_name)
    if (length(match)>0){
      df_sub<-subset(df, !window_name==win)
    }
    df<-df_sub
  }  
  exclude_windows2<-subset(df$window_name, df$constraint=="Singular")
  for (win2 in exclude_windows2){
    pa2<-paste0("^", win2, "$")
    match2<-grep(pa2, df$window_name)
    if (length(match2)>0){
      df_sub2<-subset(df, !window_name==win2)
    }
    df<-df_sub2
  }
  # exclude_windows3<-subset(df$window_name, df$SE=="NaN")
  # for (win3 in exclude_windows3){
  #   pa3<-paste0("^", win3, "$")
  #   match3<-grep(pa3, df$window_name)
  #   if (length(match3)>0){
  #     df_sub3<-subset(df, !window_name==win3)
  #   }
  #   df<-df_sub3
  # }
  df$chisq<-qchisq(df$p_value, 1, lower.tail = F)
  df$lambda<-estlambda(df$chisq, method = "regression")$estimate      #median(df$chisq)/qchisq(0.5,1)
  if (df$lambda > 1){
    df$chisq_corr<-df$chisq/df$lambda
    df$p_corrected<-pchisq(df$chisq_corr, 1 , lower.tail = F)
  } else {
    df$chisq_corr<-df$chisq
    df$p_corrected<-df$p_value
  }
  sig_df<-subset(df[,c(12:13, 17)], p_corrected< p_value)
  sig_df<-subset(sig_df, !duplicated(sig_df))
  #sig_windows_all<-rbind(sig_windows_all, sig_df)
  sig_windows_all_PCs<-rbind(sig_windows_all_PCs, sig_df)
  # if (nrow(sig_df)>0){
  #   assign(paste0(sig_df$Trait, "_reg_h2_sig_df"), df)
  # }
  assign(paste0(df$Trait, "_reg_h2_results_df"), df)
  write.table(df, file= paste0(df$Trait[1], "_reg_h2_results.txt"), sep="\t", row.names = F, quote = F)
}



#write.table(sig_windows_all, file ="reg_h2_significant_windows.txt",  sep="\t", row.names = F, quote = F)
write.table(sig_windows_all_PCs, file ="reg_h2_significant_windows_PCs.txt",  sep="\t", row.names = F, quote = F)


#antlers <- read.table("antler_data_comp.txt", header = T, stringsAsFactors = F)
antlers <- read.table("antler_resid_PCs.txt", header = T, stringsAsFactors = F)


antler_traits<-list()
k=1

 # for ( i in c(7:16)) {
 #   unit<- colnames(antlers)[i]
 #   antler_traits[[k]]<-unit
 #   k=k+1
 # }

for ( i in c(5:15)) {
  unit<- colnames(antlers)[i]
  antler_traits[[k]]<-unit
  k=k+1
}


reg_h2_results_list<-list()
l=1
for (trait in antler_traits){
  df_temp<-data.frame()
  reg_h2_results_list[[l]]<-df_temp
  l=l+1
}


n=1

for (trait in antler_traits){
  results<-list.files(pattern=paste0("^",trait, "_reg_h2_results.txt"))
  results_list<-as.list(results)
  for (r in results_list){
    assign("results_df", read.table(paste0(r), header = T))
    reg_h2_results_list[[n]]<-rbind(reg_h2_results_list[[n]], results_df)
  }
  n=n+1
} 



#save(reg_h2_results_list, file = "antler_reg_h2_qced_results_list.RData")
save(reg_h2_results_list, file = "antler_reg_h2_qced_results_list_PCs.RData")

#load("antler_reg_h2_qced_results_list.RData")
load("antler_reg_h2_qced_results_list_PCs.RData")
df<-reg_h2_results_list[[1]]

#put all results tables for each antler measure/PC into one data frame
reg_h2_results_all<-rbindlist(reg_h2_results_list, use.name=T)
#write.table(reg_h2_results_all, file="antler_reg_h2_qced_results_all.txt", sep="\t", row.names = F, quote = F)
write.table(reg_h2_results_all, file="antler_reg_h2_qced_results_PCs_all.txt", sep="\t", row.names = F, quote = F)


reg_h2_results_part<-list()
l=1
for (trait in antler_traits){
  df_temp<-data.frame()
  reg_h2_results_part[[l]]<-df_temp
  l=l+1
}

y=1

for (d in reg_h2_results_list){
  d_part<-subset(d[, c(12:13, 17)], !duplicated(d[, c(12:13, 17)]))
  d_part$window_name<-as.character(d_part$window_name)
  d_split<-strsplit(d_part$window_name, split = "_")
  link_group<-list()
  len<-length(d_split)
  x=1
  for (i in c(1:len)){
    unit<-d_split[[i]][2]
    link_group[[x]]<-unit
    x=x+1
  }
  positions<-list()
  z=1
  for (h in c(1:len)){
    element<-d_split[[h]][3]
    positions[[z]]<-element
    z=z+1
  }
  link_group_vector<-unlist(link_group)
  link_group<-as.data.frame(link_group_vector)
  colnames(link_group)[1]<-"link_group"
  positions_vector<-unlist(positions)
  positions<-as.data.frame(positions_vector)
  colnames(positions)[1]<-"positions"
  d_full_part<-cbind(d_part, link_group)
  d_full<-cbind(d_full_part, positions)
  d_full$link_group<-gsub("link", "", d_full$link_group)
  d_full$positions<-gsub("pos", "", d_full$positions)
  d_comp<-rbind(reg_h2_results_part[[y]], d_full)
  reg_h2_results_part[[y]]<-d_comp
  y=y+1
}


#save(reg_h2_results_part, file = "antler_reg_h2_results_plot_df_list.RData")
save(reg_h2_results_part, file = "antler_reg_h2_results_plot_df_list_PCs.RData")

## plotting results ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#load("antler_reg_h2_results_plot_df_list.RData")
load("antler_reg_h2_results_plot_df_list_PCs.RData")
df_plot<-reg_h2_results_part[[1]]
source("manhattan_plot_code.R")

p_value <-2.776235e-05 #- see beginning of code before filtering for significant windows


#plots using manhattan_plot code (not as nice)

for (df_plot in reg_h2_results_part){
  df_plot$positions<-as.numeric(df_plot$positions)
  df_plot$link_group<-as.numeric(df_plot$link_group)
  df_plot<-arrange(df_plot, link_group, positions)
  cum<-nrow(df_plot)
  df_plot$pos_cuml<-c(1:cum)
  
  png(filename = paste0(df_plot$Trait, "_reg_h2_manhattan_plot_lines.png"), width=800, height=600, type = "windows")
  par(mgp = c(2.5, 0.7, 0))
  manhattan_p(df_plot, chr="link_group", bp="pos_cuml", p="p_corrected", snp="window_name", 
              col=c("blue", "red"), chrlabs=c(1:33), suggestiveline = FALSE, 
              genomewideline = -log10(p_value), ylim=range(0:7), 
              cex.axis= 1.5, cex.lab=2, xlab="CEL Linkage Group", las=2)

  for (i in 1:cum){
    x<-df_plot$link_group[i]%% 2==0
    if (x==TRUE){
      colour="red"
    } else {
      colour="blue"
    }
    segments(i, 0, i, -log10(df_plot$p_corrected[i]), col=colour)
  }
  dev.off()
}


##### make facet manhattan plots for all antler measures and all PCs  #######

reg_h2_results_all_plot<-rbindlist(reg_h2_results_part, use.name=T)

reg_h2_results_all_plot$positions<-as.numeric(reg_h2_results_all_plot$positions)
reg_h2_results_all_plot$link_group<-as.numeric(reg_h2_results_all_plot$link_group)
reg_h2_results_all_plot<-arrange(reg_h2_results_all_plot, link_group, positions)
cum<-nrow(reg_h2_results_all_plot)
reg_h2_results_all_plot$pos_cuml<-c(1:cum)

reg_h2_results_all_plot<-reg_h2_results_all_plot[, c("Trait", "window_name", "p_corrected", "link_group", "pos_cuml", "positions")]


#make custom facet labels
measure.labs<-c("a. Antler Length", "b. Coronet Circumference", "c. Lower Beam Circ.", 
                "d. Upper Beam Circ.", "e. Coronet-Brow Junc.", "f. Coronet-Tray Junc.", "g. Brow Length", "h. Tray Length", 
                "i. Antler Weight", "j. Form")
names(measure.labs)<-c("Length", "CoronetCirc", "LowerBeam", "UpperBeam", "CoronetBrowJunc",
                       "CoronetTrayJunc", "BrowLength", "TrayLength", "AntlerWt", "Form" )

#force order of measures to be the same as labels above
measure_levels<-c("Length", "CoronetCirc", "LowerBeam", "UpperBeam", "CoronetBrowJunc",
                  "CoronetTrayJunc", "BrowLength", "TrayLength", "AntlerWt", "Form" )
reg_h2_results_all_plot$measure.ordered<-factor(reg_h2_results_all_plot$Trait, levels = measure_levels)

#make axis df 
axisdf_regh2 <- ddply(reg_h2_results_all_plot,.(link_group), summarize, center=(max(pos_cuml) + min(pos_cuml))/2)

label_seq<-axisdf_regh2$link_group[seq(1,length(axisdf_regh2$link_group),by=2)]
breaks_seq<-axisdf_regh2$center[seq_along(axisdf_regh2$center)%% 2 > 0]



library(forcats)
library(ggplot2)

p_regh2_wrap<-ggplot(reg_h2_results_all_plot, aes(pos_cuml, -log10(p_corrected))) +
  geom_point(aes(color=as.factor(link_group)), size=3, alpha=0.5)+
  geom_linerange(aes(x=pos_cuml, ymax=-log10(p_corrected), ymin=min(-log10(p_corrected)), color=as.factor(link_group)))+
  scale_colour_manual(values = rep(c("steelblue3", "red3"), 33))+
  scale_x_continuous(label= label_seq, breaks = breaks_seq)+
  scale_y_continuous(expand = c(0,0))+
  ylim(0,7)+ 
  geom_hline(yintercept = -log10(p_value), linetype="dashed", color="black", size=1)+
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


ggsave("manhattan_plot_regh2_wrap_ms.png", p_regh2_wrap, width = 30, height = 40, units = "cm")


#PCs wrapped manhattan plot
reg_h2_results_all_plot<-rbindlist(reg_h2_results_part, use.name=T)

reg_h2_results_all_plot$positions<-as.numeric(reg_h2_results_all_plot$positions)
reg_h2_results_all_plot$link_group<-as.numeric(reg_h2_results_all_plot$link_group)
reg_h2_results_all_plot<-arrange(reg_h2_results_all_plot, link_group, positions)
cum<-nrow(reg_h2_results_all_plot)
reg_h2_results_all_plot$pos_cuml<-c(1:cum)

reg_h2_results_all_plot<-reg_h2_results_all_plot[, c("Trait", "window_name", "p_corrected", "link_group", "pos_cuml", "positions")]


#make custom facet labels
measure.labs<-c("a. PC1", "b. PC2", "c. PC3", "d. PC4", "e. PC5", "f. PC6", "g. PC7", "h. PC8", 
                "i. PC9", "j. PC10", "k. PC11")
names(measure.labs)<-c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", 
                       "PC9", "PC10", "PC11")

#force order of measures to be the same as labels above
measure_levels<-c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", 
                  "PC9", "PC10", "PC11")
reg_h2_results_all_plot$measure.ordered<-factor(reg_h2_results_all_plot$Trait, levels = measure_levels)


#make axis df 
axisdf_regh2 <- ddply(reg_h2_results_all_plot,.(link_group), summarize, center=(max(pos_cuml) + min(pos_cuml))/2)

label_seq<-axisdf_regh2$link_group[seq(1,length(axisdf_regh2$link_group),by=2)]
breaks_seq<-axisdf_regh2$center[seq_along(axisdf_regh2$center)%% 2 > 0]


p_regh2_PCs_wrap<-ggplot(reg_h2_results_all_plot, aes(pos_cuml, -log10(p_corrected))) +
  geom_point(aes(color=as.factor(link_group)), size=3, alpha=0.5)+
  geom_linerange(aes(x=pos_cuml, ymax=-log10(p_corrected), ymin=min(-log10(p_corrected)), color=as.factor(link_group)))+
  scale_colour_manual(values = rep(c("steelblue3", "red3"), 33))+
  scale_x_continuous(label= label_seq, breaks = breaks_seq)+
  scale_y_continuous(expand = c(0,0))+
  ylim(0,7)+ 
  geom_hline(yintercept = -log10(p_value), linetype="dashed", color="black", size=1)+
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



ggsave("manhattan_plot_regh2_PCs_wrap_ms.png", p_regh2_PCs_wrap, width = 30, height = 40, units = "cm")




#make plot for manuscript - PC9 (3 regh2 hits) use ggplot instead of 'manhattan plot' function
library(plyr)
library(dplyr)
library(ggplot2)

load("antler_reg_h2_results_plot_df_list_PCs.RData")
df_plot<-reg_h2_results_part[[9]]
head(df_plot)
p_value <- 2.773156e-05 #- see beginning of code before filtering for significant windows

df_plot$positions<-as.numeric(df_plot$positions)
df_plot$link_group<-as.numeric(df_plot$link_group)
df_plot<-arrange(df_plot, link_group, positions)
cum<-nrow(df_plot)
df_plot$pos_cuml<-c(1:cum)

df_plot<-df_plot[, c("Trait", "window_name", "p_corrected", "link_group", "pos_cuml", "positions")]

#make separate x axis for chrosmosomes and add titel variable to label plot with a facet title including an expression
axisdf_regh2 <- ddply(df_plot,.(link_group), summarize, center=(max(pos_cuml) + min(pos_cuml))/2)
title.label<-paste0("PC9")
df_plot$titel<-gl(1, nrow(df_plot), labels=title.label)

label_seq<-axisdf_regh2$link_group[seq(1,length(axisdf_regh2$link_group),by=2)]
breaks_seq<-axisdf_regh2$center[seq_along(axisdf_regh2$center)%% 2 > 0]

p_PC9_regh2<-ggplot(df_plot, aes(pos_cuml, -log10(p_corrected))) +
  geom_point(aes(color=as.factor(link_group)), size=3, alpha=0.5)+
  geom_linerange(aes(x=pos_cuml, ymax=-log10(p_corrected), ymin=min(-log10(p_corrected)), color=as.factor(link_group)))+
  scale_colour_manual(values = rep(c("steelblue3", "red3"), 33))+
  scale_x_continuous(label= label_seq, breaks = breaks_seq)+ 
  scale_y_continuous(expand = c(0,0))+
  ylim(0,7)+ 
  geom_hline(yintercept = -log10(p_value), linetype="dashed", color="black", size=1)+
  xlab("CEL Linkage Group")+ylab(expression(-log[10](p)))+
  theme(legend.position = "none", axis.text.y = element_text(size=25, colour = "black"), 
        axis.text.x = element_text(size=20, colour = "black"),
        axis.title.x = element_text(size=35, margin = margin(t = 30, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size=35, margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.line = element_line(colour = "black"),
        axis.ticks=element_line(size=1), strip.text = element_text(size=40),
        strip.background = element_rect(colour="black", fill="lightgray"))+
  facet_grid(. ~ titel, labeller = "label_parsed")

ggsave("PC9_manhattan_plot_regh2_ms.png", p_PC9_regh2, width = 50, height = 30, units = "cm")


#check for over-disperion of p-values
PC9_results<-read.table("PC9_reg_h2_results.txt", header=T)

library(chopsticks)
qq.chisq(PC9_results$chisq_corr, df=1, overdisp = T)
median(PC9_results$chisq_corr)/qchisq(0.5,1)

library(GenABEL)
estlambda(PC9_results$chisq_corr, method = "regression" , plot = T)






