#------------------------------#
# Antler measure data clean-up #
#------------------------------#

library(plyr)
library(dplyr)
library(tidyr)
library(RODBC)

#create connection to deer database and import query tables directly
db<-"C:\\Users\\s1767711\\Documents\\Red_deer\\database\\System\\antler_ms_db\\RedDeer1.75.accdb"
con<-odbcConnectAccess2007(db)
queries<-sqlTables(con, tableType = "VIEW")
antler.measures<-sqlFetch(con, "AntlerMeasures_no_MuseumWeights_no_add_form")
antler.forms<-sqlFetch(con, "additional_AntlerForm_data")
antler_museumWts<-sqlFetch(con, "AntlerMuseumWeights")
odbcClose(con)

write.table(antler_museumWts, file = "Antler_museum_weights_raw.txt", sep = "\t", col.names = T,
            row.names = F, quote=F)
write.table(antler.forms, file="Antler_additional_form_data_raw.txt", sep = "\t", col.names = T,
            row.names = F, quote=F)
write.table(antler.measures, file="Antler_measures_data_raw.txt", sep = "\t", col.names = T,
            row.names = F, quote=F)

#collate data from query tables

#separate left and right antler point observations into two data frames and join by row later, so that they match format in antler measure table -
#one colum for Form1 (Tines) and one for Form2(Tops)

forms_right<-antler.forms[,c(1:2,4:5,8)]
colnames(forms_right)[c(2,3,4)]<-c("MeCaYear", "Tops", "Tines")
forms_right$Side<-"R"

forms_left<-antler.forms[,c(1:2,3,6,8)]
colnames(forms_left)[c(2,3,4)]<-c("MeCaYear", "Tops", "Tines")
forms_left$Side<-"L"
forms_joined<-rbind(forms_left,forms_right)
table(is.na(forms_joined$BirthYear))
table(is.na(antler.forms$BirthYear))
table(is.na(forms_joined$MeCaYear))
table(is.na(forms_joined$Tines))
table(is.na(forms_joined$Tops))

#combine cast and measure year into MeCaYear in antler measure table
#transform CastYear into year of antler growth by substracting 1 (year)
antler.measures$CastYear<-(antler.measures$CastYear)-1
antler.measures$MeCaYear<-ifelse(is.na(antler.measures$MeasureYear), antler.measures$CastYear, antler.measures$MeasureYear)
table(is.na(antler.measures$BirthYear))
table(is.na(antler.measures$MeCaYear))
table(is.na(antler.measures$Length))

#join additional form data to antler measure data, 
#full join because there are ID that are in the one but not the other table
antler.measures.next<-join(forms_joined, antler.measures[, c(-2,-3)], type="full")
table(is.na(antler.measures.next$Length))
antler.measures.next$Tines<-with(antler.measures.next, ifelse(is.na(Tines), Form1, Tines))
antler.measures.next$Tops<-with(antler.measures.next, ifelse(is.na(Tops), Form2, Tops))
antler.measures.next$Form1<-with(antler.measures.next, ifelse(!is.na(Form1), (Form1+Tines)/2, Tines))
antler.measures.next$Form2<-with(antler.measures.next, ifelse(!is.na(Form2), (Form2+Tops)/2, Tops))

table(is.na(antler.measures.next$BirthYear))
table(is.na(antler.measures.next$MeCaYear))
table(is.na(antler.measures.next$Form1))
table(is.na(antler.measures.next$Form2))

table(is.na(antler_museumWts$Weight))
table(is.na(antler.measures.next$RumWeight))
antler_museumWts.sub<-subset(antler_museumWts, !Code %in% antler.measures.next$Code & !CastYear %in% antler.measures.next$MeCaYear )
antler_museumWts.sub2<-subset(antler_museumWts, !Code %in% antler.measures.next$Code)

#join museum weight to antler measure data
#transform CastYear into year of antler growth by substracting 1 (year)
antler_museumWts$CastYear<-(antler_museumWts$CastYear)-1
colnames(antler_museumWts)[c(2,3)]<-c("MeCaYear","MuseumWeight")
antler.measures.full<-join(antler.measures.next[,c(-3,-4)], antler_museumWts[, c(-5)], type="full")
table(is.na(antler.measures.full$BirthYear))
table(is.na(antler.measures.full$MeCaYear))
birth_year.na<-subset(antler.measures.full, is.na(BirthYear))
antler.measures.full<-subset(antler.measures.full, !Code %in% birth_year.na$Code)

#merge Rum and museum weights and convert museum weights
antler.measures.full$AntlerWt <- with(antler.measures.full, ifelse(!is.na(RumWeight), RumWeight,
                                                                   round(16.88 + 1.0485*MuseumWeight)))

#filter out antler weight measures for broken antlers
pat<-paste0("broken")
match<-grep(pat, antler.measures.full$Comments, ignore.case = T)
antler.measures.full_sub<-antler.measures.full[match, ]

antler.measures.full$AntlerWt<-with(antler.measures.full, ifelse(Code%in%antler.measures.full_sub$Code, NA, AntlerWt))
table(is.na(antler.measures.full$AntlerWt))

#add Age column to data
antler.measures.full$Age<-antler.measures.full$MeCaYear-antler.measures.full$BirthYear
antler_data<-antler.measures.full[,c(-15,-16,-17)]
#filtering- exclude ind. younger than 3 years at measurement and entries with too few antler measures
antler_data_edit<-antler_data[rowSums(is.na(antler_data)) < 7, ]#only for imputation data set (PCA)
#antler_data_edit2<-subset(antler_data, !AntlerWt %in% NA | !CoronetCirc %in% NA | !Length %in% NA)
antler_data_edit<-subset(antler_data, Age >=3) 
antler_data_edit<-subset(antler_data_edit, Age >=3)

#add form variable
antler_data_edit$Form<-antler_data_edit$Form1+antler_data_edit$Form2
table(is.na(antler_data_edit$Length))

antler_data_avg<-ddply(antler_data_edit, c("Code", "MeCaYear"), 
                       function(x) cbind(x[1,c(1:3,16)], data.frame(t(colMeans(x[c(5:15, 17)], 
                                                                               na.rm=TRUE)))))

antler_data_avg[antler_data_avg == "NaN"]<-NA
table(is.na(antler_data_avg$Length))
table(is.na(antler_data_avg$LowerBeam))

antler_data_comp<-antler_data_avg
count<-unique(antler_data_comp$Code)
length(unique(antler_data_comp$Code))
#exclude erroneous Form enntries - if brow length is 0 the Form 1 is absent and tray length of 0 means Form1 can't be bigger than 1 (only brow present)
antler_data_comp$BrowLength[antler_data_comp$BrowLength==0 & antler_data_comp$Form1>0] <- NA 
antler_data_comp$TrayLength[antler_data_comp$TrayLength==0 & antler_data_comp$Form1>1] <- NA 
table(is.na(antler_data_comp$Form1))
table(is.na(antler_data_comp$Form2))
table(antler_data_comp$Form1)
table(antler_data_comp$Form2)
hist(antler_data_comp$Form1)
hist(antler_data_comp$Form2)
hist(antler_data_comp$Form)
hist(log(antler_data_comp$Form1))

#sort half measures of Form1 and Form2 into disrete number categories
antler_data_comp$Form1 <- ifelse(is.na(antler_data_comp$Form1), NA,
                                 ifelse(antler_data_comp$Form1 <1, 0,
                                        ifelse(antler_data_comp$Form1 <2, 1,
                                               ifelse(antler_data_comp$Form1 == 2, 2,
                                                      ifelse(antler_data_comp$Form1 < 3, 2.5, 3)))))

antler_data_comp$Form2<-ifelse(is.na(antler_data_comp$Form2), NA,
                               ifelse(antler_data_comp$Form2 <=1, 1,
                                      ifelse(antler_data_comp$Form2 <2, 1.5,
                                             ifelse(antler_data_comp$Form2 <3, 2.5,
                                                    ifelse(antler_data_comp$Form2 ==3,3,
                                                           ifelse(antler_data_comp$Form2 <4, 3.5,
                                                                  ifelse(antler_data_comp$Form2==4, 4,
                                                                         ifelse(antler_data_comp$Form2 ==4.5, 4,
                                                                                ifelse(antler_data_comp$Form2 < 5, 5, 5)))))))))

form_table<- as.data.frame(table(antler_data_comp$Form))
hist(antler_data_comp$Form)


antler_data_comp$Form<-ifelse(antler_data_comp$Form ==0, NA,
                              ifelse(antler_data_comp$Form <=1, 1,
                                     ifelse(antler_data_comp$Form <2, 1.5,
                                            ifelse(antler_data_comp$Form ==2, 2,     
                                                   ifelse(antler_data_comp$Form <2.5, 2,
                                                          ifelse(antler_data_comp$Form <3, 2.5,
                                                                 ifelse(antler_data_comp$Form ==3,3,
                                                                        ifelse(antler_data_comp$Form <3.5, 3,
                                                                               ifelse(antler_data_comp$Form <4, 3.5,
                                                                                      ifelse(antler_data_comp$Form==4, 4,
                                                                                             ifelse(antler_data_comp$Form <4.5, 4,
                                                                                                    ifelse(antler_data_comp$Form < 5, 4.5, 
                                                                                                           ifelse(antler_data_comp$Form ==5, 5,
                                                                                                                  ifelse(antler_data_comp$Form <5.5, 5,
                                                                                                                         ifelse(antler_data_comp$Form < 6, 5.5, 
                                                                                                                                ifelse(antler_data_comp$Form ==6, 6,
                                                                                                                                       ifelse(antler_data_comp$Form <6.5, 6,
                                                                                                                                              ifelse(antler_data_comp$Form < 7, 6.5, 
                                                                                                                                                     ifelse(antler_data_comp$Form ==7, 7,
                                                                                                                                                            ifelse(antler_data_comp$Form <7.5, 7,
                                                                                                                                                                   ifelse(antler_data_comp$Form < 8, 7.5, 
                                                                                                                                                                          ifelse(antler_data_comp$Form <=8, 8, 8))))))))))))))))))))))


table(antler_data_comp$Form)
hist(antler_data_comp$Form)

#test antler form data
library(lme4)
model1<-lmer(Form~Age+(1|MeCaYear)+(1|BirthYear)+(1|Code), data=antler_data_comp)
hist(residuals(model1))
model2<-lmer(Form1~Age++(1|MeCaYear)+(1|BirthYear)+(1|Code), data=antler_data_comp)
hist(residuals(model2))


#complete data set (with all Form data variables)
antler_data_comp<-antler_data_comp[rowSums(is.na(antler_data_comp)) < 12, ]
colnames(antler_data_comp)[4]<-"MeCaAge"
write.table(antler_data_comp, file = "antler_data_comp_v2.txt", sep="\t", 
            row.names = FALSE)

#data set for imputation (PCA)
antler_data_comp<-subset(antler_data_comp, !Form1 == 0) #need complete data set of antler measures

#if Form1 (tines) = 1, then either only brow lenght or only tray length can exist (part of tines)
#so if Form1 = 1 and brow length = NA then brow length should be set to 0 for imputation;
#if Form1 = 1 and tray length = NA then tray length should be set to 0 for imputation; this way
#brow and tray length won't be imputed if it is impossible for those measures to exist for that antler

antler_data_comp$BrowLength<-ifelse(antler_data_comp$Form1==1 & is.na(antler_data_comp$BrowLength)==TRUE, 0,
                                    antler_data_comp$BrowLength)

antler_data_comp$CoronetBrowJunc<-ifelse(antler_data_comp$BrowLength==0, 0, antler_data_comp$CoronetBrowJunc)

antler_data_comp$TrayLength<-ifelse(antler_data_comp$Form1==1 & is.na(antler_data_comp$TrayLength)==TRUE, 0,
                                    antler_data_comp$TrayLength)

antler_data_comp$CoronetTrayJunc<-ifelse(antler_data_comp$TrayLength==0, 0, antler_data_comp$CoronetTrayJunc)

write.table(antler_data_comp, file = "antler_data_compPCA.txt", sep="\t", row.names = FALSE)


#antler data summaries
antler_data_comp<-read.table("antler_data_comp.txt", header=T)
head(antler_data_comp)
antler_measures<-as.list(colnames(antler_data_comp)[7:16])

antler_summaries<-data.frame()

for (m in antler_measures) {
  df_sub<-antler_data_comp %>%
    select(c("Code", m)) %>%
    na.omit()
  df_temp<-data.frame(Trait=m, N_individuals = length(unique(df_sub[["Code"]])), 
                      N_observations=nrow(df_sub), Mean=round(mean(df_sub[[m]]), digits = 3),
                      SD=round(sd(df_sub[[m]]), digits = 3))
  antler_summaries<-rbind(antler_summaries, df_temp)
}

write.table(antler_summaries, file="antler_summaries_v2.txt", sep="\t", row.names=F)

