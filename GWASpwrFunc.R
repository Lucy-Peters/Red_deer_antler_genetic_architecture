                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                            # GWAS power analysis function #
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. Permute the phenotypes to break any existing correlations between
# SNPs and phenotypes in the data.
# 
# 2. Pick a SNP at random.
# 
# 3. Assign values to the three genotypes at this SNP such that the
# additive variance explained by the SNP equals the total amount of
# variance you need to explain. This will depend on the allele
# frequencies at the SNP. You could assume an additive model with the
# heterozygote exactly intermediate between the two homozygotes.
# 
# 4. Remove the SNP from the dataset. This means that you can only now
# detect the effect you have simulated in the data in step 3 via a
# linked SNP. This step is necessary because it is unlikely that there
# are causal SNPs in the real data. It implicitly assumes that QTLs segregate
# at the same frequencies as the SNPs in your data.
# 
# 5. Carry out GWAS. Do you pick up a significant  effect in the region
# of the SNP?
#   
# 6. Repeat steps 1..5 x1000.
#
# 7. Do 1...6 for a range of SNP h2 and all antler traits

getGWASpwr<-function(trait, antler_data, geno_tab, antler.gkin.sym, antler_qced.gen,
                     allele_frq, var_comp, P.threshold, h2SNPrange, nsim){
  
  PwrDf_all<-data.frame()
  GWAS.sig_all<-data.frame()
  
  antler_sub<-antler_data %>%
    dplyr::select(Code, MeCaYear, BirthYear, MeCaAge, all_of(trait))
  
  #variance of simulated causal SNP
  Va<-subset(var_comp, measure==trait & variable.names=="genetic")$component
  Vt<-sum(subset(var_comp, measure==trait)$component)
  
  
  for (h2SNP in h2SNPrange) {
    
    Vq<-h2SNP*Vt
    VarExplained<-Vq/(Va+Vq)
    #Vq<-(VarExplained*Va)/(1-VarExplained)
    SNPsidx<-sample(nrow(allele_frq), nsim)
    
    hit.count<-0
    GWAS.sig.sub<-data.frame()
    
    i<-1
    for (idx in SNPsidx) {
      
      #1.& 2. Shuffle phenotypes to break link with the ID/genotype, 
      #but keep association with age and environmental variables & pick 'causal' SNP 
      snp<-as.character(allele_frq$SNP.name[idx])
      geno<-geno_tab%>%
        dplyr::select(Code, all_of(snp))
      allele.p<-allele_frq$major.allele[idx]
      allele.q<-allele_frq$minor.allele[idx]
      hom.major<-paste0(allele_frq$major.allele[idx], allele_frq$major.allele[idx])
      hom.minor<-paste0(allele_frq$minor.allele[idx], allele_frq$minor.allele[idx])
      het<-paste0(allele_frq$minor.allele[idx], allele_frq$major.allele[idx])
      antler_data_new<-antler_sub[, c(2:ncol(antler_sub))][sample(nrow(antler_sub)), ]
      antler_data_new<-cbind(antler_sub[1], antler_data_new)
      #3. Assign additive genetic value to SNP based on selected SNP h2 (and thus Vq)
      p<-allele_frq$major.allele.freq[idx]
      q<-allele_frq$minor.allele.freq[idx]
      a<-sqrt((Vq/(2*p*q)))     
      
      #join genotype info to antler data and exclude IDs with missing genotypes 
      antler_data_new<-join(antler_data_new, geno, by="Code")
      antler_data_new<-subset(antler_data_new, !antler_data_new[, which(names(antler_data_new)==snp)]=="00")
      
      #add SNP effect to trait of interest, depending on genotype (additive, het in middle)
      antler_data_new[, which(names(antler_data_new)==trait)]<-ifelse(antler_data_new[, which(names(antler_data_new)==snp)]==hom.major,
                                                                      antler_data_new[, which(names(antler_data_new)==trait)]+2*a,
                                                                      ifelse(antler_data_new[, which(names(antler_data_new)==snp)]==hom.minor,
                                                                             antler_data_new[, which(names(antler_data_new)==trait)],
                                                                             ifelse(antler_data_new[, which(names(antler_data_new)==snp)]==het,
                                                                                    antler_data_new[, which(names(antler_data_new)==trait)]+a,
                                                                                    ifelse(is.na(antler_data_new[, which(names(antler_data_new)==trait)]), NA, NA))))
      
      
      
      # The rGLS() function that RepeatABEL uses can't deal with missing data. For 
      # each model, a dataset that is clean of NA values for that trait can be created
      # (NB but MUST contain the field "id"):
      colnames(antler_data_new)[1]<-"id"
      
      keep_ids<-data.frame(keep_ids=idnames(antler_qced.gen))
      antler_data_new<-subset(antler_data_new, id %in% keep_ids$keep_ids)
      
      #remove any NA entries from data
      antler_data_clean<-na.omit(antler_data_new)
      
      #4.remove simulated causal SNP from data
      snpkeep<-antler_qced.gen@gtdata@snpnames
      match<-grep(snp, snpkeep)
      snpkeep<-snpkeep[-match]
      antler_qced.gen_sub<-antler_qced.gen[, snpkeep]
      
      #5. run  GWAS
      
      fmla<-as.formula(paste0(trait , "~ I(MeCaAge)+ I(MeCaAge^2)"))
      
      print(paste0("GWAS ", trait, ", SNP h2: ", h2SNP, ", run", i))
      
      prefit<-preFitModel(fmla, random = ~1|id + 1|MeCaYear + 1|BirthYear, 
                          genabel.data = antler_qced.gen_sub, phenotype.data = antler_data_clean, 
                          corStruc = list(id=list("GRM","Ind"), MeCaYear=list("Ind"), BirthYear=list("Ind")),
                          GRM = antler.gkin.sym)
      
      gwas<- rGLS(fmla,
                  genabel.data = antler_qced.gen_sub,
                  phenotype.data = antler_data_clean,
                  V=prefit$V)
      
      #get results from gwas (correct for lambda inflation)
      results.measure <- results(gwas)
      print(head(results.measure))
      
      if(is.na(lambda(gwas)$estimate)){
        lambda<-estlambda(results.measure$P1df, method = "median")
        lambda.measure<-lambda$estimate
      } else {
        lambda.measure<-lambda(gwas)$estimate 
      }
      
      results.measure$chi2.1df<-qchisq(results.measure$P1df, 1, lower.tail = F)
      results.measure$chi2.1df_corrected<-results.measure$chi2.1df/lambda.measure
      results.measure$Pc1df<-pchisq(results.measure$chi2.1df_corrected, 1, lower.tail = F)
      results.ordered<-results.measure[order(results.measure$Pc1df),]
      results.measure.sig <- subset(results.measure, Pc1df < P.threshold)
      
      if(nrow(results.measure.sig)!=0){
        hit.count<-hit.count+1
        
        df.temp.sig<-data.frame(measure=rep(trait, nrow(results.measure.sig)),
                                SNP_sig=rownames(results.measure.sig), Pvalue=results.measure.sig$Pc1df,
                                sig.SNP.sim=snp, SNPh2=h2SNP, AdditiveEffect=a, SNP.prop.Va=VarExplained,
                                stringsAsFactors = F)
        
      } else {
        
        df.temp.sig<-data.frame(measure=trait, SNP_sig="NULL", Pvalue=min(na.omit(results.measure$Pc1df)), sig.SNP.sim=snp, SNPh2=h2SNP,
                                AdditiveEffect=a, SNP.prop.Va=VarExplained,
                                stringsAsFactors = F)
      }
      
      GWAS.sig.sub<-rbind(GWAS.sig.sub, df.temp.sig)
      
      i<-i+1
      
    }
    
    PwrDf<-data.frame(measure=trait, GWASpwr=hit.count/length(SNPsidx), SNPh2=h2SNP,
                      SNP.prop.Va=VarExplained)
    
    PwrDf_all<-rbind(PwrDf_all, PwrDf)
    
    GWAS.sig_all<-rbind(GWAS.sig_all, GWAS.sig.sub)
    
  }
  return(list(PwrDf_all, GWAS.sig_all))
} 



