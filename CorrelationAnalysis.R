rm(list=ls())
setwd("~/Desktop/TIG_research_project/AIBS_studies/")
source("probe2value_Marzia.R")
#load("TestExpressionData.RData")
# Flag for the non parametric approach
use.lme <- F
# Flag for the negative control via trait scrambling
random <- T
#
#load Reannotator results
annot <- read.delim("AllenInstitute_custom_Agilent_Array.txt", head=T, sep="\t", stringsAsFactors = F, strip.white = T)
# rownames(annot) <- annot[, 1] 
# for some reason it doesnt let me set the row names b/c there are duplicate probe IDs
#so remove the duplicated probe names
annot <- annot[!duplicated(annot$PROBE_ID), ]
rownames(annot) <- annot[, 1]

ngenes <- apply(annot, 1, function(x){
  length(unique(unlist( strsplit(paste(x[14]), split="[;]"))))
})
unique_probes <-  names(which(ngenes==1))

#remove intergenic probes (dist)
unique_probes <- unique_probes[grep("dist",paste(annot[unique_probes,6]),invert=T)]

unique_probes <- intersect(unique_probes, colnames(pa.dat))

##remove probes that are 'expressed' in less than 300 samples
unique_probes <- names(which(apply(pa.dat[,unique_probes],2,sum) >= 300))

unique_genes_mapped <- annot[unique_probes, c("PROBE_ID", "Gene_symbol")]

# restrict to cortical samples (slab CX)
ctx_id <- which(sample.info.new$slab_type == "CX")
expr_clean <- merged.dat[ctx_id, unique_probes]
expr_clean_withphenos <- cbind(sample.info.new[ctx_id, ], expr_clean)

rm(ctx_id, expr_clean, annot)

# Finally, fit the mixed effects models
#require(nlme)
if(use.lme){
  require(lme4)
  trait <- "ODI" # one of c("ODI","Vin","Viso","PCA")
  
  library(doParallel)
  cl<- makeCluster(detectCores())
  registerDoParallel(cl)
  
  myres <- foreach(x=1:length(unique_probes)) %dopar% {
    cat(unique_probes[x])
    # formula <- as.formula(paste(unique_probes[x], "~", trait)); # swap predictor and response!
    # mdl     <- nlme::lme(fixed = formula, random = ~1|sampleID, data = expr_clean_withphenos);
    f0 <- as.formula(paste(unique_probes[x], "~ 1 + (1 | sampleID)"))
    if (random){
      f1 <- as.formula(paste(unique_probes[x], "~ sample(", trait, ") + (1 | sampleID)"))
    } else {
      f1 <- as.formula(paste(unique_probes[x], "~", trait, "+ (1 | sampleID)"))
    }
    mdl0 <- lme4::lmer(formula = f0, data = expr_clean_withphenos, REML = F)
    mdl1 <- lme4::lmer(formula = f1, data = expr_clean_withphenos, REML = F)
    t <- summary(mdl1)$coefficients[2, 3]
    p <- anova(mdl0, mdl1)[2, 8]
    return(c(unique_probes[x], as.numeric(t), as.numeric(p)))
    # s       <- summary(mdl); 
    # return(c(unique_probes[x], as.numeric(s$tTable[2,5])))
  }; 
  
  myres <- as.data.frame(do.call(rbind, myres))
  colnames(myres) <- c("Probe", "FixedEffect.Tstat", "FixedEffect.Pvalue")
  myres$FixedEffect.Tstat  <- as.numeric(as.character(myres$FixedEffect.Tstat))
  myres$FixedEffect.Pvalue <- as.numeric(as.character(myres$FixedEffect.Pvalue))
  myres$Pvalue.fdr <- p.adjust(myres$FixedEffect.Pvalue, method = "fdr")
  myres$HGNC.symbol <- unique_genes_mapped
  
  write.table(myres, file = paste0("AIBS_", trait, "_LMEresults.txt"), quote = F, sep = "\t", row.names = F)
} else {
##########################################################################
# Non-parametric approach - with Spearman's rho

  library(dplyr)
  library(tidyr)
  
  traits <- c("NODDI_ODI","NODDI_Vin","NODDI_Viso","PCA", "Bingham_DAb", "Bingham_ODIs", "Bingham_ODIp", "Bingham_ODItot", "Bingham_Vin", "Bingham_Viso", "FA", "MD")
  subjects <- c("H0351.1009","H0351.1012","H0351.1015","H0351.1016","H0351.2001","H0351.2002")
  res <- data.frame("Trait"=factor(), "sampleID"=factor(), "Probe"=character(), "Rho"=numeric(), "P.value"=numeric())
  
  for (trait in traits){
    message(trait)
    for (subject in subjects){
      message(subject)
      tmp  <- subset(expr_clean_withphenos, subset = sampleID == subject)
      ps   <- rep(0, length(unique_probes))
      names(ps) <- unique_probes
      
      if (random){
        message("Randomising the imaging outcome...")
        rhos <- cor(x = tmp[unique_probes], y = sample(tmp[trait]), method = 'spearman', use="pairwise.complete.obs")
        
        for (probe in unique_probes){
          ps[probe]   <- cor.test(x = unlist(tmp[probe]), y = sample(unlist(tmp[trait])), method = 'spearman', exact=F)$p.value 
        }
      } else{
        rhos <- cor(x = tmp[unique_probes], y = tmp[trait], method = 'spearman', use="pairwise.complete.obs")
        
        for (probe in unique_probes){
          ps[probe]   <- cor.test(x = unlist(tmp[probe]), y = unlist(tmp[trait]), method = 'spearman', exact=F)$p.value 
        }
      }
      tmpres <- data.frame(cbind("Trait"=trait, "sampleID"=subject, "Probe"=unique_probes,"Rho"=rhos, "P.value"=ps))
      colnames(tmpres)[4] <- "Rho"
      tmpres$Rho <- as.numeric(as.character(tmpres$Rho))
      tmpres$P.value <- as.numeric(as.character(tmpres$P.value))
      res <- rbind(res, tmpres) 
    }
  }
  rm(trait, subject, tmp, rhos, ps, probe, tmpres)
  
  # the following shamelessly copied from AA
  tmp.cor <- res$Rho
  tmp.p <- res$P.value
  tpos <- tmp.p
  kkk <- tmp.cor <0
  tpos[kkk] <- 1.0 - tpos[kkk]
  tneg <- tmp.p
  tneg[!kkk] <- 1.0 - tneg[!kkk]
  res <- cbind(res, tpos, tneg)
  rm(tmp.cor, tmp.p, tpos, kkk, tneg)
  
  # combine the 6 probe-wise p-values with Stouffe's method
  res$tpos_mod = res$tpos; res$tpos_mod[res$tpos_mod==1] <- 0.999
  res$tneg_mod = res$tneg; res$tneg_mod[res$tneg_mod==1] <- 0.999
  res <- res %>% 
    group_by(.dots=c("Trait", "Probe")) %>% 
    summarise("Metaz_pos" = metap::sumz(tpos_mod)$z,
              "Metap_pos" = metap::sumz(tpos_mod)$p,
              "Metaz_neg" = -Metaz_pos,
              "Metap_neg" = metap::sumz(tneg_mod)$p) # may be replaced by (1 - Metap_pos) for speed
  
  res$Metaz <- with(res, ifelse(Metap_pos > Metap_neg, -Metaz_neg, Metaz_pos)) 
  res$Metap_min <- apply(res[, c(4,6)], 1, min) 
  
  res <- res %>% group_by(Trait) %>% mutate("Metap_min_fdr" = p.adjust(Metap_min, method = 'fdr'))
  res <- merge(res, unique_genes_mapped, by.x = "Probe", by.y = "PROBE_ID")
  colnames(res)[10] <- "HGNC.symbol"

  write.table(res, file="AIBS_AllTraits_SpearmansResults.txt", quote = F, sep = "\t", row.names = F)
}









