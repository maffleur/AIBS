rm(list = ls())
setwd("~/Desktop/TIG_research_project/AIBS_studies/")

# Libraries
require(EWCE)
require(loomR)
require(ggplot2)
require(reshape2)
require(tidyverse)
require(GO.db)

# Useful functions
gseapath <- "/Users/Marzia/Desktop/Software/GSEA_AA"
source(paste(gseapath, "/", "enrich.R", sep=""))
source("/Users/Marzia/Desktop/TIG_research_project/Scripts&Functions/ewce.plot.mod.R")
source("/Users/Marzia/Desktop/TIG_research_project/Scripts&Functions/multiplot.R")

# Flags
use.lmeres <- F
run.gsea   <- T
big.terms  <- F # if GO terms or pathways with > 500 genes are to be used for GSEA
revigo     <- F
run.ewce   <- T
use.ctx    <- T
annotationLevel <- 1
#########################################
# Prepare expression data for EWCE
download.file("https://storage.googleapis.com/linnarsson-lab-loom/l5_all.agg.loom", 
              destfile = "l5_all.agg.loom")
lfile <- connect(filename = "l5_all.agg.loom", mode = "r+", skip.validate = T)

xmat.gene.names <- lfile$get.attribute.df(MARGIN = 1, attribute.names =
                                            c("Gene"), row.names="Accession")

xmat.cell.names <- lfile$get.attribute.df(MARGIN = 2,
                                          attribute.names=c("Clusters","Class","ClusterName","Description", "Region"),
                                          col.names=c("Clusters"))

if(use.ctx){ 
  toremove = c("Enteric nervous system", "Spinal cord", "Dorsal root ganglion", "Dorsal root ganglion,Sympathetic ganglion", "Sympathetic ganglion")
  ctx.id <- !(xmat.cell.names$Region %in% toremove)
  xmat <- t(lfile[['matrix']][ctx.id,])
  rownames(xmat) <- xmat.gene.names$Gene
  l1 <- xmat.cell.names[ctx.id, "Class"]
  l2 <- xmat.cell.names[ctx.id, "Description"]
  annotLevels = list(l1 = l1, l2 = l2)
  
} else{
  
  xmat <- t(lfile[["matrix"]][, ])
  rownames(xmat) <- xmat.gene.names$Gene
  l1 <- xmat.cell.names$Class
  l2 <- xmat.cell.names$Description
  annotLevels = list(l1 = l1, l2 = l2)
  
}

ctd <- generate.celltype.data(exp = xmat, annotLevels = annotLevels, groupName = "LoomFileMouse")
lfile$close_all()

# data("mouse_to_human_homologs")
# m2h        <- unique(mouse_to_human_homologs[,c("HGNC.symbol", "MGI.symbol")])
############################################

# read the results files in
if(use.lmeres){
  resfiles         <- list.files(pattern = "LMEresults.txt")
  restables        <- lapply(resfiles, read.delim)
  names(restables) <- sub(".*AIBS_", "", sub("\\_LMEresults.txt.*", "", resfiles)) # takes care of trait names including underscores
  
  # define the gene sets for GSEA
  bg       <- paste(unique(restables[[1]]$HGNC.symbol))
  
  # PLEASE NOTE: separate lists will have to be made (and tested) for genes up- and down-regulated!
  ntopgenes <- 500 # Number of top up- and down-regulated genes to be included int the lists for GSEA
  glist.up    <- lapply(restables, FUN = function(x){
    tmp <- x %>% filter(Pvalue.fdr < 0.01) %>% top_n(ntopgenes, FixedEffect.Tstat)
    as.character(unique(tmp$HGNC.symbol))
  })
  glist.down    <- lapply(restables, FUN = function(x){
    tmp <- x %>% filter(Pvalue.fdr < 0.01) %>% top_n(-ntopgenes, FixedEffect.Tstat)
    as.character(unique(tmp$HGNC.symbol))
  })
  glist <- c(glist.up, glist.down)
  names(glist) = paste0(names(glist), c(rep(".up", length(glist.up)), rep(".down", length(glist.down))))
  rm(glist.up, glist.down)
  
  # run gsea
  if (run.gsea){
    
    kegg     <- read.csv(paste(gseapath, "/", "c2.cp.kegg.v6.1.symbols.gmt.mod", sep=""), header = F, stringsAsFactors = F)
    reactome <- read.csv(paste(gseapath, "/", "c2.cp.reactome.v6.1.symbols.gmt.mod", sep=""), header = F, stringsAsFactors = F)
    go       <- read.csv(paste(gseapath, "/", "c5.all.v6.1.symbols.gmt.mod", sep=""), header = F, stringsAsFactors = F)
    cgp      <- read.csv(paste(gseapath, "/", "c2.cgp.v6.1.symbols.gmt.mod", sep=""), header = F, stringsAsFactors = F)
    
    if (!big.terms){
      go$V4       <- apply(go, 1, FUN = function(x){length(strsplit(x[3], " ")[[1]])})
      go          <- go %>% filter(V4 > 15 & V4 < 500)
      
      kegg$V4     <- apply(kegg, 1, FUN = function(x){length(strsplit(x[3], " ")[[1]])})
      kegg        <- kegg %>% filter(V4 > 15 & V4 < 500)
      
      reactome$V4 <- apply(reactome, 1, FUN = function(x){length(strsplit(x[3], " ")[[1]])})
      reactome    <- reactome %>% filter(V4 > 15 & V4 < 500)
      
      cgp$V4      <- apply(cgp, 1, FUN = function(x){length(strsplit(x[3], " ")[[1]])})
      cgp         <- cgp %>% filter(V4 > 15 & V4 < 500)
    }
    
    message("kegg")
    glist.kegg     <- lapply(glist, function(gl){ gene_enrichment(bglist=bg, explist=gl, terms=kegg)})
    message("reactome")
    glist.reactome <- lapply(glist, function(gl){ gene_enrichment(bglist=bg, explist=gl, terms=reactome)})
    message("GO")
    glist.go       <- lapply(glist, function(gl){ gene_enrichment(bglist=bg, explist=gl, terms=go)})
    message("CPG")
    glist.cgp      <- lapply(glist, function(gl){ gene_enrichment(bglist=bg, explist=gl, terms=cgp)})
    
    # correct the p-values trait-wise, select the significant terms, and bind everything into a results dataframe
    gsea.res <- data.frame("hit"=numeric(),
                           "expected"=numeric(),
                           "OR"=numeric(),
                           "P"=numeric(),
                           "Trait"=character(),
                           "Pfdr"=numeric(),
                           "Term"=character())
    for (trait in names(glist)){
      tmp      <- data.frame(rbind(glist.kegg[[trait]], glist.go[[trait]], glist.reactome[[trait]], glist.cgp[[trait]]), "Trait"=trait)
      tmp$Pfdr <- p.adjust(tmp$P, "fdr")
      tmp$Term <- rownames(tmp)
      tmp <- filter(tmp, Pfdr < 0.05)
      gsea.res <- rbind(gsea.res, tmp)
    }
    rm(tmp)
    
    ###
    if (revigo){
      # OPTIONAL: prepare files for high-level GO visualisation with REVIGO (http://revigo.irb.hr/)
      goterms <- Term(GOTERM)
      goterms <- data.frame("Term" = as.character(goterms), "GO.ID"=as.character(names(goterms)))
      gsea.res$Term <- gsub(pattern = "GO_", replacement = "", x = gsea.res$Term)
      gsea.res$Term <- gsub(pattern = "_", replacement = " ", x = gsea.res$Term)
      gsea.res$Term <- tolower(gsea.res$Term)
      gsea.res <- merge(gsea.res, goterms, by="Term", all.x = T)

      # generate files
      for (trait in c("ODI", "Vin", "Viso", "FA", "MD")){
        write.table(gsea.res[(grepl(gsea.res$Trait, pattern = trait) & !is.na(gsea.res$GO.ID)), c("GO.ID", "Pfdr")], 
                    file = paste0("./REVIGO/GSEA_", trait, "_signif.txt"), 
                    quote=F, 
                    row.names = F)
      }
    }
  }
  # 
  ##############################################################################################################
  ##############################################################################################################
  
  # EWCE
  if (run.ewce){
    
    load(ctd)
    
    ##
    # Use the original LME results to select the top 500 up- and down-regulated genes for EWCE
    # might be useful since there are ~10k significant genes for the NODDI indices, 
    # in order to get more meaningful enrichments
    tt_results = lapply(restables, FUN = function(gl){
      ewce_expression_data(sct_data = ctd,
                           tt = gl,
                           annotLevel = annotationLevel,
                           sortBy = "FixedEffect.Tstat", 
                           ttSpecies = "human",
                           sctSpecies = "mouse",
                           thresh = 500,
                           reps = 10000)
    })
    tt_results_long <- data.frame("CellType" = factor(), 
                                  "annotLevel" = numeric(), 
                                  "p" = numeric(), 
                                  "fold_change" = numeric(),
                                  "sd_from_mean" = numeric(),
                                  "Direction" = factor(),
                                  "list" = character())
    for (trait in names(tt_results)){
      tt_results_long <- rbind(tt_results_long, 
                               data.frame(tt_results[[trait]][["joint_results"]], "list" = trait))
    }
    print(ewce.plot.mod(tt_results_long, mtc_method = 'fdr'))
  }
} else{
  
  res <- read.delim("AIBS_AllTraits_SpearmansResults.txt") 
  restables <- split(x = res, f = res$Trait)
 
  # define the gene sets for GSEA
  bg       <- paste(unique(restables[[1]]$HGNC.symbol))
  
  # PLEASE NOTE: separate lists will have to be made (and tested) for genes up- and down-regulated!
  ntopgenes <- 500 # Number of top up- and down-regulated genes to be included int the lists for GSEA
  glist.up    <- lapply(restables, FUN = function(x){
    tmp <- x %>% filter(Metap_min_fdr < 0.01) %>% top_n(ntopgenes, Metaz)
    as.character(unique(tmp$HGNC.symbol))
  })
  glist.down    <- lapply(restables, FUN = function(x){
    tmp <- x %>% filter(Metap_min_fdr < 0.01) %>% top_n(-ntopgenes, Metaz)
    as.character(unique(tmp$HGNC.symbol))
  })
  glist <- c(glist.up, glist.down)
  names(glist) = paste0(names(glist), c(rep(".up", length(glist.up)), rep(".down", length(glist.down))))
  rm(glist.up, glist.down)
  
  # run gsea
  if (run.gsea){
    
    kegg     <- read.csv(paste(gseapath, "/", "c2.cp.kegg.v6.1.symbols.gmt.mod", sep=""), header = F, stringsAsFactors = F)
    reactome <- read.csv(paste(gseapath, "/", "c2.cp.reactome.v6.1.symbols.gmt.mod", sep=""), header = F, stringsAsFactors = F)
    go       <- read.csv(paste(gseapath, "/", "c5.all.v6.1.symbols.gmt.mod", sep=""), header = F, stringsAsFactors = F)
    cgp      <- read.csv(paste(gseapath, "/", "c2.cgp.v6.1.symbols.gmt.mod", sep=""), header = F, stringsAsFactors = F)
    
    if (!big.terms){
      go$V4       <- apply(go, 1, FUN = function(x){length(strsplit(x[3], " ")[[1]])})
      go          <- go %>% filter(V4 > 15 & V4 < 500)
      
      kegg$V4     <- apply(kegg, 1, FUN = function(x){length(strsplit(x[3], " ")[[1]])})
      kegg        <- kegg %>% filter(V4 > 15 & V4 < 500)
      
      reactome$V4 <- apply(reactome, 1, FUN = function(x){length(strsplit(x[3], " ")[[1]])})
      reactome    <- reactome %>% filter(V4 > 15 & V4 < 500)
      
      cgp$V4      <- apply(cgp, 1, FUN = function(x){length(strsplit(x[3], " ")[[1]])})
      cgp         <- cgp %>% filter(V4 > 15 & V4 < 500)
    }
    
    message("kegg")
    glist.kegg     <- lapply(glist, function(gl){ gene_enrichment(bglist=bg, explist=gl, terms=kegg)})
    message("reactome")
    glist.reactome <- lapply(glist, function(gl){ gene_enrichment(bglist=bg, explist=gl, terms=reactome)})
    message("GO")
    glist.go       <- lapply(glist, function(gl){ gene_enrichment(bglist=bg, explist=gl, terms=go)})
    message("CPG")
    glist.cgp      <- lapply(glist, function(gl){ gene_enrichment(bglist=bg, explist=gl, terms=cgp)})
    
    # correct the p-values trait-wise, select the significant terms, and bind everything into a results dataframe
    gsea.res <- data.frame("hit"=numeric(),
                           "expected"=numeric(),
                           "OR"=numeric(),
                           "P"=numeric(),
                           "Trait"=character(),
                           "Pfdr"=numeric(),
                           "Term"=character())
    for (trait in names(glist)){
      tmp      <- data.frame(rbind(glist.kegg[[trait]], glist.go[[trait]], glist.reactome[[trait]], glist.cgp[[trait]]), "Trait"=trait)
      tmp$Pfdr <- p.adjust(tmp$P, "fdr")
      tmp$Term <- rownames(tmp)
      tmp <- filter(tmp, Pfdr < 0.05)
      gsea.res <- rbind(gsea.res, tmp)
    }
    rm(tmp)
    
    ###
    if (revigo){
      # OPTIONAL: prepare files for high-level GO visualisation with REVIGO (http://revigo.irb.hr/)
      goterms <- Term(GOTERM)
      goterms <- data.frame("Term" = as.character(goterms), "GO.ID"=as.character(names(goterms)))
      gsea.res$Term <- gsub(pattern = "GO_", replacement = "", x = gsea.res$Term)
      gsea.res$Term <- gsub(pattern = "_", replacement = " ", x = gsea.res$Term)
      gsea.res$Term <- tolower(gsea.res$Term)
      gsea.res <- merge(gsea.res, goterms, by="Term", all.x = T)
      
      # generate files
      for (trait in c("ODI", "Vin", "Viso", "FA", "MD")){
        write.table(gsea.res[(grepl(gsea.res$Trait, pattern = trait) & !is.na(gsea.res$GO.ID)), c("GO.ID", "Pfdr")], 
                    file = paste0("./REVIGO/GSEA_", trait, "_signif.txt"), 
                    quote=F, 
                    row.names = F)
      }
    }
  }
  
  if (run.ewce){
    load(ctd)
    tt_results = lapply(restables, FUN = function(gl){
      ewce_expression_data(sct_data = ctd,
                           tt = gl,
                           annotLevel = annotationLevel,
                           sortBy = "Metaz", 
                           ttSpecies = "human",
                           sctSpecies = "mouse",
                           thresh = 500,
                           reps = 10000)
    })
    tt_results_long <- data.frame("CellType" = factor(), 
                                    "annotLevel" = numeric(), 
                                    "p" = numeric(), 
                                    "fold_change" = numeric(),
                                    "sd_from_mean" = numeric(),
                                    "Direction" = factor(),
                                    "list" = factor())
    for (trait in names(tt_results)){
      tt_results_long <- rbind(tt_results_long, 
                                 data.frame(tt_results[[trait]][["joint_results"]], "list" = trait))
    }
    print(ewce.plot.mod(tt_results_long, mtc_method = 'fdr'))
  }

}
