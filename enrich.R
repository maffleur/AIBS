
#simplistic gene set enrichment analyis
#using fisher test
#terms contains: name, link , genes
gene_enrichment <- function(bglist, explist, terms){

  el2 <- intersect(bglist, explist) # genes in the experiment AND background
  goer <- t(apply(terms, 1, function(x){
    li <- unlist(strsplit(paste(x[3]), split=" ")[[1]]) # genes in the term
    li2 <- intersect(bglist, li) # genes in the term AND background
    
    restbg <- setdiff(bglist, li2) #genes in the background BUT NOT in the term

    inlist    <- intersect(el2, li2) # genes in the experiment AND background AND term
    notinlist <- setdiff(li2, el2) # genes in the term but not in the experiment
    onlybg    <- intersect(el2, restbg) # genes in the experiment BUT NEITHER IN TERM NOR IN BACKGROUND
    notinbg   <- setdiff(restbg, el2) # all the remaining genes in background

    vals <- c(length(inlist), length(notinlist), length(onlybg), length(notinbg))
    vals.m <- matrix(vals, nrow=2, ncol=2)
    test.m <- fisher.test(vals.m, alternative="greater")
  
    hit <- length(inlist)
    ehi <- length(li2)/length(bglist) * length(el2)

    return(c(hit, ehi, unlist(test.m$es), test.m$p))
  }))
  
  rownames(goer) <- terms[,1]
  colnames(goer) <- c("hit","expected","OR","P")
  return(goer)

}

library(ROCR)

mkUnique <- function(res){
  #column 1: gene names
  #column 2: values

  genes <- sort(unique(res[,1]))
  res <- sapply(genes, function(g){
    #idx <- res[,1] == g
    #tmp <- res[idx,2]
    tmp <- subset(res,subset=GENE==g)[,2]
    tma <- max(tmp)
    tmi <- min(tmp)
    res <- tma
    if (abs(tmi) > tma)
      res <- tmi
    return(res)
  })
  res <- data.frame(genes, res)
  colnames(res) <- c("GENE","SCORE")
  return(res)
} 

#gene set enrichment analyis using ROC curves (AUC)
#works with weighted lists
#terms contains: name, link , genes
#the list must be unique though, i.e., one value per gene
#if this is not the case run mkUnique first
gene_enrichment_auc <- function(bglist, explist, terms){

  #column 1 of explist contains the gene name
  #column 2 of explist contains the weight

  el2 <- paste(intersect(bglist, explist[,1]))

  rownames(explist) <- paste(explist[,1])
  expnames  <- explist[el2,1]
  expvalues <- explist[el2,2]

  inout <- rep(0, length(el2))
  names(inout) <- expnames

  goer <- t(apply(terms, 1, function(x){
    #print(x[1])
    li <- unlist(strsplit(paste(x[3]), split=" ")[[1]])
    li2 <- paste(intersect(bglist, li))

    inout2 <- inout
    inout2[paste(intersect(li2,expnames))] <- 1    

    len1 <- length(li2)
    len2 <- sum(inout2)
    if (len2 > 0){
      pre <- prediction(expvalues, inout2)
      auc <- unlist(performance(pre, measure="auc")@y.values)

      wil.p <- wilcox.test(expvalues[inout2==0], expvalues[inout2==1])$p.v
    } else {
      auc <- NA
      wil.p <- NA
    }
    return(c(len1, len2, auc, wil.p))
  }))
  
  rownames(goer) <- terms[,1]
  colnames(goer) <- c("List_orig","List_used","AUC","WC_P")
  return(goer)

}