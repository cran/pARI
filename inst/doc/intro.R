## ----include = FALSE, echo = FALSE, message = FALSE---------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## ----setup--------------------------------------------------------------------
#devtools::install_github("angeella/pARI")
#install.packages("pARI")
library(pARI)

## -----------------------------------------------------------------------------
datas <- simulateData(pi0 = 0.8, m = 1000, n = 30, power = 0.9, rho = 0.5,seed = 123)

## -----------------------------------------------------------------------------
out <- pARI(X = datas, ix = c(1:200), test.type = "one_sample", seed = 123)
out$TDP

## -----------------------------------------------------------------------------
out <- signTest(X = datas, B = 1000, rand = F)
P <- cbind(out$pv, out$pv_H0)
pARI(pvalues = P, ix = c(1:200),test.type = "one_sample")$TDP

## -----------------------------------------------------------------------------
ix <- sample(c(1:4), size = 1000, replace = T)
out <- pARI(pvalues = P, ix = ix,test.type = "one_sample", clusters = TRUE)$TDP
out

## ----eval = FALSE-------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE)){
#    install.packages("BiocManager")
#  }
#  
#  
#  if (!requireNamespace("dynamicTreeCut", quietly = TRUE)){
#    install.packages("dynamicTreeCut")
#  }
#  
#  
#  #BiocManager::install(c("Biobase","genefilter", "EnrichmentBrowser"))
#  
#  library(Biobase)
#  library(genefilter)
#  library(dynamicTreeCut)

## ----eval = FALSE-------------------------------------------------------------
#  load(file=url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData"))
#  
#  pdata<- pData(montpick.eset)
#  edata <- as.matrix(exprs(montpick.eset))
#  fdata <- fData(montpick.eset)
#  
#  edata <- log2(as.matrix(edata) + 1)
#  edata <- edata[rowMeans(edata) > 10, ]
#  
#  my.dist <- dist(edata)
#  my.tree <- hclust(my.dist, method="ward.D2")
#  
#  my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist),minClusterSize=10))

## ----eval = FALSE-------------------------------------------------------------
#  out <-pARI(X = edata,alpha = 0.05, test.type = "two_samples",
#             label = as.factor(pdata$population),
#             ix = my.clusters, family = "simes", clusters = TRUE)
#  out$TDP

## ----eval = FALSE-------------------------------------------------------------
#  pathways <- EnrichmentBrowser::getGenesets(org = "hsa", db = "kegg", gene.id.type = "ENSEMBL")
#  
#  out <- c()
#  for(i in seq(pathways)){
#  
#    ix <- which(rownames(edata) %in% pathways[[i]])
#    if(length(ix)!=0){
#    out[i] <-pARI(X = edata,alpha = 0.05, test.type = "two_samples",
#             label = as.factor(pdata$population),
#             ix = ix, family = "simes", clusters = TRUE)$TDP
#    }else{
#    out[i] <- NA
#    }
#  }
#  
#  db <- data.frame(TDP = out, size = sapply(seq(pathways), function(x) length(which(rownames(edata) %in% pathways[[x]]))), name = names(pathways))

## ----eval = FALSE-------------------------------------------------------------
#  library(tidyverse)
#  db <- db %>% dplyr::filter(!is.na(TDP) | TDP!=0) %>% dplyr::filter(size >10)
#  
#  db$name = factor(db$name, levels=db[order(db$TDP), "name"])
#  
#  db %>% dplyr::filter(!is.na(TDP) | TDP!=0) %>% dplyr::filter(size >10) %>% arrange(TDP)%>% ggplot(aes(x = TDP, y = name, size = size)) + geom_point() +
#    xlab(expression(bar(pi)(S[m]))) + ylab("Pathways") + labs(size = expression(paste("|", S, "|"))) + theme_classic()+
#    scale_size_continuous(breaks = c(15, 25, 35, 45))
#  

## ----eval = FALSE-------------------------------------------------------------
#  if (!requireNamespace("fMRIdata", quietly = TRUE)){
#    remotes::install_github("angeella/fMRIdata")
#  }
#  library(fMRIdata)
#  data(Auditory_clusterTH3_2)
#  data(Auditory_copes)
#  data(Auditory_mask)

## ----eval = FALSE-------------------------------------------------------------
#  auditory_out <- pARIbrain(copes = Auditory_copes, clusters = Auditory_clusterTH3_2, mask = Auditory_mask, alpha = 0.05, silent = TRUE)
#  auditory_out$out

