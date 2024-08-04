##Author: JooHee Choi, joohee.choi@wustl.edu, used by Erin Newcomer for ICU Sink project
##This script calculates pvals for genes, based on the randomly distributed mutation
##I prefer to open in RStudio and run the script line by line... haven't really tried just running the whole thing 

library(dplyr)
library(ggplot2)
library(coin)
library(FSA)



##########random mut analysis: formatting actual distribution of mutations into gene: count: snp: indel

mutfile <- read.csv('muts_perstrain.csv') ## output file of s_matrix_annotate.py

mutfile <- na.omit(mutfile) 

mutfile$gene <- gsub('hypotheticalprotein', '\\1', mutfile$gene)
mutfile$gene <- gsub('intron', '\\1', mutfile$gene)
mutfile$gene <- gsub('unknownnode', '\\1', mutfile$gene)
mutfile <- mutfile[mutfile$gene!='',]



##this part will get rid of allele numbers, if desired 
mutfile$gene <- gsub('_1','\\1',mutfile$gene)
mutfile$gene <- gsub('_2','\\1',mutfile$gene)
mutfile$gene <- gsub('_3','\\1',mutfile$gene)
mutfile$gene <- gsub('_4','\\1',mutfile$gene)
mutfile$gene <- gsub('_5','\\1',mutfile$gene)
mutfile$gene <- gsub('_6','\\1',mutfile$gene)
mutfile$gene <- gsub('_7','\\1',mutfile$gene)
mutfile$gene <- gsub('_8','\\1',mutfile$gene)

mutfile_summ <- mutfile %>% group_by(gene) %>% summarize(n()) # how many strains per gene 



genelist <- unique(mutfile$gene) 
n <- length(genelist) #3381 for random, 259 for actual , 238 if condense alleles 

counts <- c()
snps <- c()
indels <- c()



for (i in 1:n){
  temp <- as.character(genelist[i])
  m <- length(unique(mutfile[mutfile$gene==temp,]$strain)) #how many strains total? 
  counts <- c(counts, m)
  k <- length(unique(mutfile[mutfile$gene==temp & mutfile$type=='SNP',]$strain)) #which are snp ?
  l <- length(unique(mutfile[mutfile$gene==temp & mutfile$type=='indel',]$strain)) #which are indel?
  
  snps <- c(snps, k)
  indels <- c(indels, l)
}

realmuts <- data.frame(gene=genelist, count=counts, SNP=snps,indel=indels )
actual <- realmuts

#------------------permutation results-------------------------#

reference <- read.csv('randommuts.csv', header=F) ##output file of s_randommuts.py
colnames(reference) <- c('iter', 'gene','strain' )

#gonna keep alleles separate for now 

reference$gene <- gsub('_1','',reference$gene)
reference$gene <- gsub('_2','',reference$gene)
reference$gene <- gsub('_3','',reference$gene)
reference$gene <- gsub('_4','',reference$gene)
reference$gene <- gsub('_5','',reference$gene)
reference$gene <- gsub('_6','',reference$gene)
reference$gene <- gsub('_7','',reference$gene)
reference$gene <- gsub('_8','',reference$gene)



####trying to figure out true definition p-value

colnames(actual)[2] <- 'numstrain' 
genelist <- unique(actual$gene) 
k <- length(genelist)

actual$pval <- c('.')

for (i in 1:k){
  target <- as.character(genelist[i])
  real <- as.numeric(as.character(actual[actual$gene==target,]$numstrain[1]))
  if (target %in% reference$gene){
    tempdf <- reference[reference$gene==target,]
   #should have aggregated..
    tempdf$strain <- as.numeric(tempdf$strain)
    tempdf <- aggregate(tempdf$strain, by=list(iter=tempdf$iter), FUN=sum) ##this was because of the _1 and _2 genes
    colnames(tempdf) <- c('iter', 'strain')
    n <- length(rownames(tempdf))
    diff <- 1000-n
    
    iter <- rep('.', diff)
    strains <- rep(0, diff)
    
    tempdf2 <- data.frame(iter=iter, strain=strains)
    tempdf <- rbind(tempdf, tempdf2)
    
    percentile <- ecdf(tempdf$strain)
    pvalue <- 1-percentile(real)
    
    actual[actual$gene==target,]$pval <- pvalue 
  }
}

actual$pval <- unlist(actual$pval, use.names=FALSE)

write.csv(actual, 'permutation/pvals_no_alleles.csv')
