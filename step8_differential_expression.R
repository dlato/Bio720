###########################################################
#differential gene expression for the Bio 720 final project
###########################################################

####################
#setting up the data
####################
library(DESeq2)
library("RColorBrewer")
library("gplots")

setwd("C:/Users/Daniella/Documents/Bio720/SinoFinalProject/countData/")

in_dir = dir(, pattern=".txt")

counts_in <- lapply(in_dir, function(x) read.table(x, header=F, sep = "", nrows=8241))
#make matrix of counts
tot_count_matrix <- matrix(unlist(lapply(counts_in, function(x) x$V2)), nrow=8241, ncol=16)

#setting up experimental design object
parse_names <- strsplit(in_dir, split="_")

parse_names <- matrix(unlist(parse_names), nrow=16, ncol=9, byrow=T)

col_names_counts <- paste(parse_names[,1], "_", parse_names[,2], "_", parse_names[,3], "_", parse_names[,4], sep="")

colnames(tot_count_matrix) = col_names_counts # sample names as column names
# add contig names as row names
rownames(tot_count_matrix) = counts_in[[1]]$V1

dim(tot_count_matrix)

plot(log(tot_count_matrix[,1]),log(tot_count_matrix[,2]) )

head(parse_names)

experimental_design = data.frame(
  sample_names = col_names_counts,  # sample name
  individual = factor(parse_names[,1]), # each individual sample
  treatment = factor(parse_names[,6]),  # missing replicons or not (wt= all three reps, a=missing plasmid a, b=mising plasmid b, ab=missing both pa and pb)
  lane = factor(parse_names[,3])      # Which lane on the Illumina flowcell.
)

#set levels so wild type (wt) is the reference and will be compared against
experimental_design$treatment <- relevel(experimental_design$treatment, "wt")

DESeq_data <- DESeqDataSetFromMatrix(tot_count_matrix, experimental_design, 
                                     design = formula(~ treatment))

DESeq_data <- DESeq(DESeq_data)
plotDispEsts(DESeq_data, xlab="Mean of Normalized Counts", ylab="Dispersion", main="Mean Dispersion")
plotMA(DESeq_data, ylim=c(-2,2))

###########################
#comparing wt and mising a
###########################
res_wt_a <- results(DESeq_data, contrast=c("treatment","wt","a"), alpha=0.05) #comparing the different mutant strains or "treatments"
summary(res_wt_a, alpha=0.05)
head(res_wt_a)
hist(res_wt_a$padj)


wt_a_dat <- data.frame(res_wt_a@listData$baseMean, res_wt_a@listData$log2FoldChange, res_wt_a@listData$padj )
colnames(wt_a_dat) <- c("baseMean", "log2FoldChange", "padj" )

plot(y = res_wt_a$log2FoldChange, x =log2(res_wt_a$baseMean) , 
     ylab = "log2 Fold Change", xlab = "Mean Expression",
     main = "Wild Type v.s. Missing Plasmid A", pch = 20, col ="grey")
with(wt_a_dat[wt_a_dat$padj <0.05, ],
     points(y = log2FoldChange,x = log2(baseMean),
            pch=20, col = "black" ))  
abline(a=0, b=0 , col="blue")

res_wt_a_padj <- res_wt_a$padj[which(res_wt_a$padj < 0.05)] # Sorting genes based on log2 fold change
#list of differentially expressed genes
de_genes_wt_a_pos = vector()
for (j in 1:length(res_wt_a_padj)){
  de_genes_wt_a_pos <- c(de_genes_wt_a_pos,which(res_wt_a$padj==res_wt_a_padj[j]))#list of position in vector that have outliers
}
de_genes_wt_a <- row.names(res_wt_a)
de_genes_wt_a <- de_genes_wt_a[c(de_genes_wt_a_pos)]

wt_a_dat_padj<- wt_a_dat[wt_a_dat$padj < 0.05, ]#adjusted pval < 0.05
wt_a_upreg_genes<- wt_a_dat_padj[wt_a_dat_padj$log2FoldChange > 0, ]#up regulated genes
wt_a_tmp <- wt_a_upreg_genes$log2FoldChange
wt_a_upreg_vals <- wt_a_tmp[!is.na(wt_a_tmp)]#values of the upreg genes

wt_a_upreg_genes_pos = vector()
for (j in 1:length(wt_a_upreg_vals)){
  wt_a_upreg_genes_pos <- c(wt_a_upreg_genes_pos,which(wt_a_dat$log2FoldChange==wt_a_upreg_vals[j]))#list of position in vector that have log2FoldChange > 0 (upreg)
}
wt_a_upreg_genes <- row.names(res_wt_a)
wt_a_upreg_genes <- wt_a_upreg_genes[c(wt_a_upreg_genes_pos)]
sink("wt_a_upreg_genes.txt")
wt_a_upreg_genes#list of upreg genes (just names)
sink()

wt_a_downreg_genes<- wt_a_dat_padj[wt_a_dat_padj$log2FoldChange < 0, ]#genes downregulated
wt_a_tmp <- wt_a_downreg_genes$log2FoldChange
wt_a_downreg_vals <- wt_a_tmp[!is.na(wt_a_tmp)]#values of down reg genes

wt_a_downreg_genes_pos = vector()
for (j in 1:length(wt_a_downreg_vals)){
  wt_a_downreg_genes_pos <- c(wt_a_downreg_genes_pos,which(wt_a_dat$log2FoldChange==wt_a_downreg_vals[j]))#list of position in vector that have log2FoldChange > 0 (downreg)
}
wt_a_downreg_genes <- row.names(res_wt_a)
wt_a_downreg_genes <- wt_a_downreg_genes[c(wt_a_downreg_genes_pos)]
sink("wt_a_downreg_genes.txt")
wt_a_downreg_genes#list of down reg gene names
sink()

############################
#comparing wt and missing b
############################
res_wt_b <- results(DESeq_data, contrast=c("treatment","wt","b"), alpha=0.05) #comparing the different mutant strains or "treatments"
summary(res_wt_b, alpha=0.05)
head(res_wt_b)
hist(res_wt_b$pvalue)


wt_b_dat <- data.frame(res_wt_b@listData$baseMean, res_wt_b@listData$log2FoldChange, res_wt_b@listData$padj )
colnames(wt_b_dat) <- c("baseMean", "log2FoldChange", "padj" )

plot(y = res_wt_b$log2FoldChange, x =log2(res_wt_b$baseMean) , 
     ylab = "log2 Fold Change", xlab = "Mean Expression",
     main = "Wild Type v.s. Missing Plasmid B", pch = 20, col ="grey")
with(wt_b_dat[wt_b_dat$padj <0.05, ],
     points(y = log2FoldChange,x = log2(baseMean),
            pch=20, col = "black" ))  
abline(a=0, b=0 , col="blue")

res_wt_b_sorted <- res_wt_b[order(-res_wt_b$log2FoldChange),] # Sorting genes based on "significance"... yuck.
head(res_wt_b_sorted)

 res_wt_b_padj <- res_wt_b$padj[which(res_wt_b$padj < 0.05)] # Sorting genes based on log2 fold change
  
 de_genes_wt_b_pos = vector()
 for (j in 1:length(res_wt_b_padj)){
   de_genes_wt_b_pos <- c(de_genes_wt_b_pos,which(res_wt_b$padj==res_wt_b_padj[j]))#list of position in vector that have outliers
 }
 de_genes_wt_b <- row.names(res_wt_b)
 de_genes_wt_b <- de_genes_wt_b[c(de_genes_wt_b_pos)]
 head(de_genes_wt_b)
 
wt_b_dat_padj<- wt_b_dat[wt_b_dat$padj < 0.05, ]
wt_b_upreg_genes<- wt_b_dat_padj[wt_b_dat_padj$log2FoldChange > 0, ]
wt_b_tmp <- wt_b_upreg_genes$log2FoldChange
wt_b_upreg_vals <- wt_b_tmp[!is.na(wt_b_tmp)]

wt_b_upreg_genes_pos = vector()
for (j in 1:length(wt_b_upreg_vals)){
  wt_b_upreg_genes_pos <- c(wt_b_upreg_genes_pos,which(wt_b_dat$log2FoldChange==wt_b_upreg_vals[j]))#list of position in vector that have log2FoldChange > 0 (upreg)
}
wt_b_upreg_genes <- row.names(res_wt_b)
wt_b_upreg_genes <- wt_b_upreg_genes[c(wt_b_upreg_genes_pos)]
sink("wt_b_upreg_genes.txt")
wt_b_upreg_genes
sink()

wt_b_downreg_genes<- wt_b_dat_padj[wt_b_dat_padj$log2FoldChange < 0, ]
wt_b_tmp <- wt_b_downreg_genes$log2FoldChange
wt_b_downreg_vals <- wt_b_tmp[!is.na(wt_b_tmp)]

wt_b_downreg_genes_pos = vector()
for (j in 1:length(wt_b_downreg_vals)){
  wt_b_downreg_genes_pos <- c(wt_b_downreg_genes_pos,which(wt_b_dat$log2FoldChange==wt_b_downreg_vals[j]))#list of position in vector that have log2FoldChange > 0 (downreg)
}
wt_b_downreg_genes <- row.names(res_wt_b)
wt_b_downreg_genes <- wt_b_downreg_genes[c(wt_b_downreg_genes_pos)]
sink("wt_b_downreg_genes.txt")
wt_b_downreg_genes
sink()

##############################
#comparing wt and missing ab
##############################
res_wt_ab <- results(DESeq_data, contrast=c("treatment","wt","ab"), alpha=0.05) #comparing the different mutant strains or "treatments"
summary(res_wt_ab, alpha=0.05)
head(res_wt_ab)
hist(res_wt_ab$pvalue)


wt_ab_dat <- data.frame(res_wt_ab@listData$baseMean, res_wt_ab@listData$log2FoldChange, res_wt_ab@listData$padj )
colnames(wt_ab_dat) <- c("baseMean", "log2FoldChange", "padj" )

plot(y = res_wt_ab$log2FoldChange, x =log2(res_wt_ab$baseMean) , 
     ylab = "log2 Fold Change", xlab = "Mean Expression",
     main = "Wild Type v.s. Missing Both Plasmids", pch = 20, col ="grey")
with(wt_b_dat[wt_ab_dat$padj <0.05, ],
     points(y = log2FoldChange,x = log2(baseMean),
            pch=20, col = "black" ))  
abline(a=0, b=0 , col="blue")

res_wt_ab_sorted <- res_wt_ab[order(-res_wt_ab$log2FoldChange),] # Sorting genes based on "significance"... yuck.
head(res_wt_ab_sorted)

  res_wt_ab_padj <- res_wt_ab$padj[which(res_wt_ab$padj < 0.05)] # Sorting genes based on log2 fold change
   
  de_genes_wt_ab_pos = vector()
  for (j in 1:length(res_wt_ab_padj)){
    de_genes_wt_ab_pos <- c(de_genes_wt_ab_pos,which(res_wt_ab$padj==res_wt_ab_padj[j]))#list of position in vector that have outliers
  }
  de_genes_wt_ab <- row.names(res_wt_ab)
  de_genes_wt_ab <- de_genes_wt_ab[c(de_genes_wt_ab_pos)]
  head(de_genes_wt_ab)
 
wt_ab_dat_padj<- wt_ab_dat[wt_ab_dat$padj < 0.05, ]
wt_ab_upreg_genes<- wt_ab_dat_padj[wt_ab_dat_padj$log2FoldChange > 0, ]
wt_ab_tmp <- wt_ab_upreg_genes$log2FoldChange
wt_ab_upreg_vals <- wt_ab_tmp[!is.na(wt_ab_tmp)]

wt_ab_upreg_genes_pos = vector()
for (j in 1:length(wt_ab_upreg_vals)){
  wt_ab_upreg_genes_pos <- c(wt_ab_upreg_genes_pos,which(wt_ab_dat$log2FoldChange==wt_ab_upreg_vals[j]))#list of position in vector that have log2FoldChange > 0 (upreg)
}
wt_ab_upreg_genes <- row.names(res_wt_ab)
wt_ab_upreg_genes <- wt_ab_upreg_genes[c(wt_ab_upreg_genes_pos)]
sink("wt_ab_upreg_genes.txt")
wt_ab_upreg_genes
sink()

wt_ab_downreg_genes<- wt_ab_dat_padj[wt_ab_dat_padj$log2FoldChange < 0, ]
wt_ab_tmp <- wt_ab_downreg_genes$log2FoldChange
wt_ab_downreg_vals <- wt_ab_tmp[!is.na(wt_ab_tmp)]

wt_ab_downreg_genes_pos = vector()
for (j in 1:length(wt_ab_downreg_vals)){
  wt_ab_downreg_genes_pos <- c(wt_ab_downreg_genes_pos,which(wt_ab_dat$log2FoldChange==wt_ab_downreg_vals[j]))#list of position in vector that have log2FoldChange > 0 (downreg)
}
wt_ab_downreg_genes <- row.names(res_wt_ab)
wt_ab_downreg_genes <- wt_ab_downreg_genes[c(wt_ab_downreg_genes_pos)]
sink("wt_ab_downreg_genes.txt")
wt_ab_downreg_genes
sink()

####################################
  #comparing missing a and missing b
####################################
  res_a_b <- results(DESeq_data, contrast=c("treatment","a","b"), alpha=0.05) #comparing the different mutant strains or "treatments"
  summary(res_a_b, alpha=0.05)
  head(res_a_b)
  hist(res_a_b$pvalue)
   
   
  a_b_dat <- data.frame(res_a_b@listData$baseMean, res_a_b@listData$log2FoldChange, res_a_b@listData$padj )
  colnames(a_b_dat) <- c("baseMean", "log2FoldChange", "padj" )
   
  plot(y = res_a_b$log2FoldChange, x =log2(res_a_b$baseMean) , 
       ylab = "log2 Fold Change", xlab = "Mean Expression",
       main = "Missing Plasmid A v.s. Missing Plasmid B", pch = 20, col ="grey")
  with(a_b_dat[a_b_dat$padj <0.05, ],
       points(y = log2FoldChange,x = log2(baseMean),
              pch=20, col = "black" ))  
  abline(a=0, b=0 , col="blue")
   
  res_a_b_sorted <- res_a_b[order(-res_a_b$log2FoldChange),] # Sorting genes based on "significance"... yuck.
  head(res_a_b_sorted)
   
   res_a_b_padj <- res_a_b$padj[which(res_a_b$padj < 0.05)] # Sorting genes based on log2 fold change
    
   de_genes_a_b_pos = vector()
   for (j in 1:length(res_a_b_padj)){
     de_genes_a_b_pos <- c(de_genes_a_b_pos,which(res_a_b$padj==res_a_b_padj[j]))#list of position in vector that have outliers
   }
   de_genes_a_b <- row.names(res_a_b)
   de_genes_a_b <- de_genes_a_b[c(de_genes_a_b_pos)]
   head(de_genes_a_b)
   
   
#####################################   
#comparing missing both and missing a
#####################################
   res_ab_a <- results(DESeq_data, contrast=c("treatment","ab","a"), alpha=0.05) #comparing the different mutant strains or "treatments"
   summary(res_ab_a, alpha=0.05)
   head(res_ab_a)
   hist(res_ab_a$pvalue)
    
    
   ab_a_dat <- data.frame(res_ab_a@listData$baseMean, res_ab_a@listData$log2FoldChange, res_ab_a@listData$padj )
   colnames(ab_a_dat) <- c("baseMean", "log2FoldChange", "padj" )
    
   plot(y = res_ab_a$log2FoldChange, x =log2(res_ab_a$baseMean) , 
        ylab = "log2 Fold Change", xlab = "Mean Expression",
        main = "Missing Both Plasmis v.s. and Missing Plasmid A", pch = 20, col ="grey")
   with(ab_a_dat[ab_a_dat$padj <0.05, ],
        points(y = log2FoldChange,x = log2(baseMean),
               pch=20, col = "black" ))  
   abline(a=0, b=0 , col="blue")
    
   res_ab_a_sorted <- res_ab_a[order(-res_ab_a$log2FoldChange),] # Sorting genes based on "significance"... yuck.
   head(res_ab_a_sorted)
    
    res_ab_a_padj <- res_ab_a$padj[which(res_ab_a$padj < 0.05)] # Sorting genes based on log2 fold change
     
    de_genes_ab_a_pos = vector()
    for (j in 1:length(res_ab_a_padj)){
      de_genes_ab_a_pos <- c(de_genes_ab_a_pos,which(res_ab_a$padj==res_ab_a_padj[j]))#list of position in vector that have outliers
    }
    de_genes_ab_a <- row.names(res_ab_a)
    de_genes_ab_a <- de_genes_ab_a[c(de_genes_ab_a_pos)]
    head(de_genes_ab_a)
    
    
#####################################    
#comparing missing both and missing b
#####################################
    res_ab_b <- results(DESeq_data, contrast=c("treatment","ab","b"), alpha=0.05) #comparing the different mutant strains or "treatments"
    summary(res_ab_b, alpha=0.05)
    head(res_ab_b)
    hist(res_ab_b$pvalue)
     
     
    ab_b_dat <- data.frame(res_ab_b@listData$baseMean, res_ab_b@listData$log2FoldChange, res_ab_b@listData$padj )
    colnames(ab_b_dat) <- c("baseMean", "log2FoldChange", "padj" )
     
    plot(y = res_ab_b$log2FoldChange, x =log2(res_ab_b$baseMean) , 
         ylab = "log2 Fold Change", xlab = "Mean Expression",
         main = "Missing Both Plasmids and Missing Plasmid B", pch = 20, col ="grey")
    with(ab_b_dat[ab_b_dat$padj <0.05, ],
         points(y = log2FoldChange,x = log2(baseMean),
                pch=20, col = "black" ))  
    abline(a=0, b=0 , col="blue")
     
    res_ab_b_sorted <- res_ab_b[order(-res_ab_b$log2FoldChange),] # Sorting genes based on "significance"... yuck.
    head(res_ab_b_sorted)
     
     res_ab_b_padj <- res_ab_b$padj[which(res_ab_b$padj < 0.05)] # Sorting genes based on log2 fold change
      
     de_genes_ab_b_pos = vector()
     for (j in 1:length(res_ab_b_padj)){
       de_genes_ab_b_pos <- c(de_genes_ab_b_pos,which(res_ab_b$padj==res_ab_b_padj[j]))#list of position in vector that have outliers
     }
     de_genes_ab_b <- row.names(res_ab_b)
     de_genes_ab_b <- de_genes_ab_b[c(de_genes_ab_b_pos)]
     head(de_genes_ab_b)
     
ab_b_dat_padj<- ab_b_dat[ab_b_dat$padj < 0.05, ]
ab_b_upreg_genes<- ab_b_dat_padj[ab_b_dat_padj$log2FoldChange > 0, ]
ab_b_tmp <- ab_b_upreg_genes$log2FoldChange
ab_b_upreg_vals <- ab_b_tmp[!is.na(ab_b_tmp)]

ab_b_upreg_genes_pos = vector()
for (j in 1:length(ab_b_upreg_vals)){
  ab_b_upreg_genes_pos <- c(ab_b_upreg_genes_pos,which(ab_b_dat$log2FoldChange==ab_b_upreg_vals[j]))#list of position in vector that have log2FoldChange > 0 (upreg)
}
ab_b_upreg_genes <- row.names(res_ab_b)
ab_b_upreg_genes <- ab_b_upreg_genes[c(ab_b_upreg_genes_pos)]
sink("ab_b_upreg_genes.txt")
ab_b_upreg_genes
sink()

ab_b_downreg_genes<- ab_b_dat_padj[ab_b_dat_padj$log2FoldChange > 0, ]
ab_b_tmp <- ab_b_downreg_genes$log2FoldChange
ab_b_downreg_vals <- ab_b_tmp[!is.na(ab_b_tmp)]

ab_b_downreg_genes_pos = vector()
for (j in 1:length(ab_b_downreg_vals)){
  ab_b_downreg_genes_pos <- c(ab_b_downreg_genes_pos,which(ab_b_dat$log2FoldChange==ab_b_downreg_vals[j]))#list of position in vector that have log2FoldChange > 0 (downreg)
}
ab_b_downreg_genes <- row.names(res_ab_b)
ab_b_downreg_genes <- ab_b_downreg_genes[c(ab_b_downreg_genes_pos)]
sink("ab_b_downreg_genes.txt")
ab_b_downreg_genes
sink()

sink("de_genes_wt_a.txt")#puts this output into a file so I can play with it in unix instead of R
de_genes_wt_a
sink()
sink("de_genes_wt_b.txt")#puts this output into a file so I can play with it in unix instead of R
de_genes_ab_b
sink()
sink("de_genes_wt_ab.txt")#puts this output into a file so I can play with it in unix instead of R
de_genes_wt_ab
sink()
