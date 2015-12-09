#differential gene expression for the Bio 720 final project

library(DESeq2)
library("RColorBrewer")
library("gplots")

setwd("C:/Users/Daniella/Documents/Bio720/SinoFinalProject/countData/")

in_dir = dir(, pattern=".txt")
in_dir

counts_in <- lapply(in_dir, function(x) read.table(x, header=F, sep = "", nrows=8257))
head(counts_in[[1]])
head(counts_in[[3]])
#make matrix of counts
tot_count_matrix <- matrix(unlist(lapply(counts_in, function(x) x$V2)), nrow=8257, ncol=16)

head(tot_count_matrix)
#setting up experimental design object
parse_names <- strsplit(in_dir, split="_")
parse_names

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

#test if there is differential experssion based only on lane
test_lane_effects <- DESeqDataSetFromMatrix(tot_count_matrix, experimental_design, 
                                            design = formula(~ lane))

test_lane_effects2 <- DESeq(test_lane_effects) # We know fit the simple model
test_lane_effects2_results <- results(test_lane_effects2)
summary(test_lane_effects2_results) # No evidence, but this is a bit incomplete

plotDispEsts(test_lane_effects2)

for_pca <- rlog(test_lane_effects2, blind=TRUE)
plotPCA(for_pca, intgroup=c("lane")) # no obvious lane effects.

#hierarchical clustering to test lane effect
rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion
distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix
rownames(mat) <- colnames(mat) <-     with(colData(test_lane_effects2), paste(treatment, lane, sep=" : "))

hc <- hclust(distsRL)  # performs hierarchical clustering
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)  # picking our colours
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13))
