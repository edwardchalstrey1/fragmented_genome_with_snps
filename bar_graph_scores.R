random_scores <- as.vector(as.matrix(read.table("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/re_files/random_scores.txt", quote="\"")))

scores <- as.vector(as.matrix(read.table("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/re_files/scores.txt", quote="\"")))

# Standard error function:

st_err <- function(x) sd(x)/sqrt(length(x))

# So then do....

random_st <- st_err(random_scores)

serrors <- c(0, 0, random_st, 0, 0, 0, 0)

methods <- c("C0", "C1", "C2", "M1 a", "M1 b", "M2 a", "M2 b")

png("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/figures/rearrangement_methods.png")
barx <- barplot(scores, main="Comparison of Methods
                for Rearrangement of SNP Density Ordered Fragments,
                to the Original Sequenced Order", xlab="Methods", ylab="Ordinal Similarity Score", names.arg=methods)



error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

error.bar(barx, scores, serrors)

dev.off()