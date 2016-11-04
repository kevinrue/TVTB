# Define a random subset of samples
Sindices <- sample(x = seq(1, ncol(vcf)), size = 5)
Sindices
Snames <- colnames(vcf)[Sindices]
Snames

# Define a random subset of variants

Vindices <- sample(x = seq(1, nrow(vcf)), size = 10)
Vindices
Vnames <- rownames(vcf)[Vindices]
Vnames

# Subset by indices and names

dim(vcf[,Sindices])
dim(vcf[,Snames])
dim(vcf[Vindices,])
dim(vcf[Vnames,])
dim(vcf[Vindices, Snames])
