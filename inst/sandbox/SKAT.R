
library(SKAT)

tparam
genos(tparam)
str(genos(tparam))

dim(geno(vcf)[["GT"]])

g <- geno(vcf)[["GT"]]
class(g)
g[1:4,1:4]

g2 <- matrix(NA_integer_, nrow = nrow(g), ncol = ncol(g), dimnames = dimnames(g))
for (gt in ref(tparam)){
    g2[g == gt] <- 0
}
for (gt in het(tparam)){
    g2[g == gt] <- 1
}
for (gt in alt(tparam)){
    g2[g == gt] <- 2
}
g2[1:4,1:4]

geno(header(vcf)) <- rbind(
    geno(header(vcf)),
    DataFrame(
        Number=1,
        Type="Float",
        Description="Genotype (0,1,2,NA)",
        row.names = "012"
    )
)
geno(vcf)[["012"]] <- t(g2)
rm(g, g2)

VCF.rare <- subsetByFilter(vcf, vcfRules["rare"])

SKAT_Null_Model()
SKAT(Z = t(geno(vcf)[["012"]]), )
