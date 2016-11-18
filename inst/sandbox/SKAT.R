
# Run Introduction vignette code before the code below

library(SKAT)
library(TVTB)
library(VariantAnnotation)

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
geno(vcf)[["012"]] <- g2
rm(g, g2)

VCF.rare <- subsetByFilter(vcf, vcfRules["rare"])
VCF.EUR.AFR <- vcf[,colData(vcf)[,"super_pop"] %in% c("EUR", "AFR")]
colData(VCF.EUR.AFR) <- droplevels(colData(VCF.EUR.AFR))

skat.rules <- VcfInfoRules(list(
    EURorAFR = expression(super_pop_EUR_AAF > 0 | super_pop_AFR_AAF > 0)
))

summary(evalSeparately(skat.rules, VCF.EUR.AFR))

VCF.EUR.AFR <- subsetByFilter(VCF.EUR.AFR, skat.rules)

obj <- SKAT_Null_Model(
    (as.numeric(colData(VCF.EUR.AFR)[,"super_pop"]) - 1 ) ~ 1,
    out_type="D"
)

SKAT(t(geno(VCF.EUR.AFR)[["012"]]), obj, kernel = "linear.weighted")$p.value
