# tSVEParam ----

setMethod(
    f = "show",
    signature = "tSVEParam",
    definition = function(object){
        cat("Genotypes\n")
        cat(sprintf(
            "- Reference: \"%s\" {%s}\n",
            names(object@genos)[1],
            paste(object@genos[[1]], collapse = ", ")))
        cat(sprintf(
            "- Heterozygote: \"%s\" {%s}\n",
            names(object@genos)[2],
            paste(object@genos[[2]], collapse = ", ")))
        cat(sprintf(
            "- Alternative: \"%s\" {%s}\n",
            names(object@genos)[3],
            paste(object@genos[[3]], collapse = ", ")))
        cat(sprintf(
            "Alternate allele frequency: \"%s\"\n",
            object@aaf))
        cat(sprintf(
            "Minor allele frequency: \"%s\"\n",
            object@maf))
        cat(sprintf(
            "Ensembl Variant Effet Predictor: \"%s\"\n",
            object@vep))
        cat("Bioconductor parallel parameters\n")
        cat("- ")
        print(object@bp)
    }
)
