# tSVEParam ----

setMethod(
    f = "show",
    signature = "tSVEParam",
    definition = function(object){
        cat("Class: tSVEParam\n")
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
            "ranges: %i range(s) on %i sequence(s)\n",
            length(object@ranges),
            length(seqlevels(object@ranges))))
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

# VcfBasicFilter ----

setMethod(
    f = "show",
    signature = "vcfFixedFilter",
    definition = function(object){
        cat(sprintf(
            "vcfFixedFilter: %s %s %s\n",
            slot(object, "name"),
            slot(object, "condition"),
            slot(object, "value")))
    }
)

setMethod(
    f = "show",
    signature = "vcfInfoFilter",
    definition = function(object){
        cat(sprintf(
            "vcfInfoFilter: %s %s %s\n",
            slot(object, "name"),
            slot(object, "condition"),
            slot(object, "value")))
    }
)

setMethod(
    f = "show",
    signature = "vcfVepFilter",
    definition = function(object){
        cat(sprintf(
            "vcfInfoFilters: %s %s %s\n",
            slot(object, "name"),
            slot(object, "condition"),
            slot(object, "value")))
    }
)
