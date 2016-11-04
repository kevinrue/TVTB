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
    signature = "VcfFixedFilter",
    definition = function(object){
        cat(sprintf(
            "VcfFixedFilter: %s %s %s\n",
            slot(object, "name"),
            slot(object, "condition"),
            slot(object, "value")))
    }
)

setMethod(
    f = "show",
    signature = "VcfInfoFilter",
    definition = function(object){
        cat(sprintf(
            "VcfInfoFilter: %s %s %s\n",
            slot(object, "name"),
            slot(object, "condition"),
            slot(object, "value")))
    }
)

setMethod(
    f = "show",
    signature = "VcfVepFilter",
    definition = function(object){
        cat(sprintf(
            "VcfVepFilter: %s %s %s\n",
            slot(object, "name"),
            slot(object, "condition"),
            slot(object, "value")))
    }
)
