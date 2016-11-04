# TVTBparam ----

setMethod(
    f = "show",
    signature = "TVTBparam",
    definition = function(object){
        cat("Class: TVTBparam\n")
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
        cat("VcfFixedFilter:", as.character(object))
    }
)

setMethod(
    f = "show",
    signature = "VcfInfoFilter",
    definition = function(object){
        cat("VcfInfoFilter:", as.character(object))
    }
)

setMethod(
    f = "show",
    signature = "VcfVepFilter",
    definition = function(object){
        cat("VcfVepFilter:", as.character(object))
    }
)

setMethod(
    f = "show",
    signature = "VcfFilterList",
    definition = function(object){
        # TODO: "Active" column
        cat("Index\tActive\tType\tFilter\n")
        mapply(
            FUN = function(idx, active, type, filter){
                cat(sprintf(
                    "%i.\t[x]\t%s\t%s\n",
                    idx,
                    ifelse(active, "x", " "),
                    type,
                    as.character(filter)))
                },
            seq_along(filterRules(object)),
            active(object),
            sapply(X = filterRules(object), FUN = filterType),
            filterRules(object))
    }
)
