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
            "ranges: %i GRanges on %i sequence(s)\n",
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
        return(NULL)
    }
)

# VcfFixedRules ----

# setMethod(
#     f = "show",
#     signature = "VcfFixedRules",
#     definition = function(object){
#         cat(sprintf("VcfFilterRules of length %i\n", length(object)))
#         cat(sprintf(
#             "names: %s\n",
#             paste(names(object), collapse = " ")
#         ))
#         cat(sprintf(
#             "active: %s\n",
#             paste(active(object), collapse = " ")
#         ))
#     }
# )

# VcfInfoRules ----

# setMethod(
#     f = "show",
#     signature = "VcfInfoRules",
#     definition = function(object){
#         cat(sprintf("VcfFilterRules of length %i\n", length(object)))
#         cat(sprintf(
#             "names: %s\n",
#             paste(names(object), collapse = " ")
#         ))
#         cat(sprintf(
#             "active: %s\n",
#             paste(active(object), collapse = " ")
#         ))
#     }
# )

# VcfVepRules ----

# setMethod(
#     f = "show",
#     signature = "VcfVepRules",
#     definition = function(object){
#         cat(sprintf("VcfFilterRules of length %i\n", length(object)))
#         cat(sprintf(
#             "names: %s\n",
#             paste(names(object), collapse = " ")
#         ))
#         cat(sprintf(
#             "active: %s\n",
#             paste(active(object), collapse = " ")
#         ))
#     }
# )

# VcfFilterRules ----

# setMethod(
#     f = "show",
#     signature = "VcfFilterRules",
#     definition = function(object){
#         cat(sprintf("VcfFilterRules of length %i\n", length(object)))
#         cat(sprintf(
#             "names: %s\n",
#             paste(names(object), collapse = " ")
#             ))
#         cat(sprintf(
#             "type: %s\n",
#             paste(type(object), collapse = " ")
#         ))
#         cat(sprintf(
#             "active: %s\n",
#             paste(active(object), collapse = " ")
#         ))
#     }
# )
