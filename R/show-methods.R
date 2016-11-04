# Genotypes ----

setMethod(
    "show", "Genotypes",
    function(object){
        validObject(object)
        cat("Class: Genotypes\n")
        cat(sprintf(
            "  @ref (hom. ref.): \"%s\" {%s}\n",
            suffix(object)["ref"],
            paste(object@ref, collapse = ", ")))
        cat(sprintf(
            "  @het (heter.)   : \"%s\" {%s}\n",
            suffix(object)["het"],
            paste(object@het, collapse = ", ")))
        cat(sprintf(
            "  @alt (hom. alt.): \"%s\" {%s}",
            suffix(object)["alt"],
            paste(object@alt, collapse = ", ")))
        return(NULL)
    }
)

# TVTBparam ----

setMethod(
    "show", "TVTBparam",
    function(object){
        cat("class: TVTBparam\n")
        cat(  "@genos: class: Genotypes\n")
        g <- object@genos
        validObject(g)
        cat(sprintf(
            "    @ref (hom. ref.): \"%s\" {%s}\n",
            suffix(g)["ref"],
            paste(g@ref, collapse = ", ")))
        cat(sprintf(
            "    @het (heter.): \"%s\" {%s}\n",
            suffix(g)["het"],
            paste(g@het, collapse = ", ")))
        cat(sprintf(
            "    @alt (hom. alt.): \"%s\" {%s}\n",
            suffix(g)["alt"],
            paste(g@alt, collapse = ", ")))
        cat(sprintf(
            "  @ranges: %i GRanges on %i sequence(s)\n",
            length(ranges(object)),
            length(seqlevels(ranges(object)))))
        cat(sprintf(
            "  @aaf (alt. allele freq.): \"%s\"\n",
            aaf(object)))
        cat(sprintf(
            "  @maf (minor allele freq.): \"%s\"\n",
            maf(object)))
        cat(sprintf(
            "  @vep (Ensembl VEP key): \"%s\"\n",
            vep(object)))
        cat("  @svp: <ScanVcfParam object>\n")
        cat(sprintf("  @bp: <%s object>\n", class(bp(object))))
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
