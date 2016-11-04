
colnames(mcols(csq))
table(mcols(csq)[,"VARIANT_CLASS"])

vepInsertion <- VcfVepRules(exprs = list(
    isInsertion = expression(VARIANT_CLASS == "insertion")
))

evcf_vepIns <- subsetByFilter(evcf, vepInsertion)
as.data.frame(fixed(evcf_vepIns)[,c("REF", "ALT")])

insFilter <- VcfFixedRules(exprs = list(
    isInsertion = expression(
        Biostrings::width(ALT) > Biostrings::width(REF)
    )
))

evcf_fixedIns <- subsetByFilter(evcf, insFilter)
as.data.frame(fixed(evcf_fixedIns)[,c("REF", "ALT")])
