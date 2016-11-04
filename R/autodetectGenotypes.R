
setMethod(
    "autodetectGenotypes", c("VCF"),
    function(vcf){
        .autodetectGenotypes(vcf)
    }
)

# Main method ----

.autodetectGenotypes <- function(vcf){

    GT.all <- na.exclude(unique(c(geno(vcf)[["GT"]])))

    GT.OK <- grep("[[:digit:]][/|][[:digit:]]", GT.all, value=TRUE)

    GT.autoRef <- grep("(0/0)|(0\\|0)", GT.OK, value = TRUE)

    GT.split <- strsplit2(GT.OK, "[/|]")

    GT.autoHet <- NA_character_
    GT.autoAlt <- NA_character_
    if (ncol(GT.split) == 2){
        GT.autoHet <- GT.OK[GT.split[,1] != GT.split[,2]]
        GT.autoAlt <-  GT.OK[
            GT.split[,1] == GT.split[,2] &
                !GT.OK %in% GT.autoRef]
    }

    autoGenos <- Genotypes(GT.autoRef, GT.autoHet, GT.autoAlt)
    validObject(autoGenos)

    if ("TVTBparam" %in% names(metadata(vcf))){
        genos(metadata(vcf)[["TVTBparam"]]) <- autoGenos
        message("genos(metadata(vcf)[[\"TVTBparam\"]]) updated.")
    } else {
        metadata(vcf)[["TVTBparam"]] <- TVTBparam(autoGenos)
        message("metadata(vcf)[[\"TVTBparam\"]] created.")
    }

    return(vcf)
}
