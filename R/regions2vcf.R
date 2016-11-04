
# Given a chr(omosome) name [character]
# The pattern of VCF files where %s stands for the above chr [character]
# and the folder where VCF files are located
# Return the path to the VCF file corresponding to the chr
chr2file <- function(chr, pattern, folder){

    stopifnot(is.character(chr))
    stopifnot(dir.exists(folder))

    if(!grepl(pattern = "%s", x = pattern)){
        stop("\"%s\" not found in pattern")
    }

    vcfFiles <- list.files(
        path = folder,
        pattern = gsub('%s', '.*', pattern),
        ignore.case = TRUE)

    if(length(vcfFiles) == 0){
        warning("No file matching pattern in folder")
        return(character())
    }

    # Replace %s to create a pattern matching the chromosome-specific VCF
    patternVcf <- sprintf(pattern, chr)

    vcfPath <- file.path(
        folder,
        grep(
            pattern = patternVcf,
            x = vcfFiles,
            value = TRUE))

    if(length(vcfPath) < 1){
        warning(sprintf("No VCF file matching pattern: %s", patternVcf))
        return(character())
    }else if(length(vcfPath) > 1){
        stop(sprintf(
            "More than one VCF file matching pattern: %s",
            patternVcf))
    }

    return(vcfPath)
}
