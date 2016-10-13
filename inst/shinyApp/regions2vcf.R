
# Identifies VCF file(s) matching a sequence-based pattern in a folder
chr2file <- function(chr, pattern, folder){

    # Convert numeric names to character
    chr <- as.character(chr)

    # Validate input arguments
    stopifnot(length(chr) == 1)

    stopifnot(is.character(pattern))
    stopifnot(length(pattern) == 1)
    # Refuse pattern that do not contain "%s"
    stopifnot(grepl("%s", pattern))

    stopifnot(is.character(folder))
    stopifnot(dir.exists(folder))

    # Replace %s to create a pattern matching the sequence-specific VCF
    patternVcf <- sprintf(pattern, chr)

    # Full path to file(s) matching sequence-specific pattern
    matchedFiles <- file.path(
        folder,
        list.files(folder, patternVcf, ignore.case = TRUE)
    )

    # No file matched: empty vector
    if (length(matchedFiles) < 1){
        warning(sprintf("No file matching pattern: %s", patternVcf))
        return(character())
    } else if (length(matchedFiles) > 1){
        # More than one file
        warning(sprintf(
            "Multiple VCF files matching pattern: %s",
            patternVcf))
    }

    # If length(matchedFiles) == 1
    return(matchedFiles)
}
