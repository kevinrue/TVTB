
tryParsePheno <- function(file){
    message("Parsing ", file," as GRanges ...")

    rawData <- tryCatch(
        read.table(file = file, header = TRUE, row.names = 1),
        warning = function(warn){
            warning(warn)
            return(NULL)},
        error = function(err){
            warning(geterrmessage())
            return(NULL)
        }
    )

    return(rawData)
}

tryParseCsq <- function(vcf, vepKey){

    message("Parsing ", vepKey," to GRanges ...")

    rawData <- tryCatch(
        parseCSQToGRanges(
            x = vcf,
            VCFRowID = rownames(vcf),
            info.key = vepKey),
        warning = function(warn){
            warning(warn)
            return(NULL)},
        error = function(err){
            warning(geterrmessage())
            return(NULL)
        }
    )

    return(rawData)
}

tryParseBed <- function(bed){
    message("Parsing ", bed," as GRanges ...")

    rawData <- tryCatch(
        import.bed(con = bed),
        warning = function(warn){
            warning(warn)
            return(NULL)},
        error = function(err){
            warning(geterrmessage())
            return(NULL)
        }
    )

    return(rawData)
}

tryParseVcfHeader <- function(file){
    if (is.null(file))
        return(NULL)

    message("Parsing header of ", file, " ...")

    rawData <- tryCatch(
        scanVcfHeader(file = file),
        warning = function(warn){
            warning(warn)
            return(NULL)},
        error = function(err){
            warning(geterrmessage())
            return(NULL)
        }
    )

    return(rawData)
}

tryParseMultipleVcf <- function(
    folder, pattern, svp = ScanVcfParam(), yieldSize = NA_integer_,
    BPPARAM = BiocParallel::SerialParam()){

    rawData <- tryCatch(
        parseMultipleVcf(
            folder, pattern, svp = svp, yieldSize = yieldSize),
        warning = function(warn){
            warning(warn)
            return(NULL)},
        error = function(err){
            warning(geterrmessage())
            return(NULL)
        }
    )

    return(rawData)
}

tryParseSingleVcf <- function(
    file, svp = ScanVcfParam(), yieldSize = NA_integer_){

    rawData <- tryCatch(
        parseSingleVcf(
            file = file,
            svp = svp,
            yieldSize = yieldSize),
        warning = function(warn){
            warning(warn)
            return(NULL)},
        error = function(err){
            warning(geterrmessage())
            return(NULL)
        }
    )

    return(rawData)
}

parseMultipleVcf<- function(
    folder, pattern, svp = ScanVcfParam(),
    yieldSize = NA_integer_, BPPARAM = BiocParallel::SerialParam()){

    # Timing
    t1 <- Sys.time()

    # Extract gr
    gr <- vcfWhich(svp)

    # If there is at least a targeted region
    if (length(gr) > 0){

        # Identify the targeted chromosomes
        chrs <- unique(names(unlist(gr)))

        # Identify the corresponding files
        vcfFiles <- sapply(
            X = chrs,
            FUN = TVTB::chr2file,
            pattern = pattern,
            folder = folder)

        if (length(chrs) != length(vcfFiles))
            stop(sprintf(
                "length(chrs) != length(vcfFiles): %i, %i",
                length(chrs),
                length(vcfFiles)))

        if (length(vcfFiles) == 1){
            # All genomic ranges are on the same chromosome
            vcf <- parseSingleVcf(
                file = vcfFiles,
                svp = svp)
            # Timing
            t2 <- Sys.time()
            dt <- t2 - t1
            message(sprintf(
                "%i variants from %i region(s) imported in %.2f %s",
                length(vcf),
                length(vcfWhich(svp)),
                as.numeric(dt),
                units(dt)))
            return(vcf)
        }

    } else {
        # If no genomic range is given
        # Identify all VCF files matching the pattern (%s substituted to .*)
        vcfFiles <- file.path(
            folder,
            grep(
                pattern = gsub('%s', '.*', pattern),
                x = list.files(folder),
                ignore.case = TRUE,
                value = TRUE))
    }

    # The functions gets here in two situations:
    ## 1. genomic ranges span multiple chromosomes/vcfFiles
    ## 2. no genomic range was given: all vcfFiles in the folder are imported
    message("Parsing ", length(vcfFiles), " VCF files ...")
    # Import from each file in parallel (already ExpandedVCF)
    vcfs <- bplapply(
        X = vcfFiles,
        FUN = parseSingleVcf,
        svp = svp,
        yieldSize = yieldSize,
        BPPARAM = BPPARAM
    )

    message("Combining variants from multiple VCF files...")
    # Combine the imported objects
    if (length(vcfs) > 1)
        vcf <- do.call(what = rbind, args = vcfs)
    else
        vcf <- vcfs[[1]]

    # Timing
    t2 <- Sys.time()
    dt <- t2 - t1
    message(sprintf(
        "%i variants from %i region(s) in %i file(s) imported in %.2f %s",
        length(vcf),
        length(vcfWhich(svp)),
        length(vcfFiles),
        as.numeric(dt),
        units(dt)))

    return(vcf)
}

parseSingleVcf <- function(
    file, svp = ScanVcfParam(), yieldSize = NA_integer_){
    message("Parsing VCF file: ", file, " ...")
    tf <- TabixFile(file)

    # scanTabix: 'yieldSize(file)' must be 'NA_integer_' when range specified
    if (length(vcfWhich(svp)) > 0)
        yieldSize(tf) <- NA_integer_

    # Import variants
    vcf <- VariantAnnotation::readVcf(file = tf, param = svp)

    # Expand variants into bi-allelic records
    vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

    return(vcf)
}

