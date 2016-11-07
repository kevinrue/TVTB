
tryParseCsq <- function(vcf, vepKey){

  message("Parsing ", vepKey," to GRanges ...")

  rawData <- tryCatch(
    ensemblVEP::parseCSQToGRanges(
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

parseMultipleVcf<- function(
  folder, pattern, param, colData, autodetectGT, yieldSize, BPPARAM){

  # Timing
  t1 <- Sys.time()

  # Extract gr
  gr <- VariantAnnotation::vcfWhich(TVTB::svp(param))

  # If there is at least a targeted region
  if (length(gr) > 0){

    # Identify the targeted chromosomes
    chrs <- unique(names(BiocGenerics::unlist(gr)))

    # Identify the corresponding files
    vcfFiles <- sapply(
      chrs, chr2file,
      pattern = pattern,
      folder = folder)

    if (length(chrs) != length(vcfFiles))
      stop(sprintf(
        "length(chrs) != length(vcfFiles): %i, %i",
        length(chrs),
        length(vcfFiles)))

    if (length(vcfFiles) == 1){
      # All genomic ranges are on the same chromosome
      vcf <- parseSingleVcf(vcfFiles, param)
      # Timing
      t2 <- Sys.time()
      dt <- t2 - t1
      message(sprintf(
        "%i variants from %i region(s) imported in %.2f %s",
        length(vcf),
        length(VariantAnnotation::vcfWhich(param)),
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
  message("Parsing ", length(vcfFiles), " VCF file(s) ...")
  # Import from each file in parallel (already ExpandedVCF)
  vcfs <- BiocParallel::bplapply(
    vcfFiles, parseSingleVcf,
    param = param,
    colData = colData,
    autodetectGT = autodetectGT,
    yieldSize = yieldSize,
    BPPARAM = BPPARAM
  )

  # Combine the imported objects
  if (length(vcfs) > 1){
    message("Combining variants from multiple VCF files...")
    vcf <- do.call(what = rbind, args = vcfs)
  }
  else
    vcf <- vcfs[[1]]

  # Timing
  t2 <- Sys.time()
  dt <- t2 - t1
  message(sprintf(
    "%i variants from %i region(s) in %i file(s) imported in %.2f %s",
    length(vcf),
    length(VariantAnnotation::vcfWhich(TVTB::svp(param))),
    length(vcfFiles),
    as.numeric(dt),
    units(dt)))

  return(vcf)
}

parseSingleVcf <- function(
  file, param, colData, autodetectGT, yieldSize){
  message("Parsing VCF file: ", file, " ...")
  tf <- Rsamtools::TabixFile(file)

  # scanTabix: 'yieldSize(file)' must be 'NA_integer_' when range specified
  if (length(VariantAnnotation::vcfWhich(TVTB::svp(param))) > 0){
    Rsamtools::yieldSize(tf) <- NA_integer_
  }

  # Default value for the readVcf method
  if (is.null(colData)){
    vh <- VariantAnnotation::scanVcfHeader(file)
    colData <- S4Vectors::DataFrame(
      row.names = VariantAnnotation::samples(vh)
    )
  }

  # Import variants
  vcf <- TVTB::readVcf(
    tf, param = param, colData = colData, autodetectGT = autodetectGT)

  # Default value for the readVcf method
  if (is.null(colData)){
    colData <- S4Vectors::DataFrame(row.names = colnames(vcf))
  }

  # Expand variants into bi-allelic records
  vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

  return(vcf)
}
