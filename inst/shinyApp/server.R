
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

source(system.file("shinyApp", "serverRoutines.R", package = "TVTB"))

shinyServer(function(input, output, clientData, session) {

  # Reactive values ----

  Errors <- reactiveValues(
    phenotypes = NULL, # DataFrame

    EnsDbPkg = "No EnsDb package installed.", # annotation object

    BED = NULL, # GRanges
    UCSC = NULL,
    EnsDb = NULL,
    GRanges = NULL,

    VCFfolder = NULL,
    VCFfiles = NULL,
    singleVCF = .Msgs[["singleVcf"]],
    VCFheader = NULL,
    VCF = "Please import variants."
  )

  RV <- reactiveValues(

    # Phenotypes
    phenoFile = NULL, # file path (character)
    phenotypes = NULL, # DataFrame

    # VCF import
    infoKeys = NULL, # choices (character)
    genoKeys = NULL, # choices (character)

    # VCF filters
    vcfFilters = TVTB::VcfFilterRules(),
    newFilterTestResults = FALSE,
    newFilterStatus = "",

    # EnsDb
    EnsDbPkg = NULL,

    # Genomic ranges imported/analysed
    BED = GenomicRanges::GRanges(),
    UCSC = GenomicRanges::GRanges(),
    EnsDb = GenomicRanges::GRanges(),
    genomicRanges = GenomicRanges::GRanges(),

    # Single VCF file
    singleVcf = NULL, # file path (character),

    # VCF objects
    vcf = NULL,

    # TVTBparam
    GT.all = c(get("refGT", .tSVE), get("hetGT", .tSVE), get("altGT", .tSVE)),
    GT.autoRef = character(),
    GT.autoHet = character(),
    GT.autoAlt = character(),

    # Parallel settings
    parallel = BiocParallel::SerialParam()
  )

  output$Errors <- renderPrint({
    return(reactiveValuesToList(Errors))
  })

  # Import phenotype information ----
  observeEvent(
    input$selectPheno,
    {
      RV[["phenoFile"]] <- tryCatch(
        file.choose(),
        error = function(err){
          warning(geterrmessage())
          return(NULL)
        }
      )
    }
  )

  observeEvent(
    input$demoPheno,
    {
      RV[["phenoFile"]] <- system.file(
        "extdata", "integrated_samples.txt", package = "TVTB"
      )
    }
  )

  output$phenoFile <- renderText({
    phenoFile <- RV[["phenoFile"]]

    return(ifelse(is.null(phenoFile), "No file provided.", phenoFile))
  })

  # DataFrame of imported phenotypes, or NULL
  observeEvent(
    RV[["phenoFile"]],
    {
      # Start with optimism :)
      Errors[["phenotypes"]] <- NULL

      phenoFile <- RV[["phenoFile"]]

      if (is.null(phenoFile)){
        RV[["phenotypes"]] <- NULL
        return()
      }

      message("Importing phenotypes ...")
      rawData <- tryCatch(
        {
          S4Vectors::DataFrame(read.table(phenoFile, TRUE, row.names = 1))
        },
        warning = function(warn){
          Errors[["phenotypes"]] <- sprintf("Failed to parse file:\n%s", warn)
          NULL
        },
        error = function(err){
          Errors[["phenotypes"]] <- sprintf(
            "Failed to parse file:\n%s", geterrmessage()
          )
          NULL
        }
      )

      RV[["phenotypes"]] <- rawData

      if (is.null(rawData)){
        return()
      }

      if (nrow(rawData) == 0){
        Errors[["phenotypes"]] <- paste0(
          "Invalid phenotypes:\n",
          "No sample (row) detected in phenotype file."
        )
        warning(Errors[["phenotypes"]])
        return()
      }

      if (ncol(rawData) == 0){
        Errors[["phenotypes"]] <- paste0(
          "Invalid phenotypes:\n",
          "No phenotype (column) detected in phenotype file."
        )
        warning(Errors[["phenotypes"]])
        return()
      }
    },
    ignoreNULL = FALSE # update if file de-selected
  )

  # HTML summary of imported phenotypes
  output$phenoFileSummary <- renderUI({
    # Phenotype not supplied
    if (is.null(RV[["phenoFile"]])){
      return(.Msgs[["phenoFile"]])
    }

    # Invalid phenotype
    validate(need(is.null(Errors[["phenotypes"]]), Errors[["phenotypes"]]))

    phenotypes <- RV[["phenotypes"]]

    return(tagList(
      code(S4Vectors::ncol(phenotypes)),
      "phenotypes in",
      code(S4Vectors::nrow(phenotypes)),
      "samples."
    ))
  })

  # Interactive data table of raw phenotypes
  # output$phenoFileSample <- DT::renderDataTable({
  #
  #   phenotypes <- RV[["phenotypes"]]
  #
  #   validate(need(phenotypes, "No data to show."), errorClass = "optional")
  #
  #   return(
  #     DT::datatable(
  #       as.data.frame(phenotypes),
  #       options = list(
  #         pageLength = 10,
  #         searching = TRUE),
  #       filter = "top"
  #     )
  #   )
  # })

  # Column names available for selection from raw phenotypes
  output$phenoCols <- renderUI({
    # NOTE:
    # * Populate choices from the filtered phenotypes (those used)
    # * Users can easily look at their raw phenotypes separately

    validate(need(is.null(Errors[["VCF"]]), Errors[["VCF"]]))

    vcf <- RV[["filteredVcf"]]

    phenos <- SummarizedExperiment::colData(vcf)

    validate(need(
      ncol(phenos) > 0,
      ifelse(
        is.null(RV[["phenotypes"]]),
        .Msgs[["colDataEmptyNoFile"]],
        .Msgs[["colDataEmptyImported"]])
    ),
    errorClass = "optional")

    selectInput(
      "phenoCols", "Phenotypes",
      choices = colnames(phenos),
      selected = colnames(phenos)[1:5],
      multiple = TRUE
    )
  })

  # Interactive table of phenotypes attached to variants
  output$phenotypesView <- DT::renderDataTable({
    # Widget displays the error message already
    validate(need(is.null(Errors[["VCF"]]), ""))

    # Display phenotypes attached to the VCF object
    vcf <- RV[["filteredVcf"]]

    # Only display after phenotypes were attached to VCF
    validate(
      need(vcf, .Msgs[["colDataEmptyImported"]]),
      errorClass = "optional"
    )

    # Make sure the selected fields are present in data
    phenos <- SummarizedExperiment::colData(vcf)
    validate(need(
      all(input$phenoCols %in% colnames(phenos)),
      "Invalid phenotypes selected."))

    return(DT::datatable(
      as.data.frame(phenos[,input$phenoCols, drop = FALSE]),
      options = list(
        pageLength = 10,
        searching = TRUE),
      filter = "top"))
  })

  # Update choices of phenotypes when new VCF
  observeEvent(
    RV[["VCF"]],
    {
      # raw VCF
      vcf <- RV[["VCF"]]

      # phenotypes
      phenos <- SummarizedExperiment::colData(vcf)

      updateSelectInput(
        session, "phenoAddFrequencies",
        choices = colnames(phenos)
      )
    }
  )

  # Calculation of frequencies ----

  # For each phenotype
  # Create a group of checkboxes for each level
  # to select those for which frequencies should be calculated
  observeEvent(
    input$phenoAddFrequencies,
    {
      # raw VCF
      vcf <- RV[["VCF"]]

      # Without VCF, nothing to update
      req(vcf)

      # phenotypes
      phenos <- SummarizedExperiment::colData(vcf)
      phenoNames <- colnames(phenos)
      phenoLevels <- levels(phenos[,input$phenoAddFrequencies])

      if (length(phenoLevels) == 0){
        return()
      }

      ## pre-tick phenoLevels already calculated
      infoCols <- colnames(VariantAnnotation::info(vcf))
      # phenoLevels already calculated have all suffixes present
      suffixes <- TVTB::suffix(RV[["TVTBparam"]])
      levelsPresent <- sapply(
        X = phenoLevels,
        FUN = function(phenoLevel){
          all(
            paste(
              input$phenoAddFrequencies,
              phenoLevel,
              suffixes,
              sep = "_"
            ) %in% infoCols
          )
        }
      )

      # Update the checboxes
      updateCheckboxGroupInput(
        session, "phenoLevelFreqCheckboxes",
        label = paste0(
          "Levels in phenotype \"",
          input$phenoAddFrequencies,
          "\""),
        choices = phenoLevels,
        selected = names(which(levelsPresent)),
        inline = TRUE
      )
    }
  )

  # Select all phenotype levels on button click
  observeEvent(
    input$tickAllPhenoLevelsFreq,
    {
      # raw VCF
      vcf <- RV[["VCF"]]
      req(vcf)

      # phenotypes
      phenos <- SummarizedExperiment::colData(vcf)
      phenoNames <- colnames(phenos)
      choices <- levels(phenos[,input$phenoAddFrequencies])

      updateCheckboxGroupInput(
        session, "phenoLevelFreqCheckboxes",
        choices = choices,
        selected = choices,
        inline = TRUE
      )
    }
  )

  # Deselect all phenotype levels on button click
  observeEvent(
    input$untickAllPhenoLevelsFreq,
    {
      updateCheckboxGroupInput(
        session, "phenoLevelFreqCheckboxes",
        selected = character(),
        inline = TRUE
      )
    }
  )

  # Add overall frequencies on button click
  observeEvent(
    input$addOverallFrequencies,
    {

      vcf <- RV[["VCF"]]
      req(vcf)

      # collect all INFO keys
      suffixes <- TVTB::suffix(RV[["TVTBparam"]])

      withProgress(
        max = 3,
        value = 1,
        message = "Progress",
        detail = .Tracking[["preprocessing"]],{
          # Only proceed if none of the INFO keys are present
          if (!any(
            suffixes %in%
            colnames(VariantAnnotation::info(vcf)))){

            incProgress(
              amount = 1,
              detail = .Tracking[["addFreqOverall"]])

            message("Adding overall frequencies ...")
            RV[["VCF"]] <- TVTB::addOverallFrequencies(
              vcf = vcf,
              force = TRUE # watch the console for warnings
            )
            RV[["latestPhenotypeFrequency"]] <- "Overall"
            RV[["latestFrequenciesAdded"]] <- suffixes
            RV[["latestFrequenciesRemoved"]] <- character()
          }
        }
      )
    }
  )

  # Remove overall frequencies on button click
  observeEvent(
    input$removeOverallFrequencies,
    {

      vcf <- RV[["VCF"]]
      req(vcf)

      # collect all INFO keys
      suffixes <- TVTB::suffix(RV[["TVTBparam"]])

      withProgress(
        max = 3,
        value = 1,
        message = "Progress",
        detail = .Tracking[["preprocessing"]],
        {
          incProgress(
            amount = 1,
            detail = .Tracking[["rmFreqOverall"]])

          # Only proceed if all of the INFO keys are present
          if (all(
            suffixes %in% colnames(VariantAnnotation::info(vcf))
          )){

            message("Removing overall frequencies ...")
            RV[["VCF"]] <- TVTB::dropInfo(
              vcf = vcf,
              key = suffixes,
              slot = "both"
            )
            RV[["latestPhenotypeFrequency"]] <- "Overall"
            RV[["latestFrequenciesRemoved"]] <- suffixes
            RV[["latestFrequenciesAdded"]] <- character()
          }
        }
      )
    }
  )

  # Add selected & remove deselected frequencies on button click
  observeEvent(
    input$buttonFrequencies,
    {

      withProgress(
        max = 4,
        value = 1,
        message = "Progress",
        detail = .Tracking[["preprocessing"]],{
          # Remove INFO keys for unticked boxes
          # Add values for ticked boxes

          # Phenotypes selected
          selectedPhenoName <- input$phenoAddFrequencies
          selectedPhenoLevels <- input$phenoLevelFreqCheckboxes

          # Info in VCF
          vcf <- RV[["VCF"]]
          req(vcf)
          vcfInfoCols <- colnames(VariantAnnotation::info(vcf))

          # Phenotype levels available
          phenos <- SummarizedExperiment::colData(vcf)
          choices <- levels(phenos[,input$phenoAddFrequencies])

          # pre1) collect all suffixes (to check existence of fields)
          suffixes <- TVTB::suffix(RV[["TVTBparam"]])

          # pre2) identify unticked phenoLevels

          phenoLevelsUnticked <- choices[which(
            !choices %in% selectedPhenoLevels
          )]

          if (length(phenoLevelsUnticked) > 0){

            incProgress(
              amount = 1,
              detail = .Tracking[["rmFreqPhenoLevel"]])

            # Identify unticked phenoLevels present in vcf INFO
            # (to remove)
            # 1) Generate all expected key names for phenoLevels
            untickedPhenoLevelKeys <- sapply(
              X = phenoLevelsUnticked,
              FUN = function(phenoLevel){
                return(paste(
                  input$phenoAddFrequencies,
                  phenoLevel,
                  suffixes,
                  sep = "_"
                ))
              },
              simplify = FALSE
            ) # named list: names=phenoLevel, value=keys

            # 4) identify phenoLevels with all keys present
            # (in data, not header)
            areAllPhenoLevelKeysPresent <- sapply(
              X = untickedPhenoLevelKeys,
              FUN = function(keys){
                return(all(keys %in% vcfInfoCols))
              }
            ) # named boolean vector: names=phenoLevel

            phenoLevelKeysRemove <- do.call(
              c,
              untickedPhenoLevelKeys[
                which(areAllPhenoLevelKeysPresent)]
            ) # vector of keys to remove

            if (length(phenoLevelKeysRemove))
              message(
                "Removing INFO keys for levels: ",
                paste(phenoLevelsUnticked, collapse = ", "))
            RV[["VCF"]] <- TVTB::dropInfo(
              vcf = RV[["VCF"]],
              key = phenoLevelKeysRemove,
              slot = "both")

            RV[["latestFrequenciesRemoved"]] <-
              names(areAllPhenoLevelKeysPresent)[
                which(areAllPhenoLevelKeysPresent)
                ]
          } else {
            RV[["latestFrequenciesRemoved"]] <- character()
          }

          if (length(selectedPhenoLevels) > 0){

            incProgress(
              amount = 1,
              detail = .Tracking[["addFreqPhenoLevel"]])

            # Identify ticked phenoLevels absent in info slot
            # (to calculate)
            # 1) generate all keys expected for each phenoLevel
            tickedPhenoLevelKeys <- sapply(
              X = selectedPhenoLevels,
              FUN = function(phenoLevel){
                return(paste(
                  input$phenoAddFrequencies,
                  phenoLevel,
                  suffixes,
                  sep = "_"
                ))
              },
              simplify = FALSE
            ) # list: names=phenoLevels, value=keys

            # 2) identify phenoLevels with all keys absent
            # (in data, not header)
            areAllPhenoLevelKeysAbsent <- sapply(
              X = tickedPhenoLevelKeys,
              FUN = function(keys){
                return(!any(keys %in% vcfInfoCols))
              }
            ) # named boolean vector: names=phenoLevels

            phenoLevelsAdd <- names(
              areAllPhenoLevelKeysAbsent
            )[which(areAllPhenoLevelKeysAbsent)
              ] # phenoLevels to add

            # Format for addFrequencies input
            phenosAdd <- list(phenoLevelsAdd)
            names(phenosAdd)[[1]] <- selectedPhenoName

            if (length(phenoLevelsAdd) > 0)
              message(
                "Adding INFO keys for levels: ",
                paste(phenoLevelsAdd, collapse = ", "))
            RV[["VCF"]] <- TVTB::addFrequencies(
              vcf = vcf,
              phenos = phenosAdd,
              force = TRUE # watch console for warnings
            )

            RV[["latestFrequenciesAdded"]] <- phenosAdd[[1]]
          } else {
            RV[["latestFrequenciesAdded"]] <- character()
          }

          RV[["latestPhenotypeFrequency"]] <- selectedPhenoName
        }
      )
    }
  )

  # Display the frequencies added and removed in the latest update -
  output$latestFrequenciesCalculated <- renderUI({

    pheno <- RV[["latestPhenotypeFrequency"]]
    freqAdded <- RV[["latestFrequenciesAdded"]]
    freqRemoved <- RV[["latestFrequenciesRemoved"]]

    # If no calculations were done yet
    if (all(is.null(freqAdded), is.null(freqRemoved)))
      return("No frequencies calculated yet.")

    # If phenotype level frequencies were last calculated
    if (pheno != "Overall"){
      if (length(freqAdded) > 0)
        addedLevels <- paste(freqAdded, collapse = ", ")
      else
        addedLevels <- "<NA>"

      if(length(freqRemoved > 0))
        removedLevels <- paste(freqRemoved, collapse = ", ")
      else
        removedLevels <- "<NA>"

      return(tagList(
        "Phenotype", code(pheno), ":",
        tags$ul(
          tags$li("Added level(s):", code(addedLevels)),
          tags$li("Removed level(s):", code(removedLevels))
        )
      ))
    }

    if (identical(pheno, "Overall")){
      if (length(freqAdded) > 0)
        return("Overall frequencies added")
      if (length(freqRemoved) > 0)
        return("Overall frequencies removed")
    }

    return("If you see this message, please notify the maintainer !")
  })

  # Define genomic ranges ----

  # Path to BED file
  observeEvent(
    input$selectBed,
    {
      RV[["bedFile"]] <- tryCatch(
        file.choose(),
        error = function(err){
          warning(geterrmessage())
          return(NULL)
        }
      )
    }
  )

  # Demonstration input for the BED file
  observeEvent(
    input$demoBed,
    {
      RV[["bedFile"]] <- system.file(
        "extdata/SLC24A5.bed", package = "TVTB"
      )
    }
  )

  # Demonstration input for the UCSC field input
  observeEvent(
    input$demoUCSC,
    {
      updateTextInput(
        session, "ucscRanges",
        value = "15:48,413,169-48,434,869"
      )
    }
  )

  # Demonstration input for the EnsDb field input
  observeEvent(
    input$demoEnsDb,
    {
      updateTextInput(
        session, "ensDb.type",
        value = "Genename"
      )
      updateTextInput(
        session, "ensDb.condition",
        value = "="
      )
      updateTextInput(
        session, "ensDb.value",
        value = "SLC24A5"
      )
    }
  )

  # Display path to the selected BED file
  output$bedFile <- renderText({
    bedFile <- RV[["bedFile"]]

    return(ifelse(
      is.null(bedFile),
      "No file selected.",
      bedFile))
  })

  # If a BED file is de/selected, update GRangesBED
  observeEvent(
    RV[["bedFile"]],
    {
      # Start with optimism :)
      Errors[["GRanges"]] <- NULL

      # use rtracklayer::import.bed to obtain GRanges
      bedFile <- RV[["bedFile"]]

      if (is.null(bedFile)){
        RV[["BED"]] <- GenomicRanges::GRanges()
        return()
      }

      rawData <- tryCatch(
        {
          rtracklayer::import.bed(bedFile)
        }
        ,
        warning = function(warn){
          Errors[["GRanges"]] <- sprintf("import.bed failed:\n%s", warn)
          NULL
        },
        error = function(err){
          Errors[["GRanges"]] <- sprintf(
            "import.bed failed:\n%s", geterrmessage()
          )
          NULL
        }
      )

      RV[["BED"]] <- rawData
    },
    ignoreNULL = FALSE
  )

  # If UCSC text field is updated, update GRangesUCSC
  observeEvent(
    input$ucscRanges,
    {
      # Start with optimism :)
      Errors[["GRanges"]] <- NULL

      # parse the string or return NULL
      if (input$ucscRanges == ""){
        RV[["UCSC"]] <- GenomicRanges::GRanges()
        return()
      }

      # NOTE: do not trim "chr", for future UCSC support
      inputTrimmed <- gsub("[,| ]", "", input$ucscRanges)

      # Split the given string into individual regions
      inputSplit <- strsplit(inputTrimmed, ";")[[1]]

      # Ensure that all regions are UCSC-valid
      validityCheck <- sapply(
        inputSplit,
        function(x){grepl("^[[:alnum:]]+:[[:digit:]]+-[[:digit:]]+$", x)}
      )

      if (!all(validityCheck)){
        Errors[["GRanges"]] <- "Invalid UCSC text input."
        return()
      }

      rawData <- tryCatch(
        {
          chrs <- as.numeric(gsub(
            "([[:alnum:]]*):[[:digit:]]*-[[:digit:]]*",
            "\\1",
            inputSplit
          ))

          starts <- as.numeric(gsub(
            "[[:alnum:]]*:([[:digit:]]*)-[[:digit:]]*",
            "\\1",
            inputSplit
          ))

          ends <- as.numeric(gsub(
            "[[:alnum:]]*:[[:digit:]]*-([[:digit:]]*)",
            "\\1",
            inputSplit
          ))

          data.frame(
            V1 = chrs,
            V2 = starts,
            V3 = ends,
            stringsAsFactors = FALSE
          )
        },
        error = function(err){
          Errors[["GRanges"]] <- sprintf(
            "Invalid UCSC input:\n%s",
            geterrmessage()
          )
          NULL
        },
        warning = function(warn){
          Errors[["GRanges"]] <- sprintf(
            "Invalid UCSC input:\n%s",
            warn
          )
          NULL
        }
      )

      rawData <- tryCatch(
        {
          validateDataFrameGRanges(rawData, RV[["EnsDbPkg"]])
        },
        error = function(err){
          Errors[["GRanges"]] <- sprintf(
            "Invalid UCSC input:\n%s",
            geterrmessage()
          )
          NULL
        },
        warning = function(warn){
          sprintf(
            Errors[["GRanges"]] <- "Invalid UCSC input:\n%s",
            warn
          )
          NULL
        }
      )

      RV[["UCSC"]] <- rawData
    },
    ignoreNULL = FALSE
  )

  # Update EnsDb queried GRanges
  observeEvent(
    queryGenes(),
    {
      # Start with optimism :)
      Errors[["GRanges"]] <- NULL

      queryGenes <- queryGenes()

      RV[["EnsDb"]] <- queryGenes

      if (nchar(input$ensDb.value) > 0 & length(queryGenes) == 0){
        Errors[["GRanges"]] <- "EnsDb query returned empty GRanges."
        return()
      }
    }
  )

  # Update active GRanges from BED GRanges
  observeEvent(
    RV[["BED"]],
    {
      RV[["activeGRanges"]] <- RV[["BED"]]
    },
    ignoreNULL = FALSE
  )

  # Update active GRanges from UCSC GRanges
  observeEvent(
    RV[["UCSC"]],
    {
      RV[["activeGRanges"]] <- RV[["UCSC"]]
    },
    ignoreNULL = FALSE
  )

  # Update active GRanges from EnsDb GRanges
  observeEvent(
    RV[["EnsDb"]],
    {
      RV[["activeGRanges"]] <- RV[["EnsDb"]]
    },
    ignoreNULL = FALSE
  )

  # Update active GRanges based on input mode
  observeEvent(
    input$grangesInputMode,
    {
      switch (
        input$grangesInputMode,
        bed = {
          RV[["activeGRanges"]] <- RV[["BED"]]
        },
        ucsc = {
          RV[["activeGRanges"]] <- RV[["UCSC"]]
        },
        EnsDb = {
          RV[["activeGRanges"]] <- RV[["EnsDb"]]
        }
      )
    }
  )

  # How many BED records detected, show first one.
  output$rangesSummary <- renderUI({
    # Depends on genomicRanges

    genomicRanges <- RV[["activeGRanges"]]

    validate(need(is.null(Errors[["GRanges"]]), Errors[["GRanges"]]))

    if (length(genomicRanges) == 0) return(.Msgs[["noGenomicRanges"]])

    return(tagList(
      code(length(genomicRanges)), "genomic range(s)", br(),
      "[",
      as.character(head(
        x = GenomeInfoDb::seqnames(genomicRanges),
        n = 1)),
      ":",
      as.character(head(
        x = BiocGenerics::start(genomicRanges),
        n = 1)),
      "-",
      as.character(head(
        x = BiocGenerics::end(genomicRanges),
        n = 1)),
      # "{",
      # as.character(head(x = names(genomicRanges), n = 1)),
      # "}",
      " , ... ]"
    ))
  })

  # Show BED records
  # output$rangesSession <- renderPrint({
  #
  #   genomicRanges <- RV[["activeGRanges"]]
  #
  #   validate(need(
  #     !is.null(genomicRanges),
  #     .Msgs[["noGenomicRanges"]]),
  #     errorClass = "optional")
  #
  #   return(genomicRanges)
  # })

  output$rangesTableView <- DT::renderDataTable({

    genomicRanges <- RV[["activeGRanges"]]

    validate(need(is.null(Errors[["GRanges"]]), Errors[["GRanges"]]))

    validate(need(
      length(genomicRanges) > 0,
      .Msgs[["noGenomicRanges"]]),
      errorClass = "optional")

    return(DT::datatable(
      data = as.data.frame(genomicRanges),
      options = list(
        pageLength = 10,
        searching = TRUE),
      filter = "top"))
  })

  # Genome annotation ----

  observeEvent(
    input$annotationPackage,
    {
      # Start with optimism :)
      Errors[["EnsDbPkg"]] <- NULL

      ap <- input$annotationPackage

      if (identical(ap, "")){
        RV[["EnsDbPkg"]] <- NULL
        Errors[["EnsDbPkg"]] <-
          "An annotation package must be installed."
        return()
      }

      RV[["EnsDbPkg"]] <- tryCatch(
        {
          get(ap, envir = asNamespace(ap))
        },
        error = function(err){
          Errors[["EnsDbPkg"]] <- geterrmessage()
          NULL
        },
        warning = function(warn){
          Errors[["EnsDbPkg"]] <- warn
          NULL
        }
      )
    },
    ignoreNULL = FALSE
  )

  # EnsDb ----

  output$ensembl_organism <- renderUI({
    # Depends on selectedEnsDb
    edb <- RV[["EnsDbPkg"]]

    validate(need(is.null(Errors[["EnsDbPkg"]]), Errors[["EnsDbPkg"]]))

    HTML(paste(
      tags$strong("Organism:"),
      ensembldb::organism(edb))
    )

  })

  output$ensembl_version <- renderUI({
    # Depends on selectedEnsDb
    edb <- RV[["EnsDbPkg"]]

    validate(need(is.null(Errors[["EnsDbPkg"]]), Errors[["EnsDbPkg"]]))

    md <- ensembldb::metadata(edb)
    rownames(md) <- md$name

    HTML(paste(
      tags$strong("Ensembl version:"),
      md["ensembl_version", "value"])
    )
  })

  output$ensembl_genome <- renderUI({
    # Depends on selectedEnsDb
    edb <- RV[["EnsDbPkg"]]

    validate(need(is.null(Errors[["EnsDbPkg"]]), Errors[["EnsDbPkg"]]))

    md <- ensembldb::metadata(edb)
    rownames(md) <- md$name

    HTML(paste(
      tags$strong("Genome build:"),
      md["genome_build", "value"])
    )
  })

  # GRanges for the gene query ()
  queryGenes <- reactive({
    # Depends on: selectedEnsDb, ...
    #   input$ensDb.type, input$ensDb.condition, input$ensDb.value

    edb <- RV[["EnsDbPkg"]]
    validate(need(is.null(Errors[["EnsDbPkg"]]), Errors[["EnsDbPkg"]]))

    if (input$ensDb.value == ""){
      return(GenomicRanges::GRanges())
    }

    ensDbFilter <- EnsDbFilter(
      type = input$ensDb.type,
      condition = input$ensDb.condition,
      value = input$ensDb.value)

    res <- ensembldb::genes(edb, filter = ensDbFilter)

    return(res)
  })

  output$ensDb.Genes <- DT::renderDataTable({

    queryGenes <- queryGenes()

    validate(
      need(length(queryGenes) > 0, "No results."),
      errorClass = "optional"
    )

    return(DT::datatable(
      data = as.data.frame(queryGenes),
      options = list(
        pageLength = 10,
        searching = TRUE),
      filter = "top"
    ))
  })

  # output$ensDb.Transcripts <- renderDataTable({
  #   if(length(input$package) == 0)
  #     return
  #   if(!is.na(input$geneName) &
  #      length(input$geneName) > 0 &
  #      input$geneName!=""){
  #     edb <- RV[["EnsDbPkg"]]
  #     res <- transcripts(
  #       edb, filter = EnsDbFilter(input),
  #       return.type = "data.frame")
  #     # assign(".ENS_TMP_RES", res, envir = globalenv())
  #     return(res)
  #   }
  # })

  # output$ensDb.Exons <- renderDataTable({
  #   if(length(input$package) == 0)
  #     return
  #   if(!is.na(input$geneName) &
  #      length(input$geneName) > 0 &
  #      input$geneName!=""){
  #     edb <- RV[["EnsDbPkg"]]
  #     res <- exons(edb, filter=EnsDbFilter(input),
  #                  return.type="data.frame")
  #     assign(".ENS_TMP_RES", res, envir=globalenv())
  #     return(res)
  #   }
  # })

  # TVTBparam ----

  # TVTBparam object: depends on many parameters
  observe({
    # Start with optimisim :)
    Errors[["TVTBparam"]] <- NULL

    if (length(input$refSuffix) == 0){
      Errors[["TVTBparam"]] <- "Suffix for REF allele data must be defined."
    }
    else if (length(input$hetSuffix) == 0){
      Errors[["TVTBparam"]] <- "Suffix for HET allele data must be defined."
    }
    else if (length(input$altSuffix) == 0){
      Errors[["TVTBparam"]] <- "Suffix for ALT allele data must be defined."
    }
    else if (length(input$aafSuffix) == 0){
      Errors[["TVTBparam"]]<-"Suffix for ALT allele frequency must be defined."
    }
    else if (length(input$mafSuffix) == 0){
      Errors[["TVTBparam"]]<-"Suffix for minor allele frequency must be defined."
    }
    else if (length(input$vepKey) == 0){
      Errors[["TVTBparam"]] <- "INFO field of VEP prediction must be defined."
    }
    else if (!is.null(Errors[["GRanges"]])){
      Errors[["TVTBparam"]] <- Errors[["GRanges"]]
    }
    else if (!is.null(Errors[["phenotypes"]])){
      Errors[["TVTBparam"]] <- Errors[["phenotypes"]]
    }
    # If any of the above failed, clear TVTBparam and stop here
    if (!is.null(Errors[["TVTBparam"]])){
      RV[["TVTBparam"]] <- NULL
      return()
    }

    tparam <- TVTB::TVTBparam(
      genos = TVTB::Genotypes(
        ref = as.character(input$refGenotypes),
        het = as.character(input$hetGenotypes),
        alt = as.character(input$altGenotypes),
        suffix = c(
          ref = input$refSuffix,
          het = input$hetSuffix,
          alt = input$altSuffix
        )
      ),
      aaf = input$aafSuffix,
      maf = input$mafSuffix,
      vep = input$vepKey,
      bp = RV[["parallel"]],
      svp = VariantAnnotation::ScanVcfParam(
        info = c(input$vepKey, input$vcfInfoKeys),
        geno = c("GT", input$vcfFormatKeys)
      )
    )

    # Set ranges for analysis and for VCF import
    if (length(RV[["activeGRanges"]]) > 0){
      TVTB::ranges(tparam) <- GenomicRanges::GRangesList(RV[["activeGRanges"]])
      VariantAnnotation::vcfWhich(TVTB::svp(tparam)) <-
        GenomicRanges::reduce(unlist(TVTB::ranges(tparam)))
    }
    if (!is.null(rownames(RV[["phenotypes"]]))){
      VariantAnnotation::vcfSamples(TVTB::svp(tparam)) <-
        rownames(RV[["phenotypes"]])
    }

    RV[["TVTBparam"]] <- tparam
  })

  # Update TVTBparam attached to VCF
  observeEvent(
    RV[["TVTBparam"]],
    {
      if (!is.null(RV[["VCF"]])){
        S4Vectors::metadata(RV[["VCF"]])[["TVTBparam"]] <- RV[["TVTBparam"]]
      }
    }
  )

  # Session information ----

  output$TVTBparamWarning <- renderUI({
    tryCatch(
      {
        validObject(RV[["TVTBparam"]], complete = TRUE)
        return('')
      },
      warning = function(w){
        return(
          tagList(tags$span(style="color:red", as.character(w)))
        )
      }
    )
  })

  output$TVTBsettings <- renderPrint({
    validate(need(is.null(Errors[["TVTBparam"]]), Errors[["TVTBparam"]]))
    return(RV[["TVTBparam"]])
  })

  # Limited number of key settings
  output$generalSettings <- renderPrint({
    return(list(
      # Phenotypes
      phenoFile = RV[["phenoFile"]],
      # VCF file(s)
      vcfinputMode = input$vcfInputMode,
      # Genotypes
      refGenotypes = input$refGenotypes,
      hetGenotypes = input$hetGenotypes,
      altGenotypes = input$altGenotypes,
      # GRanges
      grangeInputMode = input$grangeInputMode,
      genomicRanges = RV[["activeGRanges"]],
      # VCF parsing
      vepKey = input$vepKey,
      # Genome annotation
      annotationPackage = input$annotationPackage
    ))
  })

  # More numerous settings of secondary importance
  output$advancedSettings <- renderPrint({
    return(list(
      # Parallel settings
      bpCores = input$bpCores,
      bpCores = input$bpCores,
      bpType = input$bpType,
      # VCF settings
      singleVcf = RV[["singleVCF"]],
      vcfFolder = input$vcfFolder,
      vcfPattern = input$vcfPattern,
      yieldSize = input$yieldSize,
      # GRanges
      "bedFile" = RV[["bedFile"]],
      ucscRanges = input$ucscRanges,
      ensDb.type = input$ensDb.type,
      ensDb.condition = input$ensDb.condition,
      ensDb.value = input$ensDb.value,
      # Fields displayed
      vcfCols = input$vcfCols,
      vcfInfoCols = input$vcfInfoCols,
      vepCols = input$vepCols,
      phenoCols = input$phenoCols,
      # Genotypes displayed
      genoNumRow = input$genoNumRow,
      genoFirstRow = input$genoFirstRow,
      genoNumCols = input$genoNumCols,
      genoFirstCol = input$genoFirstCol
    ))
  })

  # sessionInfo()
  output$sessionInfo <- renderPrint({
    return(sessionInfo())
  })

  output$vcfMetadata <- renderPrint({
    validate(need(is.null(Errors[["VCF"]]), Errors[["VCF"]]))
    return(S4Vectors::metadata(RV[["VCF"]]))
  })

  # VCF folder  ----

  # List of files/folders in given folder
  observeEvent(
    input$vcfFolder,
    {
      # Start with optimism :)
      Errors[["VCFfolder"]] <- NULL

      if (!dir.exists(input$vcfFolder)){
        Errors[["VCFfolder"]] <- "Invalid VCF folder."
        RV[["vcfContent"]] <- character()
        return()
      }

      RV[["vcfContent"]] <- list.files(input$vcfFolder)

      if (length(RV[["vcfContent"]]) == 0){
        Errors[["VCFfolder"]] <- "Empty folder."
        return()
      }
    }
  )

  # Display full content of selected folder
  output$vcfContent <- DT::renderDataTable({
    validate(need(is.null(Errors[["VCFfolder"]]), Errors[["VCFfolder"]]))
    return(DT::datatable(
      data.frame(Filename = RV[["vcfContent"]]),
      rownames = FALSE
    ))
  })

  # Update the list of VCF files (matching pattern) in folder
  observe({
    # Start with optimism :)
    Errors[["VCFfiles"]] <- NULL

    if (!is.null(Errors[["VCFfolder"]])){
      Errors[["VCFfiles"]] <- Errors[["VCFfolder"]]
      RV[["VCFfiles"]] <- character()
      return()
    }

    # Check pattern
    if (!grepl("%s", input$vcfPattern)){
      Errors[["VCFfiles"]] <- "VCF file pattern must contain \"%s\""
      RV[["VCFfiles"]] <- character()
      return()
    }

    RV[["VCFfiles"]] <- grep(
      gsub('%s', '.*', input$vcfPattern),
      RV[["vcfContent"]],
      ignore.case = TRUE, value = TRUE)

    if (length(RV[["VCFfiles"]]) == 0){
      Errors[["VCFfiles"]] <- "No VCF file matching pattern in folder."
      return()
    }
  })

  # Easier to read for user
  output$vcfFiles <- DT::renderDataTable({

    validate(need(is.null(Errors[["VCFfiles"]]), Errors[["VCFfiles"]]))
    vcfFiles <- RV[["VCFfiles"]]

    return(DT::datatable(
      data.frame(Filename = vcfFiles),
      rownames = FALSE
    ))
  })

  # How many files detected, show first one.
  output$vcfFolderSummary <- renderUI({

    validate(need(is.null(Errors[["VCFfiles"]]), Errors[["VCFfiles"]]))
    vcfFiles <- RV[["VCFfiles"]]

    return(tagList(
      code(length(vcfFiles)), "file(s) matching pattern in folder", br(),
      "[", basename(head(x = vcfFiles, n = 1)), " , ... ]"
    ))
  })

  # Select single VCF ----

  # Demonstration input for single VCF file
  observeEvent(
    input$demoVcf,
    {
      RV[["singleVCF"]] <- system.file(
        "extdata/chr15.phase3_integrated.vcf.gz", package = "TVTB"
      )
      Errors[["singleVCF"]] <- NULL
    })

  # Action button to select single VCF file
  # store as reactive, as file choice cancelled will not reset to NULL otherwise
  observeEvent(
    input$selectVcf,
    {

      # Save the selected file in the reactive values
      RV[["singleVCF"]] <- tryCatch(
        file.choose(),
        error = function(err){
          Errors[["singleVCF"]] <- geterrmessage()
          return(NULL)
        })

      if (is.null(RV[["singleVCF"]])){
        Errors[["singleVCF"]] <- .Msgs[["singleVcf"]]
        return()
      }

      # Check that the selected file is *vcf.gz
      if (!grepl(".*\\.vcf\\.gz$", RV[["singleVCF"]], ignore.case = TRUE)){
        Errors[["singleVCF"]] <- sprintf(
          "File is not *.vcf.gz: %s",
          RV[["singleVCF"]]
        )
        return()
      }
      # If everything went well, allow subsequent steps
      Errors[["singleVCF"]] <- NULL
    }
  )

  # Path to single VCF file selected
  output$selectedVcf <- renderText({
    # Only proceed if a file was selected
    validate(need(is.null(Errors[["singleVCF"]]), Errors[["singleVCF"]]))
    return(RV[["singleVCF"]])
  })

  # Define ScanVcfParam ----

  # Identify available fields in header of (first) VCF file
  observe({
    # Start with optimism
    Errors[["VCFheader"]] <- NULL

    # Depends on input mode and selected VCF file/folder/pattern
    vcfRead <- switch (
      input$vcfInputMode,
      SingleVcf = {
        Errors[["VCFheader"]] <- Errors[["singleVCF"]]
        RV[["singleVCF"]]
      },
      OnePerChr = {
        Errors[["VCFheader"]] <- Errors[["VCFfiles"]]
        file.path(input$vcfFolder, head(RV[["VCFfiles"]], 1))
      }
    )

    if (!is.null(Errors[["VCFheader"]])){
      rawData <- NULL
    } else {
      rawData <- tryCatch(
        VariantAnnotation::scanVcfHeader(vcfRead),
        warning = function(warn){
          Errors[["VCFheader"]] <- warn
          NULL
        },
        error = function(err){
          Errors[["VCFheader"]] <- geterrmessage()
          NULL
        }
      )
    }

    if (is.null(rawData)){
      RV[["infoKeys"]] <- RV[["genoKeys"]] <- character()
      return()
    }

    # All keys except vep_key (required)
    RV[["infoKeys"]] <- grep(
      pattern = input$vepKey,
      x = c(rownames(VariantAnnotation::info(rawData))),
      invert = TRUE,
      value = TRUE
    )

    # All keys except "GT" (required)
    RV[["genoKeys"]] <- grep(
      pattern = "GT",
      x = c(rownames(VariantAnnotation::geno(rawData))),
      invert = TRUE,
      value = TRUE
    )

  })

  observeEvent(
    RV[["infoKeys"]],
    {
      updateSelectInput(
        session, "vcfInfoKeys",
        choices = RV[["infoKeys"]], selected = RV[["infoKeys"]]
      )
    }
  )

  observeEvent(
    RV[["infoKeys"]],
    {
      updateSelectInput(
        session, "vcfFormatKeys",
        choices = RV[["genoKeys"]], selected = RV[["genoKeys"]]
      )
    }
  )

  observeEvent(
    input$tickAllInfo,
    {

      updateSelectInput(
        session, "vcfInfoKeys",
        choices = RV[["infoKeys"]], selected = RV[["infoKeys"]]
      )

    }
  )

  observeEvent(
    input$untickAllInfo,
    {

      updateSelectInput(
        session, "vcfInfoKeys",
        choices = RV[["infoKeys"]], selected = c()
      )

    }
  )

  # Import VCF information ----

  # Import and expand VCF object
  observeEvent(
    input$importVariants,
    {
      Errors[["VCF"]] <- NULL

      # Fetch the first error from pre-requisites
      Errors[["VCF"]] <- head(c(
        Errors[["VCFheader"]],
        Errors[["TVTBparam"]],
        Errors[["GRanges"]],
        Errors[["phenotypes"]] # replace by TVTBparam when possible
      ), 1)

      if (!is.null(Errors[["VCF"]])){
        RV[["VCF"]] <- NULL
        return()
      }

      withProgress(
        min = 0, max = 3, value = 1,
        message = "Progress", detail = .Tracking[["preprocessing"]],
        {
          vcf <- switch (
            input$vcfInputMode,
            SingleVcf = {

              incProgress(
                amount = 1, detail = .Tracking[["singleVcf"]])
              tryCatch(
                parseSingleVcf(
                  RV[["singleVCF"]],
                  RV[["TVTBparam"]],
                  RV[["phenotypes"]],
                  input$autodetectGTimport,
                  input$yieldSize
                ),
                error = function(err){
                  Errors[["VCF"]] <- geterrmessage()
                  return(NULL)
                }
              )
            },

            OnePerChr = {

              incProgress(
                amount = 1, detail = .Tracking[["multiVcfs"]])

              tryCatch(
                parseMultipleVcf(
                  input$vcfFolder,
                  input$vcfPattern,
                  RV[["TVTBparam"]],
                  RV[["phenotypes"]],
                  input$autodetectGTimport,
                  input$yieldSize,
                  RV[["parallel"]]
                ),
                error = function(err){
                  Errors[["VCF"]] <- geterrmessage()
                  return(NULL)
                }
              )
            }
          )

        })

      if (!is.null(Errors[["VCF"]])){
        RV[["VCF"]] <- NULL
        return()
      }

      # Clean header of INFO fields not imported
      vcf <- TVTB::dropInfo(vcf)

      # Store ExpandedVCF object in list of reactive values
      RV[["VCF"]] <- vcf

      updateActionButton(
        session, "importVariants",
        label = "Refresh variants", icon = icon("refresh")
      )

      # Update selected genotypes if required
      if (input$autodetectGTimport){
        g <- TVTB::genos(S4Vectors::metadata(RV[["VCF"]])[["TVTBparam"]])
        updateSelectInput(
          session, "refGenotypes",
          choices = TVTB::genos(g),
          selected = TVTB::ref(g)
        )
        updateSelectInput(
          session, "hetGenotypes",
          choices = TVTB::genos(g),
          selected = TVTB::het(g)
        )
        updateSelectInput(
          session, "altGenotypes",
          choices = TVTB::genos(g),
          selected = TVTB::alt(g)
        )

      }
    }
  )

  # Summary of raw variants
  output$vcfSummary <- renderUI({

    validate(need(is.null(Errors[["VCF"]]), Errors[["VCF"]]))

    vcf <- RV[["VCF"]]

    return(tagList(
      code(nrow(vcf)), "bi-allelic records and",
      code(ncol(SummarizedExperiment::colData(vcf))), "phenotypes",
      "in", code(ncol(vcf)), "samples"
    ))
  })

  # Widget to control the meta-columns of VCF shown
  output$vcfCols <- renderUI({

    validate(need(is.null(Errors[["VCF"]]), Errors[["VCF"]]))

    vcf <- RV[["VCF"]]

    colChoices <- c(colnames(S4Vectors::mcols(vcf)))

    selectInput(
      "vcfCols", "Meta-columns",
      choices = colChoices,
      selected = c(colChoices[1:min(5, length(colChoices))]),
      multiple = TRUE
    )
  })

  # Show selected VCF meta-columns
  output$vcfRowRangesView <- DT::renderDataTable({
    # Widget displays the error message already
    validate(need(is.null(Errors[["VCF"]]), ""))

    vcf <- RV[["filteredVcf"]]
    validate(need(vcf, .Msgs[["filteredVcfNULL"]]), errorClass = "optional")

    cols <- which(colnames(S4Vectors::mcols(vcf)) %in% input$vcfCols)

    displayedTable <- cbind(
      rownames = rownames(vcf),
      as.data.frame(
        SummarizedExperiment::rowRanges(vcf)[, cols],
        row.names = NULL
      )
    )

    return(DT::datatable(
      data = displayedTable,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        searching = TRUE),
      filter = "top"))
  })

  # Widget to control the INFO columns shown
  output$vcfInfoCols <- renderUI({

    validate(need(is.null(Errors[["VCF"]]), Errors[["VCF"]]))

    vcf <- RV[["VCF"]]

    validate(need(
      # VEP field are an implicite field
      ncol(VariantAnnotation::info(vcf)) > 1,
      "No INFO field available."),
      errorClass = "optional"
    )

    validate(need(
      ncol(VariantAnnotation::info(vcf)) > 0,
      .Msgs[["importVariants"]])
    )

    # All columns except the VEP predictions
    colChoices <- grep(
      pattern = input$vepKey,
      x = colnames(VariantAnnotation::info(vcf)),
      invert = TRUE,
      value = TRUE)

    return(selectInput(
      "vcfInfoCols", "Meta-columns",
      choices = colChoices,
      selected = colChoices[1:min(5, length(colChoices))],
      multiple = TRUE
    ))
  })

  # Displayed requested VCF INFO fields
  output$vcfInfoView <- DT::renderDataTable({
    # Widget displays the error message already
    validate(need(is.null(Errors[["VCF"]]), ""))

    vcf <- RV[["filteredVcf"]]

    validate(need(
      # At least one non-VEP column
      ncol(VariantAnnotation::info(RV[["VCF"]])) > 1,
      "No INFO data imported."
    ),
    errorClass = "optional")

    validate(need(
      length(input$vcfInfoCols) > 0,
      "No INFO column selected."
    ))

    cols <- which(
      colnames(VariantAnnotation::info(vcf)) %in% input$vcfInfoCols
    )

    return(DT::datatable(
      data = cbind(
        rownames = rownames(vcf),
        as.data.frame(
          VariantAnnotation::info(vcf)[, cols, drop = FALSE],
          row.names = NULL)
      ),
      rownames = FALSE,
      options = list(
        pageLength = 10,
        searching = TRUE),
      filter = "top"))
  })

  # Display ExpandedVCF summary
  output$ExpandedVCF <- renderPrint({

    validate(need(is.null(Errors[["VCF"]]), Errors[["VCF"]]))
    validate(need(RV[["VCF"]], .Msgs[["importVariants"]]))

    return(RV[["VCF"]])
  })

  # Filter variants ----

  newVcfFilter <- reactive({

    # updated when the add button in clicked
    if (input$newFilterExpression == "")
      return(NULL)

    quickFix <- gsub("[“”]", "\"", input$newFilterExpression)
    quickFix <- gsub("[‘’]", "\'", quickFix)

    return(tryCatch(
      {
        newFilter <- new(
          input$newFilterClass,
          listData = list(parse(text = quickFix, keep.source = FALSE)),
          active = input$newFilterActive)

        if (class(newFilter) == "VcfVepRules")
          TVTB::vep(newFilter) <- input$vepKey

        return(newFilter)
      },
      # warning = function(w) NULL,
      error = function(e) NULL
    ))

  })

  # Demonstration input for the Filter input field
  observeEvent(
    input$demoFilter,
    {
      updateSelectInput(
        session, "newFilterClass",
        selected = "VcfVepRules"
      )
      updateCheckboxInput(
        session, "newFilterActive",
        value = TRUE
      )
      updateTextInput(
        session, "newFilterExpression",
        value = 'grepl("missense", Consequence)'
      )
    }
  )

  # Boolean whether the rule to add evaluates successfully
  observeEvent(
    input$addNewFilter,
    {
      newFilter <- newVcfFilter()
      vcf <- head(RV[["VCF"]]) # Speed up testing!

      if (is.null(vcf)){
        RV[["newFilterStatus"]] <-
          "Variants must be imported to test the expression"
        return(FALSE)
      }

      if (is.null(newFilter)){
        RV[["newFilterStatus"]] <- "Invalid expression"
        return(FALSE)
      }

      names(newFilter) <- paste0("rule", input$addNewFilter)

      return(tryCatch(
        {
          RV[["newFilterTestResults"]] <-
            is.logical(S4Vectors::eval(newFilter, vcf))
          RV[["newFilterStatus"]] <- "Valid"
        },
        #warning = function(e) FALSE,
        error = function(e){
          RV[["newFilterTestResults"]] <- FALSE
          RV[["newFilterStatus"]] <- ge$terrmessage()
        }
      ))
    }
  )

  output$vcfFilterTest <- renderUI({
    if (RV[["newFilterTestResults"]]){
      return(
        strong(tags$span(
          style="color:green",
          RV[["newFilterStatus"]])
        )
      )
    } else {
      return(
        strong(tags$span(
          style="color:red",
          RV[["newFilterStatus"]])
        )
      )
    }
  })

  observeEvent(
    input$addNewFilter,
    {
      newFilter <- newVcfFilter()

      # Only add new filter if valid
      req(RV[["newFilterTestResults"]])

      names(newFilter) <- paste0("rule", input$addNewFilter)

      newRules <- TVTB::VcfFilterRules(
        RV[["vcfFilters"]],
        newFilter)

      RV[["vcfFilters"]] <- newRules
    }
  )

  output$vcfFilterControls <- renderUI({

    vcfFilters <- RV[["vcfFilters"]]
    countFilters <- length(vcfFilters)

    if (countFilters == 0)
      return(p("No filter.", align = "center"))

    return(lapply(
      X = 1:countFilters,
      FUN = function(filterIndex){
        fluidRow(
          # type
          shiny::column(
            width = 1,
            code(TVTB::type(vcfFilters)[filterIndex])
          ),
          # active
          shiny::column(
            width = 1,
            checkboxInput(
              inputId = gsub(
                "rule",
                "active",
                names(vcfFilters)[filterIndex]),
              label = NULL,
              value = S4Vectors::active(vcfFilters)[filterIndex])
          ),
          # expression
          shiny::column(
            width = 8,
            code(as.character(vcfFilters[filterIndex][[1]]))
          ),
          # remove
          shiny::column(
            width = 1,
            actionButton(
              inputId = gsub(
                "rule",
                "removeFilter",
                names(vcfFilters)[filterIndex]),
              label = "Remove",
              icon = icon("remove")
            )
          )
        )
      }
    ))

  })

  output$vcfRules <- renderPrint({
    return(str(RV[["vcfFilters"]]))
  })

  # Observe checkboxes defining active status
  observe({

    # If the inputs are updated
    inputs <- input

    # Obtain the status of the VCF filters
    inputNames <- names(inputs)
    activeButtonNames <- grep(
      pattern = "^active[[:digit:]]*",
      x = inputNames,
      value = TRUE)
    rulesStatus <- sapply(
      X = activeButtonNames,
      FUN = function(activeName){
        inputs[[activeName]]
      },
      simplify = TRUE)

    # Rename to match rule names
    names(rulesStatus) <- gsub("active", "rule", names(rulesStatus))

    # NOTE: cannot work on the reactiveValues themselves
    # Extract from the reactiveValues
    isolate({vcfRules <- RV[["vcfFilters"]]})
    # Match status to rule order
    idx <- match(names(vcfRules), names(rulesStatus))
    # Update the active status
    S4Vectors::active(vcfRules) <- as.logical(rulesStatus)[idx]
    # Update in the reactiveValues
    RV[["vcfFilters"]] <- vcfRules
  })

  # Observe actionButtons to remove rules
  observe({
    # If the inputs are updated
    inputs <- input

    # Obtain the status of the VCF filters
    inputNames <- names(inputs)
    removeButtonNames <- grep(
      pattern = "^removeFilter[[:digit:]]*",
      x = inputNames,
      value = TRUE)

    removeStatus <- sapply(
      X = removeButtonNames,
      FUN = function(buttonName){
        inputs[[buttonName]]
      },
      simplify = TRUE)

    # Rename to match rule names
    names(removeStatus) <- gsub(
      pattern = "removeFilter",
      replacement = "rule",
      x = names(removeStatus))

    removeNames <- names(which(removeStatus > 0))

    # NOTE: cannot work on the reactiveValues themselves
    # Extract from the reactiveValues
    isolate({vcfRules <- RV[["vcfFilters"]]})
    # Identify the idx of the rules to remove
    removeIdx <- which(names(vcfRules) %in% removeNames)

    # Remove rules
    if (length(removeIdx) > 0)
      RV[["vcfFilters"]] <- TVTB::VcfFilterRules(vcfRules[-removeIdx])
  })

  # Updates filteredVcf using raw VCF and VcfFilterRules
  # updated when button is clicked -or- new variants imported (below)
  # (makes filtered variants synonym to raw variants in absence of filters)
  observeEvent(
    input$filterVariants,
    {
      vcf <- RV[["VCF"]]
      vcfFilters <- RV[["vcfFilters"]]

      # Store the filtered VCF in the reativeValues
      RV[["filteredVcf"]] <- S4Vectors::subsetByFilter(vcf, vcfFilters)
    }
  )

  # Updates filteredVcf using raw VCF and VcfFilterRules
  # updated when new variants imported
  observeEvent(
    RV[["VCF"]],
    {
      vcf <- RV[["VCF"]]
      vcfFilters <- RV[["vcfFilters"]]

      # Store the filtered VCF in the reativeValues
      RV[["filteredVcf"]] <- S4Vectors::subsetByFilter(vcf, vcfFilters)
    }
  )

  output$filteredVcfSummary <- renderUI({
    # First make sure that raw variants were imported
    validate(need(RV[["VCF"]], .Msgs[["importVariants"]]))

    filteredVcf <- RV[["filteredVcf"]]
    validate(need(RV[["filteredVcf"]], .Msgs[["filteredVcfNULL"]]))

    return(tagList(
      code(length(filteredVcf)), "bi-allelic records and",
      code(ncol(SummarizedExperiment::colData(filteredVcf))),
      "phenotypes in",
      code(ncol(filteredVcf)), "samples filtered"
    ))
  })

  output$filtersSummary <- renderUI({
    vcfFilters <- RV[["vcfFilters"]]

    return(tagList(
      code(sum(S4Vectors::active(vcfFilters))), "active filters",
      "( of", code(length(vcfFilters)), ")"
    ))

  })

  # Parse genotypes ----

  output$genotypeEncoding <- renderUI({
    validate(need(is.null(Errors[["VCF"]]), Errors[["VCF"]]))

    vcf <- RV[["VCF"]]

    return(tagList(
      "For information and suggestion purpose only:",
      tags$ul(
        tags$li(
          "Unique genotypes codes detected:",
          tags$code(paste0(deparse(RV[["GT.all"]]), collapse = ""))
        ),
        tags$li(
          "Auto-detected", tags$em("reference homozygote"),"genotypes:",
          tags$code(paste0(deparse(RV[["GT.autoRef"]]), collapse = ""))
        ),
        tags$li(
          "Auto-detected", tags$em("heterozygote"),"genotypes:",
          tags$code(paste0(deparse(RV[["GT.autoHet"]]), collapse = ""))
        ),
        tags$li(
          "Auto-detected", tags$em("alternate homozygote"),"genotypes:",
          tags$code(paste0(deparse(RV[["GT.autoAlt"]]), collapse = ""))
        )
      )
    ))

  })

  # Widget to control the number of samples shown
  output$genoNumCols <- renderUI({
    # genoNumRows widget displays the error message already
    validate(need(is.null(Errors[["VCF"]]), ""))

    genotypes <- VariantAnnotation::geno(RV[["filteredVcf"]])[["GT"]]

    sliderInput(
      "genoNumCols", "Number of columns (samples)",
      value = min(10, ncol(genotypes)),
      min = 2,
      max = min(50, ncol(genotypes)),
      step = 1)
  })

  # Widget to control the index of the first sample shown
  output$genoFirstCol <- renderUI({
    # genoNumRows widget displays the error message already
    validate(need(is.null(Errors[["VCF"]]), ""))

    req(input$genoNumCols)

    genotypes <- VariantAnnotation::geno(RV[["filteredVcf"]])[["GT"]]

    isolate({newValue <- max(1, input$genoFirstCol, na.rm = TRUE)})

    sliderInput(
      "genoFirstCol", "First column (sample)",
      value = newValue,
      min = 1,
      max = max(c(10, ncol(genotypes) - input$genoNumCols + 1)),
      step = 1)
  })

  # Widget to control the number of variants shown
  # Contrary to samples, variants can be filtered after import
  # therefore the widget is updated using the filteredVcf reactive
  output$genoNumRows <- renderUI({
    # Widget responsible for the error message, if any
    validate(need(is.null(Errors[["VCF"]]), Errors[["VCF"]]))

    genotypes <- VariantAnnotation::geno(RV[["filteredVcf"]])[["GT"]]

    sliderInput(
      "genoNumRows", "Number of rows (variants)",
      value = min(10, nrow(genotypes)),
      min = 2,
      max = min(100, nrow(genotypes)),
      step = 1)
  })

  # Widget to control the index of the first variant shown
  # Contrary to samples, variants can be filtered after import
  # therefore the widget is updated using the filteredVcf reactive
  output$genoFirstRow <- renderUI({
    # genoNumRows widget displays the error message already
    validate(need(is.null(Errors[["VCF"]]), ""))

    req(input$genoNumRows)

    genotypes <- VariantAnnotation::geno(RV[["filteredVcf"]])[["GT"]]

    isolate({newValue <- max(1, input$genoFirstRow, na.rm = TRUE)})

    sliderInput(
      "genoFirstRow", "First row (variant)",
      value = newValue,
      min = 1,
      max = nrow(genotypes) - input$genoNumRows + 1,
      step = 1)
  })

  # Display requested genotypes as a simple table
  output$genotypesSample <- renderTable({

    # genoNumRows widget displays the error message already
    validate(need(is.null(Errors[["VCF"]]), ""))

    req(
      input$genoFirstRow, input$genoNumRows,
      input$genoFirstCol, input$genoNumCols
    )

    vcf <- RV[["filteredVcf"]]

    req(
      input$genoFirstRow, input$genoNumRows,
      input$genoFirstCol, input$genoFirstRow)

    genoSampleRanges <- list(
      rows = rep(input$genoFirstRow, 2) + c(0, input$genoNumRows - 1),
      cols = rep(input$genoFirstCol, 2) + c(0, input$genoNumCols - 1)
    )

    genotypes <- VariantAnnotation::geno(vcf)[["GT"]]

    rows <- seq(genoSampleRanges$rows[1], genoSampleRanges$rows[2])
    cols <- seq(genoSampleRanges$cols[1], genoSampleRanges$cols[2])

    req(
      max(rows) <= nrow(genotypes),
      max(cols) <= ncol(genotypes)
    )

    as.matrix(genotypes[rows,cols])
  },
  rownames = TRUE)

  # Auto-detect genotypes in VCF
  observeEvent(
    RV[["VCF"]],
    {

      RV[["GT.all"]] <- na.exclude(unique(c(
        VariantAnnotation::geno(RV[["VCF"]])[["GT"]]
      )))
      g.OK <- grep("[[:digit:]][/|][[:digit:]]", RV[["GT.all"]], value=TRUE)
      RV[["GT.autoRef"]] <- grep("(0/0)|(0\\|0)", g.OK, value = TRUE)

      g.split <- limma::strsplit2(g.OK, "[/|]")

      if (ncol(g.split) == 2){
        RV[["GT.autoHet"]] <- g.OK[g.split[,1] != g.split[,2]]
      } else {
        RV[["GT.autoHet"]] <- NA
      }

      if (ncol(g.split) == 2){
        RV[["GT.autoAlt"]] <- g.OK[
          g.split[,1] == g.split[,2] &
            !g.OK %in% RV[["GT.autoRef"]]]
      } else {
        RV[["GT.autoAlt"]] <- NA
      }

    }
  )

  observeEvent(
    input$genotypeAutofill,
    {

      updateSelectInput(
        session, "refGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(RV[["GT.autoHet"]], RV[["GT.autoAlt"]])
        ),
        selected = RV[["GT.autoRef"]]
      )

      updateSelectInput(
        session, "hetGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(RV[["GT.autoRef"]], RV[["GT.autoAlt"]])
        ),
        selected = RV[["GT.autoHet"]]
      )

      updateSelectInput(
        session, "altGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(RV[["GT.autoRef"]], RV[["GT.autoHet"]])
        ),
        selected = RV[["GT.autoAlt"]]
      )

    }
  )

  observeEvent(
    RV[["GT.all"]],
    {
      updateSelectInput(
        session, "refGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(input$hetGenotypes, input$altGenotypes)
        ),
        selected = setdiff(
          RV[["GT.all"]],
          c(input$hetGenotypes, input$altGenotypes)
        )
      )
      updateSelectInput(
        session, "hetGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(input$refGenotypes, input$altGenotypes)
        ),
        selected = setdiff(
          RV[["GT.all"]],
          c(input$refGenotypes, input$altGenotypes)
        )
      )
      updateSelectInput(
        session, "altGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(input$refGenotypes, input$hetGenotypes)
        ),
        selected = setdiff(
          RV[["GT.all"]],
          c(input$refGenotypes, input$hetGenotypes)
        )
      )
    }, ignoreNULL = FALSE
  )

  # Remove genotypes added to REF from HET and ALT
  observeEvent(
    input$refGenotypes,
    {
      # req(input$altGenotypes, input$hetGenotypes)
      updateSelectInput(
        session, "hetGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(input$refGenotypes, input$altGenotypes)
        ),
        selected = setdiff(
          input$hetGenotypes,
          c(input$refGenotypes, input$altGenotypes)
        )
      )

      updateSelectInput(
        session, "altGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(input$refGenotypes, input$hetGenotypes)
        ),
        selected = setdiff(
          input$altGenotypes,
          c(input$refGenotypes, input$hetGenotypes)
        )
      )

    }, ignoreNULL = FALSE
  )

  observeEvent(
    input$hetGenotypes,
    {
      # req(input$altGenotypes, input$refGenotypes)
      updateSelectInput(
        session, "refGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(input$hetGenotypes, input$altGenotypes)
        ),
        selected = setdiff(
          input$refGenotypes,
          c(input$hetGenotypes, input$altGenotypes)
        )
      )

      updateSelectInput(
        session, "altGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(input$hetGenotypes, input$refGenotypes)
        ),
        selected = setdiff(
          input$altGenotypes,
          c(input$hetGenotypes, input$refGenotypes)
        )
      )

    }, ignoreNULL = FALSE
  )

  observeEvent(
    input$altGenotypes,
    {
      # req(input$altGenotypes, input$refGenotypes)
      updateSelectInput(
        session, "refGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(input$altGenotypes, input$hetGenotypes)
        ),
        selected = setdiff(
          input$refGenotypes,
          c(input$altGenotypes, input$hetGenotypes)
        )
      )

      updateSelectInput(
        session, "hetGenotypes",
        choices = setdiff(
          RV[["GT.all"]],
          c(input$altGenotypes, input$refGenotypes)
        ),
        selected = setdiff(
          input$hetGenotypes,
          c(input$altGenotypes, input$refGenotypes)
        )
      )

    }, ignoreNULL = FALSE
  )

  # Parse VEP predictions ----

  # Show information about VEP field
  output$vepStructure <- renderPrint({

    vcf <- RV[["VCF"]]

    # First make sure vcf exists
    validate(need(vcf, .Msgs[["importVariants"]]))
    # If it exists, check that vepKey exist in INFO fields
    validate(need(
      input$vepKey %in% colnames(VariantAnnotation::info(vcf)),
      .Msgs[["vepKeyNotFound"]]
    ))

    csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

    validate(need(csq, .Msgs[["csq"]]))

    str(csq)
  })

  # Widget to control the VEP fields shown
  output$vepCols <- renderUI({

    validate(need(is.null(Errors[["VCF"]]), Errors[["VCF"]]))

    vcf <- RV[["VCF"]]

    # Check that vepKey exist in INFO fields # TODO, use that stored in VCF:metadata:TVTB, no need to check
    validate(need(
      input$vepKey %in% colnames(VariantAnnotation::info(vcf)),
      .Msgs[["vepKeyNotFound"]]
    ))

    csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

    vepMcols <- S4Vectors::mcols(csq)

    selectInput(
      "vepCols", "Meta-columns",
      choices = colnames(vepMcols),
      selected = colnames(vepMcols)[1:5],
      multiple = TRUE)
  })

  # Display requested VEP fields
  output$vcfVepView <- DT::renderDataTable({
    # Widget displays the error message already
    validate(need(is.null(Errors[["VCF"]]), ""))

    vcf <- RV[["filteredVcf"]]
    validate(need(RV[["filteredVcf"]], .Msgs[["filteredVcfNULL"]]))

    validate(
      need(
        input$vepKey %in% colnames(VariantAnnotation::info(vcf)),
        .Msgs[["vepKeyNotFound"]])
    )

    csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

    vepMcols <- S4Vectors::mcols(csq)

    cols <- which(colnames(vepMcols) %in% input$vepCols)

    return(DT::datatable(
      cbind(
        rownames = names(csq),
        as.data.frame(
          csq[,cols, drop = FALSE],
          row.names = NULL) # avoid duplicate rownames
      ),
      rownames = FALSE, # don't show numeric rownames
      options = list(
        pageLength = 10,
        searching = TRUE),
      filter = "top"))
  })

  # Genotype heatmap ----

  # Heatmap of genotypes (ggplot)
  output$heatmapGenotype <- renderPlot({

    if (input$doGenoHeatmap == 0)
      return()

    # Remove dependency on vcf reactive
    isolate({ vcf <- RV[["filteredVcf"]] })

    validate(need(vcf, .Msgs[["importVariants"]]))

    withProgress(min = 0, max = 3, message = "Progress", {
      incProgress(1, detail = .Tracking[["calculate"]])

      genotypes <- VariantAnnotation::geno(vcf)[["GT"]]
      # TODO: use validGenotypes to set other GT to NA
      validGenotypes <- unlist(TVTB::genos(RV[["TVTBparam"]]))
      # Currently, all genotypes found in data are considered
      # Potentially, undesired genotypes could be considered NAs
      genos.long <- reshape2::melt(genotypes, value.name = "Genotype")

      colnames(genos.long)[1:2] <- c("Variants", "Samples")

      genos.long[,"Genotype"] <- factor(
        x = genos.long[,"Genotype"],
        levels = TVTB::genos(TVTB::genos(RV[["TVTBparam"]]))
      )

      incProgress(1, detail = .Tracking[["ggplot"]])
      gg <- ggplot2::ggplot(
        data = genos.long,
        mapping = ggplot2::aes(Samples, Variants)
      ) +
        ggplot2::geom_tile(
          ggplot2::aes(
            fill = Genotype,
            colour = Genotype
          )
        ) +
        # TODO: give choice (widget)
        ggplot2::scale_fill_discrete(drop = TRUE) +
        # tie with widget above
        ggplot2::scale_colour_discrete(drop = TRUE) +
        ggplot2::theme(
          axis.text = ggplot2::element_blank(),
          axis.title = ggplot2::element_text(
            size = ggplot2::rel(1.5)),
          axis.ticks = ggplot2::element_blank()
        )

      message("Plotting heatmap of genotypes...")
      incProgress(1, detail = .Tracking[["render"]])

      return(gg)
    })
  })

  # Parallel computing ----

  # When config is changed, update choices of type & cores
  observeEvent(
    input$bpConfig, {

      switch(
        input$bpConfig,
        SerialParam = {
          updateNumericInput(
            session, "bpCores",
            value = 1, min = 1, max = 1)
          updateSelectInput(
            session, "bpType",
            choices = .PS[["choices.type"]],
            selected = .PS[["default.type"]])
        },
        MulticoreParam = {
          updateSelectInput(
            session, "bpType",
            choices = "FORK",
            selected = "FORK")
          updateNumericInput(
            session, "bpCores",
            value = BiocParallel::multicoreWorkers(),
            min = 1, max = parallel::detectCores())
        },
        SnowParam = {
          updateSelectInput(
            session, "bpType",
            choices = c("SOCK", "MPI", "FORK"))
          updateNumericInput(
            session, "bpCores",
            value = BiocParallel::snowWorkers(),
            min = 1, max = parallel::detectCores())
        }
      )
    }
  )

  observe({
    req(input$bpConfig, input$bpType, input$bpCores)
    RV[["parallel"]] <- switch(
      input$bpConfig,
      SerialParam = BiocParallel::SerialParam(),
      MulticoreParam = BiocParallel::MulticoreParam(
        workers = input$bpCores
      ),
      SnowParam = BiocParallel::SnowParam(
        workers = input$bpCores,
        type = input$bpType
      )
    )
  })

  # Report of parallel mode on various systems tested
  output$parallelReport <- DT::renderDataTable({
    reportFile <- system.file("shinyApp", "rds", "parallel.rds", package = "TVTB")
    reportTable <- readRDS(reportFile)
    return(DT::datatable(
      reportTable,
      rownames = FALSE, # don't show numeric rownames
      options = list(
        pageLength = 10,
        searching = TRUE),
      filter = "top"))
  })

  # Cleanup routine ----

  cancel.onSessionEnded <- session$onSessionEnded(function() {
    # Reset original options
    options(.originalOptions)
    # Maybe something more useful in the future
    return(TRUE)
  })
})
