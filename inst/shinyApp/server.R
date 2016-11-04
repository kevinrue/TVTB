
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

source(file.path(
    system.file(package = "TVTB"),
    "shinyApp",
    "serverRoutines.R"))

shinyServer(function(input, output, clientData, session) {

    # Reactive values ----

    RV <- reactiveValues(
        # Path to phenotype file
        phenoFile = NULL, # Path to phenotype file
        # VCF filter rules
        vcfFilters = TVTB::VcfFilterRules(), # Existing rules
        # VCF keys
        infoKeys = NULL, # choices of info key
        genoKeys = NULL, # choices of geno keys
        # GRanges
        genomicRanges = GenomicRanges::GRanges(), # active ranges
        # Single VCF file
        singleVcf = NULL, # Path to single VCF file
        singleVcfIssue = NULL # NULL if no issue
    )

    # Import phenotype information ----

    observeEvent(
        eventExpr = input$selectPheno,
        handlerExpr = {
            RV[["phenoFile"]] <- tryCatch(
                file.choose(),
                error = function(err){
                    warning(geterrmessage())
                    return(NULL)
                })
    })

    observeEvent(
        eventExpr = input$demoPheno,
        handlerExpr = {
            RV[["phenoFile"]] <- system.file(
                "extdata/integrated_samples.txt", package = "TVTB"
            )
        })

    output$phenoFile <- renderText({

        phenoFile <- RV[["phenoFile"]]

        return(ifelse(
            is.null(phenoFile),
            "No file provided.",
            phenoFile))
    })

    # DataFrame of imported phenotypes, or NULL
    observeEvent(
        RV[["phenoFile"]],
        {
            phenoFile <- RV[["phenoFile"]]

            if (is.null(phenoFile)){
                RV[["phenotypes"]] <- NULL
            } else {
                message("Importing phenotypes ...")
                rawData <- tryParsePheno(phenoFile)

                validate(need(
                    all(dim(rawData) > c(0, 0)),
                    paste(
                        "Phenotype file must have at least 1 row and 1 column",
                        "with colnames and rownames")))
                RV[["phenotypes"]] <- S4Vectors::DataFrame(rawData)
            }
        },
        ignoreNULL = FALSE
    )

    # HTML summary of imported phenotypes
    output$phenoFileSummary <- renderUI({
        # Depends on RV[["phenotypes"]]

        phenotypes <- RV[["phenotypes"]]

        if (is.null(phenotypes))
            return(Msgs[["phenotypes"]])

        return(tagList(
            code(S4Vectors::ncol(phenotypes)),
            "phenotypes in",
            code(S4Vectors::nrow(phenotypes)),
            "samples."
        ))
    })

    # Column names available for selection from phenotypes
    output$phenoCols <- renderUI({
        # NOTE:
        # phenotype view depends on the filteredVcf(), not RV[["phenotypes"]]
        # 1) users can easily look at their raw phenotypes separately
        # 2) all analyses use data stored in the VCF object only

        vcf <- RV[["vcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        phenos <- SummarizedExperiment::colData(vcf)

        validate(need(
            ncol(phenos) > 0,
            ifelse(
                is.null(RV[["phenotypes"]]),
                Msgs[["colDataEmptyOK"]],
                Msgs[["colDataEmptyImport"]])
        ),
        errorClass = "optional")

        selectInput(
            "phenoCols", "Phenotypes",
            choices = colnames(phenos),
            selected = colnames(phenos)[1:5],
            multiple = TRUE
        )
    })

    # Display structure of phenotypes attached to variants
    output$phenotypesStructure <- renderPrint({

        # Depends on RV[["vcf"]]
        vcf <- RV[["vcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        phenos <- SummarizedExperiment::colData(vcf)
        validate(need(
            ncol(phenos) > 0,
            ifelse(
                is.null(RV[["phenotypes"]]),
                Msgs[["colDataEmptyOK"]],
                Msgs[["colDataEmptyImport"]])
        ),
        errorClass = "optional")

        validate(need(
            ncol(phenos) > 0,
            ifelse(
                is.null(RV[["phenotypes"]]),
                Msgs[["colDataEmptyOK"]],
                Msgs[["colDataEmptyImport"]])
            ),
            errorClass = "optional")

        str(phenos)
    })

    # Display table of phenotypes attached to variants
    output$phenotypesSample <- DT::renderDataTable({
        # Requires: input$phenoCols, RV[["filteredVcf"]]
        # Give time to initialise columns seletion widget
        req(input$phenoCols)

        # Display phenotypes attached to the VCF object
        vcf <- RV[["filteredVcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        # Make sure the selected fields are present in data
        phenos <- SummarizedExperiment::colData(vcf)
        validate(need(
            all(input$phenoCols %in% colnames(phenos)),
            "Invalid phenotypes selected: please refresh variants"))

        DT::datatable(
            as.data.frame(phenos[,input$phenoCols, drop = FALSE]),
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top")
    })

    # Update choice of phenotypes to calculate frequencies, when new VCF
    observeEvent(
        eventExpr = RV[["vcf"]],
        handlerExpr = {

            # raw VCF
            vcf <- RV[["vcf"]]

            # phenotypes
            phenos <- SummarizedExperiment::colData(vcf)

            updateSelectInput(
                session, "phenoAddFrequencies",
                choices = colnames(phenos))
    })

    # For each phenotype
    # Create a group of checkboxes for each level
    # to select those for which frequencies should be calculated
    observeEvent(
        eventExpr = input$phenoAddFrequencies,
        handlerExpr = {

            # raw VCF
            vcf <- RV[["vcf"]]
            validate(need(vcf, Msgs[["importVariants"]]))

            # phenotypes
            phenos <- SummarizedExperiment::colData(vcf)
            phenoNames <- colnames(phenos)
            phenoLevels <- levels(phenos[,input$phenoAddFrequencies])

            ## pre-tick phenoLevels already calculated
            infoCols <- colnames(VariantAnnotation::info(vcf))
            # phenoLevels already calculated have all suffixes present
            tparam <- tparam()
            suffixes <- c(
                names(TVTB::genos(tparam)),
                TVTB::aaf(tparam),
                TVTB::maf(tparam))

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
                })


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

    })

    # Calculation of frequencies ----

    # Select all phenotype levels on button click
    observeEvent(
        eventExpr = input$tickAllPhenoLevelsFreq,
        handlerExpr = {

            # raw VCF
            vcf <- RV[["vcf"]]
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
        eventExpr = input$untickAllPhenoLevelsFreq,
        handlerExpr = {

            updateCheckboxGroupInput(
                session, "phenoLevelFreqCheckboxes",
                selected = character(),
                inline = TRUE
            )

        }
    )

    # Add overall frequencies on button click
    observeEvent(
        eventExpr = input$addOverallFrequencies,
        handlerExpr = {

            vcf <- RV[["vcf"]]
            req(vcf)

            # collect all INFO keys
            tparam <- tparam()
            suffixes <- c(
                names(TVTB::genos(tparam)),
                TVTB::aaf(tparam),
                TVTB::maf(tparam))

            withProgress(
                max = 3,
                value = 1,
                message = "Progress",
                detail = Tracking[["preprocessing"]],{
                    # Only proceed if none of the INFO keys are present
                    if (!any(
                        suffixes %in%
                        colnames(VariantAnnotation::info(vcf)))){

                        incProgress(
                            amount = 1,
                            detail = Tracking[["addFreqOverall"]])

                        RV[["vcf"]] <- TVTB::addOverallFrequencies(
                            vcf = vcf,
                            force = TRUE # watch the console for warnings
                        )
                        RV[["latestPhenotypeFrequency"]] <- "Overall"
                        RV[["latestFrequenciesAdded"]] <- suffixes
                        RV[["latestFrequenciesRemoved"]] <- character()
                    }
                }
            )



    })

    # Remove overall frequencies on button click
    observeEvent(
        eventExpr = input$removeOverallFrequencies,
        handlerExpr = {

            vcf <- RV[["vcf"]]
            req(vcf)

            # collect all INFO keys
            tparam <- tparam()
            suffixes <- c(
                names(TVTB::genos(tparam)),
                TVTB::aaf(tparam),
                TVTB::maf(tparam))

            withProgress(
                max = 3,
                value = 1,
                message = "Progress",
                detail = Tracking[["preprocessing"]],{

                    incProgress(
                        amount = 1,
                        detail = Tracking[["rmFreqOverall"]])

                    # Only proceed if all of the INFO keys are present
                    if (all(
                        suffixes %in% colnames(VariantAnnotation::info(vcf)))){

                        RV[["vcf"]] <- TVTB::dropInfo(
                            vcf = vcf,
                            key = suffixes,
                            slot = "both"
                        )
                        RV[["latestPhenotypeFrequency"]] <- "Overall"
                        RV[["latestFrequenciesRemoved"]] <- suffixes
                        RV[["latestFrequenciesAdded"]] <- character()
                    }

            })

        })

    # Add selected & remove deselected frequencies on button click
    observeEvent(
        eventExpr = input$buttonFrequencies,
        handlerExpr = {

            withProgress(
                max = 4,
                value = 1,
                message = "Progress",
                detail = Tracking[["preprocessing"]],{
                    # Remove INFO keys for unticked boxes
                    # Add values for ticked boxes

                    # Phenotypes selected
                    selectedPhenoName <- input$phenoAddFrequencies
                    selectedPhenoLevels <- input$phenoLevelFreqCheckboxes

                    # Info in VCF
                    vcf <- RV[["vcf"]]
                    req(vcf)
                    vcfInfoCols <- colnames(VariantAnnotation::info(vcf))

                    # Phenotype levels available
                    phenos <- SummarizedExperiment::colData(vcf)
                    choices <- levels(phenos[,input$phenoAddFrequencies])

                    # pre1) collect all suffixes (to check existence of fields)
                    tparam <- tparam()
                    suffixes <- c(
                        names(TVTB::genos(tparam)),
                        TVTB::aaf(tparam),
                        TVTB::maf(tparam))

                    # pre2) identify unticked phenoLevels

                    phenoLevelsUnticked <- choices[which(
                        !choices %in% selectedPhenoLevels
                    )]

                    if (length(phenoLevelsUnticked) > 0){

                        incProgress(
                            amount = 1,
                            detail = Tracking[["rmFreqPhenoLevel"]])

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
                        RV[["vcf"]] <- TVTB::dropInfo(
                            vcf = RV[["vcf"]],
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
                            detail = Tracking[["addFreqPhenoLevel"]])

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
                        RV[["vcf"]] <- TVTB::addFrequencies(
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

        return("If you see this message, contact the maintainer please!")
    })

    # Define genomic ranges ----

    # Path to BED file
    observeEvent(
        eventExpr = input$selectBed,
        handlerExpr = {
            RV[["bedFile"]] <- tryCatch(
                file.choose(),
                error = function(err){
                    warning(geterrmessage())
                    return(NULL)
            })
    })

    # Demonstration input for the BED file
    observeEvent(
        eventExpr = input$demoBed,
        handlerExpr = {
            RV[["bedFile"]] <- system.file(
                "extdata/SLC24A5.bed", package = "TVTB"
            )
        }
    )

    # Demonstration input for the UCSC field input
    observeEvent(
        eventExpr = input$demoUCSC,
        handlerExpr = {
            updateTextInput(
                session, "ucscRanges",
                value = "15:48,413,169-48,434,869"
            )
        }
    )

    # Demonstration input for the EnsDb field input
    observeEvent(
        eventExpr = input$demoEnsDb,
        handlerExpr = {
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

    output$bedFile <- renderText({
        bedFile <- RV[["bedFile"]]

        return(ifelse(
            is.null(bedFile),
            "No file provided.",
            bedFile))
    })

    # If a new BED file is selected, update GRangesBED
    observeEvent(
        eventExpr = RV[["bedFile"]],
        handlerExpr = {

            stopifnot(requireNamespace("rtracklayer"))
            # use rtracklayer::import.bed to obtain GRanges
            bedFile <- RV[["bedFile"]]

            if (is.null(bedFile)){
                RV[["GRangesBED"]] <- GenomicRanges::GRanges()
                return()
            }

            message("Importing phenotypes ...")
            rawData <- tryParseBed(bedFile)

            validate(need(rawData, "Invalid input"))

            if (length(rawData) == 0)
                RV[["GRangesBED"]] <- GenomicRanges::GRanges()
            else
                RV[["GRangesBED"]] <- rawData
        },
        ignoreNULL = FALSE
    )

    # If a new BED file is selected, update GRangesBED
    observeEvent(
        eventExpr = input$ucscRanges,
        handlerExpr = {

            # parse the string or return NULL
            if (input$ucscRanges == ""){
                RV[["GRangesUCSC"]] <- GenomicRanges::GRanges()
                return()
            }

            # NOTE: do not trim "chr", for future UCSC support
            inputTrimmed <- gsub(
                pattern = ",| ",
                replacement = "",
                x = input$ucscRanges)

            # Split the given string into individual regions
            inputSplit <- strsplit(
                x = inputTrimmed, split = ";")[[1]]

            # Ensure that all regions are UCSC-valid
            validate(need(
                all(
                    sapply(
                        X = inputSplit,
                        FUN = function(x){
                            grepl(
                                pattern =
                                    "[[:alnum:]]+:[[:digit:]]+-[[:digit:]]+",
                                x = x)
                        })
                ),
                Msgs[["invalidUcscRanges"]]),
                errorClass = "optional")

            rawData <- tryCatch({

                chrs <- as.numeric(gsub(
                    pattern =
                        "([[:alnum:]]*):[[:digit:]]*-[[:digit:]]*",
                    replacement = "\\1",
                    x = inputSplit))

                starts <- as.numeric(gsub(
                    pattern =
                        "[[:alnum:]]*:([[:digit:]]*)-[[:digit:]]*",
                    replacement = "\\1",
                    x = inputSplit))

                ends <- as.numeric(gsub(
                    pattern =
                        "[[:alnum:]]*:[[:digit:]]*-([[:digit:]]*)",
                    replacement = "\\1",
                    x = inputSplit))

                data.frame(
                    V1 = chrs,
                    V2 = starts,
                    V3 = ends,
                    stringsAsFactors = FALSE)
            },
            error = function(err){
                warning(geterrmessage())
                return(NULL)
            }
            ,
            warning = function(warn){
                warning(warn)
                return(NULL)
            })

            validate(need(rawData, "Invalid input"))
            # Check number of columns before testing content of columns
            validate(need(
                all(dim(rawData) >= c(1, 3)),
                "BED file must have at least 1 row and 3 columns"))

            rawData <- validateDataFrameGRanges(rawData, selectedEnsDb())

            RV[["GRangesUCSC"]] <- rawData
        },
        ignoreNULL = FALSE
    )

    observeEvent(
        eventExpr = queryGenes(),
        handlerExpr = {

            queryGenes <- queryGenes()

            if (is.null(queryGenes)){
                RV[["GRangesEnsDb"]] <- GenomicRanges::GRanges()
                return()
            }

            validate(need(queryGenes, "Invalid input"))

            if (length(queryGenes) == 0)
                RV[["GRangesEnsDb"]] <- GenomicRanges::GRanges()
            else
                RV[["GRangesEnsDb"]] <- queryGenes
        },
        ignoreNULL = FALSE
    )

    observeEvent(
        eventExpr = RV[["GRangesBED"]],
        handlerExpr = {
            RV[["genomicRanges"]] <- RV[["GRangesBED"]]
        },
        ignoreNULL = FALSE
    )

    observeEvent(
        eventExpr = RV[["GRangesUCSC"]],
        handlerExpr = {
            RV[["genomicRanges"]] <- RV[["GRangesUCSC"]]
        },
        ignoreNULL = FALSE
    )

    observeEvent(
        eventExpr = RV[["GRangesEnsDb"]],
        handlerExpr = {
            RV[["genomicRanges"]] <- RV[["GRangesEnsDb"]]
        },
        ignoreNULL = FALSE
    )

    observeEvent(
        eventExpr = input$grangesInputMode,
        handlerExpr = {

            switch (input$grangesInputMode,
                bed = {
                    RV[["genomicRanges"]] <- RV[["GRangesBED"]]
                },
                ucsc = {
                    RV[["genomicRanges"]] <- RV[["GRangesUCSC"]]
                },
                EnsDb = {
                    RV[["genomicRanges"]] <- RV[["GRangesEnsDb"]]
                }
            )
        }
    )

    # How many BED records detected, show first one.
    output$rangesSummary <- renderUI({
        # Depends on genomicRanges

        genomicRanges <- RV[["genomicRanges"]]

        if (is.null(genomicRanges))
            return(HTML(Msgs[["invalidGenomicRanges"]]))

        if (length(genomicRanges) == 0)
            return(HTML(Msgs[["noGenomicRanges"]]))

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
    output$rangesSession <- renderPrint({

        genomicRanges <- RV[["genomicRanges"]]

        validate(need(
            !is.null(genomicRanges),
            Msgs[["noGenomicRanges"]]),
            errorClass = "optional")

        return(genomicRanges)
    })

    output$rangesSample <- DT::renderDataTable({

        genomicRanges <- RV[["genomicRanges"]]

        validate(need(
            !is.null(genomicRanges),
            Msgs[["noGenomicRanges"]]),
            errorClass = "optional")

        return(DT::datatable(
            data = as.data.frame(genomicRanges),
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top"))
    })

    # Genome annotation ----

    # The EnsDb object
    selectedEnsDb <- reactive({
        # Depends on input$annotationPackage
        validate(need(
            input$annotationPackage,
            Msgs[["annotationPackage"]]))

        validate(need(
            requireNamespace(input$annotationPackage),
            "Failed loading annotation package."
        ))

        return(getEdb(input$annotationPackage))
    })

    genomeSeqinfo <- reactive({

        selectedEnsDb <- selectedEnsDb()

        return(ensembldb::seqinfo(selectedEnsDb))
    })

    # EnsDb ----

    output$ensembl_organism <- renderUI({
        # Depends on selectedEnsDb
        edb <- selectedEnsDb()

        validate(need(edb, label = Msgs[["edb"]]))

        HTML(paste(
            tags$strong("Organism:"),
            ensembldb::organism(edb))
        )

    })

    output$ensembl_version <- renderUI({
        # Depends on selectedEnsDb
        edb <- selectedEnsDb()

        validate(need(edb, label = Msgs[["edb"]]))

        md <- ensembldb::metadata(edb)
        rownames(md) <- md$name

        HTML(paste(
            tags$strong("Ensembl version:"),
            md["ensembl_version", "value"])
        )
    })

    output$ensembl_genome <- renderUI({
        # Depends on selectedEnsDb
        edb <- selectedEnsDb()

        validate(need(edb, label = Msgs[["edb"]]))

        md <- ensembldb::metadata(edb)
        rownames(md) <- md$name

        HTML(paste(
            tags$strong("Genome build:"),
            md["genome_build", "value"])
        )
    })

    queryGenes <- reactive({
        # Depends on: selectedEnsDb, ...
        #   input$ensDb.type, input$ensDb.condition, input$ensDb.value

        edb <- selectedEnsDb()

        if (input$ensDb.value == "")
            return(NULL)

        ensDbFilter = EnsDbFilter(
            type = input$ensDb.type,
            condition = input$ensDb.condition,
            value = input$ensDb.value)

        res <- ensembldb::genes(edb, filter = ensDbFilter)

        return(res)
    })

    output$ensDb.Genes <- DT::renderDataTable({
        # Depends on: queryGenes

        queryGenes <- queryGenes()

        validate(need(
            length(queryGenes) > 0,
            "No genomic range to show."),
            errorClass = "optional")

        DT::datatable(
            data = as.data.frame(queryGenes),
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top")
    })

    # output$ensDb.Transcripts <- renderDataTable({
    #   if(length(input$package) == 0)
    #     return
    #   if(!is.na(input$geneName) &
    #      length(input$geneName) > 0 &
    #      input$geneName!=""){
    #     edb <- selectedEnsDb()
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
    #     edb <- selectedEnsDb()
    #     res <- exons(edb, filter=EnsDbFilter(input),
    #                  return.type="data.frame")
    #     assign(".ENS_TMP_RES", res, envir=globalenv())
    #     return(res)
    #   }
    # })

    # Session information ----

    tparam <- reactive({
        # Depends on: bpParam, refGenotypes, hetGenotypes, altGenotypes, ...
        #   vepKey

        bpParam <- bpParam()

        validate(
            need(input$vepKey, label = Msgs[["vepKey"]]),
            need(input$refGenotypes, Msgs[["refGenotypes"]]),
            need(input$hetGenotypes, Msgs[["hetGenotypes"]]),
            need(input$altGenotypes, Msgs[["altGenotypes"]]),
            need(input$refSuffix, Msgs[["refSuffix"]]),
            need(input$hetSuffix, Msgs[["hetSuffix"]]),
            need(input$altSuffix, Msgs[["altSuffix"]]),
            need(input$aafSuffix, Msgs[["aafSuffix"]]),
            need(input$mafSuffix, Msgs[["mafSuffix"]])
        )

        tparam <- TVTB::TVTBparam(
            genos = TVTB::Genotypes(
                ref = input$refGenotypes,
                het = input$hetGenotypes,
                alt = input$altGenotypes,
                suffix = c(
                    ref = input$refSuffix,
                    het = input$hetSuffix,
                    alt = input$altSuffix
                )
            ),
            aaf = input$aafSuffix,
            maf = input$mafSuffix,
            vep = input$vepKey,
            bp = bpParam,
            svp = VariantAnnotation::ScanVcfParam(
                info = c(input$vepKey, input$vcfInfoKeys),
                geno = c("GT", input$vcfFormatKeys)
            )
        )

        if (!is.null(RV[["genomicRanges"]])){
            VariantAnnotation::vcfWhich(TVTB::svp(tparam)) <-
                GenomicRanges::reduce(unlist(TVTB::ranges(tparam)))
        }

        phenotypes <- RV[["phenotypes"]]

        if (!is.null(phenotypes)){
            VariantAnnotation::vcfSamples(TVTB::svp(tparam)) <-
                rownames(phenotypes)
        }

        return(tparam)
    })

    output$TVTBsettings <- renderPrint({
        return(tparam())
    })

    output$generalSettings <- renderPrint({
        return(list(
            phenoFile = RV[["phenoFile"]],
            grangeInputMode = input$grangeInputMode,
            genomicRanges = RV[["genomicRanges"]],
            "bedFile" = RV[["bedFile"]],
            ucscRanges = input$ucscRanges,
            ensDb.type = input$ensDb.type,
            ensDb.condition = input$ensDb.condition,
            ensDb.value = input$ensDb.value,
            vcfinputMode = input$vcfInputMode,
            "singleVcf" = RV[["singleVcf"]],
            vcfFolder = input$vcfFolder,
            vcfPattern = input$vcfPattern,
            vepKey = input$vepKey,
            annotationPackage = input$annotationPackage,
            vcfCols = input$vcfCols,
            vcfInfoCols = input$vcfInfoCols,
            vepCols = input$vepCols,
            phenoCols = input$phenoCols,
            genoNumRow = input$genoNumRow,
            genoFirstRow = input$genoFirstRow,
            genoNumCols = input$genoNumCols,
            genoFirstCol = input$genoFirstCol
        ))
    })

    output$advancedSettings <- renderPrint({
        return(list(
            refGenotypes = input$refGenotypes,
            hetGenotypes = input$hetGenotypes,
            altGenotypes = input$altGenotypes,
            yieldSize = input$yieldSize,
            bpCores = input$bpCores,
            bpCores = input$bpCores,
            bpType = input$bpType
        ))
    })

    output$sessionInfo <- renderPrint({
        sessionInfo()
    })

    # VCF folder  ----

    # Full content of selected folder
    vcfContent <- reactive({
        # Depends on: vcfFolder
        validate(need(input$vcfFolder, label = Msgs[["vcfFolder"]]))

        validate(need(
            dir.exists(input$vcfFolder),
            "Invalid VCF folder"))

        list.files(path = input$vcfFolder)
    })

    # Display full content of selected folder
    output$vcfContent <- DT::renderDataTable({
        # Depends on vcfContent
        vcfContent <- vcfContent()

        validate(need(vcfContent, Msgs[["vcfContent"]]))

        DT::datatable(
            data = data.frame(Filename = vcfContent),
            rownames = FALSE)
    })

    ## List VCF files detected
    vcfFiles <- reactive({
        # Depends on vcfContent, vcfPattern

        vcfContent <- vcfContent()

        validate(
            need(vcfContent, Msgs[["vcfContent"]]),
            need(input$vcfPattern, label = Msgs[["vcfPattern"]]))

        validate(need(
            grepl(pattern = "%s", x = input$vcfPattern),
            "VCF file pattern must contain \"%s\""
        ))

        matchedFiles <- grep(
            pattern = gsub('%s', '.*', input$vcfPattern), x = vcfContent,
            ignore.case = TRUE, value = TRUE)

        validate(need(
            length(matchedFiles) > 0,
            "No VCF file matching pattern in folder"))

        matchedFiles
    })

    # Easier to read for user
    output$vcfFiles <- DT::renderDataTable({

        vcfFiles <- vcfFiles()

        validate(need(vcfFiles, Msgs[["vcfFiles"]]))

        DT::datatable(
            data = data.frame(Filename = vcfFiles),
            rownames = FALSE)
    })

    # How many files detected, show first one.
    output$vcfFolderSummary <- renderUI({

        vcfFiles <- vcfFiles()

        validate(need(vcfFiles, "vcfFiles"))

        return(tagList(
            code(length(vcfFiles)), "VCF.GZ file(s) detected", br(),
            "[", basename(head(x = vcfFiles, n = 1)), " , ... ]"
        ))
    })

    # Select single VCF ----

    # Demonstration input for single VCF file
    observeEvent(
        eventExpr = input$demoVcf,
        handlerExpr = {

            RV[["singleVcf"]] <- system.file(
                "extdata/chr15.phase3_integrated.vcf.gz", package = "TVTB"
            )

    })

    # Action button to select single VCF file
    observeEvent(
        eventExpr = input$selectVcf,
        handlerExpr = {

            selected <- tryCatch(
                file.choose(),
                error = function(err){
                    warning(geterrmessage())
                    return(NULL)
                })

            # Save the selected file in the reactive values
            validate(need(selected, Msgs[["singleVcf"]]))
            RV[["singleVcf"]] <- selected

            # Check that the selected file is *vcf.gz
            if (!grepl(
                pattern = ".*\\.vcf\\.gz$",
                x = selected,
                ignore.case = TRUE)){

                RV[["singleVcfIssue"]] <- sprintf(
                    "File is not *.vcf.gz: %s",
                    selected
                )
                warning(RV[["singleVcfIssue"]])
                return(NULL)

            }

            tbiChrVcf <- paste(selected, "tbi", sep = ".")
            # Check that the selected file also has a tabix index
            if (!file.exists(tbiChrVcf)){

                RV[["singleVcfIssue"]] <- sprintf(
                    "Tabix index file does not exist: %s",
                    tbiChrVcf)
                warning(RV[["singleVcfIssue"]])
                return(NULL)

            }

            RV[["singleVcfIssue"]] <- NULL
        })

    # Path to single VCF file selected
    output$selectedVcf <- renderText({

        singleVcf <- RV[["singleVcf"]]
        singleVcfIssue <- RV[["singleVcfIssue"]]

        # Only proceed if a file was selected
        validate(need(singleVcf, Msgs[["singleVcf"]]))

        # Display warning in console as well as web-app
        if (!is.null(singleVcfIssue))
            warning(singleVcfIssue)
        validate(need(is.null(singleVcfIssue), singleVcfIssue))

        return(singleVcf)
    })

    # Define ScanVcfParam ----

    # Identify available fields in header of (first) VCF file
    observe({
        # Depends on input mode, but also change of VCF files
        vcfInputMode <- input$vcfInputMode

        vcfHeader <- switch (
            vcfInputMode,
            SingleVcf = {

                singleVcf <- RV[["singleVcf"]]
                tryParseVcfHeader(file = singleVcf)
            },
            OnePerChr = {

                vcfFiles <- vcfFiles()
                tryParseVcfHeader(file = file.path(
                    input$vcfFolder,
                    vcfFiles[[1]]))

            }
        )

        if (is.null(vcfHeader)){
            RV[["infoKeys"]] <- NA
            RV[["genoKeys"]] <- NA
        } else {
            RV[["infoKeys"]] <- grep(
                pattern = input$vepKey,
                x = c(rownames(VariantAnnotation::info(vcfHeader))),
                invert = TRUE,
                value = TRUE)
            # All keys except "GT" (required)
            RV[["genoKeys"]] <- grep(
                pattern = "GT",
                x = c(rownames(VariantAnnotation::geno(vcfHeader))),
                invert = TRUE,
                value = TRUE)
        }

        # On update, select all INFO fields by default (may change)
        updateSelectInput(
            session, "vcfInfoKeys",
            choices = RV[["infoKeys"]], selected = RV[["infoKeys"]]
        )

        # On update, select all FORMAT fields by default (may change)
        updateSelectInput(
            session, "vcfFormatKeys",
            choices = RV[["genoKeys"]], selected = RV[["genoKeys"]]
        )

    })

    observeEvent(
        eventExpr = input$tickAllInfo,
        handlerExpr = {

            updateSelectInput(
                session, "vcfInfoKeys",
                choices = RV[["infoKeys"]], selected = RV[["infoKeys"]]
            )

        }
    )

    observeEvent(
        eventExpr = input$untickAllInfo,
        handlerExpr = {

            updateSelectInput(
                session, "vcfInfoKeys",
                choices = RV[["infoKeys"]], selected = c()
            )

        }
    )

    # Import VCF information ----

    observeEvent(
        eventExpr = input$importVariants,
        handlerExpr = {

        updateActionButton(
            session, "importVariants",
            label = "Refresh variants", icon = icon("refresh"))

    })

    # Import and expand VCF object
    observeEvent(
        eventExpr = input$importVariants,
        handlerExpr = {

            vcfInputMode <- input$vcfInputMode

            tparam <- tparam()

            # TODO: is.null & RV[["..."]]
            # Message in terminal before silent validation
            if(identical(input$annotationPackage, "")){
                warning(Msgs[["annotationPackage"]])
                return(NULL)
            }
            genomeSeqinfo <- genomeSeqinfo()

            yieldSize <- input$yieldSize
            infoKeys <- unique(c(input$vcfInfoKeys, input$vepKey))
            genoKeys <- input$vcfFormatKeys

            withProgress(
                min = 0, max = 3, value = 1,
                message = "Progress", detail = Tracking[["preprocessing"]],
                {

                    # # Timing
                    # t1 <- Sys.time()

                    vcf <- switch (
                        vcfInputMode,
                        SingleVcf = {

                            isolate({singleVcf <- RV[["singleVcf"]]})
                            validate(need(singleVcf, Msgs[["singleVcf"]]))

                            incProgress(
                                amount = 1, detail = Tracking[["singleVcf"]])

                            tryParseSingleVcf(
                                file = singleVcf,
                                param = tparam,
                                yieldSize = yieldSize
                            )
                        },

                        OnePerChr = {

                            incProgress(
                                amount = 1, detail = Tracking[["multiVcfs"]])

                            isolate({
                                vcfFolder <- input$vcfFolder
                                vcfPattern <- input$vcfPattern
                            })
                            validate(need(vcfFolder, Msgs[["vcfFolder"]]))
                            validate(need(vcfPattern, Msgs[["vcfPattern"]]))

                            tryParseMultipleVcf(
                                folder = vcfFolder,
                                pattern = vcfPattern,
                                param = tparam,
                                yieldSize = yieldSize,
                                BPPARAM = TVTB::bp(tparam)
                            )
                        }
                    )

                    incProgress(
                        amount = 1, detail = Tracking[["postprocessing"]])

                    validate(need(
                        length(vcf) > 0,
                        "No variant found in BED region(s)"))

                    if (is.null(RV[["phenoFile"]])){
                        SummarizedExperiment::colData(vcf) <-
                            S4Vectors::DataFrame(row.names = colnames(vcf))
                    } else {
                        SummarizedExperiment::colData(vcf) <-
                            RV[["phenotypes"]]
                    }

                })

            # Clean header of INFO fields not imported
            vcf <- TVTB::dropInfo(vcf = vcf)

            # Store ExpandedVCF object in list of reactive values
            RV[["vcf"]] <- vcf

        }
    )

    # Summary of raw variants
    output$vcfSummary <- renderUI({

        vcf <- RV[["vcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        return(tagList(
                code(nrow(vcf)), "bi-allelic records and",
                code(ncol(SummarizedExperiment::colData(vcf))), "phenotypes",
                "in", code(ncol(vcf)), "samples"
        ))
    })

    # Widget to control the meta-columns of VCF shown
    output$vcfCols <- renderUI({

        vcf <- RV[["vcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        colChoices <- c(colnames(S4Vectors::mcols(vcf)))

        selectInput(
            "vcfCols", "Meta-columns",
            choices = colChoices,
            selected = c(colChoices[1:min(5, length(colChoices))]),
            multiple = TRUE
        )
    })

    # Show selected VCF meta-columns
    output$vcfRowRanges <- DT::renderDataTable({

        vcf <- RV[["filteredVcf"]]

        # Give time to initialise widget
        req(input$vcfCols)

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

        vcf <- RV[["vcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        validate(need(
            # VEP field are an implicite field
            ncol(VariantAnnotation::info(vcf)) > 1,
            "No INFO field available."),
            errorClass = "optional"
        )

        validate(need(
            ncol(VariantAnnotation::info(vcf)) > 0,
            Msgs[["importVariants"]])
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
    output$vcfInfo <- DT::renderDataTable({
        req(input$vcfInfoCols)
        vcf <- RV[["filteredVcf"]]

        validate(need(
            # At least one non-VEP column
            ncol(VariantAnnotation::info(RV[["vcf"]])) > 1,
            "No INFO data imported."
        ),
        errorClass = "optional")

        validate(need(
            length(input$vcfInfoCols) > 0,
            "No INFO key selected."
        ))

        cols <- which(
            colnames(VariantAnnotation::info(vcf)) %in% input$vcfInfoCols
        )

        message("returning vcfInfo")
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
    output$vcf <- renderPrint({

        vcf <- RV[["vcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        vcf
    })

    # Filter variants ----

    newVcfFilter <- reactive({

        # updated when the add button in clicked
        if (input$newFilterExpression == "")
            return(NULL)

        quickFix <- gsub("[“”]", "\"", input$newFilterExpression)
        quickFix <- gsub("[‘’]", "\'", quickFix)

        return(tryCatch({
            newFilter <- new(
                input$newFilterClass,
                listData = list(parse(text = quickFix, keep.source = FALSE)),
                active = input$newFilterActive)

            if (class(newFilter) == "VcfVepRules")
                TVTB::vep(newFilter) <- input$vepKey

            return(newFilter)},
            # warning = function(w) NULL,
            error = function(e) NULL
        ))

    })

    # Demonstration input for the Filter input field
    observeEvent(
      eventExpr = input$demoFilter,
      handlerExpr = {
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
    newFilterTestResults <- reactive({
        # Only triggered by addNewFilter (> 0)
        req(input$addNewFilter)

        isolate({
            newFilter <- newVcfFilter()
            vcf <- head(RV[["vcf"]]) # Speed up testing!
        })

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
            expr = {
                testResult <- is.logical(S4Vectors::eval(newFilter, vcf))
                RV[["newFilterStatus"]] <- "Valid"
                return(testResult)
            },
            #warning = function(e) FALSE,
            error = function(e){
                RV[["newFilterStatus"]] <- geterrmessage()
                return(FALSE)
            }
        ))
    })

    output$vcfFilterTest <- renderUI({
        testResult <- newFilterTestResults()

        if (testResult){
            return(
                strong(tags$span(
                    style="color:green",
                    RV[["newFilterStatus"]])
                ))
        } else {
            return(
                strong(tags$span(
                    style="color:red",
                    RV[["newFilterStatus"]])
                ))
        }
    })

    observeEvent(input$addNewFilter, {
        newFilter <- newVcfFilter()
        testResult <- newFilterTestResults()

        # Only add new filter if valid
        validate(need(testResult == TRUE, "Invalid VCF filter"))

        names(newFilter) <- paste0("rule", input$addNewFilter)

        newRules <- TVTB::VcfFilterRules(
            RV[["vcfFilters"]],
            newFilter)

        # names(newRules) <- paste0("rule", 1:length(newRules))

        RV[["vcfFilters"]] <- newRules
    })

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
                    column(
                        width = 1,
                        code(TVTB::type(vcfFilters)[filterIndex])
                    ),
                    # active
                    column(
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
                    column(
                        width = 8,
                        code(as.character(vcfFilters[filterIndex][[1]]))
                    ),
                    # remove
                    column(
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
    # updated when raw variants are refreshed & when button is clicked
    # (makes filtered variants synonym to raw variants in absence of filters)
    observe({
        input$filterVariants
        vcf <- RV[["vcf"]]

        # Do not update when filters change, wait for actionButton or VCF
        isolate({vcfFilters <- RV[["vcfFilters"]]})

        # Store the filtered VCF in the reativeValues
        RV[["filteredVcf"]] <- S4Vectors::subsetByFilter(vcf, vcfFilters)
    })

    output$filteredVcfSummary <- renderUI({
        filteredVcf <- RV[["filteredVcf"]]

        # Before filtered variants can be shown, raw variants must be loaded
        validate(need(RV[["vcf"]], Msgs[["importVariants"]]))
        # filtered variants are calculated as soon as raw variants are loaded
        # to avoid confusing users if they don't want to apply
        # validate(need(filteredVcf, Msgs[["filteredVcf"]]))

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

    # Structure of raw genotypes
    output$genotypeStructure <- renderPrint({

        vcf <- RV[["vcf"]]

        validate(need(vcf, Msgs[["importVariants"]]))

        genotypes <- VariantAnnotation::geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

        str(genotypes)
    })

    # Widget to control the number of samples shown
    output$genoNumCols <- renderUI({

        vcf <- RV[["vcf"]]

        validate(need(vcf, Msgs[["importVariants"]]))

        genotypes <- VariantAnnotation::geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

        sliderInput(
            "genoNumCols", "Number of columns (samples)",
            value = min(10, ncol(genotypes)),
            min = 2,
            max = min(50, ncol(genotypes)),
            step = 1)
    })

    # Widget to control the index of the first sample shown
    output$genoFirstCol <- renderUI({

        req(input$genoNumCols)

        vcf <- RV[["vcf"]]

        validate(need(vcf, Msgs[["importVariants"]]))

        genotypes <- VariantAnnotation::geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

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

        vcf <- RV[["filteredVcf"]]

        validate(need(vcf, Msgs[["importVariants"]]))

        genotypes <- VariantAnnotation::geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

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

        req(input$genoNumRows)

        vcf <- RV[["filteredVcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        genotypes <- VariantAnnotation::geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

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

        req(
            input$genoFirstRow, input$genoNumRows,
            input$genoFirstCol, input$genoNumCols
        )

        vcf <- RV[["filteredVcf"]]

        validate(need(vcf, Msgs[["filteredVcf"]]))

        req(
            input$genoFirstRow, input$genoNumRows,
            input$genoFirstCol, input$genoFirstRow)

        genoSampleRanges <- list(
            rows = rep(input$genoFirstRow, 2) + c(0, input$genoNumRows - 1),
            cols = rep(input$genoFirstCol, 2) + c(0, input$genoNumCols - 1)
        )

        genotypes <- VariantAnnotation::geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

        rows <- seq(genoSampleRanges$rows[1], genoSampleRanges$rows[2])
        cols <- seq(genoSampleRanges$cols[1], genoSampleRanges$cols[2])

        req(
            max(rows) <= nrow(genotypes),
            max(cols) <= ncol(genotypes)
        )

        as.matrix(genotypes[rows,cols])
    },
    rownames = TRUE)

    # Remove genotypes added to REF from HET and ALT
    observeEvent(input$refGenotypes, {
        req(input$altGenotypes, input$hetGenotypes)
        updateSelectInput(
            session, "hetGenotypes",
            choices = setdiff(
                GS[["all.genotypes"]],
                c(input$refGenotypes, input$altGenotypes)),
            selected = setdiff(
                input$hetGenotypes,
                c(input$refGenotypes, input$altGenotypes)))

        updateSelectInput(
            session, "altGenotypes",
            choices = setdiff(
                GS[["all.genotypes"]],
                c(input$refGenotypes, input$hetGenotypes)),
            selected = setdiff(
                input$altGenotypes,
                c(input$refGenotypes, input$hetGenotypes)))

    })

    observeEvent(input$hetGenotypes, {
        req(input$altGenotypes, input$refGenotypes)
        updateSelectInput(
            session, "refGenotypes",
            choices = setdiff(
                GS[["all.genotypes"]],
                c(input$hetGenotypes, input$altGenotypes)),
            selected = setdiff(
                input$refGenotypes,
                c(input$hetGenotypes, input$altGenotypes)))

        updateSelectInput(
            session, "altGenotypes",
            choices = setdiff(
                GS[["all.genotypes"]],
                c(input$hetGenotypes, input$refGenotypes)),
            selected = setdiff(
                input$altGenotypes,
                c(input$hetGenotypes, input$refGenotypes)))

    })

    observeEvent(input$altGenotypes, {
        req(input$altGenotypes, input$refGenotypes)
        updateSelectInput(
            session, "refGenotypes",
            choices = setdiff(
                GS[["all.genotypes"]],
                c(input$altGenotypes, input$hetGenotypes)),
            selected = setdiff(
                input$refGenotypes,
                c(input$altGenotypes, input$hetGenotypes)))

        updateSelectInput(
            session, "hetGenotypes",
            choices = setdiff(
                GS[["all.genotypes"]],
                c(input$altGenotypes, input$refGenotypes)),
            selected = setdiff(
                input$hetGenotypes,
                c(input$altGenotypes, input$refGenotypes)))

    })

    # Parse VEP predictions ----

    # Show information about VEP field
    output$vepStructure <- renderPrint({

        vcf <- RV[["vcf"]]

        # First make sure vcf exists
        validate(need(vcf, Msgs[["importVariants"]]))
        # If it exists, check that vepKey exist in INFO fields
        validate(need(
            input$vepKey %in% colnames(VariantAnnotation::info(vcf)),
            Msgs[["vepKeyNotFound"]]
        ))

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        validate(need(csq, Msgs[["csq"]]))

        str(csq)
    })

    # Widget to control the VEP fields shown
    output$vepCols <- renderUI({

        vcf <- RV[["vcf"]]
        # First make sure vcf exists
        validate(need(vcf, Msgs[["importVariants"]]))
        # If it exists, check that vepKey exist in INFO fields
        validate(need(
            input$vepKey %in% colnames(VariantAnnotation::info(vcf)),
            Msgs[["vepKeyNotFound"]]
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
    output$vcfVep <- DT::renderDataTable({
        req(input$vepCols)

        vcf <- RV[["filteredVcf"]]

        validate(
            need(vcf, Msgs[["importVariants"]]),
            need(
                input$vepKey %in% colnames(VariantAnnotation::info(vcf)),
                Msgs[["vepKeyNotFound"]])
        )

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        vepMcols <- S4Vectors::mcols(csq)

        cols <- which(colnames(vepMcols) %in% input$vepCols)

        DT::datatable(
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
            filter = "top")
    })

    # Genotype heatmap ----

    # Heatmap of genotypes (ggplot)
    output$heatmapGenotype <- renderPlot({

        if (input$doGenoHeatmap == 0)
            return()

        # Remove dependency on vcf reactive
        isolate({ vcf <- RV[["filteredVcf"]] })

        validate(
            need(vcf, Msgs[["importVariants"]]),
            need(tparam(), label = Msgs[["tparam"]])
        )

        withProgress(min = 0, max = 3, message = "Progress", {
            incProgress(1, detail = Tracking[["calculate"]])

            genotypes <- VariantAnnotation::geno(vcf)[["GT"]]
            validate(need(genotypes, Msgs[["genotypes"]]))
            # TODO: use validGenotypes to set other GT to NA
            validGenotypes <- unlist(TVTB::genos(tparam()))
            # Currently, all genotypes found in data are considered
            # Potentially, undesired genotypes could be considered NAs
            genos.long <- reshape2::melt(genotypes, value.name = "Genotype")

            colnames(genos.long)[1:2] <- c("Variants", "Samples")

            genos.long[,"Genotype"] <- factor(
                x = genos.long[,"Genotype"],
                levels = TVTB::genos(TVTB::genos(tparam()))
            )

            incProgress(1, detail = Tracking[["ggplot"]])
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
            incProgress(1, detail = Tracking[["render"]])

            return(gg)
        })
    })

    # Parallel computing ----

    # When config is changed, update choices of type & cores
    observeEvent(input$bpConfig, {

        switch(
            input$bpConfig,
            SerialParam = {
                updateNumericInput(
                    session, "bpCores",
                    value = 1, min = 1, max = 1)
                updateSelectInput(
                    session, "bpType",
                    choices = PS[["choices.type"]],
                    selected = PS[["default.type"]])
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
            })
    })

    # Update parallel configuration object
    bpParam <- reactive({
        # Give time to initialise the widget
        req(input$bpConfig)

        # req() causes hanging below
        # Instead validate allows widgets to be initialised
        validate(
            need(input$bpType, Msgs[["bpType"]]),
            need(input$bpCores, Msgs[["bpCores"]]))

        return(switch(
            input$bpConfig,
            SerialParam = BiocParallel::SerialParam(),
            MulticoreParam = BiocParallel::MulticoreParam(
                workers = input$bpCores),
            SnowParam = BiocParallel::SnowParam(
                workers = input$bpCores, type = input$bpType)
            ))
    })

    # Cleanup routine ----

    cancel.onSessionEnded <- session$onSessionEnded(function() {
        # Reset original options
        options(originalOptions)
        # Maybe something more useful in the future
        return(TRUE)
    })
})
