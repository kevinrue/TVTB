
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
source("~/Dropbox/TVTB/inst/shinyApp/fileParseFunctions.R")

shinyServer(function(input, output, clientData, session) {

    # Reactive values ----

    RV <- reactiveValues(
        # Path to phenotype file
        phenoFile = NULL,
        # VCF filter rules
        vcfFilters = VcfFilterRules(),
        # VCF keys
        infoKeys = NULL,
        genoKeys = NULL
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

    output$phenoFile <- renderText({

        phenoFile <- RV[["phenoFile"]]

        return(ifelse(
            is.null(phenoFile),
            "No file provided.",
            phenoFile))
    })

    # DataFrame of imported phenotypes, or NULL
    phenotypes <- reactive({
        # Depends on phenoFile
        phenoFile <- RV[["phenoFile"]]

        if (is.null(phenoFile))
            return(NULL)
        else {
            message("Importing phenotypes ...")
            rawData <- tryParsePheno(file = phenoFile)

            validate(need(
                all(dim(rawData) > c(0, 0)),
                paste(
                    "Phenotype file must have at least 1 row and 1 column",
                    "with colnames and rownames")))
        }

        DataFrame(rawData)
    })

    # HTML summary of imported phenotypes
    output$phenoFileSummary <- renderUI({
        # Depends on phenotypes()

        phenotypes <- phenotypes()

        if (is.null(phenotypes))
            return(Msgs[["phenotypes"]])

        return(tagList(
            code(ncol(phenotypes)),
            "phenotypes in",
            code(nrow(phenotypes)),
            "samples."
        ))
    })

    # Column names available for selection from phenotypes
    output$phenoCols <- renderUI({
        # NOTE:
        # phenotype view depends on the filteredVcf(), not phenotypes()
        # 1) users can easily look at their raw phenotypes separately
        # 2) all analyses use data stored in the VCF object only

        vcf <- RV[["filteredVcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        phenos <- colData(vcf)

        validate(need(
            ncol(phenos) > 0,
            ifelse(
                is.null(phenotypes()),
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

    # When phenotypes are updated in VCF, update the choices for plotting
    observeEvent(RV[["filteredVcf"]], {
        # Depends on RV[["filteredVcf"]] & input$phenoTVBP
        pheno.choices <- c("None", colnames(colData(RV[["filteredVcf"]])))

        # If new phenotype files also contains the current phenotype,
        # keep it active instead of resetting to the first choice
        ## Tabulate VEP by phenotype (TVBP)
        pheno.TVBPselected <- pheno.choices[
            max(1, which(pheno.choices == input$phenoTVBP))]

        updateSelectInput(
            session, "phenoTVBP",
            choices = pheno.choices,
            selected = pheno.TVBPselected)

        # Same thing for DVBP
        pheno.DVBPselected <- pheno.choices[
            max(1, which(pheno.choices == input$phenoDVBP))]

        updateSelectInput(
            session, "phenoDVBP",
            choices = pheno.choices,
            selected = pheno.DVBPselected)
    })

    # Display structure of phenotypes attached to variants
    output$phenotypesStructure <- renderPrint({
        # Depends on RV[["vcf"]]
        vcf <- RV[["vcf"]]
        phenos <- SummarizedExperiment::colData(vcf)

        validate(need(
            ncol(phenos) > 0,
            ifelse(
                is.null(phenotypes()),
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

        # Make sure the selected
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

    observeEvent(
        eventExpr = RV[["vcf"]],
        handlerExpr = {

            # raw VCF
            vcf <- RV[["vcf"]]

            # phenotypes
            phenos <- colData(vcf)

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
            phenos <- colData(vcf)
            phenoNames <- colnames(phenos)
            phenoLevels <- levels(phenos[,input$phenoAddFrequencies])

            ## pre-tick phenoLevels already calculated
            infoCols <- colnames(info(vcf))
            # phenoLevels already calculated have all suffixes present
            tparam <- tparam()
            suffixes <- c(
                names(genos(tparam)),
                aaf(tparam),
                maf(tparam))

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

    observeEvent(
        eventExpr = input$tickAllPhenoLevelsFreq,
        handlerExpr = {

            # raw VCF
            vcf <- RV[["vcf"]]

            # phenotypes
            phenos <- colData(vcf)
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

    observeEvent(
        eventExpr = input$buttonFrequencies,
        handlerExpr = {

            # Remove INFO keys for unticked boxes
            # Add values for ticked boxes

            # Phenotypes selected
            selectedPhenoName <- input$phenoAddFrequencies
            selectedPhenoLevels <- input$phenoLevelFreqCheckboxes

            # Info in VCF
            vcf <- RV[["vcf"]]
            vcfInfoCols <- colnames(info(vcf))

            # Phenotype levels available
            phenos <- colData(vcf)
            choices <- levels(phenos[,input$phenoAddFrequencies])

            # pre1) collect all suffixes (to check existence of fields)
            tparam <- tparam()
            suffixes <- c(
                names(genos(tparam)),
                aaf(tparam),
                maf(tparam))

            # pre2) identify unticked phenoLevels

            phenoLevelsUnticked <- choices[which(
                !choices %in% selectedPhenoLevels
            )]

            if (length(phenoLevelsUnticked) > 0){

                # Identify unticked phenoLevels present in vcf INFO (to remove)
                # 1) Generate all expected key names for those phenoLevels
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
                    untickedPhenoLevelKeys[which(areAllPhenoLevelKeysPresent)]
                ) # vector of keys to remove

                if (length(phenoLevelKeysRemove))
                    RV[["vcf"]] <- dropInfo(
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

            # Identify ticked phenoLevels absent in info slot (to calculate)

            if (length(selectedPhenoLevels) > 0){

                # 1) generate all keys expected for each phenoLevel (named list)
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

                phenoLevelsAdd <- names(areAllPhenoLevelKeysAbsent)[which(
                    areAllPhenoLevelKeysAbsent
                )] # phenoLevels to add

                # Format for addFrequencies input
                phenosAdd <- list(phenoLevelsAdd)
                names(phenosAdd)[[1]] <- selectedPhenoName

                if (length(phenoLevelsAdd) > 0)
                    RV[["vcf"]] <- addFrequencies(
                        vcf = vcf,
                        phenos = phenosAdd,
                        param = tparam(),
                        force = TRUE # keep an eye on the console for warnings
                    )

                RV[["latestFrequenciesAdded"]] <- phenosAdd[[1]]
            } else {
                RV[["latestFrequenciesAdded"]] <- "NA"
            }

            RV[["latestPhenotypeFrequency"]] <- selectedPhenoName

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
                addedLevels <- "NA"

            if(length(freqRemoved > 0))
                removedLevels <- paste(freqRemoved, collapse = ", ")
            else
                removedLevels <- "NA"

            return(tagList(
                "Phenotype", code(pheno), ":",
                tags$ul(
                    tags$li("Added level(s):", code(addedLevels)),
                    tags$li("Removed level(s):", code(removedLevels))
                )
            ))
        }

        # TODO: overall frequencies

    })

    # Define genomic ranges ----

    # Path to BED file
    bedFile <- reactive({
        # Triggered by input$selectBed
        if (input$selectBed > 0){
            selected <- tryCatch(
                file.choose(),
                error = function(err){
                    warning(geterrmessage())
                    return(NULL)
                })
        } else {
            return(NULL)
        }

        selected
    })

    output$bedFile <- renderText({
        bedFile <- bedFile()

        return(ifelse(
            is.null(bedFile),
            "No file provided.",
            bedFile))
    })

    # Import BED records, format as GRanges
    genomicRanges <- reactive({
        # Depends on input$regionInputMode
        # Either...
        # Requires: bedFile()
        # Requires: input$ucscRegions
        # Requires: input$ensDb.type, input$ensDb.condition, input$ensDb.value
        switch (input$regionInputMode,
                bed = {
                    # use rtracklayer::import.bed to obtain GRanges
                    bedFile <- bedFile()

                    if (is.null(bedFile))
                        return(NULL)
                    else {
                        message("Importing phenotypes ...")
                        rawData <- tryParseBed(bedFile)
                    }
                },
                ucsc = {
                    # parse the string or return NULL
                    if (input$ucscRegions == "")
                        return(NULL)

                    # NOTE: do not trim "chr", for future UCSC support
                    inputTrimmed <- gsub(
                        pattern = ",| ",
                        replacement = "",
                        x = input$ucscRegions)

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
                                        pattern = "[[:alnum:]]+:[[:digit:]]+-[[:digit:]]+",
                                        x = x)
                                })
                        ),
                        Msgs[["invalidUcscRegions"]]),
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
                        return(NULL)}
                    ,
                    warning = function(warn){
                        warning(warn)
                        return(NULL)
                    })
                },
                EnsDb = {
                    # use ensembldb query to obtain GRanges
                    queryGenes <- queryGenes()

                    if (is.null(queryGenes))
                        return(NULL)

                    rawData <- queryGenes # Re-use existing value!
                }
        )
        # UCSC is converted to a data.frame which requires a few extra checks
        validate(need(
            class(rawData) %in% c("data.frame", "GRanges"),
            "Invalid input"))

        switch (class(rawData),
                data.frame = {
                    # Test the number of columns before trying to access them
                    validate(need(
                        all(dim(rawData) >= c(1, 3)),
                        "BED file must have at least 1 row and 3 columns"))
                    validate(
                        need(
                            all(
                                rawData[,1] %in%
                                    seqlevels(selectedPackage())),
                            "First column contains invalid chromosome names"),
                        need(
                            is.numeric(rawData[,2]),
                            "Second column is not numeric"),
                        need(
                            is.numeric(rawData[,3]),
                            "Third column is not numeric"))

                    return(
                        GRanges(
                            seqnames = rawData[,1],
                            ranges = IRanges(
                                start = rawData[,2],
                                end = rawData[,3]))
                    )
                },
                GRanges = {
                    if (length(rawData) == 0)
                        return(NULL)

                    return(rawData)
                }
        )
    })

    # How many BED records detected, show first one.
    output$rangesSummary <- renderUI({
        # Depends on genomicRanges

        genomicRanges <- genomicRanges()

        if (is.null(genomicRanges))
            return(HTML(Msgs[["genomicRanges"]]))

        return(tagList(
            code(length(genomicRanges)), "genomic range(s)", br(),
            "[",
            as.character(head(x = seqnames(genomicRanges), n = 1)),
            ":",
            as.character(head(x = start(genomicRanges), n = 1)),
            "-",
            as.character(head(x = end(genomicRanges), n = 1)),
            # "{",
            # as.character(head(x = names(genomicRanges), n = 1)),
            # "}",
            " , ... ]"
        ))
    })

    # Show BED records
    output$rangesStructure <- renderPrint({

        genomicRanges <- genomicRanges()

        validate(need(
            genomicRanges,
            Msgs[["genomicRanges"]]),
            errorClass = "optional")

        str(genomicRanges)
    })

    output$rangesSample <- DT::renderDataTable({

        genomicRanges <- genomicRanges()

        validate(need(
            genomicRanges,
            Msgs[["genomicRanges"]]),
            errorClass = "optional")

        DT::datatable(
            data = as.data.frame(genomicRanges),
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top")
    })

    # Genome annotation ----

    # The EnsDb object
    selectedPackage <- reactive({
        # Depends on input$annotationPackage
        validate(need(
            input$annotationPackage,
            Msgs[["annotationPackage"]]))

        validate(need(
            require(input$annotationPackage, character.only = TRUE),
            "Failed loading annotation package."
        ))

        return(getEdb(input$annotationPackage))
    })

    genomeSeqinfo <- reactive({

        selectedPackage <- selectedPackage()

        seqinfo(selectedPackage)
    })

    # EnsDb ----

    output$ensembl_organism <- renderUI({
        # Depends on selectedPackage
        edb <- selectedPackage()

        validate(need(edb, label = Msgs[["edb"]]))

        HTML(paste(
            tags$strong("Organism:"),
            organism(edb))
        )

    })

    output$ensembl_version <- renderUI({
        # Depends on selectedPackage
        edb <- selectedPackage()

        validate(need(edb, label = Msgs[["edb"]]))

        md <- metadata(edb)
        rownames(md) <- md$name

        HTML(paste(
            tags$strong("Ensembl version:"),
            md["ensembl_version", "value"])
        )
    })

    output$ensembl_genome <- renderUI({
        # Depends on selectedPackage
        edb <- selectedPackage()

        validate(need(edb, label = Msgs[["edb"]]))

        md <- metadata(edb)
        rownames(md) <- md$name

        HTML(paste(
            tags$strong("Genome build:"),
            md["genome_build", "value"])
        )
    })

    queryGenes <- reactive({
        # Depends on: selectedPackage, ...
        #   input$ensDb.type, input$ensDb.condition, input$ensDb.value

        edb <- selectedPackage()

        if (input$ensDb.value == "")
            return(NULL)

        ensDbFilter = EnsDbFilter(
            type = input$ensDb.type,
            condition = input$ensDb.condition,
            value = input$ensDb.value)

        res <- genes(
            edb,
            filter = ensDbFilter)

        return(res)
    })

    output$ensDb.Genes <- DT::renderDataTable({
        # Depends on: queryGenes

        queryGenes <- queryGenes()

        validate(need(
            length(queryGenes) > 0,
            "No genomic region to show."),
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
    #     edb <- selectedPackage()
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
    #     edb <- selectedPackage()
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
            need(input$altGenotypes, Msgs[["altGenotypes"]])
        )

        return(new(
            Class = "TVTBparam",
            genos = list(
                REF = input$refGenotypes,
                HET = input$hetGenotypes,
                ALT = input$altGenotypes),
            aaf = "AAF",
            maf = "MAF",
            vep = input$vepKey,
            bp = bpParam
        ))
    })

    output$TVTBsettings <- renderPrint({
        return(tparam())
    })

    output$generalSettings <- renderPrint({
        return(list(
            phenoFile = RV[["phenoFile"]],
            regionInputMode = input$regionInputMode,
            "bedFile()" = bedFile(),
            ucscRegions = input$ucscRegions,
            ensDb.type = input$ensDb.type,
            ensDb.condition = input$ensDb.condition,
            ensDb.value = input$ensDb.value,
            vcfinputMode = input$vcfInputMode,
            "singleVcf()" = singleVcf(),
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
            genoFirstCol = input$genoFirstCol,
            vepTVBP = input$vepTVBP,
            phenoTVBP = input$phenoTVBP,
            unique2phenoTVBP = input$unique2phenoTVBP,
            vepFacetKeyTVBP = input$vepFacetKeyTVBP,
            vepFacetsTVBP = input$vepFacetsTVBP,
            stackedPercentageTVBP = input$stackedPercentageTVBP
        ))
    })

    output$advancedSettings <- renderPrint({
        return(list(
            legendTextSizeTVBP = input$legendTextSizeTVBP,
            xAxisAngleTVBP = input$xAxisAngleTVBP,
            xAxisSizeTVBP = input$xAxisSizeTVBP,
            xAxisHjustTVBP = input$xAxisHjustTVBP,
            xAxisVjustTVBP = input$xAxisVjustTVBP,
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

    singleVcf <- reactive({

        if (input$selectVcf > 0){
            selected <- tryCatch(
                file.choose(),
                error = function(err){
                    warning(geterrmessage())
                    return(NULL)
                })
        } else {
            return(NULL)
        }

        validate(
            need(
                grepl(
                    pattern = ".*\\.vcf(\\.gz)?$",
                    x = selected,
                    ignore.case = TRUE),
                "File is not *.vcf.gz"))

        tbiChrVcf <- paste(selected, "tbi", sep = ".")
        validate(
            need(
                file.exists(tbiChrVcf),
                sprintf("Tabix index file does not exist: %s", tbiChrVcf)))

        return(selected)
    })

    output$selectedVcf <- renderText({

        singleVcf <- singleVcf()

        validate(need(singleVcf, Msgs[["singleVcf"]]))

        singleVcf
    })

    # Define ScanVcfParam ----

    # Identify available fields in header of (first) VCF file
    observe({
        # Depends on input mode, but also change of VCF files
        vcfInputMode <- input$vcfInputMode

        vcfHeader <- switch (
            vcfInputMode,
            SingleVcf = {

                singleVcf <- singleVcf()
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
                x = c(rownames(info(vcfHeader))),
                invert = TRUE,
                value = TRUE)
            # All keys except "GT" (required)
            RV[["genoKeys"]] <- grep(
                pattern = "GT",
                x = c(rownames(geno(vcfHeader))),
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
            genomicRanges <- genomicRanges()
            tparam <- tparam()
            vepKey <- input$vepKey
            genomeSeqinfo <- genomeSeqinfo()
            yieldSize <- input$yieldSize
            infoKeys <- c(input$vcfInfoKeys, input$vepKey)
            genoKeys <- input$vcfFormatKeys
            if (length(genoKeys) == 0)
                genoKeys <- "GT" # Mandatory

            if (is.null(RV[["phenoFile"]])) {
                phenotypes <- DataFrame()
            } else {
                phenotypes <- phenotypes()
            }

            validate(
                need(vepKey, label = Msgs[["vepKey"]]),
                need(genomeSeqinfo, Msgs[["genomeSeqinfo"]])
            )

            shiny::withProgress(
                min = 0, max = 3, value = 1,
                message = "Progress", detail = Tracking[["preprocessing"]],
                {

                    # NOTE: ALL FIXED fields imported
                    svp <- ScanVcfParam(
                        info = infoKeys,
                        geno = genoKeys
                    )

                    # Only import samples with phenotype information
                    if (nrow(phenotypes) > 0)
                        vcfSamples(svp) <- rownames(phenotypes)

                    # Only import variants in targeted regions
                    if (length(genomicRanges) > 0)
                        vcfWhich(svp) <- reduce(genomicRanges) # optimise import


                    # # Timing
                    # t1 <- Sys.time()

                    vcf <- switch (
                        vcfInputMode,
                        SingleVcf = {

                            isolate({singleVcf <- singleVcf()})
                            validate(need(singleVcf, Msgs[["singleVcf"]]))

                            shiny::incProgress(1, detail = Tracking[["singleVcf"]])

                            tryParseSingleVcf(
                                file = singleVcf,
                                svp = svp,
                                yieldSize = yieldSize
                            )
                        },

                        OnePerChr = {

                            shiny::incProgress(1, detail = Tracking[["multiVcfs"]])

                            isolate({
                                vcfFolder <- input$vcfFolder
                                vcfPattern <- input$vcfPattern
                            })
                            validate(need(singleVcf, Msgs[["singleVcf"]]))

                            tryParseMultipleVcf(
                                folder = vcfFolder,
                                pattern = vcfPattern,
                                svp = svp,
                                yieldSize = yieldSize,
                                BPPARAM = bp(tparam)
                            )
                        }
                    )

                    shiny::incProgress(1, detail = Tracking[["postprocessing"]])

                    validate(need(
                        length(vcf) > 0,
                        "No variant found in BED region(s)"))

                    if (nrow(phenotypes) == 0)
                        phenotypes <- DataFrame(row.names = colnames(vcf))

                    colData(vcf) <- phenotypes
                })

            RV[["vcf"]] <- vcf

        }
    )

    output$vcfSummary <- renderUI({
        # Summary of raw variants
        vcf <- RV[["vcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        return(tagList(
                code(nrow(vcf)), "bi-allelic records and",
                code(ncol(colData(vcf))), "phenotypes",
                "in", code(ncol(vcf)), "samples"
        ))
    })

    output$vcfCols <- renderUI({

        vcf <- RV[["filteredVcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        colChoices <- c(colnames(mcols(vcf)))

        selectInput(
            "vcfCols", "Meta-columns",
            choices = colChoices,
            selected = c(colChoices[1:min(5, length(colChoices))]),
            multiple = TRUE
        )
    })

    output$vcfRowRanges <- DT::renderDataTable({

        vcf <- RV[["filteredVcf"]]

        # Give time to initialise widget
        req(input$vcfCols)

        cols <- which(colnames(mcols(vcf)) %in% input$vcfCols)

        displayedTable <- cbind(
            rownames = rownames(vcf),
            as.data.frame(rowRanges(vcf)[, cols], row.names = NULL)
        )

        DT::datatable(
            data = displayedTable,
            rownames = FALSE,
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top")
    })

    output$vcfInfoCols <- renderUI({

        vcf <- RV[["filteredVcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        validate(need(
            ncol(info(vcf)) > 1, # VEP field are an implicite field
            "No INFO field available."),
            errorClass = "optional"
        )

        validate(need(ncol(info(vcf)) > 0, Msgs[["importVariants"]])) # TODO

        # All columns except the VEP predictions
        colChoices <- grep(
            pattern = input$vepKey,
            x = colnames(info(vcf)),
            invert = TRUE,
            value = TRUE)

        selectInput(
            "vcfInfoCols", "Meta-columns",
            choices = colChoices,
            selected = colChoices[1:min(5, length(colChoices))],
            multiple = TRUE
        )
    })

    output$vcfInfo <- DT::renderDataTable({
        req(input$vcfInfoCols)
        vcf <- RV[["filteredVcf"]]

        validate(need(
            ncol(info(RV[["vcf"]])) > 1, # At least one non-VEP column
            "No INFO data imported."
        ),
        errorClass = "optional")

        validate(need(
            length(input$vcfInfoCols) > 0,
            "No INFO key selected."
        ))

        cols <- which(colnames(info(vcf)) %in% input$vcfInfoCols)

        DT::datatable(
            data = cbind(
                rownames = rownames(vcf),
                as.data.frame(
                    info(vcf)[, cols, drop = FALSE],
                    row.names = NULL)
            ),
            rownames = FALSE,
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top")
    })

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
                Class = input$newFilterClass,
                exprs = list(quickFix),
                active = input$newFilterActive)

            if (class(newFilter) == "VcfVepRules")
                vep(newFilter) <- input$vepKey

            return(newFilter)},
            # warning = function(w) NULL,
            error = function(e) NULL
        ))

    })

    newFilterTestResults <- reactive({
        # Only triggered by addNewFilter (> 0)
        req(input$addNewFilter)

        isolate({
            newFilter <- newVcfFilter()
            vcf <- head(RV[["vcf"]]) # Speed up testing!
        })

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

        if (testResult)
            return()
        else
            return(
                strong(tags$span(style="color:red", RV[["newFilterStatus"]])))
    })

    observeEvent(input$addNewFilter, {
        newFilter <- newVcfFilter()
        testResult <- newFilterTestResults()

        # Only add new filter if valid
        validate(need(testResult == TRUE, "Invalid VCF filter"))

        names(newFilter) <- paste0("rule", input$addNewFilter)

        newRules <- VcfFilterRules(
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
                    shiny::column(
                        width = 1,
                        code(type(vcfFilters)[filterIndex])
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
                            value = active(vcfFilters)[filterIndex])
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

    # Observe checboxes defining active status
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
        active(vcfRules) <- as.logical(rulesStatus)[idx]
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
            RV[["vcfFilters"]] <- VcfFilterRules(vcfRules[-removeIdx])
    })

    observe({
        # Depend on the actionButton *and* the vcf object
        # makes filtered variants synonym to raw variants in absence of filters
        input$filterVariants
        vcf <- RV[["vcf"]]

        # Do not update when filters change, wait for actionButton/vcf
        isolate({vcfFilters <- RV[["vcfFilters"]]})

        # Store the filtered VCF in the reativeValues
        RV[["filteredVcf"]] <- subsetByFilter(x = vcf, filter = vcfFilters)
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
            code(ncol(colData(filteredVcf))), "phenotypes in",
            code(ncol(filteredVcf)), "samples filtered"
        ))
    })

    output$filtersSummary <- renderUI({
        vcfFilters <- RV[["vcfFilters"]]

        return(tagList(
            code(sum(active(vcfFilters))), "active filters",
            "( of", code(length(vcfFilters)), ")"
            ))

    })

    # Parse genotypes ----

    output$genotypeStructure <- renderPrint({

        vcf <- RV[["vcf"]]

        validate(need(vcf, Msgs[["importVariants"]]))

        genotypes <- geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

        str(genotypes)
    })

    output$genoNumCols <- renderUI({

        vcf <- RV[["filteredVcf"]]

        validate(need(vcf, Msgs[["importVariants"]]))

        genotypes <- geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

        sliderInput(
            "genoNumCols", "Number of columns (samples)",
            value = min(10, ncol(genotypes)),
            min = 2,
            max = min(50, ncol(genotypes)),
            step = 1)
    })

    output$genoFirstCol <- renderUI({

        req(input$genoNumCols)

        vcf <- RV[["filteredVcf"]]

        validate(need(vcf, Msgs[["importVariants"]]))

        genotypes <- geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

        isolate({newValue <- max(1, input$genoFirstCol, na.rm = TRUE)})

        sliderInput(
            "genoFirstCol", "First column (sample)",
            value = newValue,
            min = 1,
            max = max(c(10, ncol(genotypes) - input$genoNumCols + 1)),
            step = 1)
    })

    output$genoNumRows <- renderUI({

        vcf <- RV[["filteredVcf"]]

        validate(need(vcf, Msgs[["importVariants"]]))

        genotypes <- geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

        sliderInput(
            "genoNumRows", "Number of rows (variants)",
            value = min(10, nrow(genotypes)),
            min = 2,
            max = min(100, nrow(genotypes)),
            step = 1)
    })

    output$genoFirstRow <- renderUI({

        req(input$genoNumRows)

        vcf <- RV[["filteredVcf"]]
        validate(need(vcf, Msgs[["importVariants"]]))

        genotypes <- geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

        isolate({newValue <- max(1, input$genoFirstRow, na.rm = TRUE)})

        sliderInput(
            "genoFirstRow", "First row (variant)",
            value = newValue,
            min = 1,
            max = nrow(genotypes) - input$genoNumRows + 1,
            step = 1)
    })

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

        genotypes <- geno(vcf)[["GT"]]
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

    # Show information about consequence field
    output$vepStructure <- renderPrint({

        vcf <- RV[["vcf"]]

        validate(
            need(vcf, Msgs[["importVariants"]]),
            need(
                input$vepKey %in% colnames(info(vcf)),
                Msgs[["vepKeyNotFound"]])
        )

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        validate(need(csq, Msgs[["csq"]]))

        str(csq)
    })

    output$vepCols <- renderUI({

        vcf <- RV[["filteredVcf"]]
        # First make sure vcf exists
        validate(need(vcf, Msgs[["importVariants"]]))
        # If it exists, check that vepKey exist in INFO fields
        validate(need(
            input$vepKey %in% colnames(info(vcf)),
            Msgs[["vepKeyNotFound"]]
        ))

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        vepMcols <- mcols(csq)

        selectInput(
            "vepCols", "Meta-columns",
            choices = colnames(vepMcols),
            selected = colnames(vepMcols)[1:5],
            multiple = TRUE)
    })

    output$vcfVep <- DT::renderDataTable({
        req(input$vepCols)

        vcf <- RV[["filteredVcf"]]

        validate(
            need(vcf, Msgs[["importVariants"]]),
            need(
                input$vepKey %in% colnames(info(vcf)),
                Msgs[["vepKeyNotFound"]])
        )

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        vepMcols <- mcols(csq)

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

    # When variants are updated, update the choices for plotting
    observeEvent(RV[["filteredVcf"]], {

        vcf <- RV[["filteredVcf"]]

        validate(
            need(vcf, Msgs[["importVariants"]]),
            need(
                input$vepKey %in% colnames(info(vcf)),
                Msgs[["vepKeyNotFound"]])
        )

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        validate(need(csq, Msgs[["csq"]]))

        vepMcols <- mcols(csq)

        vepKey.choices <- colnames(vepMcols)

        # Initialise by order of preference if present:
        # Current selection > Consequence > Impact > First field
        ## TVBP (Tabulate VEP by phenotype)
        vep.TVBPdefault <- max(
            1,
            head(x = na.omit(
                match(
                    x = c(input$vepTVBP, "Consequence","IMPACT"),
                    table = vepKey.choices
                    )
                ), n = 1)
            )

        updateSelectInput(
            session, "vepTVBP",
            choices = vepKey.choices,
            selected = vepKey.choices[vep.TVBPdefault])

        ## DVBP (Density VEP by phenotype)
        vep.DVBPdefault <- max(
            1,
            head(x = na.omit(
                match(
                    x = c(input$vepTVBP, "CADD_PHRED","Protein_position"),
                    table = vepKey.choices
                )
            ), n = 1)
        )

        updateSelectInput(
            session, "vepDVBP",
            choices = vepKey.choices,
            selected = vepKey.choices[vep.TVBPdefault])

        vepKey.choices <- c("None", vepKey.choices)

        # Initialise by order of preference if present:
        # Current selection > None
        ## TVBP
        vepFacetKey.TVBPselected <- max(
            1,
            head(x = na.omit(
                match(
                    x = c(input$vepFacetKeyTVBP),
                    table = vepKey.choices
                )
            ), n = 1)
        )

        updateSelectInput(
            session, "vepFacetKeyTVBP",
            choices = vepKey.choices,
            selected = vepKey.choices[vepFacetKey.TVBPselected])

        ## DVBP
        vepFacetKey.DVBPselected <- max(
            1,
            head(x = na.omit(
                match(
                    x = c(input$vepFacetKeyDVBP),
                    table = vepKey.choices
                )
            ), n = 1)
        )

        updateSelectInput(
            session, "vepFacetKeyDVBP",
            choices = vepKey.choices,
            selected = vepKey.choices[vepFacetKey.DVBPselected])
    })

    # Tabulate VEP by phenotype (TVBP) ----

    # Phenotype to summarise VEP: or NULL if "None"
    varVepPlotPhenoTVBP <- reactive({

        phenotype <- input$phenoTVBP
        # Give time to initialise the widget
        req(phenotype)

        if (phenotype == "None"){
            return("Phenotype")
        } else {
            return(phenotype)
        }
    })

    # Update list of VEP facets if the faceting variable changes
    observeEvent(eventExpr = input$vepFacetKeyTVBP, handlerExpr = {
        # Raw variants must exist
        validate(need(RV[["vcf"]], Msgs[["importVariants"]]))

        vcf <- RV[["filteredVcf"]]
        # Give time to filter variants
        req(vcf)

        validate(need(
            input$vepKey %in% colnames(info(vcf)),
            Msgs[["vepKeyNotFound"]]))
        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        validate(need(csq, Msgs[["csq"]]))

        vepMcols <- mcols(csq)

        # Initialise to avoid crash
        facetsSelected <- NULL
        if (input$vepFacetKeyTVBP != "None"){
            vepFacets.choices <- unique(vepMcols[,input$vepFacetKeyTVBP])
        } else {
            vepFacets.choices <- NA_character_
        }

        updateSelectInput(
            session, "vepFacetsTVBP",
            choices = vepFacets.choices,
            selected = vepFacets.choices)
    })

    # Base of the ggplot summarising counts of VEP
    TVBPggplot <- reactive({

        # Update only on button click
        validate(need(
            input$buttonTVBP,
            "Please click \"Apply\" to apply parameters."),
            errorClass = "optional")

        isolate({
            vepTVBP <- input$vepTVBP
            vepFacetKey <- input$vepFacetKeyTVBP
            vepFacets <- input$vepFacetsTVBP
            stackedPercentage <- input$stackedPercentageTVBP
            unique2pheno <- input$unique2phenoTVBP
        })

        withProgress(min = 0, max = 2, message = "Progress", {
            incProgress(1, detail = Tracking[["calculate"]])

            if (vepFacetKey != "None"){
                # req(vepFacets) # Give time to observeEvent
                validate(need(
                    length(vepFacets) > 0,
                    "The plot must contain at least one facet"))
            }

            # App does not display "NULL" choices
            if (vepFacetKey == "None"){
                vepFacetKey <- NULL
            }

            validate(need(RV[["vcf"]], Msgs[["importVariants"]]))

            vcf <- RV[["filteredVcf"]]

            validate(need(
                nrow(vcf) > 0,
                Msgs[["filterVcfEmpty"]]),
                errorClass = "optional")

            # Special case of phenotype "None"
            isolate({phenoTVBP <- input$phenoTVBP})
            if (phenoTVBP == "None"){
                colData(vcf)[,"Phenotype"] <- factor("All")
                plotPhenotype <- "Phenotype"
            } else {
                plotPhenotype <- phenoTVBP
            }

            isolate({varVepPlotPheno <- varVepPlotPhenoTVBP()})
            validate(
                need(vepTVBP, Msgs[["vepTVBP"]]),
                need(
                    is.logical(stackedPercentage),
                    Msgs[["stackedPercentageTVBP"]])
                )

            incProgress(1, detail = Tracking[["ggplot"]])

            gg <- tabulateVepByPhenotype(
                vcf = vcf,
                phenoCol = plotPhenotype,
                vepCol = vepTVBP,
                param = tparam(),
                unique = unique2pheno,
                facet = vepFacetKey,
                plot = TRUE,
                percentage = stackedPercentage
            )

            if (!is.null(vepFacetKey)){
                gg$data <- gg$data[
                    gg$data[,vepFacetKey] %in% vepFacets,]
                gg$data <- droplevels(gg$data)
            }

            return(gg)
        })

    })

    # Stacked bars classifying consequence of variants
    output$vepCountBarplot <- renderPlot({
        withProgress(min = 0, max = 2, message = "Progress", {
            incProgress(1, detail = Tracking[["calculate"]])

            # Update only on button click
            gg <- TVBPggplot()

            isolate({
                legendTextSize <- input$legendTextSizeTVBP
                xAxisAngle <- input$xAxisAngleTVBP
                xAxisHjust <- input$xAxisHjustTVBP
                xAxisVjust <- input$xAxisVjustTVBP
                xAxisSize <- input$xAxisSizeTVBP
                legend <- input$legendTVBP
            })

            message("Plotting predictions...")
            incProgress(1, detail = Tracking[["render"]])

            gg <- gg +
                theme(
                    axis.text = element_text(size = rel(1.5)),
                    axis.title = element_text(size = rel(1.5)),
                    legend.text = element_text(size = rel(legendTextSize)),
                    legend.title = element_text(size = rel(legendTextSize)),
                    axis.text.x = element_text(
                        angle = xAxisAngle,
                        hjust = xAxisHjust,
                        vjust = xAxisVjust,
                        size = rel(xAxisSize))
                ) +
                guides(fill = c("none", "legend")[legend + 1])

            return(gg)
        })
    })

    vepTableCount <- reactive({
        # Depeds on: vepTabulated

        vepTabulated <- TVBPggplot()$data

        vepTable <- as.data.frame(table(vepTabulated))

        return(vepTable)
    })

    # Table of VEP counts by phenotype level by VEP facet
    output$vepTableDecreasing <- DT::renderDataTable({

        vepTableCount <- vepTableCount()

        vepTableDecreasing <- arrange(vepTableCount, desc(Freq))

        DT::datatable(
            vepTableDecreasing,
            options = list(
                pageLength = 20,
                searching = TRUE),
            filter = "top",
            rownames = FALSE)
    })

    # Print the count of consequences in the area hovered.
    output$varVepCount <- renderUI({

        req(input$plotVarClass_hover)

        vepTableCount <- vepTableCount()

        hover <- input$plotVarClass_hover
        vepTVBP <- input$vepTVBP

        # If not hovered, show empty values
        if (is.null(hover)){
            x_lvl <- NULL
            y_lvl <- NULL
            countVarLevel <- NULL
        }
        else{
            # Identify the level of the phenotype (x axis)
            x_lvls <- levels(vepTableCount[,hover$mapping$x])
            x_lvl <- x_lvls[round(hover$x)]
            # Identify the level of the variant category (y axis)
            filters <- data.frame(
                pheno = vepTableCount[,hover$mapping$x] == x_lvl
            )
            # Identify data for the panel (vepFacetKey) hovered
            if (input$vepFacetKeyTVBP != "None"){

                filters$vepFacetKeyTVBP <- vepTableCount[
                    ,hover$mapping$panelvar1] ==
                    hover$panelvar1
            }
            # Extract VEP counts for the panel hovered:
            # matches x (phenotype level) and y (VEP level)
            y_lvls <- vepTableCount[
                apply(X = filters, MARGIN = 1, FUN = all),
                c(vepTVBP, "Freq")]

            # If in percentage mode, rescale the values
            if (input$stackedPercentageTVBP){
                y_lvls[,"Freq"] <- y_lvls[,"Freq"] / sum(y_lvls[,"Freq"])
            }

            # If out of range, show empty values
            if (
                input$plotVarClass_hover$y < 0 ||
                input$plotVarClass_hover$y > sum(y_lvls[,"Freq"])){
                y_lvl <- NULL
                countVarLevel <- NULL
            } else {
                y_lvls[,"cumsum"] <- cumsum(y_lvls[,"Freq"])
                y_lvls <- y_lvls[
                    y_lvls[,"cumsum"] >= input$plotVarClass_hover$y,]
                y_lvl <- as.character(y_lvls[1, vepTVBP])
                countVarLevel <- y_lvls[y_lvls[,vepTVBP] == y_lvl, "Freq"]
                countVarLevel <- ifelse(
                    input$stackedPercentageTVBP,
                    yes = sprintf(
                        "%.1f%%", 100 * countVarLevel),
                    no = as.character(countVarLevel))
            }
        }

        html1Start <- "<ul>"
        html2Facet <- ifelse(
            input$vepFacetKeyTVBP == "None",
            yes = "",
            no = sprintf(
                "<li>%s : <code>%s</code></li>",
                hover$mapping$panelvar1,
                hover$panelvar1
            ))
        html3Pheno <- sprintf(
            "<li>%s : <code>%s</code></li>",
            hover$mapping$x,
            x_lvl
        )

        if(is.null(y_lvl)){
            html4Vep <- ""
            html5Value <- ""
        } else {
            html4Vep <- sprintf(
                "<li>%s : <code>%s</code></li>",
                vepTVBP,
                y_lvl
            )
            html5Value <- ifelse(
                input$stackedPercentageTVBP,
                yes = sprintf(
                    "<li>Percentage : <code>%s</code></li>",
                    countVarLevel,
                    y_lvl),
                no = sprintf(
                    "<li>Count : <code>%s</code></li>",
                    countVarLevel)
            )
        }

        html6End <- "</ul>"

        HTML(paste(
            html1Start, html2Facet, html3Pheno, html4Vep, html5Value, html6End,
            sep = ""))
    })

    # Density VEP by phenotype (DVBP) ----

    # Phenotype to summarise VEP density: or NULL if "None"
    varVepPlotPhenoDVBP <- reactive({

        phenotype <- input$phenoDVBP
        # Give time to initialise the widget
        req(phenotype)

        if (phenotype == "None"){
            return("Phenotype")
        } else {
            return(phenotype)
        }
    })

    # Update list of VEP facets if the faceting variable changes
    observeEvent(eventExpr = input$vepFacetKeyDVBP, handlerExpr = {
        # Raw variants must exist
        validate(need(RV[["vcf"]], Msgs[["importVariants"]]))

        vcf <- RV[["filteredVcf"]]
        # Give time to filter variants
        req(vcf)

        validate(need(
            input$vepKey %in% colnames(info(vcf)),
            Msgs[["vepKeyNotFound"]]))
        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        validate(need(csq, Msgs[["csq"]]))

        vepMcols <- mcols(csq)

        # Select a subset of VEP facets
        if (input$vepFacetKeyDVBP != "None"){
            # List all choices of facets
            vepFacets.choices <- unique(vepMcols[,input$vepFacetKeyDVBP])
        } else {
            vepFacets.choices <- NA_character_
        }

        updateSelectInput(
            session, "vepFacetsDVBP",
            choices = vepFacets.choices,
            selected = vepFacets.choices)
    })

    # Base of the ggplot summarising counts of VEP
    DVBPggplot <- reactive({

        # Update only on button click
        validate(need(
            input$buttonDVBP,
            "Please click \"Apply\" to apply parameters."),
            errorClass = "optional")

        isolate({
            vepDVBP <- input$vepDVBP
            vepFacetKey <- input$vepFacetKeyDVBP
            vepFacets <- input$vepFacetsDVBP
            patternDVBP <- input$patternDVBP
            unique2pheno <- input$unique2phenoDVBP
            layerDVBP <- input$layerDVBP
        })

        withProgress(min = 0, max = 2, message = "Progress", {
            incProgress(1, detail = Tracking[["calculate"]])

            if (vepFacetKey != "None"){
                # req(vepFacets) # Give time to observeEvent
                validate(need(
                    length(vepFacets) > 0,
                    "The plot must contain at least one facet"))
            }

            # App does not display "NULL" choices
            if (vepFacetKey == "None"){
                vepFacetKey <- NULL
            }

            validate(need(RV[["vcf"]], Msgs[["importVariants"]]))

            vcf <- RV[["filteredVcf"]]

            validate(need(
                nrow(vcf) > 0,
                Msgs[["filterVcfEmpty"]]),
                errorClass = "optional")

            # Special case of phenotype "None"
            isolate({phenoDVBP <- input$phenoDVBP})
            if (input$phenoDVBP == "None"){
                colData(vcf)[,"Phenotype"] <- factor("All")
                plotPhenotype <- "Phenotype"
            } else {
                plotPhenotype <- phenoDVBP
            }

            isolate({varVepPlotPheno <- varVepPlotPhenoDVBP()})
            validate(need(vepDVBP, Msgs[["vepDVBP"]]))

            # Convert to default value for the method
            if (patternDVBP == "")
                patternDVBP <- NULL

            incProgress(1, detail = Tracking[["ggplot"]])

            gg <- densityVepByPhenotype(
                vcf = vcf,
                phenoCol = plotPhenotype,
                vepCol = vepDVBP,
                param = tparam(),
                unique = unique2pheno,
                facet = vepFacetKey,
                plot = TRUE,
                pattern = patternDVBP,
                layer = layerDVBP
            )

            if (!is.null(vepFacetKey)){
                gg$data <- gg$data[
                    gg$data[,vepFacetKey] %in% vepFacets,]
                gg$data <- droplevels(gg$data)
            }

            return(gg)
        })

    })

    # Stacked bars classifying consequence of variants
    output$vepDensityPlot <- renderPlot({
        withProgress(min = 0, max = 2, message = "Progress", {
            incProgress(1, detail = Tracking[["calculate"]])

            # Update only on button click
            gg <- DVBPggplot()

            isolate({
                legendTextSize <- input$legendTextSizeDVBP
                xAxisSize <- input$xAxisSizeDVBP
                yAxisSize <- input$yAxisSizeDVBP
                legend <- input$legendTVBP
            })

            message("Plotting predictions...")
            incProgress(1, detail = Tracking[["render"]])

            gg <- gg +
                theme(
                    axis.text = element_text(size = rel(1.5)),
                    axis.title = element_text(size = rel(1.5)),
                    legend.text = element_text(size = rel(legendTextSize)),
                    legend.title = element_text(size = rel(legendTextSize)),
                    axis.text.x = element_text(size = rel(xAxisSize)),
                    axis.text.y = element_text(size = rel(yAxisSize))
                ) +
                guides(fill = c("none", "legend")[legend + 1])

            return(gg)
        })
    })

    # Genotype heatmap ----

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

            genotypes <- geno(vcf)[["GT"]]
            validate(need(genotypes, Msgs[["genotypes"]]))

            validGenotypes <- unlist(genos(tparam()))

            genos.long <- melt(genotypes, value.name = "Genotype")

            colnames(genos.long)[1:2] <- c("Variants", "Samples")

            genos.long[,"Genotype"] <- factor(
                x = genos.long[,"Genotype"],
                levels = unlist(genos(tparam()))
            )

            incProgress(1, detail = Tracking[["ggplot"]])
            gg <- ggplot(
                data = genos.long,
                mapping = aes(Samples, Variants)
            ) +
                geom_tile(aes(fill = Genotype, colour = Genotype)) +
                scale_fill_discrete() +
                theme(
                    axis.text = element_blank(),
                    axis.title = element_text(size = rel(1.5)),
                    axis.ticks = element_blank()
                )

            message("Plotting heatmap of genotypes...")
            incProgress(1, detail = Tracking[["render"]])
            gg
        })
    })

    # Parallel computing ----

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
                    value = multicoreWorkers(),
                    min = 1, max = detectCores())
            },
            SnowParam = {
                updateSelectInput(
                    session, "bpType",
                    choices = c("SOCK", "MPI", "FORK"))
                updateNumericInput(
                    session, "bpCores",
                    value = snowWorkers(),
                    min = 1, max = detectCores())
            })
    })

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
            SerialParam = SerialParam(),
            MulticoreParam = MulticoreParam(workers = input$bpCores),
            SnowParam = SnowParam(
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
