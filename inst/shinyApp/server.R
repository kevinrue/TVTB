
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

source(file.path(
    system.file(package = "tSVE"),
    "shinyApp",
    "serverRoutines.R"))

shinyServer(function(input, output, clientData, session) {

    # Genome annotation ----

    selectedPackage <- reactive({

        validate(need(
            input$annotationPackage,
            label = Msgs[["annotationPackage"]]))
        validate(
            need(
                require(input$annotationPackage, character.only = TRUE),
                message = "Cannot load EnsDb package. Invalid?"
            ))

        return (getEdb(input$annotationPackage))
    })

    # EnsDb ----

    output$ensembl_organism <- renderUI({
        edb <- selectedPackage()

        validate(need(edb, label = Msgs[["edb"]]))

        HTML(paste(
            tags$strong("Organism:"),
            organism(edb))
        )

    })

    output$ensembl_version <- renderUI({
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

        edb <- selectedPackage()

        validate(
            need(input$ensDb.type, label = Msgs[["ensDb.type"]]),
            need(input$ensDb.condition, label = Msgs[["ensDb.condition"]]),
            need(input$ensDb.values, label = Msgs[["ensDb.values"]])
        )

        ensDbFilter = EnsDbFilter(
            type = input$ensDb.type,
            condition = input$ensDb.condition,
            value = input$ensDb.values)

        validate(need(ensDbFilter, Msgs[["ensDbFilter"]]))

        res <- genes(
            edb,
            filter = ensDbFilter,
            return.type = "data.frame")

        return(res)
    })

    output$ensDb.Genes <- DT::renderDataTable({

        queryGenes <- queryGenes()

        validate(need(queryGenes, Msgs[["queryGenes"]]))

        DT::datatable(
            data = queryGenes,
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

    genomeSeqinfo <- reactive({

        edb <- getEdb(input$annotationPackage)
        validate(need(edb, label = "annotationPackage"))

        seqinfo(edb)
    })

    # Session information ----

    tParam <- reactive({

        bpParam <- bpParam()

        validate(
            need(input$refGenotypes, label = Msgs[["refGenotypes"]]),
            need(input$hetGenotypes, label = Msgs[["hetGenotypes"]]),
            need(input$altGenotypes, label = Msgs[["altGenotypes"]]),
            need(bpParam, Msgs[["bpParam"]])
        )

        new(
            Class = "tSVEParam",
            genos = list(
                REF = input$refGenotypes,
                HET = input$hetGenotypes,
                ALT = input$altGenotypes),
            aaf = "AAF", maf = "MAF",
            vep = input$csqField,
            bp = bpParam
        )
    })

    output$tSVESettings <- renderPrint({
        tParam()
    })

    output$advancedSettings <- renderPrint({
        list(
            csqPhenoUnique = input$csqPhenoUnique,
            annotationPackage = input$annotationPackage,
            legendText = input$legendText,
            regionInputMode = input$regionInputMode,
            stackedPercentage = input$stackedPercentage,
            vcfInputMode = input$vcfInputMode,
            xAxisHjust = input$xAxisHjust,
            xAxisVjust = input$xAxisVjust,
            xAxisAngle = input$xAxisAngle,
            xAxisSize = input$xAxisSize
        )
    })

    output$sessionInfo <- renderPrint({
        sessionInfo()
    })

    # VCF folder  ----

    # List of files detected
    vcfContent <- reactive({

        validate(need(input$vcfFolder, label = Msgs[["vcfFolder"]]))

        validate(need(
            dir.exists(input$vcfFolder),
            "Invalid VCF folder"))

        list.files(path = input$vcfFolder)
    })

    # Easier to read for user
    output$vcfContent <- DT::renderDataTable({

        vcfContent <- vcfContent()

        validate(need(vcfContent, Msgs[["vcfContent"]]))

        DT::datatable(
            data = data.frame(Filename = vcfContent),
            rownames = FALSE)
    })

    ## List VCF files detected
    vcfFiles <- reactive({

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

        HTML(paste(
            length(vcfFiles),
            "VCF.GZ file(s) detected.</br>",
            "[",
            basename(head(x = vcfFiles, n = 1)),
            " , ... ]"
        ))
    })

    # BED file import ----

    bedFile <- reactive({

        validate(need(input$selectVcf, Msgs[["selectBed"]]))

        selected <- tryCatch(
            file.choose(),
            error = function(err){
                warning(geterrmessage())
            })

        # Catch the error message if user cancelled file seletion pop-up
        validate(need(
            !grepl(pattern = selected, x = Msgs[["fileChooseCancelled"]]),
            Msgs[["fileChooseCancelled"]]))

        selected
    })

    # Import BED records, format as GRanges
    bedRecords <- reactive({

        switch (input$regionInputMode,
                bed = {
                    validate(need(input$bedFile, label = Msgs[["bedFile"]]))

                    message("Importing BED records ...")
                    rawData <- tryParseBed(input$bedFile$datapath)
                },
                ucsc = {

                    validate(need(
                        input$ucscRegions,
                        label = Msgs[["ucscRegions"]]))
                    inputTrimmed <- gsub(
                        pattern = "chr|,| ",
                        replacement = "",
                        x = input$ucscRegions)

                    validate(need(
                        grepl(
                            pattern =
                                "([[:alnum:]]*:[[:digit:]]*-[[:digit:]]*;?)",
                            x = inputTrimmed),
                        label = Msgs[["ucscRegions"]]))
                    rawData <- tryCatch({
                        inputSplit <- strsplit(
                            x = inputTrimmed, split = ";")[[1]]

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

                        data.frame(V1 = chrs, V2 = starts, V3 = ends)
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
                    validate(need(
                        input$ensDb.values,
                        label = Msgs[["ensDb.values"]]))

                    edb <- selectedPackage()

                    ensDbFilter = EnsDbFilter(
                        type = input$ensDb.type,
                        condition = input$ensDb.condition,
                        value = input$ensDb.values)

                    validate(need(ensDbFilter, Msgs[["ensDbFilter"]]))

                    rawData <- genes(
                        edb,
                        filter = ensDbFilter,
                        return.type = "GRanges")
                }
        )

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
                                    seqlevels(EnsDb.Hsapiens.v75)),
                            "First column contains invalid chromosomes"),
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
                    validate(need(
                        length(rawData) > 0,
                        "At least one region must be specified"
                    ))

                    return(rawData)
                }
        )
    })

    # How many BED records detected, show first one.
    output$bedFileSummary <- renderUI({

        bedRecords <- bedRecords()

        validate(need(bedRecords, Msgs[["bedRecords"]]))

        HTML(paste(
            length(bedRecords),
            "BED record(s) detected.</br>",
            "[",
            as.character(head(x = seqnames(bedRecords), n = 1)),
            ":",
            as.character(head(x = start(bedRecords), n = 1)),
            "-",
            as.character(head(x = end(bedRecords), n = 1)),
            # "{",
            # as.character(head(x = names(bedRecords), n = 1)),
            # "}",
            " , ... ]"
        ))
    })

    # Show BED records
    output$bedStructure <- renderPrint({

        bedRecords <- bedRecords()

        validate(need(bedRecords, Msgs[["bedRecords"]]))

        str(bedRecords)
    })

    output$bedSample <- DT::renderDataTable({

        bedRecords <- bedRecords()

        validate(need(bedRecords, Msgs[["bedRecords"]]))

        DT::datatable(
            data = as.data.frame(bedRecords),
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top")
    })

    # Select single VCF ----

    singleVcf <- reactive({

        validate(need(input$selectVcf, Msgs[["selectVcf"]]))

        selected <- tryCatch(
            file.choose(),
            error = function(err){
                warning(geterrmessage())
            })

        # Catch the error message if user cancelled file seletion pop-up
        validate(need(
            !grepl(pattern = selected, x = Msgs[["fileChooseCancelled"]]),
            Msgs[["fileChooseCancelled"]]))

        validate(
            need(
                grepl(
                    pattern = ".*\\.vcf\\.gz$",
                    x = selected,
                    ignore.case = TRUE),
                "File is not *.vcf.gz"))

        tbiChrVcf <- paste(selected, "tbi", sep = ".")
        validate(
            need(
                file.exists(tbiChrVcf),
                sprintf("Tabixed VCF file not found: %s", tbiChrVcf)))

        selected
    })

    output$selectedVcf <- renderText({

        singleVcf <- singleVcf()

        validate(need(singleVcf, Msgs[["singleVcf"]]))

        singleVcf
    })

    # Import VCF information ----

    mafRange <- reactive({

        validate(
            need(input$maf.min, label = Msgs[["maf.min"]]),
            need(input$maf.max, label = Msgs[["maf.max"]])
        )
        validate(need(
            input$maf.min < input$maf.max,
            "MAF: Minimum must be lower than maximum!")
        )

        range(input$maf.min, input$maf.max)

    })

    output$mafRange <- renderUI({

        mafRange <- mafRange()

        validate(need(mafRange, Msgs[["mafRange"]]))

        HTML(paste(
            "Minor allele frequencies:</br>",
            "[ ",
            min(mafRange()),
            " - ",
            max(mafRange()),
            " ]"
        ))
    })

    # Import Vcf object
    vcf <- reactive({

        bedRecords <- bedRecords()
        validate(need(bedRecords, Msgs[["bedRecords"]]))

        # Only proceed if the variants should be refreshed
        validate(need(
            input$refreshVariants == "Active",
            Msgs[["refreshVariants"]]))

        if (all(mafRange() == c(0, 1))){
            mafRange <- NULL
        } else {
            mafRange <- mafRange()
        }

        withProgress(
            min = 0, max = 3, value = 0, message = "Progress", {

                isolate({bpParam <- bpParam()})

                csqField <- input$csqField
                genomeSeqinfo <- genomeSeqinfo()
                validate(
                    need(csqField, label = Msgs[["csqField"]]),
                    need(genomeSeqinfo, Msgs[["genomeSeqinfo"]])
                )

                incProgress(1, detail = Tracking[["checkPhenotypes"]])
                if (is.null(input$phenoFile)) {
                    phenotypes <- DataFrame()
                } else {
                    phenotypes <- phenotypes()
                }
                # Timing
                t1 <- Sys.time()

                switch (
                    input$vcfInputMode,
                    SingleVcf = {

                        singleVcf <- singleVcf()
                        validate(need(singleVcf, Msgs[["singleVcf"]]))
                        incProgress(2, detail = Tracking[["singleVcf"]])

                        vcf <- preprocessVariants(
                            file = TabixFile(file = singleVcf),
                            regions = bedRecords,
                            param = tParam(),
                            phenos = phenotypes)
                    },
                    OnePerChr = {
                        bedChrs <- levels(seqnames(bedRecords))
                        vcfFiles <- vcfFiles()

                        validate(
                            need(bedChrs, Msgs[["bedChrs"]]),
                            need(
                                input$vcfPattern,
                                label = Msgs[["vcfPattern"]]))
                        validate(need(
                            grepl(pattern = "%s", x = input$vcfPattern),
                            "VCF file pattern must contain \"%s\""
                        ))

                        incProgress(1, detail = Tracking[["multiPaths"]])
                        vcfPaths <- tryCatch(
                            sapply(
                                X = bedChrs,
                                FUN = chr2file,
                                pattern = input$vcfPattern,
                                folder = input$vcfFolder),
                            error = function(err){
                                warning(geterrmessage())
                                return(NULL)
                            },
                            warning = function(warn){
                                warning(warn)
                                return(NULL)
                            })

                        validate(need(vcfPaths, Msgs[["vcfPaths"]]))

                        # Wrap in tryCatch if not tabixed
                        tfl <- tryCatch({TabixFileList(vcfPaths)},
                                        error = function(err){
                                            warning(geterrmessage())
                                            return(NULL)
                                        },
                                        warning = function(warn){
                                            warning(warn)
                                            return(NULL)
                                        })

                        if (length(bedChrs) == 1){
                            incProgress(2, detail = Tracking[["singleVcf"]])

                            vcf <- preprocessVariants(
                                file = tfl[[1]],
                                regions = bedRecords,
                                param = tParam(),
                                phenos = phenotypes)

                        } else {
                            incProgress(1, detail = Tracking[["multiVcfs"]])
                            vcfs <- tryCatch({
                                bpmapply(
                                    FUN = preprocessVariants,
                                    chr = bedChrs,
                                    file = tfl,
                                    MoreArgs = list(
                                        regions = bedRecords,
                                        param = tParam(),
                                        phenos = phenotypes),
                                    BPPARAM = bpParam)
                            },
                            error = function(err){
                                warning(geterrmessage())
                                return(NULL)
                            },
                            warning = function(warn){
                                warning(warn)
                                return(NULL)
                            })
                            validate(need(vcfs, Msgs[["vcfs"]]))
                            incProgress(1, detail = Tracking[["mergeVcfs"]])
                            vcf <- do.call(rbind, vcfs)
                        }

                    }
                )

                validate(need(vcf, Msgs[["vcf"]]))

                validate(need(
                    length(vcf) > 0,
                    "No variant found in BED region(s)"))
            })
        t2 <- Sys.time()
        dt <- t2 - t1
        message(sprintf(
            "%i variants from %i region(s) imported in %.2f %s",
            length(vcf),
            length(bedRecords),
            as.numeric(dt),
            units(dt)))

        vcf
    })

    output$vcfCols <- renderUI({

        vcf <- vcf()

        validate(need(vcf, label = Msgs[["vcf"]]))

        colChoices <- c(colnames(mcols(vcf)))

        selectInput(
            "vcfCols", "Meta-columns",
            choices = colChoices,
            selected = c(colChoices[1:min(4, length(colChoices))]),
            multiple = TRUE
        )
    })

    output$vcfRowRanges <- DT::renderDataTable({

        vcf <- vcf()

        req(vcf, input$vcfCols)

        cols <- which(colnames(mcols(vcf)) %in% input$vcfCols)

        DT::datatable(
            data = as.data.frame(rowRanges(vcf)[, cols]),
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top")
    })

    output$vcfInfoCols <- renderUI({

        vcf <- vcf()

        validate(need(vcf, label = Msgs[["vcf"]]))

        # All columns except the VEP predictions
        colChoices <- grep(
            pattern = input$csqField,
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

        vcf <- vcf()

        validate(need(
            length(input$vcfInfoCols) > 0,
            "No data to show"
        ))

        cols <- which(colnames(info(vcf)) %in% input$vcfInfoCols)

        DT::datatable(
            data = as.data.frame(info(vcf)[, cols]),
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top")
    })

    output$vcf <- renderPrint({

        vcf <- vcf()

        validate(need(vcf, Msgs[["vcf"]]))

        vcf
    })

    observeEvent(input$vcfInputMode, {

        validate(need(input$bedFile, label = Msgs[["bedFile"]]))

        updateTabsetPanel(session, "tabset.settings", selected = "BED")
        updateTabsetPanel(session, "tabset.bed", selected = "Sample")
    })

    # Filter variants post-import ----

    # TODO" at the moment variants are still filtered during import
    # vcfFilters <- reactive({
    #
    #     .vcfFilters(
    #         mafMin = input$maf.min,
    #         mafMax = input$maf.max)
    #
    # })

    # Parse genotypes ----

    output$genotypeStructure <- renderPrint({

        vcf <- vcf()

        validate(need(vcf, Msgs[["vcf"]]))

        genotypes <- geno(vcf)[["GT"]]
        validate(need(genotypes, Msgs[["genotypes"]]))

        str(genotypes)
    })

    output$genoNumCols <- renderUI({

        vcf <- vcf()

        validate(need(vcf, Msgs[["vcf"]]))

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

        vcf <- vcf()

        validate(need(vcf, Msgs[["vcf"]]))

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

        vcf <- vcf()

        validate(need(vcf, Msgs[["vcf"]]))

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

        vcf <- vcf()

        validate(need(vcf, Msgs[["vcf"]]))

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

        vcf <- vcf()

        validate(need(vcf, Msgs[["vcf"]]))

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
    output$csqStructure <- renderPrint({

        vcf <- vcf()

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(input$csqField, label = Msgs[["csqField"]]))

        csq <- tryParseCsq(vcf = vcf, csqField = input$csqField)

        validate(need(csq, Msgs[["csq"]]))

        str(csq)
    })

    output$csqCols <- renderUI({

        vcf <- vcf()

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(input$csqField, Msgs[["csqField"]]))

        csq <- tryParseCsq(vcf = vcf, csqField = input$csqField)

        validate(need(csq, Msgs[["csq"]]))

        csqMcols <- mcols(csq)

        validate(need(csqMcols, Msgs[["csqMcols"]]))

        selectInput(
            "csqCols", "Meta-columns",
            choices = colnames(csqMcols),
            selected = colnames(csqMcols)[1:5],
            multiple = TRUE)
    })

    output$csqSample <- DT::renderDataTable({

        req(input$csqCols)

        vcf <- vcf()

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(input$csqField, Msgs[["csqField"]]))

        csq <- tryParseCsq(vcf = vcf, csqField = input$csqField)

        validate(need(csq, Msgs[["csq"]]))

        csqMcols <- mcols(csq)

        validate(need(csqMcols, Msgs[["csqMcols"]]))

        cols <- which(colnames(csqMcols) %in% input$csqCols)

        DT::datatable(
            cbind(
                Variant = names(csq),
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

    # When predictions are updated, update the choices for plotting
    observeEvent(vcf(), {

        vcf <- vcf()

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(input$csqField, label = Msgs[["csqField"]]))

        csq <- tryParseCsq(vcf = vcf, csqField = input$csqField)

        validate(need(csq, Msgs[["csq"]]))

        csqMcols <- mcols(csq)

        csq.choices <- colnames(csqMcols)

        csq.default <- max(
            1,
            head(
                x = match(c("consequence","impact"), tolower(csq.choices)),
                n = 1),
            na.rm = TRUE)

        updateSelectInput(
            session, "variantCsq",
            choices = csq.choices,
            selected = csq.choices[csq.default])

        facet.choices <- c("None", csq.choices)

        updateSelectInput(
            session, "facet",
            choices = facet.choices)
    })

    # Import phenotype information ----

    phenotypes <- reactive({

        validate(need(input$phenoFile, Msgs[["phenoFile"]]))

        message("Importing phenotypes ...")

        rawData <- read.table(
            file = input$phenoFile$datapath,
            header = TRUE, row.names = 1)

        validate(need(
            all(dim(rawData) > c(0, 0)),
            paste(
                "Phenotype file must have at least 1 row and 1 column",
                "with colnames and rownames")))

        DataFrame(rawData)
    })

    output$phenoCols <- renderUI({

        phenos <- phenotypes()

        validate(need(
            nrow(phenos) > 0,
            "No phenotype information")
        )

        selectInput(
            "phenoCols", "Phenotypes",
            choices = colnames(phenos),
            selected = colnames(phenos)[1:5],
            multiple = TRUE
        )
    })

    # When phenotypes are updated, update the choices for plotting
    observeEvent(vcf(), {

        pheno.choices <- c("None", colnames(phenotypes()))

        # If new phenotype files also contains the current phenotype,
        # keep it active instead of resetting to the first choice
        pheno.selected <- pheno.choices[
            max(1, which(pheno.choices == input$phenotype))]

        updateSelectInput(
            session, "phenotype",
            choices = pheno.choices,
            selected = pheno.selected)
    })

    output$phenoFileSummary <- renderUI({

        phenos <- phenotypes()

        validate(need(phenos, "phenos"))

        HTML(paste(
            ncol(phenos),
            "phenotypes in",
            nrow(phenos),
            "samples detected."
        ))
    })

    output$phenotypeStructure <- renderPrint({

        phenos <- phenotypes()

        validate(need(phenos, Msgs[["phenos"]]))

        str(phenos)
    })

    output$phenotypesSample <- DT::renderDataTable({

        req(input$phenoCols)

        phenos <- phenotypes()

        # req gives time to update data when new phenotype file provided
        validate(need(phenos, Msgs[["phenos"]]))

        req(all(input$phenoCols %in% colnames(phenos)))

        DT::datatable(
            as.data.frame(phenos[,input$phenoCols, drop = FALSE]),
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top")
    })

    # Summarise consequences ----

    csqTable <- reactive({

        vcf <- vcf()

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(input$csqField, label = Msgs[["csqField"]]),
            need(input$variantCsq, label = Msgs[["variantCsq"]]),
            need(input$facet, label = Msgs[["facet"]]))

        if (input$facet != "None"){
            req(input$csqFacets) # Give time to observeEvent
            validate(need(
                length(input$csqFacets) > 0,
                "The plot must contain at least one facet"))
        }

        csq <- tryParseCsq(vcf = vcf, csqField = input$csqField)

        validate(need(csq, Msgs[["csq"]]))

        csqMcols <- mcols(csq)

        validate(need(csqMcols, Msgs[["csqMcols"]]))

        if (input$facet == "None"){
            facet <- NULL
        } else {
            facet <- input$facet
        }

        # If no Phenotype file given or None are explicitely required
        if (is.null(input$phenoFile) || input$phenotype == "None"){
            message("Summarising predictions across all samples...")
            # Note that facet may have been updated to NULL above
            rawTable <- data.frame(
                Prediction = csqMcols[,input$variantCsq],
                facet = csqMcols[,facet],
                Phenotype = "All")
            # Do not rely on position as facet column is absent if NULL
            colnames(rawTable) <- gsub(
                "Prediction", input$variantCsq, colnames(rawTable))
            colnames(rawTable) <- gsub(
                "facet", input$facet, colnames(rawTable))
        }
        # If Phenotype file given AND one phenotype is selected
        else {
            message("Summarising prediction by phenotype: ", input$phenotype)
            rawTable <- tabulateCsqByPhenotype(
                vcf = vcf,
                phenoCol = input$phenotype,
                csqCol = input$variantCsq,
                param = tParam(),
                unique = input$csqPhenoUnique,
                facet = facet,
                plot = FALSE,
                percentage = FALSE
            )
        }

        if (input$facet != "None"){
            rawTable <- rawTable[rawTable[,input$facet] %in% input$csqFacets,]
            rawTable <- droplevels(rawTable)
        }

        validate(need(
            nrow(rawTable) > 0,
            "No data to show"))

        rawTable
    })

    output$csqTableDecreasing <- DT::renderDataTable({

        csqTable <- csqTable()

        validate(need(csqTable, label = Msgs[["csqTable"]]))

        csqTableTable <- as.data.frame(table(csqTable))

        validate(need(csqTableTable, Msgs[["csqTableTable"]]))

        csqTableTableDecreasing <- csqTableTable[
            order(csqTableTable[,"Freq"], decreasing = TRUE),]

        validate(need(
            csqTableTableDecreasing,
            Msgs[["csqTableTableDecreasing"]]))

        DT::datatable(
            csqTableTableDecreasing,
            options = list(
                pageLength = 20,
                searching = TRUE),
            filter = "top",
            rownames = FALSE)
    })

    # Used by 2+ reactives, might as well have as a small reactive
    varCsqPlotPheno <- reactive({

        phenotype <- input$phenotype

        validate(need(phenotype, label = Msgs[["phenotype"]]))

        if (phenotype == "None"){
            return("Phenotype")
        } else {
            return(phenotype)
        }
    })

    # Update list of facets if the faceting variable changes
    observeEvent(input$facet, {

        vcf <- vcf()

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(input$csqField, label = Msgs[["csqField"]]))

        csq <- tryParseCsq(vcf = vcf, csqField = input$csqField)

        validate(need(csq, Msgs[["csq"]]))

        csqMcols <- mcols(csq)

        if (input$facet != "None"){
            facet.choices <- unique(csqMcols[,input$facet])
        } else {
            facet.choices <- c()
        }

            updateSelectInput(
            session, "csqFacets",
            choices = facet.choices,
            selected = facet.choices)
    })

    # Plot consequence barplot ----

    # Stacked bars classifying consequence of variants
    output$variantsCsqBarplot <- renderPlot({
        withProgress(min = 0, max = 3, message = "Progress", {
            incProgress(1, detail = Tracking[["calculate"]])

            csqTable <- csqTable()
            varCsqPlotPheno <- varCsqPlotPheno()

            validate(
                need(csqTable, Msgs[["csqTable"]]),
                need(varCsqPlotPheno, Msgs[["varCsqPlotPheno"]]),
                need(
                    input$variantCsq != "",
                    label = Msgs[["variantCsq"]]),
                need(input$facet, label = Msgs[["facet"]]),
                need(input$legendText, label = Msgs[["legendText"]]),
                need(input$xAxisAngle, label = Msgs[["xAxisAngle"]]),
                need(input$xAxisHjust, label = Msgs[["xAxisHjust"]]),
                need(input$xAxisVjust, label = Msgs[["xAxisVjust"]]),
                need(input$xAxisSize, label = Msgs[["xAxisSize"]]),
                need(
                    is.logical(input$stackedPercentage),
                    label = Msgs[["stackedPercentage"]]))

            incProgress(1, detail = Tracking[["ggplot"]])

            gg <- ggplot(
                data = csqTable,
                mapping = aes_string(
                    x = varCsqPlotPheno,
                    fill = input$variantCsq)) +
                theme(
                    axis.text = element_text(size = rel(1.5)),
                    axis.title = element_text(size = rel(1.5)),
                    legend.text = element_text(size = rel(input$legendText)),
                    legend.title = element_text(size = rel(input$legendText)),
                    axis.text.x = element_text(
                        angle = input$xAxisAngle,
                        hjust = input$xAxisHjust,
                        vjust = input$xAxisVjust,
                        size = rel(input$xAxisSize))
                ) +
                guides(
                    fill = c("none", "legend")[input$legend + 1]
                ) +
                scale_x_discrete(drop = FALSE)

            if (input$stackedPercentage){
                gg <- gg +
                    geom_bar(position = "fill")
            } else {
                gg <- gg + geom_bar()
            }

            if (input$facet != "None"){
                gg <- gg + facet_wrap(facets = input$facet)
            }

            message("Plotting predictions...")
            incProgress(1, detail = Tracking[["render"]])
            gg
        })
    })

    # Print the count of consequences in the area hovered.
    output$varCsqCount <- renderUI({

        csqTable <- csqTable()
        req(csqTable, input$plotVarClass_hover)

        hover <- input$plotVarClass_hover

        validate(
            need(
                input$variantCsq != "",
                label = Msgs[["variantCsq"]]),
            need(csqTable, label = Msgs[["csqTable"]]))
        if (input$facet == "None"){
            csqTableTable <- as.data.frame(table(
                csqTable[,c(
                    input$variantCsq,
                    hover$mapping$x)]))
        } else {
            csqTableTable <- as.data.frame(table(
                csqTable[,c(
                    input$variantCsq,
                    hover$mapping$x,
                    hover$mapping$panelvar1)]))
        }

        validate(need(csqTableTable, Msgs[["csqTableTable"]]))

        variantCsq <- input$variantCsq

        # If not hovered, show empty values
        if (is.null(hover)){
            x_lvl <- NULL
            y_lvl <- NULL
            countVarLevel <- NULL
        }
        else{
            # Identify the level of the phenotype (x axis)
            x_lvls <- levels(csqTableTable[,hover$mapping$x])
            x_lvl <- x_lvls[round(hover$x)]
            # Identify the level of the variant category (y axis)
            filters <- data.frame(
                pheno = csqTableTable[,hover$mapping$x] == x_lvl
            )
            if (input$facet != "None"){
                filters$facet <- csqTableTable[
                    ,hover$mapping$panelvar1] == hover$panelvar1
            }

            y_lvls <- csqTableTable[
                apply(X = filters, MARGIN = 1, FUN = all),
                c(variantCsq, "Freq")]

            # If in percentage mode, rescale the values
            if (input$stackedPercentage){
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
                y_lvl <- as.character(y_lvls[1, variantCsq])
                countVarLevel <- y_lvls[y_lvls[,variantCsq] == y_lvl, "Freq"]
                countVarLevel <- ifelse(
                    input$stackedPercentage,
                    yes = sprintf(
                        "%.1f%%", 100 * countVarLevel),
                    no = as.character(countVarLevel))
            }
        }

        html1Start <- "<ul>"
        html2Facet <- ifelse(
            input$facet == "None",
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
            html4Csq <- ""
            html5Value <- ""
        } else {
            html4Csq <- sprintf(
                "<li>%s : <code>%s</code></li>",
                variantCsq,
                y_lvl
            )
            html5Value <- ifelse(
                input$stackedPercentage,
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
            html1Start, html2Facet, html3Pheno, html4Csq, html5Value, html6End,
            sep = ""))
    })

    # Genotype heatmap ----

    output$heatmapGenotype <- renderPlot({

        if (input$doGenoHeatmap == 0)
            return()
        validate(need(input$doGenoHeatmap > 0, Msgs[["doGenoHeatmap"]]))

        # Remove dependency on vcf reactive
        isolate({ vcf <- vcf() })

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(tParam(), label = Msgs[["tParam"]])
        )

        withProgress(min = 0, max = 3, message = "Progress", {
            incProgress(1, detail = Tracking[["calculate"]])

            genotypes <- geno(vcf)[["GT"]]
            validate(need(genotypes, Msgs[["genotypes"]]))

            validGenotypes <- unlist(genos(tParam()))

            genos.long <- melt(genotypes, value.name = "Genotype")

            colnames(genos.long)[1:2] <- c("Variants", "Samples")

            genos.long[,"Genotype"] <- factor(
                x = genos.long[,"Genotype"],
                levels = unlist(genos(tParam()))
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

        validate(need(input$bpConfig, label = Msgs[["bpConfig"]]))

        switch(
            input$bpConfig,
            SerialParam = SerialParam(),
            MulticoreParam = MulticoreParam(workers = input$bpCores),
            SnowParam = SnowParam(
                workers = input$bpCores, type = input$bpType))
    })

    # Cleanup routine ----

    cancel.onSessionEnded <- session$onSessionEnded(function() {
        return(TRUE)
    })
})
