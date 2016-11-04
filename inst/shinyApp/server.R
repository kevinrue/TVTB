
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

    # Import phenotype information ----

    # Path to phenotype file, or NULL
    phenoFile <- reactive({
        # Triggered by input$tryCatch
        if (input$selectPheno > 0){
            selected <- tryCatch(
                file.choose(),
                error = function(err){
                    warning(geterrmessage())
                    return(NULL)
                })
        } else {
            return(NULL)
        }

        return(selected)
    })

    # DataFrame of imported phenotypes, or NULL
    phenotypes <- reactive({
        # Depends on phenoFile()
        phenoFile <- phenoFile()

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

        return(HTML(paste(
            ncol(phenotypes),
            "phenotypes in",
            nrow(phenotypes),
            "samples detected."
        )))
    })

    # Column names available for selection from phenotypes
    output$phenoCols <- renderUI({
        # NOTE:
        # phenotype view depends on the vcf(), not phenotypes()
        # 1) users can easily look at their raw phenotypes separately
        # 2) all analyses use data stored in the VCF object only

        vcf <- vcf()
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
    observeEvent(vcf(), {
        # Depends on vcf() & input$phenoAnalysed
        pheno.choices <- c("None", colnames(colData(vcf())))

        # If new phenotype files also contains the current phenotype,
        # keep it active instead of resetting to the first choice
        pheno.selected <- pheno.choices[
            max(1, which(pheno.choices == input$phenoAnalysed))]

        updateSelectInput(
            session, "phenotype",
            choices = pheno.choices,
            selected = pheno.selected)
    })

    # Display structure of phenotypes attached to variants
    output$phenotypesStructure <- renderPrint({
        # Depends on vcf()
        vcf <- vcf()
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
        # Requires: input$phenoCols, vcf()
        # Give time to initialise columns seletion widget
        req(input$phenoCols)

        # Display phenotypes attached to the VCF object
        vcf <- vcf()
        validate(need(vcf, Msgs[["vcf"]]))

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

        return(HTML(paste(
            length(genomicRanges),
            "BED record(s) detected.</br>",
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
        )))
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
            need(input$vepKey, Msgs[["vepKey"]]),
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
            "phenoFile()" = phenoFile(),
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
            vepAnalysed = input$vepAnalysed,
            phenoAnalysed = input$phenoAnalysed,
            unique2pheno = input$unique2pheno,
            vepFacetKey = input$vepFacetKey,
            vepFacets = input$vepFacets,
            stackedPercentage = input$stackedPercentage
        ))
    })

    output$advancedSettings <- renderPrint({
        return(list(
            legendTextSize = input$legendTextSize,
            xAxisAngle = input$xAxisAngle,
            xAxisSize = input$xAxisSize,
            xAxisHjust = input$xAxisHjust,
            xAxisVjust = input$xAxisVjust,
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

        HTML(paste(
            length(vcfFiles),
            "VCF.GZ file(s) detected.</br>",
            "[",
            basename(head(x = vcfFiles, n = 1)),
            " , ... ]"
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

    # Import VCF information ----

    observeEvent(input$importVariants, {

        updateActionButton(
            session, "importVariants",
            label = "Refresh variants", icon = icon("refresh"))

    })

    # Import Vcf object
    vcf <- reactive({

        # Only proceed if the variants should be refreshed
        validate(need(
            input$importVariants > 0,
            Msgs[["importVariants"]]))

        isolate({
            vcfInputMode <- input$vcfInputMode
            genomicRanges <- genomicRanges()
            tparam <- tparam()
            vepKey <- input$vepKey
            genomeSeqinfo <- genomeSeqinfo()
            yieldSize <- input$yieldSize

            if (is.null(phenoFile())) {
                phenotypes <- DataFrame()
            } else {
                phenotypes <- phenotypes()
            }
        })

        validate(
            need(vepKey, label = Msgs[["vepKey"]]),
            need(genomeSeqinfo, Msgs[["genomeSeqinfo"]])
        )

        shiny::withProgress(
            min = 0, max = 3, value = 1,
            message = "Progress", detail = Tracking[["preprocessing"]],
            {

                # TODO: let user choose fixed/info/format keys to import
                # for speed and memory
                svp <- ScanVcfParam()

                # Only import samples with phenotype information
                if (nrow(phenotypes) > 0)
                    vcfSamples(svp) <- rownames(phenotypes)

                # Only import variants in targeted regions
                if (length(genomicRanges) > 0)
                    vcfWhich(svp) <- reduce(genomicRanges) # optimise import

                # # Timing
                # t1 <- Sys.time()

                switch (
                    vcfInputMode,
                    SingleVcf = {

                        isolate({singleVcf <- singleVcf()})
                        validate(need(singleVcf, Msgs[["singleVcf"]]))

                        shiny::incProgress(1, detail = Tracking[["singleVcf"]])

                        vcf <- tryParseSingleVcf(
                            file = singleVcf,
                            svp = svp,
                            yieldSize = yieldSize)
                    },

                    OnePerChr = {

                        shiny::incProgress(1, detail = Tracking[["multiVcfs"]])

                        isolate({
                            vcfFolder <- input$vcfFolder
                            vcfPattern <- input$vcfPattern
                        })
                        validate(need(singleVcf, Msgs[["singleVcf"]]))

                        vcf <- tryParseMultipleVcf(
                            folder = vcfFolder,
                            pattern = vcfPattern,
                            svp = svp,
                            yieldSize = yieldSize,
                            BPPARAM = bp(tparam))
                    }
                )

                shiny::incProgress(1, detail = Tracking[["postprocessing"]])

                validate(need(
                    length(vcf) > 0,
                    "No variant found in BED region(s)"))

                if (nrow(phenotypes) == 0)
                    phenotypes <- DataFrame(row.names = colnames(vcf))
                SummarizedExperiment::colData(vcf) <- phenotypes
            })

        return(vcf)
    })

    output$vcfSummary <- renderUI({

        vcf <- vcf()

        # validate(need(vcf, "vcf"))

        HTML(paste(
            nrow(vcf),
            "bi-allelic records and",
            ncol(colData(vcf)),
            "phenotypes in",
            ncol(vcf),
            "samples imported"
        ))
    })

    output$vcfCols <- renderUI({

        vcf <- vcf()

        # validate(need(vcf, label = Msgs[["vcf"]]))

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

        # Give time to initialise widget
        req(input$vcfCols)

        cols <- which(colnames(mcols(vcf)) %in% input$vcfCols)

        displayedTable <- cbind(
            rownames = rownames(vcf),
            as.data.frame(rowRanges(vcf)[, cols], row.names = NULL)
        )

        DT::datatable(
            data = displayedTable,
            options = list(
                pageLength = 10,
                searching = TRUE),
            filter = "top")
    })

    output$vcfInfoCols <- renderUI({

        vcf <- vcf()

        # validate(need(vcf, label = Msgs[["vcf"]]))

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

        vcf <- vcf()

        validate(need(
            length(input$vcfInfoCols) > 0,
            "No INFO key selected."
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

    # Filter variants post-import ----

    # TODO

    # Filter variants ----

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
    output$vepStructure <- renderPrint({

        vcf <- vcf()

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(input$vepKey, label = Msgs[["vepKey"]]))

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        validate(need(csq, Msgs[["csq"]]))

        str(csq)
    })

    output$vepCols <- renderUI({

        vcf <- vcf()

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(input$vepKey, Msgs[["vepKey"]]))

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        vepMcols <- mcols(csq)

        selectInput(
            "vepCols", "Meta-columns",
            choices = colnames(vepMcols),
            selected = colnames(vepMcols)[1:5],
            multiple = TRUE)
    })

    output$vepSample <- DT::renderDataTable({

        req(input$vepCols)

        vcf <- vcf()

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(input$vepKey, Msgs[["vepKey"]])
            )

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        vepMcols <- mcols(csq)

        cols <- which(colnames(vepMcols) %in% input$vepCols)

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
            need(input$vepKey, label = Msgs[["vepKey"]]))

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        validate(need(csq, Msgs[["csq"]]))

        vepMcols <- mcols(csq)

        vepKey.choices <- colnames(vepMcols)

        # Initialise by order of preference if present:
        # Consequence > Impact > First field
        vepAnalysed.default <- max(
            1,
            head(x = na.omit(
                match(
                    x = c("consequence","impact"),
                    table = tolower(vepKey.choices)
                    )
                ), n = 1)
            )

        updateSelectInput(
            session, "vepAnalysed",
            choices = vepKey.choices,
            selected = vepKey.choices[vepAnalysed.default])

        vepKey.choices <- c("None", vepKey.choices)

        updateSelectInput(
            session, "vepFacetKey",
            choices = vepKey.choices)
    })

    # Summarise consequences ----

    # Phenotype to summarise VEP: or NULL if "None"
    varVepPlotPheno <- reactive({

        phenotype <- input$phenoAnalysed
        # Give time to initialise the widget
        req(phenotype)

        if (phenotype == "None"){
            return("Phenotype")
        } else {
            return(phenotype)
        }
    })

    # Base of the ggplot summarising counts of VEP
    vepTabulated <- reactive({

        withProgress(min = 0, max = 2, message = "Progress", {
            incProgress(1, detail = Tracking[["calculate"]])

            if (input$vepFacetKey != "None"){
                req(input$vepFacets) # Give time to observeEvent
                validate(need(
                    length(input$vepFacets) > 0,
                    "The plot must contain at least one facet"))
            }

            if (input$vepFacetKey == "None"){
                vepFacetKey <- NULL
            } else {
                vepFacetKey <- input$vepFacetKey
            }

            vcf <- vcf()

            if (input$phenoAnalysed == "None"){
                colData(vcf)[,"Phenotype"] <- factor("All")
                plotPhenotype <- "Phenotype"
            } else {
                plotPhenotype <- input$phenoAnalysed
            }

            varVepPlotPheno <- varVepPlotPheno()

            # Give time to initialise the widgets

            # req() causes hanging below
            # Instead validate allows widgets to be initialised
            validate(
                need(input$vepAnalysed, Msgs[["vepAnalysed"]]),
                need(input$vepFacetKey, Msgs[["vepFacetKey"]]),
                need(input$legendTextSize, Msgs[["legendTextSize"]]),
                need(input$xAxisAngle, Msgs[["xAxisAngle"]]),
                need(input$xAxisHjust, Msgs[["xAxisHjust"]]),
                need(input$xAxisVjust, Msgs[["xAxisVjust"]]),
                need(input$xAxisSize, Msgs[["xAxisSize"]]),
                need(input$stackedPercentage, Msgs[["stackedPercentage"]])
                )

            incProgress(1, detail = Tracking[["ggplot"]])

            gg <- tabulateVepByPhenotype(
                vcf = vcf,
                phenoCol = plotPhenotype,
                vepCol = input$vepAnalysed,
                param = tparam(),
                unique = input$unique2pheno,
                facet = vepFacetKey,
                plot = TRUE,
                percentage = input$stackedPercentage
            )

            if (input$vepFacetKey != "None"){
                gg$data <- gg$data[
                    gg$data[,input$vepFacetKey] %in% input$vepFacets,]
                gg$data <- droplevels(gg$data)
            }

            return(gg)
        })

    })

    # Stacked bars classifying consequence of variants
    output$vepCountBarplot <- renderPlot({
        withProgress(min = 0, max = 2, message = "Progress", {
            incProgress(1, detail = Tracking[["calculate"]])

            gg <- vepTabulated()

            message("Plotting predictions...")
            incProgress(1, detail = Tracking[["render"]])

            gg <- gg +
                theme(
                    axis.text = element_text(size = rel(1.5)),
                    axis.title = element_text(size = rel(1.5)),
                    legend.text = element_text(size = rel(input$legendTextSize)),
                    legend.title = element_text(size = rel(input$legendTextSize)),
                    axis.text.x = element_text(
                        angle = input$xAxisAngle,
                        hjust = input$xAxisHjust,
                        vjust = input$xAxisVjust,
                        size = rel(input$xAxisSize))
                ) +
                guides(fill = c("none", "legend")[input$legend + 1])

            return(gg)
        })
    })

    vepTableCount <- reactive({
        # Depeds on: vepTabulated, vepAnalysed

        vepTabulated <- vepTabulated()$data

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

    # Update list of VEP facets if the faceting variable changes
    observeEvent(input$vepFacetKey, {

        vcf <- vcf()

        validate(
            need(vcf, Msgs[["vcf"]]),
            need(input$vepKey, label = Msgs[["vepKey"]]))

        csq <- tryParseCsq(vcf = vcf, vepKey = input$vepKey)

        validate(need(csq, Msgs[["csq"]]))

        vepMcols <- mcols(csq)

        if (input$vepFacetKey != "None"){
            vepFacets.choices <- unique(vepMcols[,input$vepFacetKey])
        } else {
            vepFacets.choices <- c()
        }

            updateSelectInput(
            session, "vepFacets",
            choices = vepFacets.choices,
            selected = vepFacets.choices)
    })

    # Print the count of consequences in the area hovered.
    output$varVepCount <- renderUI({

        vepTableCount <- vepTableCount()

        req(vepTableCount, input$plotVarClass_hover)
        hover <- input$plotVarClass_hover
        vepAnalysed <- input$vepAnalysed

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
            if (input$vepFacetKey != "None"){
                filters$vepFacetKey <- vepTableCount[
                    ,hover$mapping$panelvar1] ==
                    hover$panelvar1
            }
            # Extract VEP counts for the panel hovered:
            # matches x (phenotype level) and y (VEP level)
            y_lvls <- vepTableCount[
                apply(X = filters, MARGIN = 1, FUN = all),
                c(vepAnalysed, "Freq")]

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
                y_lvl <- as.character(y_lvls[1, vepAnalysed])
                countVarLevel <- y_lvls[y_lvls[,vepAnalysed] == y_lvl, "Freq"]
                countVarLevel <- ifelse(
                    input$stackedPercentage,
                    yes = sprintf(
                        "%.1f%%", 100 * countVarLevel),
                    no = as.character(countVarLevel))
            }
        }

        html1Start <- "<ul>"
        html2Facet <- ifelse(
            input$vepFacetKey == "None",
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
                vepAnalysed,
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
            html1Start, html2Facet, html3Pheno, html4Vep, html5Value, html6End,
            sep = ""))
    })

    # Genotype heatmap ----

    output$heatmapGenotype <- renderPlot({

        if (input$doGenoHeatmap == 0)
            return()

        # Remove dependency on vcf reactive
        isolate({ vcf <- vcf() })

        validate(
            need(vcf, Msgs[["vcf"]]),
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
