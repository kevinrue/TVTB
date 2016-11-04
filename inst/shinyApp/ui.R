
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(navbarPage(theme = "bootstrap.css",

    title = "tSVE",
    id = "tSVE",
    navbarMenu(
        title = "Input",

        tabPanel(
            title = "Phenotypes",

            # Phenotype options ----

            wellPanel(
                h4("Phenotypes"),
                hr(),

                fluidRow(

                    column(
                        width = 2,
                        p(strong("Phenotype file")),
                        actionButton(
                            "selectPheno", "Browse",
                            icon = icon("file"), width = '100%')
                    ),
                    column(
                        width = 6,
                        strong("File"), br(),
                        uiOutput("phenoFile"),
                        # Wrap long file names
                        tags$head(tags$style(
                            "#phenoFile{
                            display:block;
                            word-wrap:break-word;
                            }"
                            ))
                    ),
                    column(
                        width = 4,
                        strong("Summary"), br(),
                        htmlOutput("phenoFileSummary")
                    )

                ),
                fluidRow(
                    column(
                        width = 2,
                        br(),
                        actionButton(
                            "demoPheno", "Sample file",
                            icon = icon("file-text"), width = '100%')
                    )
                )

            )

        ),

        tabPanel(
            title = "Genomic ranges",

            # GRanges options ----
            wellPanel(
                h4("Genomic ranges"),
                hr(),

                fluidRow(

                    column(
                        width = 2,
                        selectInput(
                            "grangesInputMode", "Input type",
                            choices = GS[["choices.grangesInputMode"]],
                            selected = GS[["default.grangesInputMode"]],
                            width = '100%')
                    ),

                    column(
                        width = 4, offset = 6,
                        strong("Summary"),
                        htmlOutput("rangesSummary")
                    )

                ),

                fluidRow(

                    conditionalPanel(
                        condition = "input.grangesInputMode == 'bed'",
                        fluidRow(
                            column(
                                width = 2,
                                br(),
                                actionButton(
                                    "selectBed", "Browse",
                                    icon = icon("file"), width = '100%')
                            ),
                            column(
                                width = 6,
                                strong("File"), br(),
                                uiOutput("bedFile")
                            ),
                            # Wrap long file names
                            tags$head(tags$style(
                                "#bedFile{
                            display:block;
                           scanVcfHeader word-wrap:break-word;
                            }"
                            ))
                        ),
                        fluidRow(
                            column(
                                width = 2,
                                br(),
                                actionButton(
                                    "demoBed", "Sample file",
                                    icon = icon("file-text"), width = '100%')
                            )
                        )
                    ),
                    conditionalPanel(
                        condition = "input.grangesInputMode == 'ucsc'",
                        column(
                            width = 8,
                            textInput(
                                "ucscRanges", "UCSC-type genomic ranges",
                                value = "",
                                placeholder = paste(
                                    "chr21:33,031,597-33,041,570",
                                    "chr2:2,031,597-2,041,570",
                                    "...",
                                    sep = " ; "),
                                width = "100%")
                        ),
                        column(
                            width = 2,
                            br(),
                            actionButton(
                                "demoUCSC", "Sample input",
                                icon = icon("font"), width = '100%')
                        )
                    ),
                    conditionalPanel(
                        condition = "input.grangesInputMode == 'EnsDb'",
                        column(
                            width = 2,
                            selectInput(
                                "ensDb.type", NA,
                                choices = GS[["choices.ensDbType"]],
                                selected = GS[["default.ensDbType"]])
                        ),
                        column(
                            width = 1,
                            selectInput(
                                "ensDb.condition", NA,
                                choices = GS[["choices.ensDbFilters"]],
                                selected = GS[["default.ensDbFilters"]])
                        ),
                        column(
                            2,
                            textInput(
                                "ensDb.value", NA,
                                value = "",
                                placeholder = "SLC24A5,IL17A,...")
                        ),
                        column(
                            width = 2,
                            actionButton(
                                "demoEnsDb", "Sample input",
                                icon = icon("font"), width = '100%')
                        ),
                        column(
                            width = 4, offset = 1,
                            strong("Note"),
                            p(
                                "For the ", code("like"), "filter,",
                                "use ", code("%"), "as wildcard."
                            )
                        ),
                        tabsetPanel(
                            id = "ensDb.resultTab",
                            selected = "Genes",

                            tabPanel(
                                title = 'Genes',
                                DT::dataTableOutput("ensDb.Genes")
                            )#, # TODO
                            # tabPanel('Transcripts',
                            #          dataTableOutput("Transcripts")
                            # ),
                            # tabPanel('Exons',
                            #          dataTableOutput("Exons")
                            # )
                        )

                    )

                )

            )

        ),

        tabPanel(
            title = "Variants",

            # VCF options ----

            wellPanel(
                h4("VCF file(s)"),
                hr(),

                fluidRow(

                    column(
                        width = 2,
                        selectInput(
                            "vcfInputMode", "VCF input type",
                            choices = GS[["choices.vcfInputMode"]],
                            selected = GS[["default.vcfInputMode"]],
                            width = '100%')
                    ),
                    conditionalPanel(
                        condition = "input.vcfInputMode == 'SingleVcf'",
                        column(
                            width = 2,
                            br(),
                            actionButton(
                                "selectVcf", "Browse",
                                icon = icon("file"), width = '100%')
                        ),
                        column(
                            width = 2,
                            br(),
                            actionButton(
                                "demoVcf", "Sample file",
                                icon = icon("file-text"), width = '100%')
                        ),
                        column(
                            width = 4, offset = 2,
                            strong("Summary"),
                            br(),
                            textOutput("selectedVcf"),
                            # Wrap long file names
                            tags$head(tags$style(
                                "#selectedVcf{
                                display:block;
                                word-wrap:break-word;
                                }"
                            ))
                        )
                    ),
                    conditionalPanel(
                        condition = "input.vcfInputMode == 'OnePerChr'",
                        column(
                            width = 6,
                            textInput(
                                "vcfFolder", "Folder of VCF files",
                                value = GS[["default.vcfFolder"]],
                                width = '100%',
                                placeholder = "./extdata")
                        ),
                        column(
                            width = 4,
                            strong("Summary"),
                            htmlOutput("vcfFolderSummary")
                        )
                    )

                ),

                fluidRow(

                    conditionalPanel(
                        condition = "input.vcfInputMode == 'OnePerChr'",
                        fluidRow(

                            column(
                                width = 6, offset = 2,
                                textInput(
                                    "vcfPattern",
                                    paste(
                                        "Pattern of VCF files",
                                        "(%s : chromosome placeholder)"
                                    ),
                                    value = GS[["default.vcfPattern"]],
                                    width = '100%',
                                    placeholder = "^chr%s_.*\\.vcf\\.gz$")
                            )

                        )

                    )

                )

            ),

            wellPanel(
                h4("VCF scan parameters"),
                hr(),

                fluidRow(

                    # ScanVcfParam ----

                    # INFO fields (except VEP)
                    column(
                        width = 8,
                        fluidRow(
                            column(
                                width = 12,
                            selectInput(
                                "vcfInfoKeys", "INFO fields",
                                choices = c(), selected = c(),
                                multiple = TRUE)
                            )
                        ),
                        fluidRow(
                            column(
                                width = 2,
                                actionButton(
                                    "tickAllInfo", "Select all",
                                    icon = icon("check-square-o"))
                            ),
                            column(
                                width = 2,
                                actionButton(
                                    "untickAllInfo", "Deselect all",
                                    icon = icon("square-o"))
                            ),
                            column(
                                width = 7, offset = 1,
                                strong("Note:"),
                                "VEP field implicitely required"
                            )
                        )

                    ),


                    # VEP prediction INFO field
                    column(
                        width = 2,
                        textInput(
                            "vepKey", "VEP field (INFO)",
                            value = GS[["default.vep"]],
                            placeholder = 'CSQ, ANN, ...')
                    ),

                    # FORMAT fields
                    column(
                        width = 2,
                        selectInput(
                            "vcfFormatKeys", "FORMAT fields",
                            choices = c(), selected = c(), # TODO ScanVcfParam
                            multiple = TRUE),
                        strong("Note:"), "\"GT\" implicitely required"
                    )

                )
            ),

            wellPanel(
                fluidRow(

                    # VEP prediction INFO field ----
                    column(
                        width = 2, offset = 5,
                        br(),
                        actionButton(
                            "importVariants", "Import variants",
                            icon = icon("open", lib = "glyphicon")
                        )
                    ),
                    column(
                        width = 4, offset = 1,
                        strong("Summary"),
                        htmlOutput("vcfSummary")
                    )

                )
            ),

            hr(),

            fluidRow(
                column(
                    width = 6,
                    h4("Content of folder"),
                    hr(),
                    DT::dataTableOutput("vcfContent")
                ),
                column(
                    width = 6,
                    h4("VCF file(s) matching pattern"),
                    hr(),
                    DT::dataTableOutput("vcfFiles")
                )
            )

        ),

        tabPanel(
            title = "Annotations",

            # Genome annotation package ----

            wellPanel(
                h4("Annotations"),
                hr(),
                fluidRow(

                    column(
                        width = 3,
                        selectInput(
                            "annotationPackage",
                            "Select installed EnsDb package",
                            choices = as.list(EnsDbPacks),
                            width = '100%')
                    ),

                    column(
                        3,
                        strong("EnsDb annotation"),
                        htmlOutput("ensembl_organism"),
                        htmlOutput("ensembl_version"),
                        htmlOutput("ensembl_genome")
                    ),

                    column(
                        width = 6,
                        strong("Note"),
                        p(
                            "Only",
                            tags$code("EnsDb"),
                            "annotation packages supported for starters.",
                            "Ideally, ",
                            tags$code("TxDb"),
                            "and",
                            tags$code("OrganismDb"),
                            "packages supported soon.")
                    )

                )

            )

        )

    ),

    # Calculate frequencies ----

    tabPanel(
        title = "Frequencies", icon = icon("calculator "),

        wellPanel(
            fluidRow(
                column(
                    width = 2,
                    strong("Latest changes:")
                ),
                column(
                    width = 10,
                    uiOutput("latestFrequenciesCalculated")
                )
            )
        ),

        wellPanel(
            fluidRow(
                h4("Overall frequencies"), hr(),
                column(
                    width = 1, offset = 1,
                    actionButton(
                        "addOverallFrequencies", "Add",
                        icon = icon("plus")

                    )
                ),
                column(
                    width = 1,
                    actionButton(
                        "removeOverallFrequencies", "Remove",
                        icon = icon("minus")
                    )
                )
            )
        ),

        wellPanel(
            fluidRow(
                h4("Frequencies in phenotype levels"), hr(),
                column(
                    width = 2,
                    selectInput(
                        "phenoAddFrequencies",
                        "Phenotype",
                        choices = c()
                    )
                ),
                column(
                    width = 2, offset = 1,
                    actionButton(
                        "tickAllPhenoLevelsFreq",
                        "Select all",
                        icon = icon("check-square-o"),
                        width = "100%"
                    ), br(),
                    actionButton(
                        "untickAllPhenoLevelsFreq",
                        "Deselect all",
                        icon = icon("square-o"),
                        width = "100%"
                    )
                ),
                column(
                    width = 2, offset = 1,
                    br(),
                    actionButton(
                        "buttonFrequencies", "Refresh",
                        icon = icon("refresh"), width = "100%"
                    )
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    checkboxGroupInput(
                        "phenoLevelFreqCheckboxes", "Phenotype levels",
                        choices = c(), inline = TRUE
                    )
                )
            )
        )

    ),

    # VCF filter Rules ----

    tabPanel(
        title = "Filters", icon = icon("filter"),

        wellPanel(
            h4("Add filter"),
            fluidRow(
                column(
                    width = 1,
                    br(),
                    actionButton(
                        "addNewFilter", "Add filter",
                        icon = icon("plus")
                    )
                ),
                column(
                    width = 1,
                    selectInput(
                        "newFilterClass", "Type",
                        choices = GS[["vcfFilterClass.choices"]],
                        selected = GS[["vcfFilterClass.default"]]
                    )
                ),
                column(
                    width = 1,
                    br(),
                    checkboxInput(
                        "newFilterActive", "Active?",
                        value = TRUE
                    )
                ),
                column(
                    width = 8,
                    textInput(
                        "newFilterExpression", "Expression",
                        placeholder = paste(
                            "grepl(\"pass\", tolower(FILTER))",
                            "ALT + HET > 0",
                            "IMPACT %in% c(\"HIGH\", \"MODERATE\")",
                            sep = " - or - "
                        )
                    )
                ),
                column(
                  width = 1,
                  br(),
                  actionButton(
                    "demoFilter", "Sample input",
                    icon = icon("font"), width = '100%')
                )
            ),
            fluidRow(
                column(
                    width = 12,
                    uiOutput("vcfFilterTest")
                )
            ),
            fluidRow(
                br(),
                p(strong("Notes:")),
                tags$ul(
                    tags$li(
                        "Filters are tested against variants to ensure the",
                        "validity of filters. Therefore, variants must be",
                        "loaded", em("before"), "filters can be created."
                    ),
                    tags$li(
                        "Currently, filters are not re-tested if variants are",
                        "updated. If variants are refreshed, users should",
                        "ensure filters remain valid, or remove filters",
                        "manually."
                    ),
                    tags$li(
                        "Users may ignore auto-correction of quotes in the",
                        strong("Expression"), "field. The application",
                        "automatically substitutes",
                        "curly quotes (single and double) by their",
                        "corresponding regular quotes (",
                        em("i.e."), code("\""), "and", code("'"), ")"
                    )
                )
            )
        ),
        wellPanel(
            fluidRow(
                column(
                    width = 4, offset = 1,
                    strong("Summary"), br(),
                    uiOutput("filtersSummary")
                ),
                column(
                    width = 2,
                    actionButton(
                        "filterVariants", "Apply filters",
                        icon = icon("filter"), width = "100%"
                    )
                ),
                column(
                    width = 4,
                    strong("Summary"), br(),
                    uiOutput("filteredVcfSummary")
                )
            )
        ),
        wellPanel(
            fluidRow(
                column(
                    width = 1,
                    strong("Class")
                ),
                column(
                    width = 1,
                    strong("Active?")
                ),
                column(
                    width = 8,
                    strong("Expression")
                )
            ),
            br(),
            uiOutput("vcfFilterControls")
        ),
        wellPanel(
            fluidRow(
                column(
                    width = 12,
                    verbatimTextOutput("vcfRules")
                )
            )
        )
    ),

    tabPanel(
        title = "Views", icon = icon("picture-o"),

        tabsetPanel(
            id = "tabset.views",

            # Genomic ranges view ----
            tabPanel(
                title = "Genomic ranges",

                fluidRow(
                    column(
                        width = 12,
                        DT::dataTableOutput("rangesSample")
                    )
                )

            ),

            # Variants view ----
            tabPanel(
                title = "Variants",

                fluidRow(
                    column(
                        width = 12,
                        wellPanel(
                            uiOutput("vcfCols")
                        )
                    )
                ),

                fluidRow(
                    column(
                        width = 12,
                        DT::dataTableOutput("vcfRowRanges")
                    )
                )

            ),

            # Variants INFO view ----
            tabPanel(
                title = "INFO",

                fluidRow(
                    column(
                        width = 12,
                        wellPanel(
                            uiOutput("vcfInfoCols")
                        )
                    )
                ),

                fluidRow(
                    column(
                        width = 12,
                        DT::dataTableOutput("vcfInfo")
                    )
                )

            ),

            # VEP predictions view ----
            tabPanel(
                title = "VEP",

                fluidRow(
                    column(
                        width = 12,
                        wellPanel(
                            uiOutput("vepCols")
                        )
                    )
                ),

                fluidRow(
                    column(
                        width = 12,
                        DT::dataTableOutput("vcfVep")
                    )
                )

            ),

            # Phenotypes view ----
            tabPanel(
                title = "Phenotypes",

                fluidRow(
                    column(
                        width = 12,
                        wellPanel(
                            uiOutput("phenoCols")
                        )
                    )
                ),

                fluidRow(
                    column(
                        width = 12,
                        DT::dataTableOutput("phenotypesSample")
                    )
                )

            ),

            # Genotypes view ----
            tabPanel(
                title = "Genotypes",

                tabsetPanel(
                    tabPanel(
                        title = "Matrix",
                        fluidRow(
                            wellPanel(
                                column(
                                    width = 6,
                                    uiOutput("genoNumRows")
                                ),
                                column(
                                    width = 6,
                                    uiOutput("genoFirstRow")
                                )
                            )
                        ),
                        fluidRow(
                            wellPanel(
                                column(
                                    width = 6,
                                    uiOutput("genoNumCols")
                                ),
                                column(
                                    width = 6,
                                    uiOutput("genoFirstCol")
                                )
                            )
                        ),
                        fluidRow(
                            column(
                                width = 12,
                                tableOutput("genotypesSample")
                            )
                        )
                    ),

                    tabPanel(
                        title = "Heatmap",

                        p(
                            "Click the button after loading variants",
                            "to generate/update the figure",
                            actionButton(
                                "doGenoHeatmap", "Go!",
                                icon = icon("time")
                            )
                        ),

                        fluidRow(
                            column(
                                width = 12,
                                plotOutput("heatmapGenotype")
                            )
                        ),

                        p(
                            "Notes",
                            tags$ul(
                                tags$li(
                                    "This may take some time to plot.",
                                    em(
                                        "(~15s for 218 variants & 5844",
                                        "samples)"
                                    )
                                ),
                                tags$li(
                                    "Only genotypes codes found in the data",
                                    "are listed in the legend, irrespective",
                                    "of those defined in the",
                                    tags$strong("Advanced settings"), "."
                                )
                            )
                        )

                    )

                )

            )

        )

    ),

    navbarMenu(
        title = "Plots", icon = icon("pie-chart"),

        # VEP count barplot ----

        tabPanel(
            title = "VEP counts",

            sidebarLayout(

                # Sidebar with a slider input
                sidebarPanel(
                    width = 3,

                    actionButton(
                        "buttonTVBP", "Apply",
                        icon = icon("picture"), width = "100%"
                    ), hr(),

                    selectInput(
                        "vepTVBP",
                        "Variant effect prediction",
                        choices = c()
                    ),

                    selectInput(
                        "phenoTVBP",
                        "Phenotype field",
                        choices = c("None"),
                        selected = "None"
                    ),

                    conditionalPanel(
                        condition = "input.phenoTVBP != 'None'",
                        checkboxInput(
                            "unique2phenoTVBP",
                            "Unique to phenotype?",
                            value = FALSE
                        )
                    ),

                    selectInput(
                        "vepFacetKeyTVBP",
                        "VEP faceting key",
                        choices = c("None"),
                        selected = "None"
                    ),

                    conditionalPanel(
                        condition = "input.facetTVBP != 'None'",
                        selectInput(
                            "vepFacetsTVBP",
                            "Facets",
                            choices = c(),
                            selected = c(),
                            multiple = TRUE
                        )
                    ),

                    checkboxInput(
                        "stackedPercentageTVBP",
                        "Show as percentage?",
                        value = FALSE
                    ),

                    checkboxInput(
                        "advancedTVBP",
                        "Advanced controls",
                        value = FALSE),

                    conditionalPanel(
                        condition = "input.advancedTVBP == true",
                        checkboxInput(
                            "legendTVBP",
                            "Show legend",
                            value = TRUE)
                    ),

                    conditionalPanel(
                        condition = paste(
                            "input.advancedTVBP == true",
                            "input.legendTVBP == true",
                            sep = " && "
                        ),
                        sliderInput(
                            "legendTextSizeTVBP", "Legend font size",
                            value = 1, min = 0.1, max = 2, step = 0.1)
                    ),

                    conditionalPanel(
                        condition = "input.advancedTVBP == true",
                        sliderInput(
                            "xAxisAngleTVBP",
                            "Angle of X labels",
                            min = 0, max = 90, value = 0, step = 5
                        ),

                        sliderInput(
                            "xAxisSizeTVBP",
                            "Relative size of X text",
                            value = 1, min = 0.1, max = 2, step = .1
                        ),

                        sliderInput(
                            "xAxisVjustTVBP",
                            "Vertical just. of X labels",
                            min = 0, max = 1, value = 0.5, step = 0.1
                        ),

                        sliderInput(
                            "xAxisHjustTVBP",
                            "Horiz. just. of X labels",
                            min = 0, max = 1, value = 0.5, step = 0.1
                        )
                    )

                ),

                mainPanel(

                    tabsetPanel(
                        id = "vepCount",

                        tabPanel(
                            title = "Barplot",

                            fluidRow(
                                column(
                                    width = 12,
                                    plotOutput(
                                        "vepCountBarplot",
                                        height = vepCountBarplotHeight,
                                        hover = hoverOpts(
                                            "plotVarClass_hover",
                                            delayType = "debounce")
                                    )
                                )
                            ),

                            fluidRow(
                                column(
                                    width = 12,
                                    htmlOutput("varVepCount")
                                )
                            )
                        ),

                        tabPanel(
                            title = "Table",
                            DT::dataTableOutput("vepTableDecreasing")
                        )

                    )

                ) # main panel

            ) # sidebar layout

        ),

        # VEP density plot ----

        tabPanel(
            title = "VEP density",

            sidebarLayout(

                # Sidebar with a slider input
                sidebarPanel( # width = 4
                    width = 3,

                    actionButton(
                        "buttonDVBP", "Apply",
                        icon = icon("picture"), width = "100%"
                    ),
                    hr(),

                    selectInput(
                        "vepDVBP",
                        "Variant effect prediction",
                        choices = c()
                    ),

                    selectInput(
                        "phenoDVBP",
                        "Phenotype field",
                        choices = c("None"),
                        selected = "None"
                    ),

                    selectInput(
                        "layerDVBP",
                        "Layer(s)",
                        choices = list(
                            "Density & Dot plot" = "density+dotplot",
                            "Density" = "density",
                            "Dot plot" = "dotplot"
                        ),
                        selected = "density+dotplot"
                    ),

                    conditionalPanel(
                        condition = "input.phenoDVBP != 'None'",
                        checkboxInput(
                            "unique2phenoDVBP",
                            "Unique to phenotype?",
                            value = FALSE
                        )
                    ),

                    selectInput(
                        "vepFacetKeyDVBP",
                        "VEP faceting key",
                        choices = c("None"),
                        selected = "None"
                    ),

                    conditionalPanel(
                        condition = "input.facetDVBP != 'None'",
                        selectInput(
                            "vepFacetsDVBP",
                            "Facets",
                            choices = c(),
                            selected = c(),
                            multiple = TRUE
                        )
                    ),

                    textInput(
                        "patternDVBP", "Pattern",
                        value = "", placeholder = ".*:(.*)"
                    ),

                    checkboxInput(
                        "advancedDVBP",
                        "Advanced controls",
                        value = FALSE),

                    conditionalPanel(
                        condition = "input.advancedDVBP == true",
                        checkboxInput(
                            "legendDVBP",
                            "Show legend",
                            value = TRUE)
                    ),

                    conditionalPanel(
                        condition = paste(
                            "input.advancedDVBP == true",
                            "input.legendDVBP == true",
                            sep = " && "
                        ),
                        sliderInput(
                            "legendTextSizeDVBP", "Legend font size",
                            value = 1, min = 0.1, max = 2, step = 0.1)
                    ),

                    conditionalPanel(
                        condition = "input.advancedDVBP == true",

                        sliderInput(
                            "xAxisSizeDVBP",
                            "Relative size of X text",
                            value = 1, min = 0.1, max = 2, step = .1
                        ),

                        sliderInput(
                            "yAxisSizeDVBP",
                            "Relative size of Y text",
                            value = 1, min = 0.1, max = 2, step = .1
                        )
                    )

                ),

                mainPanel( # width = 8
                    plotOutput(
                        outputId = "vepDensityPlot",
                        height = vepDensityplotHeight,
                        )
                )

            ) #sideBarLayout

        ) # tabPanel

    ),

    navbarMenu(
        title = "Settings", icon = icon("wrench"),

        # Advanced settings ----

        tabPanel(
            title = "Advanced",

            wellPanel(
                h4("Genotypes"),
                hr(),

                fluidRow(
                    column(
                        width = 4,
                        selectInput(
                            "refGenotypes", "Reference homozygote genotype(s)",
                            choices = setdiff(
                                GS[["all.genotypes"]],
                                c(
                                    GS[["default.hetGenotypes"]],
                                    GS[["default.altGenotypes"]])),
                            selected = GS[["default.refGenotypes"]],
                            multiple = TRUE
                        )
                    ),

                    column(
                        width = 4,
                        selectInput(
                            "hetGenotypes", "Heterozygote genotype(s)",
                            choices = setdiff(
                                GS[["all.genotypes"]],
                                c(
                                    GS[["default.refGenotypes"]],
                                    GS[["default.altGenotypes"]])),
                            selected = GS[["default.hetGenotypes"]],
                            multiple = TRUE
                        )
                    ),

                    column(
                        width = 4,
                        selectInput(
                            "altGenotypes", "Alternate homozygote genotype(s)",
                            choices = setdiff(
                                GS[["all.genotypes"]],
                                c(
                                    GS[["default.refGenotypes"]],
                                    GS[["default.hetGenotypes"]])),
                            selected = GS[["default.altGenotypes"]],
                            multiple = TRUE
                        )
                    )

                ),

                fluidRow(
                    column(
                        width = 1,
                        textInput(
                            "refSuffix", "Suffix",
                            value = GS[["default.refSuffix"]],
                            placeholder = GS[["default.refSuffix"]]
                        )
                    ),
                    column(
                        width = 1, offset = 3,
                        textInput(
                            "hetSuffix", "Suffix",
                            value = GS[["default.hetSuffix"]],
                            placeholder = GS[["default.hetSuffix"]]
                        )
                    ),
                    column(
                        width = 1, offset = 3,
                        textInput(
                            "altSuffix", "Suffix",
                            value = GS[["default.altSuffix"]],
                            placeholder = GS[["default.altSuffix"]]
                        )
                    )
                )

            ),

            wellPanel(
                h4("INFO suffixes"),
                hr(),

                fluidRow(
                    column(
                        width = 3,
                        textInput(
                            "aafSuffix", "ALT allele freq.",
                            value = GS[["default.aafSuffix"]],
                            placeholder = GS[["default.aafSuffix"]]
                        )
                    ),
                    column(
                        width = 3,
                        textInput(
                            "mafSuffix", "Minor allele freq.",
                            value = GS[["default.mafSuffix"]],
                            placeholder = GS[["default.mafSuffix"]]
                        )
                    )
                )
            ),

            wellPanel(
                h4("VCF file(s)"),
                hr(),

                fluidRow(

                    column(
                        width = 2,
                        numericInput(
                            "yieldSize", "VCF yield size (100-100^3)",
                            min = 100, max = 100E3,
                            value = 4E3,
                            step = 1E3
                        )
                    )

                )
            )
        ),

        tabPanel(
            title = "Parallel",

            tags$span(
                style="color:red",
                tags$em("Experimental")),
            tags$sup("1,2"),

            wellPanel(
                h4("Parallel settings"),
                hr(),
                fluidRow(

                    column(
                        width = 3,
                        numericInput(
                            "bpCores", "Cores",
                            value = PS[["default.bpCores"]],
                            min = 1, max = PS[["default.bpCores"]], step = 1)
                    ),

                    column(
                        width = 3,
                        selectInput(
                            "bpConfig", "Cluster configuration",
                            choices = structure(
                                PS[["choices.bpClass"]],
                                names = gsub(
                                    "Param", "", PS[["choices.bpClass"]])),
                            selected = PS[["default.bpClass"]])
                    ),

                    conditionalPanel(
                        condition = "input.bpConfig != 'SerialParam'",
                        column(
                            width = 3,
                            selectInput(
                                "bpType", "Cluster type",
                                choices = structure(
                                    PS[["choices.bpType"]],
                                    names = gsub(
                                        "Param", "", PS[["choices.bpType"]])),
                                selected = PS[["default.bpType"]])
                        )
                    )

                )

            ),

            p(
                tags$sup("1"),
                "Parallel workers succesfully tested on",
                tags$code("Ubuntu"),
                "and",
                tags$code("Scientific Linux"),
                "systems."
            ),
            p(
                tags$sup("2"),
                "Known issue on", code("Mac OS X El Capitan"), ":",
                "Web-app hangs while CPUs work infinitely at full capacity."
            )
        )

    ),

    # Session settings view ----

    tabPanel(
        title = "Session",

        tabsetPanel(
            id = "tabset.session",

            tabPanel(
                title = "Session info",
                verbatimTextOutput("sessionInfo")
            ),

            tabPanel(
                title = "TVTB settings",
                verbatimTextOutput("TVTBsettings")
            ),

            tabPanel(
                title = "General settings",
                verbatimTextOutput("generalSettings")
            ),

            tabPanel(
                title = "Advanced settings",
                verbatimTextOutput("advancedSettings")
            ),

            tabPanel(
                title = "BED",
                verbatimTextOutput("rangesStructure")
            ),

            tabPanel(
                title = "Variants",
                verbatimTextOutput("vcf")
            ),

            tabPanel(
                title = "VEP",
                verbatimTextOutput("vepStructure")
            ),

            tabPanel(
                title = "Phenotypes",
                verbatimTextOutput("phenotypesStructure")
            ),

            tabPanel(
                title = "Genotypes",
                verbatimTextOutput("genotypeStructure")
            )
        )

    )

))
