
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(navbarPage(

    title = "tSVE",
    id = "tSVE",
    navbarMenu(
        title = "Input",

        tabPanel(
            title = "Data",

            ## Keep for a possible header message
            # fluidRow(
            #     shiny::column(
            #         width = 6, offset = 3,
            #         p(
            #             "Welcome. Please define general settings here.",
            #             style = "font-family: 'Source Sans Pro';"
            #         )
            #     )
            # ),
            # hr(),

            # VCF options ----

            wellPanel(
                h4("Variants"),
                hr(),

                fluidRow(

                    shiny::column(
                        width = 3,
                        selectInput(
                            "vcfInputMode", "VCF input type",
                            choices = GS[["choices.vcfInputMode"]],
                            selected = GS[["default.vcfInputMode"]],
                            width = '100%')
                    ),
                    conditionalPanel(
                        condition = "input.vcfInputMode == 'SingleVcf'",
                        shiny::column(
                            width = 2,
                            br(),
                            actionButton(
                                "selectVcf", "Browse",
                                icon = icon("file"), width = '100%')
                        ),
                        strong("Summary"),
                        br(),
                        shiny::column(
                            width = 7,
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
                        shiny::column(
                            width = 5,
                            textInput(
                                "vcfFolder", "Folder of VCF files",
                                value = GS[["default.vcfFolder"]],
                                width = '100%',
                                placeholder = "./extdata")
                        ),
                        shiny::column(
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

                            shiny::column(
                                width = 5, offset = 3,
                                textInput(
                                    "vcfPattern",
                                    "Pattern of VCF files (%s : chromosome)",
                                    value = GS[["default.vcfPattern"]],
                                    width = '100%',
                                    placeholder = "^chr%s_.*\\.vcf\\.gz$")
                            )

                        )

                    )

                ),

                fluidRow(

                    # VEP prediction INFO field ----
                    shiny::column(
                        width = 2,
                        textInput(
                            "csqField", "VEP field",
                            value = GS[["default.csq"]],
                            placeholder = 'CSQ, ANN, ...'))

                )

            ),

            # Phenotype options ----

            wellPanel(
                h4("Phenotypes"),
                hr(),

                fluidRow(

                    shiny::column(
                        width = 3,
                        fileInput(
                            "phenoFile", "Phenotype file",
                            accept = c("txt"))
                    ),
                    shiny::column(
                        width = 4, offset = 5,
                        strong("Summary"),
                        htmlOutput("phenoFileSummary")
                    )

                )

            ),

            hr(),

            fluidRow(
                shiny::column(
                    width = 6,
                    h4("Contents of folder"),
                    hr(),
                    DT::dataTableOutput("vcfContent")
                ),
                shiny::column(
                    width = 6,
                    h4("VCF files in folder"),
                    hr(),
                    DT::dataTableOutput("vcfFiles")
                )
            )

        ),

        tabPanel(
            title = "Region(s)",

            fluidRow(

                shiny::column(
                    width = 2, offset = 5,
                    radioButtons(
                        "refreshVariants", "Refresh variants",
                        choices = c("Paused", "Active"),
                        selected = "Active", inline = TRUE)
                )

            ),

            # Genome annotation package ----

            wellPanel(

                fluidRow(

                    shiny::column(
                        width = 3,
                        selectInput(
                            "annotationPackage",
                            "Select installed EnsDb package",
                            choices = as.list(EnsDbPacks),
                            width = '100%')
                    ),

                    shiny::column(
                        3,
                        strong("EnsDb annotation"),
                        htmlOutput("ensembl_organism"),
                        htmlOutput("ensembl_version"),
                        htmlOutput("ensembl_genome")
                    ),

                    shiny::column(
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

            ),

            # GRanges options ----
            wellPanel(

                fluidRow(

                    shiny::column(
                        width = 2,
                        selectInput(
                            "regionInputMode", "Region(s) input type",
                            choices = GS[["choices.regionInputMode"]],
                            selected = GS[["default.regionInputMode"]],
                            width = '100%')
                    ),
                    shiny::column(
                        width = 4, offset = 6,
                        strong("Summary"),
                        htmlOutput("bedFileSummary")
                    )

                ),

                fluidRow(

                    conditionalPanel(
                        condition = "input.regionInputMode == 'bed'",
                        shiny::column(
                            width = 3,
                            fileInput(
                                "bedFile", "BED file",
                                accept = c("bed"))
                        )
                    ),
                    conditionalPanel(
                        condition = "input.regionInputMode == 'ucsc'",
                        shiny::column(
                            width = 7,
                            textInput(
                                "ucscRegions", "UCSC-type region(s)",
                                value = "",
                                placeholder = paste(
                                    "chr21:33,031,597-33,041,570",
                                    "chr2:2,031,597-2,041,570",
                                    "...",
                                    sep = " ; "),
                                width = "100%")
                        )
                    ),
                    conditionalPanel(
                        condition = "input.regionInputMode == 'EnsDb'",
                        fluidRow(
                            shiny::column(
                                width = 2,
                                selectInput(
                                    "ensDb.type", NA,
                                    choices = GS[["choices.ensDbType"]],
                                    selected = GS[["default.ensDbType"]])
                            ),
                            shiny::column(
                                width = 1,
                                selectInput(
                                    "ensDb.condition", NA,
                                    choices = GS[["choices.ensDbFilters"]],
                                    selected = GS[["default.ensDbFilters"]])
                            ),
                            shiny::column(
                                2,
                                textInput(
                                    "ensDb.values", NA,
                                    value = "",
                                    placeholder = "IL17,IL12B,...")
                            )
                        ),
                        tabsetPanel(
                            id = "ensDb.resultTab",
                            selected = "Genes",

                            tabPanel(
                                title = 'Genes',
                                DT::dataTableOutput("ensDb.Genes")
                            )#,
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
            title = "Filters",

            wellPanel(
                fluidRow(

                    shiny::column(
                        width = 2,
                        numericInput(
                            "maf.min", "Minimum MAF",
                            value = 0 , min = 0, max = 0.5,
                            step = 0.01)
                    ),
                    shiny::column(
                        width = 2,
                        numericInput(
                            "maf.max", "Maximum MAF (e.g. 1e-6)",
                            value = 0.5 , min = 0, max = 0.5,
                            step = 0.01)
                    ),
                    shiny::column(
                        width = 4
                    ),
                    shiny::column(
                        width = 4,
                        strong("Summary"),
                        htmlOutput("mafRange")
                    )

                )
            )

        )

    ),

    tabPanel(
        title = "Views",

        tabsetPanel(
            id = "tabset.views",

            # BED view ----
            tabPanel(
                title = "BED",

                fluidRow(
                    shiny::column(
                        width = 12,
                        DT::dataTableOutput("bedSample")
                    )
                )

            ),

            # Variants view ----
            tabPanel(
                title = "Variants",

                fluidRow(
                    shiny::column(
                        width = 12,
                        wellPanel(
                            uiOutput("vcfCols")
                        )
                    )
                ),

                fluidRow(
                    shiny::column(
                        width = 12,
                        DT::dataTableOutput("vcfRowRanges")
                    )
                )

            ),

            # Variants INFO view ----
            tabPanel(
                title = "INFO",

                fluidRow(
                    shiny::column(
                        width = 12,
                        wellPanel(
                            uiOutput("vcfInfoCols")
                        )
                    )
                ),

                fluidRow(
                    shiny::column(
                        width = 12,
                        DT::dataTableOutput("vcfInfo")
                    )
                )

            ),

            # VEP predictions view ----
            tabPanel(
                title = "VEP",

                fluidRow(
                    shiny::column(
                        width = 12,
                        wellPanel(
                            uiOutput("csqCols")
                        )
                    )
                ),

                fluidRow(
                    shiny::column(
                        width = 12,
                        DT::dataTableOutput("csqSample")
                    )
                )

            ),

            # Phenotypes view ----
            tabPanel(
                title = "Phenotypes",

                fluidRow(
                    shiny::column(
                        width = 12,
                        wellPanel(
                            uiOutput("phenoCols")
                        )
                    )
                ),

                fluidRow(
                    shiny::column(
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
                                shiny::column(
                                    width = 6,
                                    uiOutput("genoNumRows")
                                ),
                                shiny::column(
                                    width = 6,
                                    uiOutput("genoFirstRow")
                                )
                            )
                        ),
                        fluidRow(
                            wellPanel(
                                shiny::column(
                                    width = 6,
                                    uiOutput("genoNumCols")
                                ),
                                shiny::column(
                                    width = 6,
                                    uiOutput("genoFirstCol")
                                )
                            )
                        ),
                        fluidRow(
                            shiny::column(
                                width = 12,
                                tableOutput("genotypesSample")
                            )
                        )
                    ),

                    tabPanel(
                        title = "Heatmap",

                        p(
                        strong("This may take some time to plot"),
                        em("(10-20 seconds for 218 variants in 5844 sample)"),
                        "."),

                        p(
                            "Please click the button after loading variants",
                            "to generate/update the figure",
                            actionButton(
                                "doGenoHeatmap", "Go!",
                                icon = icon("time")
                            )
                        ),

                        fluidRow(
                            shiny::column(
                                width = 12,
                                plotOutput("heatmapGenotype")
                            )
                        )

                    )

                )

            )

        )

    ),

    tabPanel(
        title = "Barplot",

        # Variants predictions barplot ----

        sidebarLayout(

            # Sidebar with a slider input
            sidebarPanel(
                width = 3,

                selectInput(
                    "variantCsq",
                    "Variant effect prediction",
                    choices = c()
                ),

                selectInput(
                    "phenotype",
                    "Phenotype field",
                    choices = c("None"),
                    selected = "None"
                ),

                conditionalPanel(
                    condition = "input.phenotype != 'None'",
                    checkboxInput(
                        "csqPhenoUnique",
                        "Unique to phenotype?",
                        value = FALSE
                    )
                ),

                selectInput(
                    "facet",
                    "Faceting field",
                    choices = c("None"),
                    selected = "None"
                ),

                conditionalPanel(
                    condition = "input.facet != 'None'",
                    selectInput(
                        "csqFacets",
                        "Facets",
                        choices = c(),
                        selected = c(),
                        multiple = TRUE
                    )
                ),

                checkboxInput(
                    "stackedPercentage",
                    "Show as percentage?",
                    value = FALSE
                ),

                checkboxInput("advanced", "Advanced controls", value = FALSE),

                conditionalPanel(
                    condition = "input.advanced == true",
                    checkboxInput("legend", "Show legend", value = TRUE)
                ),

                conditionalPanel(
                    condition = paste(
                        "input.advanced == true",
                        "input.legend == true",
                        sep = " && "
                    ),
                    sliderInput(
                        "legendText", "Legend font size",
                        value = 1.5, min = 0.1, max = 2, step = 0.1)
                ),

                conditionalPanel(
                    condition = "input.advanced == true",
                    sliderInput(
                        "xAxisAngle",
                        "Angle of X labels",
                        min = 0, max = 90, value = 0, step = 5
                    ),

                    sliderInput(
                        "xAxisSize",
                        "Relative size of X text",
                        value = 1, min = 0.1, max = 2, step = .1
                    ),

                    sliderInput(
                        "xAxisVjust",
                        "Vertical just. of X labels",
                        min = 0, max = 1, value = 0.5, step = 0.1
                    ),

                    sliderInput(
                        "xAxisHjust",
                        "Horiz. just. of X labels",
                        min = 0, max = 1, value = 0.5, step = 0.1
                    )
                )

            ),

            mainPanel(
                # fluidRow(
                #     shiny::column(
                #         width = 6, offset = 3,
                #         p(
                #             "Stacked bar chart.",
                #             style = "font-family: 'Source Sans Pro';"
                #         )
                #     )
                # ),

                tabsetPanel(
                    id = "csqCount",

                    tabPanel(
                        title = "Barplot",

                        fluidRow(
                            shiny::column(
                                width = 12,
                                plotOutput(
                                    "variantsCsqBarplot",
                                    height = csqCountHeight,
                                    hover = hoverOpts(
                                        "plotVarClass_hover",
                                        delayType = "debounce")
                                )
                            )
                        ),

                        fluidRow(
                            shiny::column(
                                width = 12,
                                htmlOutput("varCsqCount")
                            )
                        )
                    ),

                    tabPanel(
                        title = "Table",
                        DT::dataTableOutput("csqTableDecreasing")
                    )

                )

            )
        )
    ),

    navbarMenu(
        title = "Settings",

        # Advances settings view ----

        tabPanel(
            title = "Advanced",

            wellPanel(
                h4("Genotypes"),
                hr(),

                fluidRow(
                    shiny::column(
                        width = 4,
                        selectInput(
                            "refGenotypes", "REF genotype(s)",
                            choices = setdiff(
                                GS[["all.genotypes"]],
                                c(
                                    GS[["default.hetGenotypes"]],
                                    GS[["default.altGenotypes"]])),
                            selected = GS[["default.refGenotypes"]],
                            multiple = TRUE
                        )
                    ),

                    shiny::column(
                        width = 4,
                        selectInput(
                            "hetGenotypes", "ALT heterozygote genotype(s)",
                            choices = setdiff(
                                GS[["all.genotypes"]],
                                c(
                                    GS[["default.refGenotypes"]],
                                    GS[["default.altGenotypes"]])),
                            selected = GS[["default.hetGenotypes"]],
                            multiple = TRUE
                        )
                    ),

                    shiny::column(
                        width = 4,
                        selectInput(
                            "altGenotypes", "ALT homozygote genotype(s)",
                            choices = setdiff(
                                GS[["all.genotypes"]],
                                c(
                                    GS[["default.refGenotypes"]],
                                    GS[["default.hetGenotypes"]])),
                            selected = GS[["default.altGenotypes"]],
                            multiple = TRUE
                        )
                    )

                )

            ),

            wellPanel(
                h4("VCF file(s)"),
                hr(),

                fluidRow(

                    shiny::column(
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
                "Experimental"),
            tags$sup("1,2"),

            wellPanel(

                fluidRow(

                    shiny::column(
                        width = 3,
                        numericInput(
                            "bpCores", "Cores",
                            value = PS[["default.bpCores"]],
                            min = 1, max = PS[["default.bpCores"]], step = 1)
                    ),

                    shiny::column(
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
                        shiny::column(
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
                "Parallel workers validated on",
                tags$code("Ubuntu"),
                "and",
                tags$code("Scientific Linux"),
                "systems."
            ),
            p(
                tags$sup("2"),
                "Known issue on Mac OS X El Capitan:",
                "App hangs while CPUs infinitely working at full speed."
            )
        )

    ),

    # Session view ----

    tabPanel(
        title = "Session",

        tabsetPanel(
            id = "tabset.session",

            tabPanel(
                title = "Session info",
                verbatimTextOutput("sessionInfo")
            ),

            tabPanel(
                title = "tSVE settings",
                verbatimTextOutput("tSVESettings")
            ),

            tabPanel(
                title = "Advanced settings",
                verbatimTextOutput("advancedSettings")
            ),

            tabPanel(
                title = "BED",
                verbatimTextOutput("bedStructure")
            ),

            tabPanel(
                title = "Variants",
                verbatimTextOutput("vcf")
            ),

            tabPanel(
                title = "VEP",
                verbatimTextOutput("csqStructure")
            ),

            tabPanel(
                title = "Phenotypes",
                verbatimTextOutput("phenotypeStructure")
            ),

            tabPanel(
                title = "Genotypes",
                verbatimTextOutput("genotypeStructure")
            )
        )

    )

))
