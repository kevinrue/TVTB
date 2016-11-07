
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(navbarPage(
  theme = "bootstrap.css",

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
            strong("File"),
            br(),
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

      ),

      fluidRow(
        column(
          width = 12,
          tags$h4("Notes"),
          "The phenotype file must be formatted as follows:",
          tags$ul(
            tags$li(
              "Fields separator must be 'white space', the default",
              tags$a(
                href="http://stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html",
                tags$code("read.table")
              ),
              "field separator."
            ),
            tags$li(
              "First row must be phenotype names."
            ),
            tags$li(
              "First column must be samples identifiers matching those in",
              "the VCF file(s)."
            )
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
              choices = list(
                "BED file" = "bed",
                "UCSC browser" = "ucsc",
                "EnsDb package" = "EnsDb"
              ),
              selected = "bed",
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
                "ensDb.type", "Type",
                choices = list("Gene name" = "Genename"),
                selected = "Genename")
            ),
            column(
              width = 1,
              selectInput(
                "ensDb.condition", "Condition",
                choices = c("=", "!=", "like", "in"),
                selected = "=")
            ),
            column(
              2,
              textInput(
                "ensDb.value", "Value",
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
              choices = list(
                "Single VCF" = "SingleVcf",
                "One per chromosome" = "OnePerChr"
              ),
              selected = "OnePerChr",
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
                value = system.file("extdata", package = "TVTB"),
                width = '100%',
                placeholder = "/path/to/VCF/folder")
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
                  value = "^chr%s\\..*\\.vcf\\.gz$",
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

          # INFO fields (except VEP) ----

          column(
            width = 8,
            fluidRow(
              column(
                width = 12,
                selectInput(
                  "vcfInfoKeys", "INFO fields",
                  choices = character(),
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

          # VEP prediction INFO field ----

          column(
            width = 2,
            textInput(
              "vepKey", "VEP field (INFO)",
              value = get("vepKey", .tSVE),
              placeholder = 'CSQ, ANN, ...')
          ),

          # FORMAT fields ----

          column(
            width = 2,
            selectInput(
              "vcfFormatKeys", "FORMAT fields",
              choices = character(),
              multiple = TRUE),
            strong("Note:"), "\"GT\" implicitely required"
          )

        )
      ),

      wellPanel(
        fluidRow(

          # VCF import button! ----

          column(
            width = 2,
            checkboxInput(
              "autodetectGTimport", "Autodetect genotypes",
              value = get("autodetectGTimport", .tSVE)
            )
          ),

          column(
            width = 2, offset = 3,
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
              choices = as.list(.EnsDbPacks),
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

    uiOutput("TVTBparamWarning"),

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
            choices = character()
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
            choices = list(
              "fixed" = "VcfFixedRules",
              "info" = "VcfInfoRules",
              "VEP" = "VcfVepRules"
            ),
            selected = "VcfFixedRules"
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
          width = 7,
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
          width = 2,
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
            DT::dataTableOutput("rangesTableView")
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
            DT::dataTableOutput("vcfRowRangesView")
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
              uiOutput("vcfInfoCols"),
              br(),
              p(strong("Notes:")),
              tags$ul(
                tags$li(
                  "Fields that contain more than one value",
                  "(", tags$em("e.g."), "confidence intervals)",
                  "may not display properly."
                )
              )
            )
          )
        ),

        fluidRow(
          column(
            width = 12,
            DT::dataTableOutput("vcfInfoView")
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
            DT::dataTableOutput("vcfVepView")
          )
        )

      ),

      # Phenotypes view ----
      tabPanel(
        title = "Phenotypes",

        fluidRow(
          column(
            width = 12,
            "This panel displays phenotype information attached to",
            "the imported VCF object.",
            wellPanel(
              uiOutput("phenoCols")
            )
          )
        ),

        fluidRow(
          column(
            width = 12,
            DT::dataTableOutput("phenotypesView")
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
                plotOutput(
                  "heatmapGenotype",
                  height = get("genoHeatmap.height", .tSVE)
                )
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

          ),

          tabPanel(
            title = "Info",
            shiny::h4("Encoding"),
            uiOutput("genotypeEncoding")
          )

        )

      )

    )

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
            width = 1,
            br(),
            actionButton(
              "genotypeAutofill", "Autofill", icon("magic")
            )
          ),
          column(
            width = 3,
            selectInput(
              "refGenotypes", "Reference homozygote genotype(s)",
              choices = c(get("refGT", .tSVE), get("hetGT", .tSVE), get("altGT", .tSVE)),
              selected = get("refGT", .tSVE),
              multiple = TRUE
            )
          ),

          column(
            width = 4,
            selectInput(
              "hetGenotypes", "Heterozygote genotype(s)",
              choices = c(get("refGT", .tSVE), get("hetGT", .tSVE), get("altGT", .tSVE)),
              selected = get("hetGT", .tSVE),
              multiple = TRUE
            )
          ),

          column(
            width = 4,
            selectInput(
              "altGenotypes", "Alternate homozygote genotype(s)",
              choices = c(get("refGT", .tSVE), get("hetGT", .tSVE), get("altGT", .tSVE)),
              selected = get("altGT", .tSVE),
              multiple = TRUE
            )
          )

        ),

        fluidRow(
          column(
            width = 1,
            textInput(
              "refSuffix", "Suffix",
              value = get("refSuffix", .tSVE),
              placeholder = get("refSuffix", .tSVE)
            )
          ),
          column(
            width = 1, offset = 3,
            textInput(
              "hetSuffix", "Suffix",
              value = get("hetSuffix", .tSVE),
              placeholder = get("hetSuffix", .tSVE)
            )
          ),
          column(
            width = 1, offset = 3,
            textInput(
              "altSuffix", "Suffix",
              value = get("altSuffix", .tSVE),
              placeholder = get("altSuffix", .tSVE)
            )
          )
        ),

        fluidRow(
          column(
            width = 12,
            tags$strong("Notes:"), br(),
            tags$ul(
              tags$li(
                "The",tags$strong("choices"),"of genotypes are updated when",
                "new variants are imported."
              ),
              tags$li(
                "The",tags$strong("selected"),"genotypes may be automatically",
                "updated immediately after import using the",
                tags$strong("Autodetect genotypes"), "checkbox in the",
                tags$strong("Input"), "panel, or manually after import using",
                "the", tags$strong("Autofill"), "button in this panel."
              ),
              tags$li(
                "Selected genotypes are not allowed to overlap.",
                "Selecting a genotype removes it from the choices",
                "available in the other widgets. As a consequence, genotypes",
                "must first be unselected from a widget before it can be",
                "selected in another one."
              )
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
              value = get("aafSuffix", .tSVE),
              placeholder = get("aafSuffix", .tSVE)
            )
          ),
          column(
            width = 3,
            textInput(
              "mafSuffix", "Minor allele freq.",
              value = get("mafSuffix", .tSVE),
              placeholder = get("mafSuffix", .tSVE)
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

      wellPanel(
        h4("Parallel settings"),
        hr(),
        fluidRow(

          column(
            width = 3,
            numericInput(
              "bpCores", "Cores",
              value = .PS[["default.bpCores"]],
              min = 1, max = .PS[["default.bpCores"]], step = 1)
          ),

          column(
            width = 3,
            selectInput(
              "bpConfig", "Cluster configuration",
              choices = structure(
                .PS[["choices.bpClass"]],
                names = gsub(
                  "Param", "", .PS[["choices.bpClass"]])),
              selected = .PS[["default.bpClass"]])
          ),

          conditionalPanel(
            condition = "input.bpConfig != 'SerialParam'",
            column(
              width = 3,
              selectInput(
                "bpType", "Cluster type",
                choices = structure(
                  .PS[["choices.bpType"]],
                  names = gsub(
                    "Param", "", .PS[["choices.bpType"]])),
                selected = .PS[["default.bpType"]])
            )
          )

        ) # fluidRow

      ), # wellPanel
      wellPanel(
        fluidRow(
          column(
            width = 12,
            h1("Platforms tested"),
            DT::dataTableOutput("parallelReport")
          )
        )
      ),
      tags$h4(
        "Notes",
        tags$ul(
          tags$li(
            "Report"
          ), br(),
          tags$ul(
            tags$li(
              tags$strong("Hang:"),
              "Application hangs while CPUs work infinitely at full capacity."
            )
          )
        )

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
        title = "ExpandedVCF",
        "This panel displays the structure of the imported",
        tags$code("ExpandedVCF"), "object:",
        verbatimTextOutput("ExpandedVCF"),
        "and the attached", tags$code("metadata"), ":",
        verbatimTextOutput("vcfMetadata")
      ),

      tabPanel(
        title = "VEP",
        verbatimTextOutput("vepStructure")
      ),

      tabPanel(
        title = "Errors",
        verbatimTextOutput("Errors")
      )
    )

  )

))
