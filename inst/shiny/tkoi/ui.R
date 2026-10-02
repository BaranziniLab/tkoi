ui = fluidPage(
  theme = shinytheme("journal"),
  tags$head(
    tags$title(glue::glue("tKOI | tkoi {tkoi_version}")),
    tags$link(rel = "icon", type = "image/x-icon", href = "tkoi_favicon.ico"),
    tags$style(HTML("
      body, label, input, button, select {
        font-family: Arial, sans-serif !important;
      }
    "))
  ),

  div(
    style = "display: flex; align-items: center; gap: 15px; padding-bottom: 2px; padding-top: 10px;",
    tags$img(src = "tkoi_logo.png", height = "100px"),
    div(tags$h4("Transcriptomic Knowledge-graph-driven Omics Integration", style = "margin: 0; color: #555;"),
        tags$p(glue::glue("Powered by tkoi {tkoi_version}"), style = "margin: 4px 0 0; color: #555;"))
  ),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      fileInput("upload_data",
                "Upload Expression Data",
                accept = c(".csv", ".tsv", ".xlsx"),
                placeholder = "Use Example Data",
                multiple = FALSE),
      helpText(paste(
        "Upload a data file with columns for gene identifiers (gene_name), log fold changes",
        "(logfc), and p-values (pvalue)."
      )),
      downloadButton("download_example", "Download Example"),
      helpText("Download the example input data. Run tKOI to compute fresh results."),
      numericInput("pvalue_threshold", "P-value Threshold",
        value = tkoi_defaults$pvalue_threshold, min = 0, max = 1, step = 0.001),
      helpText("Genes with p-values at or below this threshold are included."),
      numericInput("logfc_threshold", "Log Fold Change Threshold",
        value = tkoi_defaults$logfc_threshold, min = 0, step = 0.01),
      helpText(paste(
        "Threshold for filtering genes based on absolute log fold change. Genes below this",
        "threshold will be excluded."
      )),
      actionButton("run_tkoi_analysis", label = "Run tKOI"),
      tags$hr(),
      h4("Advanced Options"),
      numericInput("indirect_link_threshold", "Minimum Seed Support (Ranking)",
        value = tkoi_defaults$indirect_link_threshold, min = 0, max = Inf, step = 1),
      helpText(paste(
        "Prioritizes nodes connected to at least this many seed genes within two hops. Results are",
        "ranked, not filtered."
      )),
      numericInput("topology_similarity", "Topology Similarity",
        value = tkoi_defaults$topology_similarity, min = 0, max = 1, step = 0.01),
      helpText(paste(
        "Similarity threshold (0-1) used for selecting substitute genes with similar network",
        "topology during permutation."
      )),
      numericInput("n_permutation", "Number of Permutations", value = tkoi_defaults$n_permutation, min = 2, step = 1),
      helpText("Number of permutations to run for the network enrichment test."),
      numericInput("damping_factor", "Damping Factor",
        value = tkoi_defaults$damping_factor, min = 0, max = 1, step = 0.01),
      helpText("Damping factor used in the personalized PageRank algorithm."),
      numericInput("maximum_iteration", "Maximum Iterations",
        value = tkoi_defaults$maximum_iteration, min = 1, step = 1),
      helpText("Maximum number of iterations allowed for PageRank convergence."),
      numericInput("n_cores", "CPU Cores", value = 4, min = 1, step = 1),
      helpText("Number of CPU cores used by the tkoi library.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Welcome",
                 br(),
                 h3("Welcome to tKOI"),
                 p(paste(
                   "Transcriptomic Knowledge-graph-driven Omics Integration (tKOI) is a computational",
                   "framework designed to enhance the interpretation of transcriptomic data by leveraging a",
                   "large-scale, human-specific heterogeneous biological knowledge graph."
                 )),

                 h4("What is tKOI?"),
                 p(paste(
                   "tKOI stands for Transcriptomic Knowledge-graph-driven Omics Integration. It integrates",
                   "RNA-seq expression data with a large, curated biological knowledge graph. By applying",
                   "network propagation and statistical permutation, tKOI can highlight biologically relevant",
                   "features, even when genes are not statistically significant by conventional criteria."
                 )),

                 h4("What is tKOI used for?"),
                 p(paste(
                   "tKOI is used to identify pathways, cell types, disease associations, and other biological",
                   "entities enriched in transcriptomic data. It can rescue biologically meaningful but",
                   "statistically underpowered genes and contextualize expression data across 18 types of",
                   "biological concepts."
                 )),

                 h4("What are the components of the tKOI knowledge graph?"),
                 p(paste(
                   "The knowledge graph includes 18 node types representing major biological concepts, all",
                   "derived from the SPOKE database and related resources:"
                 )),
                 tableOutput("graph_components_table"),

                 h4("Analysis method"),
                 p(glue::glue(
                   "Network propagation, permutation testing, statistical inference, and annotation are ",
                   "performed by the installed tkoi {tkoi_version} library."
                 )),
                 p(tags$a("Read the tkoi method documentation",
                          href = "https://baranzinilab.github.io/tkoi/",
                          target = "_blank", rel = "noopener")),

                 h4("What does the input data look like and how should it be interpreted?"),
                 p("Input data is a differential gene expression table with three columns:"),
                 tags$ul(
                   tags$li("gene_name – Ensembl gene identifier"),
                   tags$li("logfc – log2 fold change between conditions"),
                   tags$li("pvalue – p-value from statistical test")
                 ),
                 p(paste(
                   "tKOI uses logfc and pvalue thresholds to determine input genes. It then projects these",
                   "genes into the network and ranks all biological entities based on how strongly they are",
                   "connected to the expression signal."
                 )),

                 h4("About the Author"),
                 p("This tool was developed by Wanjun Gu, a member of the ",
                   tags$a(href = "https://baranzinilab.ucsf.edu/", "Baranzini Lab", target = "_blank"),
                   " at the University for California, San Francisco (UCSF)."),
                 p("For questions, feedback, or collaborations, please contact Wanjun at ",
                   tags$a(href = "mailto:wanjun.gu@ucsf.edu", "wanjun.gu@ucsf.edu"), ".")

        ),

        tabPanel("Expression Data",
                 br(),
                 h4("Volcano Plot"),
                 uiOutput("volcano_ui"),
                 helpText("The volcano plot shows gene-level changes:
                          log fold change on the x-axis and –log(p-value) on the y-axis.
                          Red and blue points represent significantly up- and downregulated genes,
                          respectively, based on user-defined thresholds."),

                 tags$hr(),
                 h4("Significance Summary"),
                 tableOutput("significance_summary"),
                 helpText("Summary of total genes analyzed and those significantly
                          up- or downregulated based on selected thresholds."),
                 tags$hr(),
                 h4("Expression Table"),
                 DT::dataTableOutput("expression_table"),
                 helpText("Table of genes mapped onto the gene network, ranked by p-value.
                          Only genes present in the knowledge graph are included.")

        ),

        tabPanel("tKOI Result",
                 br(),
                 h4("tKOI Result Output"),
                 uiOutput("tkoi_result_ui")
        ),

        tabPanel("Visualize Network",
                 br(),
                 h4("tKOI Result Lookup Table"),
                 helpText("Please find the target node to visualize and copy the node_id to the visualizer below."),
                 uiOutput("tkoi_network_visualization_lookup"),
                 h4("Network Visualization"),
                 helpText("In a browser window, right click the generated image to copy or save."),
                 uiOutput("tkoi_network_visualization")
        )
      )
    )
  )
)
