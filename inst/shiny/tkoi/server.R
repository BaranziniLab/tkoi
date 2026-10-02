# Tests can exercise the real tkoi engine on a small injected graph.
# The running app uses the graph from the installed tkoi library.
create_server = function(subnetwork = tkoi::tkoi_net) {
  force(subnetwork)
  function(input, output, session) {
    expression_data = reactiveVal(example_data)
    tkoi_output = reactiveVal(NULL)
    analysis_parameters = reactiveVal(NULL)
    output$graph_components_table = renderTable({
      data.frame(
        `NodeType` = c("Anatomy", "BiologicalProcess", "CellType", "CellularComponent", "ClinicalLab",
                       "Complex", "Compound", "Disease", "EC", "Gene", "MiRNA", "MolecularFunction",
                       "Pathway", "Protein", "ProteinDomain", "ProteinFamily", "PwGroup", "Reaction"),
        Description = c("Anatomical structures (e.g., heart, liver)",
                        "GO Biological Process terms",
                        "Cell types (e.g., T cell, neuron)",
                        "GO Cellular Component terms",
                        "Clinical lab tests (e.g., LDL, CRP)",
                        "Protein complexes",
                        "Metabolites and small molecules",
                        "Disease terms (e.g., rheumatoid arthritis)",
                        "Enzyme Commission categories",
                        "Human genes",
                        "MicroRNAs",
                        "GO Molecular Function terms",
                        "Biological pathways (Reactome, KEGG)",
                        "Proteins",
                        "Protein domains (e.g., SH3)",
                        "Protein families",
                        "Pathway groupings",
                        "Biochemical reactions")
      )
    })


    output$tkoi_result_ui = renderUI({
      req(tkoi_output())
      tagList(
        selectInput("result_category", "Select Node Type:",
          choices = names(tkoi_output()@network_summary_statistics),
          selected = "BiologicalProcess"
        ),
        DT::dataTableOutput("tkoi_result_table"),
        tags$hr(),
        fluidRow(
          column(width = 4, downloadButton("download_rda", "Download .rda")),
          column(width = 4, downloadButton("download_excel", "Download .xlsx")),
          column(width = 4, downloadButton("download_excel_significant", "Download tKOIAgent Context"))
        ),
        br()
      )
    })
    output$tkoi_network_visualization_lookup = renderUI({
      req(tkoi_output())
      tagList(
        selectInput("lookup_category", "Select Node Type:",
          choices = names(tkoi_output()@network_summary_statistics),
          selected = "BiologicalProcess"
        ),
        DT::dataTableOutput("tkoi_result_table_lookup"),
        br()
      )
    })
    output$tkoi_network_visualization = renderUI({
      req(tkoi_output())
      tagList(
        textInput("target_node_id", "Node ID", placeholder = "Paste a node_id from the lookup table"),
        numericInput("degree_expansion", "Degree Expansion", min = 1, max = 5, value = 2),
        numericInput("resolution", "Resolution", min = 80, max = 500, value = 150),
        selectInput("layout_type", "Network Layout",
          choices = c("Kamada-Kawai" = "kk", "Fruchterman-Reingold" = "fr",
            "GEM" = "gem", "Graphopt" = "graphopt", "LGL" = "lgl", "MDS" = "mds"),
          selected = "kk"
        ),
        actionButton("plot_network", "Plot Network"),
        plotOutput("network_plot")
      )
    })
    observeEvent(input$plot_network, {
      result = tkoi_output()
      req(result, input$target_node_id)
      target = input$target_node_id
      expansion = input$degree_expansion
      layout = input$layout_type
      output$network_plot = renderPlot({
        validate(need(target %in% igraph::V(subnetwork)$name, "Enter a node_id from the lookup table."))
        tkoi::plot_network(
          tkoi_result = result,
          target_node_id = target,
          degree_expansion = expansion,
          network_layout_type = layout,
          subnetwork = subnetwork
        )
      }, height = 600, res = input$resolution)
    })
    output$download_rda = downloadHandler(
      filename = function() "tkoi_result.rda",
      content = function(file) {
        req(tkoi_output())
        tkoi_result = tkoi_output()
        tkoi_engine_version = tkoi_version
        save(tkoi_result, tkoi_engine_version, file = file)
      },
      contentType = "application/octet-stream"
    )
    export_excel = function(file, significant_only = FALSE) {
      req(tkoi_output())
      workbook = openxlsx::createWorkbook()
      summary = tkoi_output()@network_summary_statistics
      for (name in names(summary)) {
        data = summary[[name]]
        if (significant_only) {
          data = dplyr::filter(data, fdr <= 0.05)
        }
        sheet = gsub("[\\/:*?\"<>|]", "_", substr(name, 1, 31))
        write_result_worksheet(workbook, sheet, data)
      }
      openxlsx::saveWorkbook(workbook, file, overwrite = TRUE)
    }
    output$download_excel = downloadHandler(
      filename = function() "tkoi_result.xlsx",
      content = function(file) export_excel(file)
    )
    output$download_excel_significant = downloadHandler(
      filename = function() "tKOIAgent Context.xlsx",
      content = function(file) export_excel(file, significant_only = TRUE)
    )
    output$tkoi_result_table = DT::renderDataTable({
      req(tkoi_output(), input$result_category)
      data = tkoi_output()@network_summary_statistics[[input$result_category]]
      req(data)
      probability_datatable(data, page_length = 10)
    })
    output$tkoi_result_table_lookup = DT::renderDataTable({
      req(tkoi_output(), input$lookup_category)
      data = tkoi_output()@network_summary_statistics[[input$lookup_category]]
      req(data)
      probability_datatable(data, page_length = 5)
    })
    observeEvent(input$run_tkoi_analysis, {
      tkoi_output(NULL)
      analysis_parameters(NULL)
      start_time = Sys.time()
      result = tryCatch({
        # Read the current upload before analysis; examples and uploads both
        # run through the same installed package.
        data = read_expression_data(input$upload_data, example = example_data)
        parameters = list(
          pvalue_threshold = input$pvalue_threshold,
          logfc_threshold = input$logfc_threshold,
          indirect_link_threshold = input$indirect_link_threshold,
          topology_similarity = input$topology_similarity,
          n_permutation = input$n_permutation,
          damping_factor = input$damping_factor,
          maximum_iteration = input$maximum_iteration,
          n_cores = input$n_cores
        )
        expression_data(data)
        analysis_parameters(parameters)
        withProgress(message = "Running tKOI analysis", value = 0.1, {
          setProgress(0.1, detail = paste("Using tkoi", tkoi_version))
          result = tkoi::run_tkoi(
            expression_data = data,
            subnetwork = subnetwork,
            pvalue_threshold = parameters$pvalue_threshold,
            logfc_threshold = parameters$logfc_threshold,
            indirect_link_threshold = parameters$indirect_link_threshold,
            topology_similarity = parameters$topology_similarity,
            n_permutation = parameters$n_permutation,
            damping_factor = parameters$damping_factor,
            maximum_iteration = parameters$maximum_iteration,
            n_cores = parameters$n_cores,
            keep_permutations = FALSE,
            verbose = FALSE
          )
          setProgress(1, detail = "Analysis complete")
          result
        })
      }, error = function(error) {
        message("tKOI analysis failed: ", conditionMessage(error))
        analysis_parameters(NULL)
        showNotification(conditionMessage(error), type = "error", duration = NULL)
        NULL
      })
      if (is.null(result)) {
        return(invisible(NULL))
      }
      tkoi_output(result)
      elapsed = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      showModal(modalDialog(
        title = "tKOI analysis complete.",
        glue::glue(
          "Computed with tkoi {tkoi_version} in {round(elapsed, 2)} seconds ",
          "using {result@n_permutation} permutations."
        ),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
    output$download_example = downloadHandler(
      filename = function() "example_data.csv",
      content = function(file) {
        data = example_data
        data$pvalue = formatC(data$pvalue, format = "e", digits = 16)
        write.csv(data, file, row.names = FALSE)
      }
    )
    processed_data = reactive({
      parameters = analysis_parameters()
      req(parameters)
      dplyr::mutate(expression_data(),
        nlogp = -log10(pvalue),
        regulation = dplyr::case_when(
          pvalue <= parameters$pvalue_threshold & logfc >= parameters$logfc_threshold ~ "Up",
          pvalue <= parameters$pvalue_threshold & logfc <= -parameters$logfc_threshold ~ "Down",
          TRUE ~ "Not Significant"
        )
      )
    })
    output$expression_table = DT::renderDataTable({
      req(processed_data())
      data = processed_data() |>
        dplyr::arrange(pvalue) |>
        dplyr::inner_join(tkoi::genes, by = dplyr::join_by("gene_name" == "ensembl")) |>
        dplyr::mutate(
          logfc = formatC(logfc, format = "e", digits = 3)
        ) |>
        dplyr::select(gene_name, symbol = name, logfc, pvalue)
      probability_datatable(data, page_length = 10)
    })
    output$significance_summary = renderTable({
      req(processed_data())
      data = processed_data()
      data.frame(
        TotalGenes = nrow(data),
        Significant = sum(data$regulation != "Not Significant", na.rm = TRUE),
        Upregulated = sum(data$regulation == "Up", na.rm = TRUE),
        Downregulated = sum(data$regulation == "Down", na.rm = TRUE)
      )
    })
    output$volcano_ui = renderUI({
      req(processed_data())
      plotly::plotlyOutput("volcano_plot", height = "500px")
    })
    output$volcano_plot = plotly::renderPlotly({
      req(processed_data())
      parameters = analysis_parameters()
      colors = c("Up" = "#bb2a2d", "Down" = "#4e5d71", "Not Significant" = "gray")
      plot = ggplot(processed_data(), aes(x = logfc, y = nlogp, color = regulation, text = gene_name)) +
        geom_point(alpha = 0.7) +
        scale_color_manual(values = colors) +
        geom_hline(yintercept = -log10(parameters$pvalue_threshold), linetype = "dashed") +
        geom_vline(xintercept = c(-parameters$logfc_threshold, parameters$logfc_threshold), linetype = "dashed") +
        labs(x = "Log Fold Change", y = "-log10(P-value)", color = "") +
        theme_minimal() +
        theme(text = element_text(size = 15))
      plotly::ggplotly(plot, tooltip = c("x", "y", "text", "color"))
    })
  }
}
server = create_server()
