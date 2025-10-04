library(tkoi)

# Load example data file path from the tkoi package
file_path = system.file("extdata", "example_data.csv", package = "tkoi")

# Read gene expression data from the CSV file using data.table for fast processing
expression_data = data.table::fread(file_path)

# Perform tKOI analysis on the loaded gene expression data
# This calculates PageRank scores, performs permutation tests, and integrates the data with the network
tkoi_result = run_tkoi(
  expression_data = expression_data,
  subnetwork = tkoi::tkoi_net,  # Predefined network in the tkoi package
  pvalue_threshold = 0.05,      # Threshold for filtering significant genes based on p-value
  logfc_threshold = 0.25,       # Threshold for log fold change significance
  indirect_link_threshold = 3,  # A numeric value specifying the threshold for the least amount of genes to be 2-hop-connected to the propagated node.
  topology_similarity = 0.9,    # Degree of similarity for substitute genes during permutation
  n_permutation = 100,          # Number of permutations for statistical tests
  damping_factor = 0.85,        # Damping factor for PageRank computation
  maximum_iteration = 500       # Maximum iterations for PageRank convergence
)

# Perform Gene Ontology enrichment analysis and integrate results with tKOI data
# The results include enriched terms, FDR values, and scatter plots comparing enrichment and tKOI statistics
tkoi_result = run_gene_enrichment(tkoi_result)

# Access the first scatter plot visualizing tKOI FDR vs enrichment q-values
tkoi_result@gene_enrichment_comparison$comparison_scatter1

# Access the second scatter plot visualizing the same relationship with facets for namespaces
tkoi_result@gene_enrichment_comparison$comparison_scatter2

# Uncomment the following line to save the tKOI analysis results for future use
# save(tkoi_result, file = "tkoi_result.rda")

# Generate a scatter plot to explore gene-level data based on experimental p-values and tKOI-adjusted FDR
plt1 = make_gene_exploration_plot(
  tkoi_list = tkoi_result,       # tKOI analysis results
  sig_color = "#F39B7FB2",       # Color for significant genes
  non_sig_color = "gray"         # Color for non-significant genes
)

# Export detailed gene exploration data as a data frame
gene_data = export_gene_exploration_data(tkoi_result)

# Visualize the top 25 Biological Processes with the highest network enrichment statistics
plt2 = visualize_topn(
  tkoi_list = tkoi_result,        # tKOI analysis results
  category = "Gene", # Category to visualize (e.g., BiologicalProcess, Gene)
  top_n = 25,                     # Number of top terms to display
  high_color = "#FF5733",         # Color representing high FDR significance
  low_color = "#154360"           # Color representing low FDR significance
)

