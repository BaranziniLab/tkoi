###########################
# Author: Wanjun Gu
# Email: wanjun.gu@ucsf.edu
# Date: 2026-04-20
###########################

library(tkoi)

# Load example data file path from the tkoi package
file_path = system.file("extdata", "example_data.csv", package = "tkoi")

# Read gene expression data from the CSV file using data.table for fast processing
expression_data = data.table::fread(file_path)

# Perform tKOI analysis on the loaded gene expression data
tkoi_result = run_tkoi(
  expression_data = expression_data,
  subnetwork = tkoi::tkoi_net,
  pvalue_threshold = 0.05,
  logfc_threshold = 0.25,
  indirect_link_threshold = 3,
  topology_similarity = 0.9,
  n_permutation = 100,
  damping_factor = 0.85,
  maximum_iteration = 500
)

# Perform Gene Ontology enrichment analysis and integrate results with tKOI data
tkoi_result = run_gene_enrichment(tkoi_result)

# Access the first scatter plot visualizing tKOI FDR vs enrichment q-values
tkoi_result@gene_enrichment_comparison$comparison_scatter1

# Access the second scatter plot visualizing the same relationship with facets for namespaces
tkoi_result@gene_enrichment_comparison$comparison_scatter2

# save(tkoi_result, file = "tkoi_result.rda")

# Generate a scatter plot to explore gene-level data based on experimental p-values and tKOI-adjusted FDR
plt1 = make_gene_exploration_plot(
  tkoi_list = tkoi_result,
  sig_color = "#F39B7FB2",
  non_sig_color = "gray"
)

# Export detailed gene exploration data as a data frame
gene_data = export_gene_exploration_data(tkoi_result)

# Visualize the top 25 Biological Processes with the highest network enrichment statistics
plt2 = visualize_topn(
  tkoi_list = tkoi_result,
  category = "Gene",
  top_n = 25,
  high_color = "#FF5733",
  low_color = "#154360"
)
