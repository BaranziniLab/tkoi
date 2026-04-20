# Run tKOI Analysis

Performs tKOI (Transcriptomic Knowledge-graph-driven Omics Integration)
analysis on a gene expression dataset. The analysis integrates
transcriptomic measurements with a biological knowledge graph to
identify and interpret biologically enriched subnetworks. The pipeline
includes the following steps:

## Usage

``` r
run_tkoi(
  expression_data,
  subnetwork = tkoi::tkoi_net,
  pvalue_threshold = 0.05,
  logfc_threshold = 0.25,
  indirect_link_threshold = 3,
  topology_similarity = 0.9,
  n_permutation = 100,
  damping_factor = 0.85,
  maximum_iteration = 500,
  n_cores = 1
)
```

## Arguments

- expression_data:

  A data frame with columns `gene_name` (Ensembl IDs), `logfc`, and
  `pvalue`.

- subnetwork:

  An igraph object representing the biological knowledge graph. Default
  is [`tkoi::tkoi_net`](tkoi_net.md).

- pvalue_threshold:

  Numeric threshold for filtering genes by p-value. Default `0.05`.

- logfc_threshold:

  Numeric threshold for filtering genes by absolute log fold change.
  Default `0.25`.

- indirect_link_threshold:

  Minimum number of seed genes that must be 2-hop-connected to a node
  for it to be prioritized. Default `3`.

- topology_similarity:

  Numeric in (0, 1) controlling how similar in degree a substitute gene
  must be to the original during permutation sampling. Default `0.9`.

- n_permutation:

  Integer number of permutations for the null distribution. Default
  `100`.

- damping_factor:

  Damping factor for personalized PageRank. Default `0.85`.

- maximum_iteration:

  Maximum iterations for PageRank convergence. Default `500`.

- n_cores:

  Number of cores for parallel permutation computation. Default `1`
  (sequential). Values `> 1` use
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) and
  require macOS or Linux; note that parallel mode produces statistically
  equivalent but numerically different results due to independent
  per-worker RNG streams.

## Value

An S4 object of class `tKOIList` with slots:

- `expression_data`:

  The input gene expression data frame.

- `pagerank_data`:

  A data frame of observed and permutation PageRank scores for every
  node in the subnetwork.

- `network_summary_statistics`:

  A named list of data frames, one per node type, containing beta,
  p_value, fdr, and biological annotations.

## Details

1.  **Gene filtering:** Significant genes are selected based on
    user-defined thresholds for p-value and log fold change.

2.  **Network mapping:** Filtered genes are mapped onto the biological
    subnetwork using Ensembl gene identifiers.

3.  **Personalized PageRank:** A PageRank propagation is performed using
    log fold change-weighted probabilities, quantifying each node's
    influence.

4.  **Permutation testing:** The propagation is repeated `n_permutation`
    times using substitute genes with similar topological properties to
    build a null distribution.

5.  **Network enrichment scoring:** Each node's z-score (beta) and
    p-value are computed by comparing observed PageRank to the
    permutation null.

6.  **Node annotation:** Nodes are annotated with domain-specific
    metadata (GO terms, diseases, cell types, compounds) from curated
    knowledge resources.

7.  **Prioritization:** Results are organized by node type with FDR
    adjustment and connectivity statistics for biological
    prioritization.

Probability vectors and metapath filters are built with vectorized
operations ([`match()`](https://rdrr.io/r/base/match.html), batch
[`igraph::ego()`](https://r.igraph.org/reference/ego.html)) rather than
per-gene dplyr operations. Candidate gene pools for permutation sampling
are pre-computed once before the loop. Enrichment z-scores are computed
with `rowMeans`/`sweep`/`rowSums` across the full permutation matrix
rather than row-by-row.

## See also

[`page_rank`](https://r.igraph.org/reference/page_rank.html),
[`compute_network_enrichment`](compute_network_enrichment.md),
[`run_gene_enrichment`](run_gene_enrichment.md)

## Examples

``` r
if (FALSE) { # \dontrun{
expression_data = data.frame(
  gene_name = c("ENSG00000000003", "ENSG00000000005"),
  logfc = c(1.5, -0.8),
  pvalue = c(0.01, 0.02)
)

result = run_tkoi(
  expression_data = expression_data,
  subnetwork = tkoi::tkoi_net,
  pvalue_threshold = 0.05,
  logfc_threshold = 0.25,
  n_permutation = 100
)

result@pagerank_data
result@network_summary_statistics
} # }
```
