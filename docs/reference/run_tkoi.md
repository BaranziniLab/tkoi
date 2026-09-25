# Run tKOI Analysis

Performs tKOI (Transcriptomic Knowledge-graph-driven Omics Integration)
analysis on a differential expression table. Transcriptomic signals are
propagated through a biological knowledge graph with personalized
PageRank, and a degree-matched permutation null identifies enriched
concepts. The pipeline has these steps:

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
  n_cores = NULL,
  keep_permutations = TRUE,
  tolerance = 1e-14,
  verbose = TRUE
)
```

## Arguments

- expression_data:

  A data frame with columns `gene_name` (Ensembl gene IDs), `logfc`, and
  `pvalue`. Rows with a missing or blank `gene_name` are ignored, and
  only the first row of each gene is used.

- subnetwork:

  An igraph knowledge graph. Default [`tkoi::tkoi_net`](tkoi_net.md).

- pvalue_threshold:

  Keep genes with `pvalue` at or below this value. Default `0.05`.

- logfc_threshold:

  Keep genes with `abs(logfc)` at or above this value. Default `0.25`.

- indirect_link_threshold:

  Nodes within two hops of at least this many seed genes are ranked
  first in each result table. This only orders the tables; no node is
  excluded. Default `3`.

- topology_similarity:

  Number in `[0, 1]`. A replacement gene's degree must lie within
  `[topology_similarity * d, (2 - topology_similarity) * d]`, where `d`
  is the seed gene's degree in [`tkoi::genes`](genes.md) (its degree in
  `tkoi_net`, also for custom networks). Default `0.9`.

- n_permutation:

  Number of permutations in the null (at least 2). Default `100`.

- damping_factor:

  PageRank damping factor in `(0, 1)`. Default `0.85`.

- maximum_iteration:

  Maximum solver iterations per PageRank vector. The solver needs about
  55 iterations on `tkoi_net` at the default `tolerance`. Default `500`.

- n_cores:

  Number of CPU threads. `NULL` (the default) uses every available core,
  or `getOption("tkoi.n_cores")` when that is set. Requests above the
  number of available cores (which respects Slurm and Linux cgroup CPU
  limits) are capped. Results do not depend on this value.

- keep_permutations:

  If `TRUE` (the default), `pagerank_data` holds every permutation's
  PageRank vector (`perm.1`, `perm.2`, ...). If `FALSE`, it holds only
  the null mean and standard deviation (`null_mean`, `null_sd`), which
  saves memory for large `n_permutation`.

- tolerance:

  Relative residual at which the PageRank solver stops. A vector is
  accepted when its residual is below `tolerance` both in the symmetric
  system the solver works on and, in 1-norm, in the original PageRank
  system, which bounds the 1-norm error of the PageRank vector by about
  `2 * tolerance`. The default `1e-14` is at least as accurate as
  [`igraph::page_rank()`](https://r.igraph.org/reference/page_rank.html)
  on `tkoi_net`, per node as well as overall.

- verbose:

  If `TRUE` (the default), print progress messages.

## Value

An S4 object of class `tKOIList` with slots:

- `expression_data`:

  The input data frame.

- `pagerank_data`:

  A data frame with `node_id`, observed `pagerank`, and either every
  permutation (`perm.1`, ...) or `null_mean` and `null_sd`.

- `network_summary_statistics`:

  A named list of tibbles, one per node type, with `node_id`,
  `node_type`, `identifier`, `pagerank`, `beta`, `p_value`, `fdr`,
  `direct_links`, `indirect_links`, `valency`, and annotation columns
  (missing for nodes without a curated annotation). Every node of the
  network is listed once, except Compound nodes that are not human
  metabolites.

## Details

1.  **Gene filtering:** genes are kept when `pvalue <= pvalue_threshold`
    and `abs(logfc) >= logfc_threshold`.

2.  **Network mapping:** kept genes are mapped onto the network through
    their Ensembl IDs (see [`genes`](genes.md)).

3.  **Personalized PageRank:** PageRank is propagated from the seed
    genes, each weighted by `abs(logfc)`.

4.  **Permutation null:** every seed gene is replaced by a random gene
    of similar degree, `n_permutation` times, and PageRank is recomputed
    for each replacement set.

5.  **Network enrichment scoring:** each node's z-score (`beta`) and
    one-sided p-value compare its observed PageRank with the null.

6.  **Node annotation:** nodes are annotated with curated metadata (GO
    terms, diseases, cell types, compounds, and more). Every node is
    reported; nodes without a curated annotation have missing annotation
    columns. Among compounds, only human metabolites
    ([`human_metabolites`](human_metabolites.md)) are reported.

7.  **Prioritization:** results are split by node type, adjusted for
    multiple testing (Benjamini-Hochberg FDR over the reported nodes of
    each type), and ranked.

**Speed.** The PageRank solver is written in C++. It solves the
symmetric form of the personalized PageRank system with conjugate
gradient, up to 16 PageRank vectors per pass over the network, on
`n_cores` threads. Two-hop neighborhood counts use bitsets, and the
permutation null is drawn in C++.

**Reproducibility.** The null seed sets are drawn in one sequence from
R's random number generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) fixes the result, and
the same seed gives bit-identical statistics for any `n_cores`. On
`tkoi_net`, the same seed also draws the same null seed sets as tkoi
1.0.0 run with `n_cores = 1`, and PageRank values, `beta`, and `p_value`
agree with 1.0.0 to numerical precision. Result tables differ from 1.0.0
by design: unannotated nodes are kept, the Compound FDR family is the
reported human metabolites, and nodes with a constant null are
untestable (see below).

**Custom networks.** Seed genes that are not vertices of `subnetwork`
are dropped before the null is drawn, so the observed run and every null
run use the same seeds with the same weights. Replacement genes are
drawn only from [`genes`](genes.md) that are vertices of `subnetwork`.
Nodes with a missing type are grouped as `"Unknown"`, and nodes of types
without a curated annotation are reported without annotation columns.

**Memory.** Beyond `tkoi_net` itself (about 0.5 GB), the network matrix
takes about 0.3 GB and the solver about 0.6 GB. Keeping all permutations
takes `8 * vcount(subnetwork) * (n_permutation + 1)` bytes (about 0.75
GB for 100 permutations on `tkoi_net`). Before any PageRank is computed,
the run checks this against the available memory (physical memory, or a
Linux container limit) and stops with advice if it would not fit.

**Untestable nodes.** When a node's null PageRank is the same in every
permutation (standard deviation 0), its z-score is undefined, so its
`beta`, `p_value`, and `fdr` are `NaN`; such nodes are left out of the
FDR adjustment and listed last. This covers nodes that no seed gene or
replacement gene reaches, and seed genes in small disconnected
components that no replacement gene reaches (which would otherwise get
an infinite `beta`). `direct_links` and `indirect_links` are `NA` for
nodes that have no seed gene within one or two hops. A gene whose first
row in `expression_data` has a missing `pvalue` or `logfc` is not used
as a seed.

## See also

[`visualize_topn`](visualize_topn.md),
[`export_gene_exploration_data`](export_gene_exploration_data.md),
[`run_gene_enrichment`](run_gene_enrichment.md)

## Examples

``` r
if (FALSE) { # \dontrun{
expression_data = data.table::fread(
  system.file("extdata", "example_data.csv", package = "tkoi")
)

set.seed(1)
result = run_tkoi(
  expression_data = expression_data,
  n_permutation = 100
)

result@network_summary_statistics$BiologicalProcess
} # }
```
