# Plot a Local Network Around a Target Node from tKOI Results

Draws the part of the knowledge graph that links a target node (for
example an enriched GO term) to the significant genes near it. Genes are
colored by log fold change (blue down, white zero, red up), the target
is orange, other nodes are gray, and node size follows the tKOI effect
size (`beta`).

## Usage

``` r
plot_network(
  tkoi_result,
  target_node_id,
  degree_expansion = 2,
  network_layout_type = c("kk", "fr", "gem", "graphopt", "lgl", "mds"),
  subnetwork = tkoi::tkoi_net
)
```

## Arguments

- tkoi_result:

  A `tKOIList` returned by [`run_tkoi`](run_tkoi.md).

- target_node_id:

  Node ID (vertex name) of the node to center on.

- degree_expansion:

  Maximum number of hops between the target and a gene. Default `2`.

- network_layout_type:

  Layout algorithm: `"kk"` (Kamada-Kawai, the default), `"fr"`
  (Fruchterman-Reingold), `"gem"`, `"graphopt"`, `"lgl"`, or `"mds"`.

- subnetwork:

  The igraph network used for the analysis. Default
  [`tkoi::tkoi_net`](tkoi_net.md).

## Value

The plotted igraph subgraph, invisibly.

## Details

Significant genes pass the p-value and log fold change thresholds stored
in `tkoi_result`. The plot shows every node on a simple path of at most
`degree_expansion` edges between the target and one of these genes. For
`degree_expansion <= 2` these nodes are found directly from neighbor
sets; longer paths use
[`igraph::all_simple_paths()`](https://r.igraph.org/reference/all_simple_paths.html),
which can be slow around highly connected nodes.

## Examples

``` r
if (FALSE) { # \dontrun{
top_term = tkoi_result@network_summary_statistics$BiologicalProcess$node_id[1]
plot_network(tkoi_result, target_node_id = top_term, degree_expansion = 2)
} # }
```
