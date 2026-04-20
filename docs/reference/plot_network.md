# Plot a Local Network Around a Target Node from tKOI Results

This function visualizes a subnetwork of biologically relevant nodes
around a specified target node from a tKOI result. It includes all
significant genes (based on p-value and logFC thresholds) within a
user-defined network neighborhood, as well as all connecting paths
between the target and these genes. The resulting network is colored by
log fold change (logFC), sized by beta value, and plotted with an igraph
layout.

## Usage

``` r
plot_network(
  tkoi_result,
  target_node_id,
  degree_expansion = 2,
  network_layout_type = "kk"
)
```

## Arguments

- tkoi_result:

  An object of class `tkoi_result` containing differential expression
  results, thresholds, and network information from a tKOI analysis.

- target_node_id:

  A character string specifying the node ID (e.g., GO term) to center
  the network around.

- degree_expansion:

  Integer specifying the maximum graph distance from the target node to
  explore. Defaults to 2.

- network_layout_type:

  A character string indicating the layout algorithm to use for
  plotting. Options include: `"kk"` (Kamada-Kawai), `"fr"`
  (Fruchterman-Reingold), `"gem"` (Graph Embedding), `"graphopt"` (Graph
  Optimization), `"lgl"` (Large Graph Layout), and `"mds"`
  (Multidimensional Scaling).

## Value

A network plot is rendered directly using `plot.igraph()`.

## Details

The function first identifies differentially expressed genes that meet
user-defined thresholds and are within the specified graph neighborhood
of the target node. It then computes all simple paths between the target
node and each significant gene, compiles a list of nodes involved, and
builds a subgraph.

Each node is colored by logFC (gradient), or gray if not a gene. The
target node is highlighted in orange. Node sizes are scaled to beta
values. The layout algorithm used can be selected with
`network_layout_type`.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_network(tkoi_result = my_tkoi_result,
             target_node_id = "4:c77f6410-bc08-43ba-a172-0503ab1c93db:1234567",
             degree_expansion = 2,
             network_layout_type = "kk")
} # }
```
