# Retrieve Neighboring Nodes in a Knowledge Graph

Returns every node within `degree_expansion` hops of a node in the tKOI
knowledge graph, including the node itself.

## Usage

``` r
get_neighboring_nodes(
  gene_node_id,
  degree_expansion,
  subnetwork = tkoi::tkoi_net
)
```

## Arguments

- gene_node_id:

  A node ID (a vertex name of `subnetwork`), for example a gene's `id`
  in [`genes`](genes.md).

- degree_expansion:

  Number of hops to expand. `1` returns the node and its direct
  neighbors; `2` also returns their neighbors.

- subnetwork:

  An igraph object to search. Default [`tkoi::tkoi_net`](tkoi_net.md).

## Value

A character vector of node IDs.

## Details

The search uses
[`igraph::ego()`](https://r.igraph.org/reference/ego.html) and ignores
edge direction.

## See also

[`ego`](https://r.igraph.org/reference/ego.html),
[`plot_network`](plot_network.md)

## Examples

``` r
if (FALSE) { # \dontrun{
gene_id = tkoi::genes$id[tkoi::genes$name == "TP53"]

# The gene and its direct neighbors
get_neighboring_nodes(gene_id, degree_expansion = 1)

# Everything within two hops
get_neighboring_nodes(gene_id, degree_expansion = 2)
} # }
```
