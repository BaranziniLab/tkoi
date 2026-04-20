# Retrieve Neighboring Nodes in a Knowledge Graph

This function returns all neighboring nodes for a given gene node ID in
the tKOI knowledge graph. The neighborhood is defined by the number of
hops (`degree_expansion`) from the input node.

## Usage

``` r
get_neighboring_nodes(gene_node_id, degree_expansion, subnetwork)
```

## Arguments

- gene_node_id:

  A character string representing the unique node ID of a gene in the
  [`tkoi::tkoi_net`](tkoi_net.md) igraph object.

- degree_expansion:

  An integer indicating the radius (in number of hops) to expand the
  search for neighbors. For example, `degree_expansion = 1` returns
  first-degree (direct) neighbors, while `degree_expansion = 2` includes
  both direct and second-degree (indirect) neighbors.

- subnetwork:

  An `igraph` object (typically a subgraph of `tkoi_net`) in which the
  neighborhood search will be performed.

## Value

A character vector of node IDs representing the neighbors of the
specified gene node.

## Details

The function leverages the
[`igraph::ego()`](https://r.igraph.org/reference/ego.html) function to
extract the subgraph containing all nodes within a specified path length
from the target node. It then returns the node IDs (as character vector)
of these neighbors. This is useful in the tKOI pipeline for assessing
local network topology around significantly altered genes.

## See also

[`ego`](https://r.igraph.org/reference/ego.html),
[`run_tkoi`](run_tkoi.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Get direct neighbors of a gene node with ID "ENSG00000123456"
get_neighboring_nodes("ENSG00000123456", degree_expansion = 1)

# Get both direct and indirect (2-hop) neighbors
get_neighboring_nodes("ENSG00000123456", degree_expansion = 2)
} # }
```
