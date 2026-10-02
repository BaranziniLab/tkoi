# Retrieve the Graph Used for a tKOI Analysis

Returns the exact igraph object supplied to
[`run_tkoi()`](https://baranzinilab.github.io/tkoi/reference/run_tkoi.md),
including its vertex and edge attributes. Save the result with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) to retain this
network for later contextualization or agent traversal.

## Usage

``` r
get_analysis_graph(tkoi_result)
```

## Arguments

- tkoi_result:

  A `tKOIList` returned by
  [`run_tkoi()`](https://baranzinilab.github.io/tkoi/reference/run_tkoi.md).

## Value

The stored igraph object.

## Details

Results created before tkoi 1.3.0 may not contain an analysis graph.
This function refuses those results instead of substituting the
currently installed `tkoi_net`, which may differ from the network
actually analyzed. Rerun the analysis with its original graph to create
a complete result.

## Examples

``` r
if (FALSE) { # \dontrun{
result = run_tkoi(expression_data, subnetwork = tkoi::tkoi_net)
saveRDS(result, "analysis.rds")
graph = get_analysis_graph(readRDS("analysis.rds"))
get_neighboring_nodes(result@pagerank_data$node_id[1], 1, subnetwork = graph)
} # }
```
