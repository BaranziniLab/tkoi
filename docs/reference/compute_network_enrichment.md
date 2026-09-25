# Compute Z-Score and P-Value for Node Pagerank

This function calculates the z-score and p-value for a node's pagerank
by comparing it to a set of permutation values. It returns a data frame
containing the z-score (`beta`) and p-value (`p_value`).

## Usage

``` r
compute_network_enrichment(node)
```

## Arguments

- node:

  Either a named list / one-row data frame with `pagerank` and
  permutation values whose names start with `"perm"` (e.g. `perm.1`,
  `perm.2`), or a data frame with one row per node in that layout, such
  as the `pagerank_data` slot of a [`run_tkoi()`](run_tkoi.md) result
  with `keep_permutations = TRUE`. A list of equal-length vectors is
  treated as such a data frame. `node` needs a `pagerank` value and at
  least two `perm*` values, or an error is raised.

## Value

A data frame with one row per node and columns:

- `beta`: The computed z-score.

- `p_value`: The one-tailed p-value derived from the z-score.

## Details

The function calculates the z-score as: \$\$z = \frac{\text{pagerank} -
\text{mean}(\text{perm_values})}{\text{sd}(\text{perm_values})}\$\$
(`NaN` when the permutation values have no spread, as in
[`run_tkoi`](run_tkoi.md)), and the p-value is calculated as the
survival function of the z-score: \$\$p = 1 - \Phi(z)\$\$ where \\\Phi\\
is the cumulative distribution function of the standard normal
distribution.

## Examples

``` r
compute_network_enrichment(list(pagerank = 0.3, perm.1 = 0.1, perm.2 = 0.2, perm.3 = 0.15))
#>   beta     p_value
#> 1    3 0.001349898

nodes = data.frame(pagerank = c(0.3, 0.1), perm.1 = c(0.1, 0.1), perm.2 = c(0.2, 0.12))
compute_network_enrichment(nodes)
#>         beta    p_value
#> 1  2.1213203 0.01694743
#> 2 -0.7071068 0.76024994
```
