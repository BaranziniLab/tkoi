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

  A named list or data structure containing:

  - `pagerank`: A numeric value representing the pagerank of the node.

  - Permutation values: One or more numeric values with names matching
    the pattern `"^perm.*"` (e.g., `perm1`, `perm2`), used for
    calculating the z-score.

## Value

A data frame with the following columns:

- `beta`: The computed z-score.

- `p_value`: The one-tailed p-value derived from the z-score.

## Details

The function calculates the z-score as: \$\$z = \frac{\text{pagerank} -
\text{mean}(\text{perm_values})}{\text{sd}(\text{perm_values})}\$\$ and
the p-value is calculated as the survival function of the z-score: \$\$p
= 1 - \Phi(z)\$\$ where \\\Phi\\ is the cumulative distribution function
of the standard normal distribution.
