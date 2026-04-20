# Visualize Top Network Enrichment Statistics

This function generates a bar plot to visualize the top `n` network
enrichment results for a specified category. The plot highlights the
effect size (`beta`) and the false discovery rate (`-log10(FDR)`),
enabling quick assessment of the most enriched terms or entities in the
network.

## Usage

``` r
visualize_topn(
  tkoi_list,
  category = "BiologicalProcess",
  top_n = 25,
  ranknorm = FALSE,
  lognorm = TRUE,
  high_color = "#FF5733",
  low_color = "#154360"
)
```

## Arguments

- tkoi_list:

  A `tKOIList` object containing network summary statistics.

- category:

  A character string specifying the category to visualize. Default is
  `"BiologicalProcess"`. Accepted values include:

  - `"Anatomy"`

  - `"BiologicalProcess"`

  - `"CellType"`

  - `"CellularComponent"`

  - `"ClinicalLab"`

  - `"Complex"`

  - `"Compound"`

  - `"Disease"`

  - `"EC"`

  - `"Gene"`

  - `"MiRNA"`

  - `"MolecularFunction"`

  - `"Pathway"`

  - `"Protein"`

  - `"ProteinDomain"`

  - `"ProteinFamily"`

  - `"PwGroup"`

  - `"Reaction"`.

- top_n:

  An integer specifying the number of top results to display. Default is
  25.

- ranknorm:

  A logical value indicating whether to apply inverse rank-based
  normalization to the `beta` values. Default is `FALSE`.

- lognorm:

  A logical value indicating whether to apply log base-2 transformation
  to the `beta` values. Default is `TRUE`.

- high_color:

  A character string specifying the high value color for the gradient
  scale representing `-log10(FDR)`. Default is `"#FF5733"`.

- low_color:

  A character string specifying the low value color for the gradient
  scale representing `-log10(FDR)`. Default is `"#154360"`.

## Value

A `ggplot` object representing the bar plot of top network enrichment
statistics.

## Details

The function performs the following steps:

1.  Extracts the top `n` entries from the specified category within the
    `network_summary_statistics` slot of the `tKOIList` object.

2.  Filters out entries with missing `beta` values.

3.  Optionally applies transformations to the `beta` values:

    - If `ranknorm` is set to `TRUE`, applies inverse rank-based
      normalization.

    - If `lognorm` is set to `TRUE`, applies log base-2 transformation.

    - If both transformations are enabled, `ranknorm` is overridden and
      set to `FALSE`.

4.  Adjusts very small `fdr` values to avoid extreme values in plotting.

5.  If the category is `"Gene"`, the function joins additional metadata
    from the [`tkoi::genes`](genes.md) table to obtain gene names. For
    other categories, the `name` column is used as the identifier.

6.  Filters out entries with missing identifiers and ensures consistent
    ordering of identifiers for plotting.

7.  Creates a horizontal bar plot where:

    - The x-axis represents the (possibly transformed) effect size
      (`beta`).

    - The y-axis represents the identifiers (e.g., terms or entities).

    - The color gradient of the bars represents `-log10(FDR)`, with
      customizable low and high colors.

## See also

[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html),
[`arrange`](https://dplyr.tidyverse.org/reference/arrange.html),
[`filter`](https://dplyr.tidyverse.org/reference/filter.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Visualize the top 10 Biological Processes with default transformations and colors
plt = visualize_topn(tkoi_list, category = "BiologicalProcess", top_n = 10)
print(plt)

# Visualize the top 20 Genes with custom color gradient
plt = visualize_topn(tkoi_list, category = "Gene", top_n = 20, high_color = "#E74C3C", low_color = "#3498DB")
print(plt)

# Visualize the top 15 Pathways without rank-based normalization or log transformation
plt = visualize_topn(tkoi_list, category = "Pathway", top_n = 15, ranknorm = FALSE, lognorm = FALSE)
print(plt)

# Visualize the top 5 Diseases with rank-based normalization enabled
plt = visualize_topn(tkoi_list, category = "Disease", top_n = 5, ranknorm = TRUE, lognorm = FALSE)
print(plt)
} # }
```
