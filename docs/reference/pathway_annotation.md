# Pathway Annotations (WikiPathways)

A curated annotation table mapping biological pathway identifiers to
pathway names. These annotations are derived primarily from WikiPathways
and describe metabolic, signaling, and disease-related pathways relevant
to cellular and physiological processes. This dataset is used in the
tKOI framework to annotate nodes of type "Pathway" in biological
knowledge graphs.

## Usage

``` r
pathway_annotation
```

## Format

A data frame with 4831 rows and 2 columns:

- identifier:

  A character string representing the pathway ID, typically from
  WikiPathways (e.g., "WP2875").

- name:

  The human-readable name of the pathway (e.g., "Constitutive androstane
  receptor pathway").

## Source

WikiPathways <https://www.wikipathways.org/>

## Details

Pathways describe sets of molecular interactions and regulatory events
that carry out defined biological functions. This dataset enables
biological interpretation and functional enrichment of graph nodes
associated with pathway membership or influence.

The `identifier` column is used to map pathway node IDs from graph
analyses (e.g., [`run_tkoi`](run_tkoi.md)), while the `name` column
provides context for labeling and visualization.

## See also

[`reaction_annotation`](reaction_annotation.md),
[`disease_annotation`](disease_annotation.md), [`run_tkoi`](run_tkoi.md)

## Examples

``` r
data(pathway_annotation)
subset(pathway_annotation, grepl("thyroid", name, ignore.case = TRUE))
#>      identifier
#> 5        WP1981
#> 384      WP3859
#> 586    PWY-6261
#> 2270   PWY-6260
#> 2531     WP3872
#> 2864   PWY-6241
#> 3885     WP4746
#> 4260     WP4928
#> 4380     WP2032
#>                                                                           name
#> 5                                         Thyroxine thyroid hormone production
#> 384  TGF beta signaling in thyroid cells for epithelial mesenchymal transition
#> 586         thyroid hormone metabolism II (via conjugation and/or degradation)
#> 2270                           thyroid hormone metabolism I (via deiodination)
#> 2531            Regulation of apoptosis by parathyroid hormone related protein
#> 2864                                              thyroid hormone biosynthesis
#> 3885   Thyroid hormones production and peripheral downstream signaling effects
#> 4260                                 MAPK pathway in congenital thyroid cancer
#> 4380                         Thyroid stimulating hormone TSH signaling pathway
```
